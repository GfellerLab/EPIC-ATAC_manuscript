README
================
Aurélie Gabriel
2023-10-11

This readme describes how the code available in the scripts/ folder was
used. In the following sections we will refer to a data/ folder
containing files that were too large to be hosted on the github
repository. The complete data/ folder can be retrieved on zenodo
(<https://zenodo.org/record/8431792>, additional_data.zip).

## Preprocessing of the bulk ATAC-Seq data

The following command line was used to process the ATAC-Seq samples
listed in *data/samples/Ref/SRA_samples_metadata.txt*. Note that the
git_repo_path variable corresponds to the path to this github
repository. The genome_folder_path variable corresponds to the folder
containing the refgenie genome assets which are PEPATAC prerequisites
for ATAC-Seq data processing, please refer to the PEPATAC documentation
to download these data (<http://pepatac.databio.org/en/latest/assets/>).
sra_output_folder corresponds to the folder where SRA data will be
downloaded.

``` bash
samples_metadata=data/samples/Ref/SRA_samples_metadata.txt  
genome_folder_path=user_defined_path 
pepatac_output_folder=user_defined_path 
sra_output_folder=user_defined_path 

nextflow run scripts/bulk_preprocessing_scripts/ATAC_processing.nf \
-c scripts/bulk_preprocessing_scripts/nextflow.config -profile singularity \
--samples_metadata ${samples_metadata} \
--sra_folder sra_output_folder/ \   
--output_folder ${pepatac_output_folder} \
--genome_version hg38 \
--genome_folder ${genome_folder_path} \
--multiqc_config scripts/bulk_preprocessing_scripts/multiqc_config.yaml \
--blacklist_file scripts/bulk_preprocessing_scripts/hg38-blacklist.v2.bed \
--original_config scripts/bulk_preprocessing_scripts/pepatac_original.yaml \
--samples_origin SRA --adapters scripts/bulk_preprocessing_scripts/all_adapters_PE.fa \
--preprocess true --getcounts false -bg -resume
```

The same pipeline was used to process the ENCODE samples listed in
*data/samples/ENCODE/SRA_samples_metadata.txt*

## Definition of the set of consensus peaks and counts extraction

The same pipeline can then be used to identify a set of consensus peaks
across studies and cell-types and extract the raw counts for each peak
and sample with the following options `--preprocess false` and
`--getcounts true`.

``` bash
metadata_file=data/samples/Ref/SRA_samples_metadata.txt  
raw_counts_output_folder=user_defined_path 

nextflow run scripts/bulk_preprocessing_scripts/ATAC_processing.nf \
-c scripts/bulk_preprocessing_scripts/nextflow.config -profile singularity \
--samples_metadata ${metadata_file} \
--output_folder ${pepatac_output_folder} \
--genome_version hg38 --peak_score_thr 2 \
--genome_folder ${genome_folder_path} \
--output_folder_counts ${raw_counts_output_folder} \
--pepatac_config scripts/bulk_preprocessing_scripts/pepatac_config.yaml \
--preprocess false --getcounts true -bg -resume
```

Raw counts for each peaks and sample and the metadata associated to each
sample are provided as output in *raw_counts_output_folder*
(raw_counts.txt and metadata.txt). These two files are used later in
“Computing reference profiles and identifying marker peaks” and are also
provided in the data/ folder.

We extracted the raw counts from the ENCODE data for the peaks (listed
in the .saf file) identified in the reference samples using the
following command line. The resulting counts and the encode metadata are
located in *data/encode_counts.txt* and *data/encode_metadata.txt*.

``` bash
nextflow run scripts/bulk_preprocessing_scripts/extract_peaks_counts.nf \
-c scripts/bulk_preprocessing_scripts/nextflow.config \
-profile singularity \
--bam_folder bam_files_ENCODE/ \
--output_folder_counts ${raw_counts_output_folder}ENCODE/ \
--saf_file ${raw_counts_output_folder}all_consensusPeaks.saf -bg -resume

# bam_files_ENCODE/ contains all bam files obtained from the preprocessing of the ENCODE ATAC-Seq data
```

## Computing reference profiles and identifying marker peaks

Markers are identified in 10 subsets of the reference samples (list of
samples in each subset located in *data/samples/Ref/*) and a consensus
of these markers are retrieved. The following command line will generate
reference profiles and identify cell-type specific marker peaks for
EPIC-ATAC. It will also use CIBERSORTx to build reference profiles. To
do so a token which can be retrieved on CIBERSORTx website should be
obtained and provided to the cibersortx_token parameter.

``` bash
# Without considering Tcell subtypes  
ref_output_path=user_defined_path 

nextflow run scripts/reference_profiles_scripts/build_references.nf \
-c scripts/reference_profiles_scripts/nextflow.config -profile singularity \
--counts_path ${raw_counts_output_folder}raw_counts.txt \
--metadata_path ${raw_counts_output_folder}metadata.txt \
--output_folder ${ref_output_path} \
--cibersortx_token CIBERSORTx_token \
--cross_validation_files data/samples/Ref/subsample \
--encode_count_file data/encode_counts.txt \
--encode_metadata_file data/encode_metadata.txt \
--TCGA_path data/TCGA_data.rds \
--HA_path data/Human_Atlas_peaks.txt \
-bg -resume

# Considering Tcell subtypes
ref_output_path=user_defined_path 

nextflow run scripts/reference_profiles_scripts/build_references.nf \
-c scripts/reference_profiles_scripts/nextflow.config -profile singularity \
--counts_path ${raw_counts_output_folder}raw_counts.txt \
--metadata_path ${raw_counts_output_folder}metadata.txt \
--output_folder ${ref_output_path} \
--cibersortx_token CIBERSORTx_token \
--cross_validation_files data/samples/Ref/subsample \
--encode_count_file data/encode_counts.txt \
--encode_metadata_file data/encode_metadata.txt \
--TCGA_path data/TCGA_data.rds \
--HA_path data/Human_Atlas_peaks.txt \
--major_groups FALSE \
-bg -resume
```
