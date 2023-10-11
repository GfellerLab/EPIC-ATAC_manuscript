#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

params.help = null

log.info ""
log.info "-----------------------------------------------------------------------------------"
log.info "                     ATAC-Seq bulk preprocessing using PEPATAC                     "
log.info "-----------------------------------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run ATAC_processing.nf [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments if internal data are preprocessed:"
    log.info '    --input_files          STRING          Path to the folder containing the fastq files to process'
    log.info '    --adapters   STRING          Path to the file listing adapters for the adapter trimming step'
    log.info '    --genome_folder   STRING          Path to the refgenie genome folder'
    log.info '    --multiqc_config   STRING          Path to the multiQC configuration yaml file'
    log.info '    --blacklist_file   STRING          Path to the file listing the blacklist regions'
    log.info '    --original_config   STRING          Path to the pepatac config file'
    log.info "Mandatory arguments if SRA data are preprocessed:"
    log.info '    --samples_metadata          STRING          Path to the file listing SRR and/or EGA IDs'
    log.info '    --sra_folder       STRING          Path to the ncbi folder, where the sra tools will download the SRA files'
    log.info '    --adapters   STRING          Path to the file listing adapters for the adapter trimming step'
    log.info '    --genome_folder   STRING          Path to the refgenie genome folder'
    log.info '    --multiqc_config   STRING          Path to the multiQC configuration yaml file'
    log.info '    --blacklist_file   STRING          Path to the file listing the blacklist regions'
    log.info '    --original_config   STRING          Path to the pepatac config file'
    log.info "Mandatory arguments if EGA data are preprocessed:"
    log.info '    --ega_config     STRING          Path to the EGA identification file, this parameter is need when samples_origin is set to EGA'
    log.info '    --EGA_folder     STRING          Path to the folder where the EGA files will be downloaded, this parameter is need when samples_origin is set to EGA'
    log.info '    --adapters   STRING          Path to the file listing adapters for the adapter trimming step'
    log.info '    --genome_folder   STRING          Path to the refgenie genome folder'
    log.info '    --multiqc_config   STRING          Path to the multiQC configuration yaml file'
    log.info '    --blacklist_file   STRING          Path to the file listing the blacklist regions'
    log.info '    --original_config   STRING          Path to the pepatac config file'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --samples_origin    STRING            String spefifying whether the samples listed in samples_metadata.txt needs to be retrieved from SRA (default: SRA, other values: EGA or internal).'
    log.info '    --genome_version    STRING            Genome version (default: hg38).'
    log.info '    --output_folder    STRING            Output folder (default: .).'
    log.info '    --organism    STRING            Organism (default: human).'
    log.info ''
    exit 0
}


params.samples_metadata = null
params.samples_origin = "SRA"
params.ega_config = null
params.sra_folder = null
params.EGA_folder = null
params.input_files = null
params.output_folder = null
params.genome_version = "hg38"
params.organism = "human"
params.adapters=null
params.genome_folder = null
params.multiqc_config = null
params.blacklist_file = null
params.original_config = null
params.get_fragments_files = "false"

atac_wdl = file("/atac-seq-pipeline-2.0.0/atac.wdl")
params.chrom_sizes = params.genome_folder + "alias/"+ params.genome_version +"/fasta/default/"+ params.genome_version +".chrom.sizes"
params.TSS_bed = params.genome_folder + "alias/"+ params.genome_version +"/refgene_anno/default/"+ params.genome_version +"_TSS.bed"
if(params.genome_version == "hg38"){
  rCRSd = Channel.fromPath( params.genome_folder+'alias/rCRSd/bowtie2_index/default/rCRSd*' ).collect()
  human_repeats = Channel.fromPath( params.genome_folder+'alias/human_repeats/bowtie2_index/default/human_repeats*' ).collect()
}else if(params.genome_version == "hg19"){
  rCRSd = Channel.fromPath( params.genome_folder+'alias/rCRSd/bowtie2_index/default/rCRSd*' ).collect()
  human_repeats = Channel.fromPath( params.genome_folder+'alias/human_repeats/bowtie2_index/default/human_repeats*' ).collect()
}else{
  mouse_chrM2x = Channel.fromPath( params.genome_folder+'alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x*' ).collect()
}
genome_index = Channel.fromPath( params.genome_folder+'alias/'+ params.genome_version +'/bowtie2_index/default/*' ).collect()


params.pepatac_config = null
params.peak_score_thr = 0
params.output_folder_counts = null
params.preprocess = true
params.getcounts = false

if(params.preprocess){
  if(params.samples_origin == "SRA" ){
    process processing_samples_metadata {
      input:
      file samples_metadata from Channel.fromPath(params.samples_metadata)

    	output:
      file 'SRR*'  optional true into download_files,download_files2 mode flatten
      file 'EGA*'  optional true into download_ega_files mode flatten
    	shell:
        '''
        Rscript !{baseDir}/bin/generate_download_files.r !{params.samples_metadata} !{params.output_folder}
        '''
    }

    process SRA_download {
      cache 'lenient'
      time '15h'
    	input:
    	file download_files

    	output:
      file '*.end_file.txt' optional true into downloaded_samples mode flatten

    	shell:
        '''
        line=$(cat !{download_files})
        echo ${line}

        ls !{params.sra_folder}

        #FILE=!{params.sra_folder}${line}_1.fastq
        FILE=!{params.output_folder}${line}/PEPATAC_commands.sh
        echo $FILE
        if [ -f "$FILE" ]; then
            echo "already downloaded"
        else
            prefetch -v ${line} -O !{params.sra_folder} --max-size 50000000000
            fasterq-dump -f -O !{params.sra_folder} !{params.sra_folder}/${line}/${line}.sra
            vdb-validate -x !{params.sra_folder}/${line}/${line}.sra
            echo ${line} > ${line}.end_file.txt
        fi

        '''
    }

  } else if (params.samples_origin == "EGA" ) {
    process processing_EGA_table {
      input:
      file samples_metadata from Channel.fromPath(params.samples_metadata)

    	output:
      file 'SRR*'  optional true into download_files,download_files2 mode flatten
      file 'EGA*'  optional true into download_ega_files mode flatten
    	shell:
        '''
        Rscript !{baseDir}/bin/generate_download_json.r !{params.samples_metadata} !{params.output_folder}
        '''
    }
  }else {
    fastq1_files=Channel.fromPath( params.input_files+'/*R1*fastq.gz' ).collect()
    fastq2_files=Channel.fromPath( params.input_files+'/*R2*fastq.gz' ).collect()
    fastqPairs = Channel.fromFilePairs( params.input_files+'/*R{1,2}.fastq.gz' )
  }

  if ( params.samples_origin == "SRA" ) {
    process run_PEPATAC_workflow_sra {
      cache 'lenient'
      // config infos
      cpus 10
      memory '100 GB'
      time '30h'
    	input:
      file downloaded_samples
      file 'genome_index/*' from genome_index
      file config from Channel.fromPath(params.original_config).collect()

    	output:
       path "$sampleID/fastqc/*zip" into fastqc_res
       path "$sampleID/aligned_${params.genome_version}/*_sort.bam.{flagstat,idxstats,stats}" into filtered_bam_flagstat_mqc
       path "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam.{flagstat,idxstats,stats}" into dedup_bam_flagstat_mqc
       path "$sampleID/QC_${params.genome_version}/*_preseq_yield.txt" into preseq_mqc
       path "$sampleID/peak_calling_${params.genome_version}/*_peaks.xls" into macs2_mqc
       path "$sampleID/peak_calling_${params.genome_version}/*_peaks_normalized.narrowPeak" into macs2_res
       path "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam*" into bam_files,bam_files2
       path "$sampleID/prealignments/*zip" into fq_filtered_files
       path "$sampleID/PEPATAC*" into pepatac_log
       path "$sampleID/*.tsv" into tsv_files

       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/QC_${params.genome_version}/*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/PEPATAC*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/*.tsv"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/peak_calling_${params.genome_version}/*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam*"

    	shell:
        sampleID = downloaded_samples.getSimpleName()
        '''
        sampleID=$(cat !{downloaded_samples})

        echo ${sampleID}
        pathAdapters=!{params.adapters}
        sed 's|adapters: null|adapters: '"${pathAdapters}"'|' pepatac_original.yaml  > pepatac.yaml
        path=$PWD
        genome=!{params.genome_version}
        if [ $genome = "hg38" -o $genome = "hg19" ]; then
          prealignment_params="--prealignment-index rCRSd=!{params.genome_folder}alias/rCRSd/bowtie2_index/default/rCRSd human_repeats=!{params.genome_folder}alias/human_repeats/bowtie2_index/default/human_repeats"
        fi
        if [ $genome = "mm10" ]; then
          prealignment_params="--prealignment-index mouse_chrM2x=!{params.genome_folder}alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x"
        fi
        pepatac.py --single-or-paired paired -C $PWD/pepatac.yaml \
        ${prealignment_params} \
        --genome !{params.genome_version} --genome-index genome_index/. \
        --chrom-sizes !{params.chrom_sizes} \
        --sample-name ${sampleID} --input !{params.sra_folder}${sampleID}_1.fastq --input2 !{params.sra_folder}${sampleID}_2.fastq  -O . --peak-caller macs2 \
        --trimmer trimmomatic \
        --aligner bowtie2 \
        --deduplicator samtools \
        -P 10 \
        --TSS-name !{params.TSS_bed} \
        --blacklist !{params.blacklist_file}

        prefix=${sampleID}/aligned_!{params.genome_version}/${sampleID}_sort
        samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
        samtools idxstats ${prefix}.bam > ${prefix}.bam.idxstats
        samtools stats ${prefix}.bam > ${prefix}.bam.stats

        prefix=${sampleID}/aligned_!{params.genome_version}/${sampleID}_sort_dedup
        samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
        samtools idxstats ${prefix}.bam > ${prefix}.bam.idxstats
        samtools stats ${prefix}.bam > ${prefix}.bam.stats

        prefix=${sampleID}/prealignments/
        if [ $genome = "hg38" -o $genome = "hg19" ]; then
          fastqc --noextract --outdir ${prefix} ${prefix}*human_repeats_unmap_R1.fq.gz
          fastqc --noextract --outdir ${prefix} ${prefix}*human_repeats_unmap_R2.fq.gz
          rm -r ${sampleID}/prealignments/rCRSd_bt2
        fi
        if [ $genome = "mm10" ]; then
          fastqc --noextract --outdir ${prefix} ${prefix}*mouse_chrM2x_unmap_R1.fq.gz
          fastqc --noextract --outdir ${prefix} ${prefix}*mouse_chrM2x_unmap_R2.fq.gz
          rm -r ${sampleID}/prealignments/mouse_chrM2x_bt2
        fi

        frag_param=!{params.get_fragments_files}
        if [ $frag_param = "true" ]; then
          # get fragments files from the bam file:
          prefix=${sampleID}/aligned_!{params.genome_version}/${sampleID}_sort_dedup
          sinto fragments -b ${prefix}.bam -f ${sampleID}_fragments.tsv -p 10 -t RG
          sort -k 1,1 -k2,2n ${sampleID}_fragments.tsv > !{params.output_folder}fragments_files/${sampleID}_sorted_fragments.tsv
          bgzip !{params.output_folder}fragments_files/${sampleID}_sorted_fragments.tsv
          tabix -p bed !{params.output_folder}fragments_files/${sampleID}_sorted_fragments.tsv.gz
          rm ${sampleID}_fragments.tsv
        fi


        mkdir -p !{params.output_folder}/${sampleID}/peak_calling_all/
        cp -r ${sampleID}/peak_calling_!{params.genome_version}/* !{params.output_folder}/${sampleID}/peak_calling_all/
        rm -r ${sampleID}/prealignments/*.fq.gz
        rm -r ${sampleID}/fastq ${sampleID}/raw
        rm -r !{params.sra_folder}${sampleID}*
        '''
    }
  } else if (params.samples_origin == "EGA" ) {
    process run_PEPATAC_workflow_ega {
      // config infos
      cpus 10
      memory '50 GB'
      time '7h'
      input:
      file download_ega_files
      file 'genome_index/*' from genome_index
      file config from Channel.fromPath(params.original_config).collect()

      output:
       path "$sampleID/fastqc/*zip" into fastqc_res
       path "$sampleID/aligned_${params.genome_version}/*_sort.bam.{flagstat,idxstats,stats}" into filtered_bam_flagstat_mqc
       path "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam.{flagstat,idxstats,stats}" into dedup_bam_flagstat_mqc
       path "$sampleID/QC_${params.genome_version}/*_preseq_yield.txt" into preseq_mqc
       path "$sampleID/peak_calling_${params.genome_version}/*_peaks.xls" into macs2_mqc
       path "$sampleID/peak_calling_${params.genome_version}/*_peaks_normalized.narrowPeak" into macs2_res
       path "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam*" into bam_files,bam_files2
       path "$sampleID/prealignments/*zip" into fq_filtered_files
       path "$sampleID/PEPATAC*" into pepatac_log
       path "$sampleID/*.tsv" into tsv_files

       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/QC_${params.genome_version}/*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/fastqc/*zip"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/PEPATAC*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/*.tsv"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/peak_calling_${params.genome_version}/*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam*"

      shell:
        sampleID = download_ega_files.baseName
        '''
        index=0
        while read line
        do
           myarr[$index]=$line
           index=$(($index+1))
        done < !{download_ega_files}

        ID1="${myarr[0]}"
        ID2="${myarr[1]}"
        file=!{download_ega_files}
        sampleID="${file%.*}"
        echo ${sampleID}
        path=$PWD

        cd !{params.EGA_folder}

        pyega3 -cf !{params.ega_config} -c 10 -ms 1000000000 fetch ${ID1} --max-retries 3000
        pyega3 -cf !{params.ega_config} -c 10 -ms 1000000000 fetch ${ID2} --max-retries 3000

        cd ${path}

        pathAdapters=!{params.adapters}
        sed 's|adapters: null|adapters: '"${pathAdapters}"'|' pepatac_original.yaml  > pepatac.yaml

        genome=!{params.genome_version}
        if [ $genome = "hg38" -o $genome = "hg19" ]; then
          prealignment_params="--prealignment-index rCRSd=!{params.genome_folder}alias/rCRSd/bowtie2_index/default/rCRSd human_repeats=!{params.genome_folder}alias/human_repeats/bowtie2_index/default/human_repeats"
        fi
        if [ $genome = "mm10" ]; then
          prealignment_params="--prealignment-index mouse_chrM2x=!{params.genome_folder}alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x"
        fi

        pepatac.py --single-or-paired paired -C $PWD/pepatac.yaml \
        ${prealignment_params} \
        --genome !{params.genome_version} --genome-index genome_index/. \
        --chrom-sizes !{params.chrom_sizes} \
        --sample-name ${sampleID} --input !{params.EGA_folder}${ID1}/*_R1.fastq.gz --input2 !{params.EGA_folder}${ID2}/*_R2.fastq.gz  -O . --peak-caller macs2 \
        --trimmer trimmomatic \
        --aligner bowtie2 \
        --deduplicator samtools \
        -P 10 \
        --TSS-name !{params.TSS_bed} \
        --blacklist !{params.blacklist_file}

        prefix=${sampleID}/aligned_!{params.genome_version}/${sampleID}_sort
        samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
        samtools idxstats ${prefix}.bam > ${prefix}.bam.idxstats
        samtools stats ${prefix}.bam > ${prefix}.bam.stats

        prefix=${sampleID}/aligned_!{params.genome_version}/${sampleID}_sort_dedup
        samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
        samtools idxstats ${prefix}.bam > ${prefix}.bam.idxstats
        samtools stats ${prefix}.bam > ${prefix}.bam.stats

        prefix=${sampleID}/prealignments/
        if [ $genome = "hg38" -o $genome = "hg19" ]; then
          fastqc --noextract --outdir ${prefix} ${prefix}*human_repeats_unmap_R1.fq.gz
          fastqc --noextract --outdir ${prefix} ${prefix}*human_repeats_unmap_R2.fq.gz
          rm -r ${sampleID}/prealignments/rCRSd_bt2
        fi
        if [ $genome = "mm10" ]; then
          fastqc --noextract --outdir ${prefix} ${prefix}*mouse_chrM2x_unmap_R1.fq.gz
          fastqc --noextract --outdir ${prefix} ${prefix}*mouse_chrM2x_unmap_R2.fq.gz
          rm -r ${sampleID}/prealignments/mouse_chrM2x_bt2
        fi

        frag_param=!{params.get_fragments_files}
        if [ $frag_param = "true" ]; then
          # get fragments files from the bam file:
          prefix=${sampleID}/aligned_!{params.genome_version}/${sampleID}_sort_dedup
          sinto fragments -b ${prefix}.bam -f ${sampleID}_fragments.tsv -p 10 -t RG
          sort -k 1,1 -k2,2n ${sampleID}_fragments.tsv > !{params.output_folder}fragments_files/${sampleID}_sorted_fragments.tsv
          bgzip !{params.output_folder}fragments_files/${sampleID}_sorted_fragments.tsv
          tabix -p bed !{params.output_folder}fragments_files/${sampleID}_sorted_fragments.tsv.gz
          rm ${sampleID}_fragments.tsv
        fi

        mkdir -p !{params.output_folder}/${sampleID}/peak_calling_all/
        cp -r ${sampleID}/peak_calling_!{params.genome_version}/* !{params.output_folder}/${sampleID}/peak_calling_all/
        rm -r ${sampleID}/prealignments/*.fq.gz
        rm -r ${sampleID}/fastq ${sampleID}/raw
        '''
    }
  } else {
    process run_PEPATAC_workflow_internalFastq {
      // config infos
      cpus 30
      memory '50 GB'
      time '7h'
    	input:
      set sampleID, file(fastq_file) from fastqPairs
      file 'genome_index/*' from genome_index
      file config from Channel.fromPath(params.original_config).collect()

    	output:
    	 // file ("$sampleID/*") into sample_level_output
       path "$sampleID/fastqc/*zip" into fastqc_res
       path "$sampleID/aligned_${params.genome_version}/*_sort.bam.{flagstat,idxstats,stats}" into filtered_bam_flagstat_mqc
       path "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam.{flagstat,idxstats,stats}" into dedup_bam_flagstat_mqc
       path "$sampleID/QC_${params.genome_version}/*_preseq_yield.txt" into preseq_mqc
       path "$sampleID/peak_calling_${params.genome_version}/*_peaks.xls" into macs2_mqc
       path "$sampleID/peak_calling_${params.genome_version}/*_peaks_normalized.narrowPeak" into macs2_res
       path "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam*" into bam_files,bam_files2
       path "$sampleID/prealignments/*zip" into fq_filtered_files
       path "$sampleID/PEPATAC*" into pepatac_log
       path "$sampleID/*.tsv" into tsv_files

       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/QC_${params.genome_version}/*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/PEPATAC*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/*.tsv"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/peak_calling_${params.genome_version}/*"
       publishDir "${params.output_folder}/", mode: 'copy', pattern: "$sampleID/aligned_${params.genome_version}/*_sort_dedup.bam*"

    	shell:
        '''
        echo !{sampleID}
        ls
        pathAdapters=!{params.adapters}
        sed 's|adapters: null|adapters: '"${pathAdapters}"'|' pepatac_original.yaml  > pepatac.yaml
        genome=!{params.genome_version}
        if [ $genome = "hg38" -o $genome = "hg19" ]
        then
          prealignment_params="--prealignment-index rCRSd=!{params.genome_folder}alias/rCRSd/bowtie2_index/default/rCRSd human_repeats=!{params.genome_folder}alias/human_repeats/bowtie2_index/default/human_repeats"
        elif [ $genome = 'mm10' ]
        then
          prealignment_params="--prealignment-index mouse_chrM2x=!{params.genome_folder}alias/mouse_chrM2x/bowtie2_index/default/mouse_chrM2x"
        fi
        pepatac.py --single-or-paired paired  -C $PWD/pepatac.yaml \
        ${prealignment_params} \
        --genome !{params.genome_version} --genome-index genome_index/. \
        --chrom-sizes !{params.chrom_sizes} \
        --sample-name !{sampleID} --input !{sampleID}_R1.fastq.gz --input2 !{sampleID}_R2.fastq.gz  -O . --peak-caller macs2 \
        --trimmer trimmomatic \
        --aligner bowtie2 \
        --deduplicator samtools \
        -P 30 \
        --TSS-name !{params.TSS_bed} \
        --blacklist !{params.blacklist_file}

        prefix=!{sampleID}/aligned_!{params.genome_version}/!{sampleID}_sort
        samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
        samtools idxstats ${prefix}.bam > ${prefix}.bam.idxstats
        samtools stats ${prefix}.bam > ${prefix}.bam.stats

        prefix=!{sampleID}/aligned_!{params.genome_version}/!{sampleID}_sort_dedup
        samtools flagstat ${prefix}.bam > ${prefix}.bam.flagstat
        samtools idxstats ${prefix}.bam > ${prefix}.bam.idxstats
        samtools stats ${prefix}.bam > ${prefix}.bam.stats

        prefix=!{sampleID}/peak_calling_!{params.genome_version}/!{sampleID}
        #cat ${prefix}_peaks.narrowPeak | wc -l | awk -v OFS='\t' '{ print "!{sampleID}", \$1 }' | cat $mlib_peak_count_header - > ${prefix}_peaks.count_mqc.tsv

        prefix=!{sampleID}/prealignments/
        if [ $genome = "hg38" -o $genome = "hg19" ]
        then
          fastqc --noextract --outdir ${prefix} ${prefix}*human_repeats_unmap_R1.fq.gz
          fastqc --noextract --outdir ${prefix} ${prefix}*human_repeats_unmap_R2.fq.gz
          rm -r !{sampleID}/prealignments/rCRSd_bt2
        elif [ $genome = 'mm10' ]
        then
          fastqc --noextract --outdir ${prefix} ${prefix}*mouse_chrM2x_unmap_R1.fq.gz
          fastqc --noextract --outdir ${prefix} ${prefix}*mouse_chrM2x_unmap_R2.fq.gz
          rm -r !{sampleID}/prealignments/mouse_chrM2x_bt2
        fi

        mkdir -p !{params.output_folder}/!{sampleID}/peak_calling_all/
        cp -r !{sampleID}/peak_calling_!{params.genome_version}/* !{params.output_folder}/!{sampleID}/peak_calling_all/
        rm -r !{sampleID}/prealignments/*.fq.gz
        rm -r !{sampleID}/fastq !{sampleID}/raw

        '''
    }
  }

  process run_multiQC {
    time '3h'
  	input:
      path ('fastqc/*') from fastqc_res.collect()
      path ('alignment/filtered/*') from filtered_bam_flagstat_mqc.collect()
      path ('alignment/dedup/*') from dedup_bam_flagstat_mqc.collect()
      path ('preseq/*') from preseq_mqc.collect()
      path ('macs2/*') from macs2_mqc.collect()
      path ('fastqc_filtered/*') from fq_filtered_files.collect()
      file mqc_config from Channel.fromPath(params.multiqc_config).collect()
    output:
      file("multiqc_report.html") into multiqc_post
      file("multiqc_report_data") into multiqc_post_data

       publishDir "${params.output_folder}/QC", mode: 'copy'

    shell:
    '''
    multiqc -v . -n multiqc_report.html --comment "ATAC-Seq QC report" --config multiqc_config.yaml
    '''
  }
}


if(params.getcounts){
  process get_consensus_peaks {
    time '2h'
    memory '50 GB'
  	input:
      file data from Channel.fromPath(params.output_folder).collect()
      file SRA_info from Channel.fromPath(params.samples_metadata)
      file chr_sizes from Channel.fromPath(params.chrom_sizes)
      file pepatac_project_config from Channel.fromPath(params.pepatac_config)
      path ('macs2/*') from Channel.fromPath(params.output_folder+"*/peak_calling_"+ params.genome_version +"/*_peaks_normalized.narrowPeak").collect()
    output:
      path "summary/*_consensusPeaks.narrowPeak" into consensus_peaks_per_group
      path "summary_stats.txt" into qc_metrics_summary

      publishDir "${params.output_folder_counts}/", mode: 'copy'

    shell:
    '''
    Rscript !{baseDir}/bin/gather_stats.r results/ "SRR|EGA|Sample" !{params.samples_metadata}
    Rscript !{baseDir}/bin/get_consensus_peaks.r !{params.samples_metadata} pepatac_config.yaml !{params.genome_version}.chrom.sizes 2 summary_stats.txt !{params.peak_score_thr} !{params.organism}
    '''
  }
  process merge_consensus_peaks {
  	input:
      path peaks_files from consensus_peaks_per_group.collect()
    output:
      path "all_consensusPeaks.*" into merged_peaks

      publishDir "${params.output_folder_counts}/", mode: 'copy'

    shell:
    peaks = peaks_files.collect{it.toString()}.sort().join(' ')
    '''
    sort -T '.' -k1,1 -k2,2n !{peaks} \\
            | mergeBed > all_consensusPeaks.txt

    echo -e "GeneID\tChr\tStart\tEnd\tStrand" > all_consensusPeaks.saf
    awk -v FS='\t' -v OFS='\t' '{ print $1":"$2"-"$3, \$1, \$2, \$3,  "." }' all_consensusPeaks.txt >> all_consensusPeaks.saf
    '''
  }
  process get_counts {
    time '10h'
    cpus 5
  	input:
      // path bams from bam_files.collect()
      path bams from Channel.fromPath(params.output_folder+"*/aligned_"+ params.genome_version +"/*_sort_dedup.bam").collect()
      path saf from merged_peaks.collect()
    output:
      path 'featureCounts.txt' into ch_mlib_macs_consensus_counts
      path 'featureCounts.txt.summary' into ch_mlib_macs_consensus_counts_mqc

      publishDir "${params.output_folder_counts}/", mode: 'copy'
    shell:
    bam_files = bams.findAll{it.toString().endsWith('.bam')}.sort().join(' ')
    '''
    echo !{bam_files}
    featureCounts -F SAF -O --fracOverlap 0.2 -T 5 -p -a all_consensusPeaks.saf -o featureCounts.txt !{bam_files}
    # Warning -p only for paired end, to count fragments rather than reads
    '''
  }
  process get_filtered_count_matrix{
    container 'docker://agabriel/deconvolution_work:latest'
    input:
      path counts from ch_mlib_macs_consensus_counts.collect()
      path stats from qc_metrics_summary.collect()
    output:
      file 'raw_counts.txt' into matrix_counts
      file 'tpm.txt' into tmp
      file 'metadata.txt' into metadata

      publishDir "${params.output_folder_counts}/", mode: 'copy'
    shell:
    '''
    Rscript !{baseDir}/bin/generate_count_matrix.r featureCounts.txt !{params.samples_metadata} summary_stats.txt
    '''
  }
}
