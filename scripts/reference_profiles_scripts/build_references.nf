#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

params.help = null

log.info ""
log.info "-----------------------------------------------------------------------------------"
log.info "                      Building profiles for deconvolution                          "
log.info "-----------------------------------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run build_references.nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --metadata_path       STRING   Metadata associated to the reference samples.'
    log.info '    --counts_path       STRING     Raw counts associated to each peak and sample.'
    log.info '    --cibersortx_token       STRING     CIBERSORTx token to be able to generate CIBERSORTx profiles.'
    log.info '    --CIBERSORTx_singularity       STRING     CIBERSORTx singularity image to be able to run CIBERSORTx.'
    log.info '    --CIBERSORTx_username       STRING     CIBERSORTx username associated to the CIBERSORTx token.'
    log.info '    --encode_count_file       STRING     Raw counts associated to each peak and ENCODE sample.'
    log.info '    --encode_metadata_file       STRING     Metadata associated to the ENCODE samples.'
    log.info '    --TCGA_path       STRING     Path to TCGA data for markers filtering.'
    log.info '    --HA_path       STRING     Path to human atlas module regions for markers filtering.'
    log.info ""
    log.info "Optional arguments:"
    log.info '    --output_folder    STRING            Output folder (default: .).'
    log.info '    --cross_validation_files    STRING           Folder name pattern containing samples subsets is provided'
    log.info '    --major_groups    STRING            Bolean indicaing whether T cell subtypes should be considered'

    log.info ''
    exit 0
}

params.metadata_path = null
params.counts_path = null

params.output_folder = "./"
params.cibersortx_token = null
params.CIBERSORTx_singularity = null
params.CIBERSORTx_username = null
params.cross_validation_files = "none"
params.major_groups = "TRUE"
params.encode_count_file = null
params.encode_metadata_file = null
params.TCGA_path = null
params.HA_path = null

if ( params.cross_validation_files == "none" ) {
  samples_selection = Channel.from(0)
}else{
  // to perform markers selection on multiple folds of the collected data
  samples_selection = Channel.from(1..10)
}

process perform_DA {
  time '3h'
  memory '50G'
  publishDir params.output_folder + '/DA_res/', mode: 'copy'

	input:
  each samples_selection_id from samples_selection

	output:
  file 'Peaks_pairwise_*' into DA_res

	shell:
    '''
    if [ !{samples_selection_id} -eq "0" ]; then
        cv_file="none"
    else
        cv_file=!{params.cross_validation_files}!{samples_selection_id}/SRA_samples_metadata.txt
    fi
  	Rscript !{baseDir}/bin/1_get_markers.r !{params.counts_path} !{params.metadata_path} !{params.major_groups} !{params.encode_count_file} !{params.encode_metadata_file} ${cv_file} !{samples_selection_id}
    '''
}

if ( params.cross_validation_files == "none" ) {
  process select_markers_all {
    publishDir params.output_folder + "/markers/", mode: 'copy'

    time '2h'
    memory = { task.exitStatus == 137 ? 10.GB * task.attempt : 20.GB }
  	input:
    each DA_tables from DA_res.flatten()

  	output:
    file '*markers*' into markers_selection

  	shell:
      '''
    	Rscript !{baseDir}/bin/2_process_markers.r !{DA_tables}
      '''
  }
}else{
  process select_markers_cv {
    publishDir params.output_folder + "/markers/", mode: 'copy'

    time '2h'
    memory = { task.exitStatus == 137 ? 10.GB * task.attempt : 20.GB }
  	input:
    each DA_tables from DA_res.flatten()

  	output:
    file '*markers*' into markers_selection_per_cv

  	shell:
      '''
    	Rscript !{baseDir}/bin/2_process_markers.r !{DA_tables}
      '''
  }

  process get_consensus_markers {
    publishDir params.output_folder + "/markers/", mode: 'copy'

  	input:
    file marker_tables from markers_selection_per_cv.collect()

  	output:
    file '*consensus*' into markers_selection

  	shell:
      '''
    	Rscript !{baseDir}/bin/2_consensus_markers.r
      '''
  }
}


process save_profiles {
  publishDir params.output_folder + "/profiles/", mode: 'copy'

  input:
  each markers_list from markers_selection.flatten()

	output:
  file 'profile*' into profile_matrices

	shell:
    '''
  	Rscript !{baseDir}/bin/3_save_profiles.r !{params.counts_path} !{params.metadata_path} !{markers_list} !{params.major_groups}
    '''
}

process markers_filtering {
  time '1h'
  publishDir params.output_folder + "/profiles/", mode: 'copy'

  input:
  each profile from profile_matrices.flatten()

	output:
  file '*profile.Rdata' into filtered_profile
	shell:
    '''
  	Rscript !{baseDir}/bin/4_markers_filtering.r !{profile} !{params.TCGA_path} !{params.HA_path}
    '''
}

process generate_external_softwares_input {
  memory '20 GB'
  publishDir params.output_folder + "/profiles/", mode: 'copy'

  output:
  file 'deconpeaker_input_ref.txt' into deconpeaker_input_ref
  file 'deconpeaker_input_pheno.txt' into deconpeaker_input_pheno
  file 'cibersortx_input_ref.txt' into cibersortx_input_ref
  file 'cibersortx_input_pheno.txt' into cibersortx_input_pheno
  file 'deconpeaker_input_PBMCref.txt' into deconpeaker_input_PBMCref
  file 'deconpeaker_input_PBMCpheno.txt' into deconpeaker_input_PBMCpheno
  file 'cibersortx_input_PBMCref.txt' into cibersortx_input_PBMCref
  file 'cibersortx_input_PBMCpheno.txt' into cibersortx_input_PBMCpheno

	shell:
    '''
  	Rscript !{baseDir}/bin/cibersort_deconpeaker_inputs.r !{params.counts_path} !{params.metadata_path} !{params.major_groups}
    '''
}

process save_Deconpeaker_profiles {
  publishDir params.output_folder + '/deconpeaker/', mode: 'copy'
  cpus 5
  memory '50 GB'
  time '3h'
  input:
  file pheno_file from deconpeaker_input_pheno.collect()
  file ref_file from deconpeaker_input_ref.collect()
  file PBMCpheno_file from deconpeaker_input_PBMCpheno.collect()
  file PBMCref_file from deconpeaker_input_PBMCref.collect()
	output:
  file 'findctsps/*'

	shell:
    '''
    deconPeaker.py findctsps --lib-strategy=ATAC-Seq --profile=deconpeaker_input_ref.txt --phenotype=deconpeaker_input_pheno.txt --outdir=./
    deconPeaker.py findctsps --lib-strategy=ATAC-Seq --profile=deconpeaker_input_PBMCref.txt --phenotype=deconpeaker_input_PBMCpheno.txt --outdir=./
    '''
}

process save_CIBERSORT_profiles {
  publishDir params.output_folder, mode: 'copy'
  time '3h'
  input:
  file pheno_file from cibersortx_input_ref.collect()
  file ref_file from cibersortx_input_pheno.collect()
  file PBMCpheno_file from cibersortx_input_PBMCref.collect()
  file PBMCref_file from cibersortx_input_PBMCpheno.collect()

	output:
  file 'cibersortx/*'
  file 'cibersortx_PBMC/*'

	shell:
    '''
    token=$(cat !{params.cibersortx_token})
    mkdir cibersortx
    singularity exec --no-home -c -B ./:/src/data -B cibersortx/:/src/outdir !{params.CIBERSORTx_singularity} /src/CIBERSORTxFractions --username !{params.CIBERSORTx_username} --token ${token} --refsample /src/data/cibersortx_input_ref.txt --phenoclasses /src/data/cibersortx_input_pheno.txt --outdir cibersortx
    mkdir cibersortx_PBMC
    singularity exec --no-home -c -B ./:/src/data -B cibersortx_PBMC/:/src/outdir !{params.CIBERSORTx_singularity} /src/CIBERSORTxFractions --username !{params.CIBERSORTx_username} --token ${token} --refsample /src/data/cibersortx_input_PBMCref.txt --phenoclasses /src/data/cibersortx_input_PBMCpheno.txt --outdir cibersortx_PBMC
    '''
}
