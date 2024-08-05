

process Gather_predictions {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  publishDir path: params.benchmarking_output_path, mode: 'copy'

  input:
    path benchmarking_res_epic
    path benchmarking_res_deconpeaker_ourMarkers
    path benchmarking_res_deconpeaker_ourRef
    path benchmarking_res_cibersortx_newRef
    path benchmarking_res_cibersortx_ourRef
    path benchmarking_res_ga
	output:
    path "*txt", emit: predictions_summary

	script:
    """
    Rscript ${projectDir}/bin/figures/gather_atac_deconvolution_results.r ${params.benchmarking_output_path} ${projectDir}/bin/figures/cell_types_matching.xlsx FALSE
    Rscript ${projectDir}/bin/figures/gather_GA_results.r ${params.benchmarking_output_path} ${params.benchmarking_output_path}/GA_deconvolution/ ${projectDir}/bin/figures/cell_types_matching.xlsx
    """
}

process Gather_predictions_withSubtypes {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  publishDir path: params.benchmarking_output_path, mode: 'copy'

  input:
    path benchmarking_res_epic
    path benchmarking_res_deconpeaker_ourMarkers
    path benchmarking_res_deconpeaker_ourRef
    path benchmarking_res_cibersortx_newRef
    path benchmarking_res_cibersortx_ourRef
	output:
    path "*txt", emit: predictions_summary

	script:
    """
    Rscript ${projectDir}/bin/figures/gather_atac_deconvolution_results.r ${params.benchmarking_output_path} ${projectDir}/bin/figures/cell_types_matching_withSubtypes.xlsx TRUE
    """
}

process Benchmarking_figures {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  time '1h'
  input:
  file predictions

	script:
    """
    withSubtypes=${params.withSubtypes}
    if [ \$withSubtypes = false ]; then
      Rscript ${projectDir}/bin/figures/benchmarking_figures_revisions.r ${params.benchmarking_output_path} ${params.bulk_path} ${params.output_path}manuscript_figures/
      Rscript ${projectDir}/bin/figures/RNA_vs_ATAC_figures.r ${params.benchmarking_output_path} ${params.output_path}manuscript_figures/
    fi
    if [ \$withSubtypes = true ]; then
      Rscript ${projectDir}/bin/figures/T_cell_subtypes.r ${params.benchmarking_output_path} ${params.output_path}manuscript_figures/
    fi
    """
}

process ref_PBMC_Tumor_profiles_figures {
  memory '20G'
  time '1h'
  container 'docker://agabriel/epicatac_main_analyses:v1.0'

  script:
    """
  	Rscript ${projectDir}/bin/figures/reference_profiles_plots.r ${params.withSubtypes} ${params.data_path} ${params.output_path}manuscript_figures/
    """
}

process annotation_figures {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  input:
  file annotation_files

  script:
    """
    Rscript ${projectDir}/bin/figures/annotation_figures.r ${params.output_path}manuscript_figures/ ${params.output_path}markers_annotation/
    """
}

process true_bulk_applications {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'

	script:
    """
    Rscript ${projectDir}/bin/figures/true_bulk_application_breastCancer.r ${params.data_path}/true_bulk_data/Kumegawa_data/ ${params.output_path}manuscript_figures/
    Rscript ${projectDir}/bin/figures/true_bulk_application_PBMC.r ${params.data_path}/true_bulk_data/Morandini_data/ ${params.rna_gtf_file} ${params.output_path}manuscript_figures/
    """
}
