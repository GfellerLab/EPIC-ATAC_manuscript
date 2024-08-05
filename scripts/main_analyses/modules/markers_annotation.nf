
process Markers_annotation {
  memory '50G'
  time '3h'
  container 'docker://agabriel/peak_annotation:latest'
  publishDir path: params.output_path + "markers_annotation/", mode: 'copy'
  publishDir path: params.benchmarking_output_path+ '/manuscript_figures/', mode: 'copy', pattern: ".pdf"

  input:
  path reference

	output:
  path 'PBMC_noSubtypes/*.txt', emit: PBMC_annotations_noSubtypes, optional: true
  path 'TME_noSubtypes/*.txt', emit: TME_annotations_noSubtypes, optional: true
  path 'PBMC_withSubtypes/*.txt', emit: PBMC_annotations_withSubtypes, optional: true
  path 'TME_withSubtypes/*.txt', emit: TME_annotations_withSubtypes, optional: true
  path '*pdf', emit: annot_figure

	script:
    """
    export BFC_CACHE=/tmp/R_cache/BiocFileCache
  	Rscript ${projectDir}/bin/markers_annotation/markers_annotation.r ${reference} ${params.remap2022} ${params.jaspar2022}
    """
}
