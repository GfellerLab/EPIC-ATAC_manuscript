process Inputs_processing {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  time '2h'
  memory = { task.exitStatus == 137 ? 10.GB * task.attempt : 20.GB }
  publishDir path: params.benchmarking_output_path + '/intermediate_files/', mode: 'copy'

	input:
  tuple val(key), val(bulk_data), val(bulk_name)

	output:
  path "${bulk_name}/${bulk_name}*_epic_bulk_input.xlsx", emit: bulk_input
	path "${bulk_name}/${bulk_name}_deconPeaker_bulk_input.txt", emit: deconPeaker_bulk
  path "${bulk_name}/${bulk_name}_deconPeaker_ref_input.txt", emit: deconPeaker_ref
  path "${bulk_name}/${bulk_name}_cibersort_mixture.txt", emit: cibersort_bulk
  path "${bulk_name}/${bulk_name}_CIBERSORTx_ref_input.txt", emit: cibersortx_ref_ourMarkers
  path "${bulk_name}/*proportions.txt", emit: bulk_proportions
  path "${bulk_name}/*OM_deconPeaker_bulk_input.txt", emit: DeconPeaker_OM_inputs
  path "*.rda"

  tag "${bulk_name}"
	script:
    """
    Rscript ${projectDir}/bin/benchmarking/save_EPIC_ATAC_ref.r ${params.profile_path}
  	Rscript ${projectDir}/bin/benchmarking/process_bulks.r ${bulk_data} ${bulk_name} ${params.deconPeaker_signature} ${params.withSubtypes}
    mkdir -p ${bulk_name}/
    mv *OM_deconPeaker_bulk_input.txt *proportions.txt ${bulk_name}/
    cp ${bulk_name}_epic_bulk_input.xlsx ${bulk_name}_deconPeaker_* ${bulk_name}_cibersort_mixture.txt ${bulk_name}_CIBERSORTx_ref_input.txt ${bulk_name}/
    """
}

// Run R deconvolution methods
process R_methods {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  time '2h'
  memory = { task.exitStatus == 137 ? 50.GB * task.attempt : 50.GB }
  publishDir path: params.benchmarking_output_path+ '/R_methods/', mode: 'copy', pattern: "*_predictions*"
  publishDir path: params.benchmarking_output_path+ '/time/', mode: 'copy', pattern: "*-time.txt"


	input:
  tuple val(bulk_name), val(tool)
  path bulk_data
  path bulk_data_initial

	output:
  path "*_predictions*", emit: epic_res
  path "*-time.txt", emit: benchmarking_r

  tag "${bulk_name}-${tool}"
	script:
    """
    Rscript ${projectDir}/bin/benchmarking/save_EPIC_ATAC_ref.r ${params.profile_path}
    withSubtypes=${params.withSubtypes}
    if [ \${withSubtypes} = true ]; then
      (time Rscript ${projectDir}/bin/benchmarking/run_r_methods.r ${bulk_name}_epic_bulk_input.xlsx ${bulk_name} TRUE ${tool} ) 2>  ${bulk_name}-${tool}-time.txt
    else
      (time Rscript ${projectDir}/bin/benchmarking/run_r_methods.r ${bulk_name}_epic_bulk_input.xlsx ${bulk_name} FALSE ${tool} ) 2>  ${bulk_name}-${tool}-time.txt
    fi
    """
}

// Run of DeconPeaker on the peak matrix with original markers
process DeconPeaker_original_marker {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  publishDir path: params.benchmarking_output_path+ '/DeconPeaker_original_markers/', mode: 'copy'
  publishDir path: params.benchmarking_output_path+ '/time/', mode: 'copy', pattern: "*-time.txt"

	input:
  val bulk_name
  path deconpeaker_data

	output:
  path "deconvolution/*deconPeaker-Results*", emit: deconPeaker_originalMarker_res
  path "*-time.txt", emit: benchmarking_deconpeaker_OM

  tag "${bulk_name}-deconpeaker_OM"

	script:
    """
    out_path=\$(pwd)
    cd /DeconPeaker
  	(time deconPeaker.py deconvolution --lib-strategy=ATAC-Seq --mixture=${params.benchmarking_output_path}/intermediate_files/${bulk_name}/${bulk_name}_OM_deconPeaker_bulk_input.txt --pure=${params.deconPeaker_signature} --format=TABLE --pvalue=FALSE --outdir=\${out_path} ) 2>  \${out_path}/${bulk_name}-deconpeaker_OM-time.txt
    mv \${out_path}/deconvolution/deconPeaker-Results.xls \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_original_markers.xls
    mv \${out_path}/deconvolution/deconPeaker-Results.png \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_original_markers.png
    """
}

// Run DeconPeaker using our markers
process DeconPeaker_ourMarkers {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  time '2h'
  memory = { task.exitStatus == 137 ? 10.GB * task.attempt : 20.GB }
  publishDir path: params.benchmarking_output_path + '/DeconPeaker/', mode: 'copy'
  publishDir path: params.benchmarking_output_path+ '/time/', mode: 'copy', pattern: "*-time.txt"

	input:
  val bulk_name
  path ref
  path bulk

	output:
	path "deconvolution/*deconPeaker-Results*", emit: deconPeaker_res
  path "*-time.txt", emit: benchmarking_DeconPeaker_ourMarkers

  tag "${bulk_name}-DeconPeaker_ourMarkers"
	script:
    """
    out_path=\$(pwd)
    cd /DeconPeaker
  	(time deconPeaker.py deconvolution --lib-strategy=ATAC-Seq --mixture=\${out_path}/${bulk_name}_deconPeaker_bulk_input.txt \\
     --pure=\${out_path}/${bulk_name}_deconPeaker_ref_input.txt --format=TABLE --pvalue=FALSE --outdir=\${out_path} ) 2>  \${out_path}/${bulk_name}-DeconPeaker_ourMarkers-time.txt
     mv \${out_path}/deconvolution/deconPeaker-Results.xls \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_epic_markers.xls
     mv \${out_path}/deconvolution/deconPeaker-Results.png \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_epic_markers.png
    """
}

// Run DeconPeaker using the markers identified with DeconPeaker from the pure ATAC-Seq reference samples collected in this work
process DeconPeaker_ourRef {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  time '2h'
  memory = { task.exitStatus == 137 ? 10.GB * task.attempt : 20.GB }
  publishDir path: params.benchmarking_output_path+ '/DeconPeaker_newRef/', mode: 'copy'
  publishDir path: params.benchmarking_output_path+ '/time/', mode: 'copy', pattern: "*-time.txt"

	input:
  val bulk_name
  path bulk

	output:
	path "deconvolution/*deconPeaker-Results*", emit: deconPeaker_newRef_res
  path "*-time.txt", emit: benchmarking_DeconPeaker_ourRef

  tag "${bulk_name}-DeconPeaker_ourRef"

	script:
    """
    out_path=\$(pwd)
    cd /DeconPeaker

    if [[ ${bulk_name} = *"PBMC"* ]]; then
      (time deconPeaker.py deconvolution --lib-strategy=ATAC-Seq --mixture=\${out_path}/${bulk_name}_deconPeaker_bulk_input.txt \\
        --pure=${params.deconPeaker_PBMCsignature_ourRef} --format=TABLE --pvalue=FALSE --outdir=\${out_path} ) 2>  \${out_path}/${bulk_name}-DeconPeaker_ourRef-time.txt
        mv \${out_path}/deconvolution/deconPeaker-Results.xls \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_newRef.xls
        mv \${out_path}/deconvolution/deconPeaker-Results.png \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_newRef.png
    else
      (time deconPeaker.py deconvolution --lib-strategy=ATAC-Seq --mixture=\${out_path}/${bulk_name}_deconPeaker_bulk_input.txt \\
       --pure=${params.deconPeaker_signature_ourRef} --format=TABLE --pvalue=FALSE --outdir=\${out_path} ) 2>  \${out_path}/${bulk_name}-DeconPeaker_ourRef-time.txt
       mv \${out_path}/deconvolution/deconPeaker-Results.xls \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_newRef.xls
       mv \${out_path}/deconvolution/deconPeaker-Results.png \${out_path}/deconvolution/${bulk_name}_deconPeaker-Results_newRef.png
    fi

    """
}

// Run CIBERSORTx using the markers identified with CIBERSORTx from the pure ATAC-Seq reference samples collected in this work
process CIBERSORTx_ourRef_runs {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  time '1h'
  publishDir path: params.benchmarking_output_path + 'CIBERSORTx_newRef/', mode: 'copy'
  publishDir path: params.benchmarking_output_path + '/time/', mode: 'copy', pattern: "*-time.txt"

	input:
  val bulk_name
  path ref_signature
  path PBMCref_signature
  path mixture

  output:
  path "CIBERSORTx_newRef/*", emit: cibersortx_newRef_res
  path "*-time.txt", emit: benchmarking_CIBERSORTx_ourRef

  tag "${bulk_name}-CIBERSORTx_ourRef"
	script:
    """
    token=\$(cat ${params.cibersortx_token})
    mkdir CIBERSORTx_newRef

    withSubtypes=${params.withSubtypes}
    if [ \${withSubtypes} = true ]; then
      ref_suffix="ref_withSubtypes.txt"
    else
      ref_suffix="ref.txt"
    fi

    if [[ ${bulk_name} = *"PBMC"* ]]; then
      (time singularity exec --no-home -c -B ./:/src/data -B CIBERSORTx_newRef/:/src/outdir ${params.CIBERSORTx_singularity} /src/CIBERSORTxFractions --username ${params.CIBERSORTx_username} --token \${token} --sigmatrix /src/data/CIBERSORTx_PBMC_\${ref_suffix} --mixture ${bulk_name}_cibersort_mixture.txt --outdir CIBERSORTx_newRef ) 2>  ${bulk_name}-CIBERSORTx_ourRef-time.txt
      mv CIBERSORTx_newRef/CIBERSORTx_Results.txt CIBERSORTx_newRef/${bulk_name}_CIBERSORTx_noAbs_Results_newRef.txt
    else
      (time singularity exec --no-home -c -B ./:/src/data -B CIBERSORTx_newRef/:/src/outdir ${params.CIBERSORTx_singularity} /src/CIBERSORTxFractions --username ${params.CIBERSORTx_username} --token \${token} --sigmatrix /src/data/CIBERSORTx_TME_\${ref_suffix} --mixture ${bulk_name}_cibersort_mixture.txt --outdir CIBERSORTx_newRef ) 2>  ${bulk_name}-CIBERSORTx_ourRef-time.txt
      mv CIBERSORTx_newRef/CIBERSORTx_Results.txt CIBERSORTx_newRef/${bulk_name}_CIBERSORTx_noAbs_Results_newRef.txt
    fi
    """
}

// Run CIBERSORTx using our markers
process CIBERSORTx_ourMarkers_runs {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  time '1h'
  publishDir path: params.benchmarking_output_path + 'CIBERSORTx_ourMarkers/', mode: 'copy'
  publishDir path: params.benchmarking_output_path+ '/time/', mode: 'copy', pattern: "*-time.txt"

	input:
  val bulk_name
  path ref
  path mixture

	output:
  path "CIBERSORTx_ourMarkers/*", emit: cibersortx_ourMarkers_res
  path "*-time.txt", emit: benchmarking_CIBERSORTx_ourMarkers

  tag "${bulk_name}-CIBERSORTx_ourMarkers"

	script:
    """
    token=\$(cat ${params.cibersortx_token})
    mkdir CIBERSORTx_ourMarkers

    (time singularity exec --no-home -c -B ./:/src/data -B CIBERSORTx_ourMarkers/:/src/outdir ${params.CIBERSORTx_singularity} /src/CIBERSORTxFractions --username ${params.CIBERSORTx_username} --token \${token} --sigmatrix /src/data/${bulk_name}_CIBERSORTx_ref_input.txt --mixture ${bulk_name}_cibersort_mixture.txt --outdir CIBERSORTx_ourMarkers ) 2>  ${bulk_name}-CIBERSORTx_ourMarkers-time.txt
    mv CIBERSORTx_ourMarkers/CIBERSORTx_Results.txt CIBERSORTx_ourMarkers/${bulk_name}_CIBERSORTx_noAbs_Results_ourMarkers.txt
    """
}

// Run the deconvolution on the gene activity matrix
process GA_deconvolution {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  publishDir path: params.benchmarking_output_path + '/GA_deconvolution/', mode: 'copy'

  input:
  val bulk_data

  output:
  path "*GA_predictions*", emit: epic_ga_res

	script:
    """
    Rscript ${projectDir}/bin/benchmarking/run_ga_deconvolution.R ${bulk_data}
    """
}

// Run the deconvolution on the gene expression matrix
process RNA_deconvolution {
  container 'docker://agabriel/epicatac_main_analyses:v1.0'
  publishDir path: params.benchmarking_output_path + '/RNA_deconvolution/', mode: 'copy'

  input:
  val bulk_data

  output:
  path "*_RNA_deconvolution_summary.txt", emit: epic_rna_res

	script:
    """
    Rscript ${projectDir}/bin/benchmarking/run_rna_deconvolution.r ${bulk_data}
    """
}
