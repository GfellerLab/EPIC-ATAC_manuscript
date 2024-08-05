#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Pseudobulks parameters
ATAC_dataset_names = Channel.from(["PBMC_Satpathy","BCC_Satpathy","PBMC_Granja","Gynecological_cancers","PBMC_multiome","PBMC_experiment"])
multiome_dataset_names = Channel.from(["PBMC_multiome_simulated","HNSCC","BRCA","CRC","CESC","SKCM","PDAC", "OV"])
all_datasets = ATAC_dataset_names.concat(multiome_dataset_names)

// Benchmarking parameters:
params.data_path = null
params.output_path = null
params.cibersortx_token = null
params.CIBERSORTx_singularity = null
params.CIBERSORTx_username = null
params.bulk_path = params.data_path + "pseudobulks/"
params.withSubtypes = null
params.benchmarking_output_path = params.withSubtypes.toBoolean() ? params.output_path + "benchmarking_withSubtypes/" : params.output_path + "benchmarking/"

// References:
params.profile_path = params.data_path + "reference_profiles/"
params.deconPeaker_signature = params.profile_path + "GSE74912_Corces_MR_pure_readcounts_signature_matrix.xls"
params.deconPeaker_signature_ourRef = params.withSubtypes.toBoolean() ? params.profile_path + "deconpeaker_TME_ref_withSubtypes.xls" : params.profile_path + "deconpeaker_TME_ref.xls"
params.deconPeaker_PBMCsignature_ourRef = params.withSubtypes.toBoolean() ? params.profile_path + "deconpeaker_PBMC_ref_withSubtypes.xls" : params.profile_path + "deconpeaker_PBMC_ref.xls"
params.CIBERSORTx_signature_ourRef = params.withSubtypes.toBoolean() ? params.profile_path + "CIBERSORTx_TME_ref_withSubtypes.txt" : params.profile_path + "CIBERSORTx_TME_ref.txt"
params.CIBERSORTx_PBMCsignature_ourRef = params.withSubtypes.toBoolean() ? params.profile_path + "CIBERSORTx_PBMC_ref_withSubtypes.txt" : params.profile_path + "CIBERSORTx_PBMC_ref.txt"

// Annotation files:
params.remap2022 = params.data_path + "annotation_files/remap2022_nr_macs2_hg38_v1_0.bed.gz"
params.jaspar2022 = params.data_path + "annotation_files/JASPAR2022.sqlite"
params.rna_gtf_file = params.data_path + "annotation_files/gencode.v32.primary_assembly.annotation.gtf.gz"

//  Include modules to run benchmarking
include { Inputs_processing } from './modules/benchmarking'
include { R_methods } from './modules/benchmarking'
include { DeconPeaker_original_marker } from './modules/benchmarking'
include { DeconPeaker_ourMarkers } from './modules/benchmarking'
include { DeconPeaker_ourRef } from './modules/benchmarking'
include { CIBERSORTx_ourRef_runs } from './modules/benchmarking'
include { CIBERSORTx_ourMarkers_runs } from './modules/benchmarking'
include { GA_deconvolution } from './modules/benchmarking'
include { RNA_deconvolution } from './modules/benchmarking'

//  Include modules to generate the Figures
include { Gather_predictions } from './modules/figures'
include { Gather_predictions_withSubtypes } from './modules/figures'
include { Benchmarking_figures } from './modules/figures'
include { ref_PBMC_Tumor_profiles_figures } from './modules/figures'
include { true_bulk_applications } from './modules/figures'
include { annotation_figures } from './modules/figures'

// Include the modules related to the markers annotation
include { Markers_annotation } from './modules/markers_annotation'

// Workflow definition
workflow {

    // Markers annotations
    profiles_file = Channel.fromPath(params.profile_path + '*_profile_*.Rdata')
    // profiles_file.view { "profiles_files: $it" }
    annotation_res = Markers_annotation(profiles_file)
    all_annotations_ch = Channel
    .empty()
    .concat(annotation_res.TME_annotations_withSubtypes.ifEmpty { Channel.empty() })
    .concat(annotation_res.PBMC_annotations_withSubtypes.ifEmpty { Channel.empty() })
    .concat(annotation_res.TME_annotations_noSubtypes.ifEmpty { Channel.empty() })
    .concat(annotation_res.PBMC_annotations_noSubtypes.ifEmpty { Channel.empty() })
    annotation_figures(all_annotations_ch.first().collect())

    all_bulks = Channel.fromPath(params.bulk_path + '*_bulk.Rdata')

    // Benchmarking
    all_bulks_key = all_bulks.map { file ->
      def key = file.getName().toString().replaceAll("_bulk.Rdata", "")
      return [key, file]
    }
    dataset_list = all_datasets.map {[it, it]}
    matched_channels = all_bulks_key.join(dataset_list, by: 0)
    // matched_channels.view { "Dataset: $it" }

    deconv_inputs = Inputs_processing(matched_channels)

    res_r_methods = R_methods(all_datasets.combine(Channel.from(["EPIC", "quantiSeq", "MPC_counter", "ABIS"])),
                              deconv_inputs.bulk_input.collect(),
                              all_bulks.collect())

    res_deconpeaker_OM = DeconPeaker_original_marker(all_datasets, deconv_inputs.DeconPeaker_OM_inputs.collect())
    res_deconpeaker_ourMarkers = DeconPeaker_ourMarkers(all_datasets, deconv_inputs.deconPeaker_ref.collect(), deconv_inputs.deconPeaker_bulk.collect())
    res_deconpeaker_ourRef = DeconPeaker_ourRef(all_datasets, deconv_inputs.deconPeaker_bulk.collect())

    res_CIBERSORTx_ourRef = CIBERSORTx_ourRef_runs(all_datasets,
                                                   Channel.fromPath(params.CIBERSORTx_signature_ourRef).collect(),
                                                   Channel.fromPath(params.CIBERSORTx_PBMCsignature_ourRef).collect(),
                                                   deconv_inputs.cibersort_bulk.collect())
    res_CIBERSORTx_ourMarkers = CIBERSORTx_ourMarkers_runs(all_datasets,
                                                           deconv_inputs.cibersortx_ref_ourMarkers.collect(),
                                                           deconv_inputs.cibersort_bulk.collect())

    if (!params.withSubtypes) {
      // Run Deconvolution on RNA and GA samples
      all_ga_bulks = Channel.fromPath(params.bulk_path + '*_bulk_GA.rdata')
      all_rna_bulks = Channel.fromPath(params.bulk_path + '*_bulk_RNA.rdata')
      res_GA_deconv = GA_deconvolution(all_ga_bulks)
      res_RNA_deconv = RNA_deconvolution(all_rna_bulks)

      // Generate figures
      all_predictions = Gather_predictions(res_r_methods.epic_res.collect(),
                                           res_deconpeaker_ourMarkers.deconPeaker_res.collect(),
                                           res_deconpeaker_ourRef.deconPeaker_newRef_res.collect(),
                                           res_CIBERSORTx_ourRef.cibersortx_newRef_res.collect(),
                                           res_CIBERSORTx_ourMarkers.cibersortx_ourMarkers_res.collect(),
                                           res_GA_deconv.collect())

      // Figures related to the evaluation of the method
      Benchmarking_figures(all_predictions.predictions_summary.collect())

      // Figure 2 panels
      ref_PBMC_Tumor_profiles_figures()

      // Run the application of EPIC-ATAC on the breast cancer and on the PBMC true bulk datasets
      true_bulk_applications()
    } else {
      // Generate figures
      all_predictions = Gather_predictions_withSubtypes(res_r_methods.epic_res.collect(),
                                                        res_deconpeaker_ourMarkers.deconPeaker_res.collect(),
                                                        res_deconpeaker_ourRef.deconPeaker_newRef_res.collect(),
                                                        res_CIBERSORTx_ourRef.cibersortx_newRef_res.collect(),
                                                        res_CIBERSORTx_ourMarkers.cibersortx_ourMarkers_res.collect())

      // Figures related to the evaluation of the method
      Benchmarking_figures(all_predictions.predictions_summary.collect())

      // refeference profiles heatmaps
      ref_PBMC_Tumor_profiles_figures()
    }

}
