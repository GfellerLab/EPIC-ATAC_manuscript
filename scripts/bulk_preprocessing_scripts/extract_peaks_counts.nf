#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

params.help = null

log.info ""
log.info "-----------------------------------------------------------------------------------"
log.info "               Extract counts from bam files for a set of peaks                    "
log.info "-----------------------------------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run get_consensus_peaks.nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info '    --output_folder_counts    STRING            Output folder in which the counts file will be saved (default: .).'
    log.info '    --saf_file    STRING          Path to saf files containing the list of regions in which the number of counts should be obtained. If null, a saf file is generated from the PEPATAC output files.'
    log.info '    --bam_folder    STRING          Path to the folder containing the bam files. Should be provided only if saf file is not null, otherwise the bam files will be found in the PEPATAC outpu folder (processed_data).'
    log.info ''
    exit 0
}


params.output_folder_counts = null
params.saf_file = null
params.bam_folder = null

process get_counts {
  time '4h'
  cpus 20
  input:
    path bams from Channel.fromPath(params.bam_folder+"/*_sort_dedup.bam").collect()
  output:
    path 'featureCounts.txt' into macs_consensus_counts
    path 'featureCounts.txt.summary' into macs_consensus_counts_mqc

    publishDir "${params.output_folder_counts}/", mode: 'copy'
  shell:
  bam_files = bams.findAll{it.toString().endsWith('.bam')}.sort().join(' ')
  '''
  featureCounts -F SAF -O --fracOverlap 0.2 -T 20 -p -a !{params.saf_file} -o featureCounts.txt !{bam_files}
  # Warning -p only for paired end, to count fragments rather than reads
  '''
}
