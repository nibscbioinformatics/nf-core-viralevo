// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TRIMLOG {
    tag 'Summarising Cutadapt logs'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'tsv2vcf', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path trim_logs

    output:
    path 'trimming-summary.csv', emit: summary

    script: // This script is bundled with the pipeline, in nf-core-viralevo/bin.. provide directory as input
    """
    python $baseDir/bin/logger.py . trimming-summary.csv cutadapt
    """
}
