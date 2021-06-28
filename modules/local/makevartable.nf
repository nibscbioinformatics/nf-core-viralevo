include { saveFiles } from './functions'

params.options = [:]

process MAKEVARTABLE {
    tag "$vcf"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }


    input:
    path vcf
    val alt_depth_threshold
    val vaf_threshold

    output:
    path "varianttable.csv", emit: csv

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    python $projectDir/bin/tablefromvcf.py $vcf varianttable.csv ${alt_depth_threshold} ${vaf_threshold} 
    """
}
