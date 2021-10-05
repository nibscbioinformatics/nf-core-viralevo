include { saveFiles } from './functions'

params.options = [:]

process MAKEVARTABLE {
    tag "Creating Variant Table for R.."
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
    path "variant-table.csv", emit: csv
    path ("*_filtered.vcf"), emit: filteredvars

    script:
    """
    python $projectDir/bin/tablefromvcf.py . variant-table.csv ${alt_depth_threshold} ${vaf_threshold}
    """
}
