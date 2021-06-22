include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNPEFF_ANN {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::snpeff=5.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/quay.io/biocontainers/snpeff:5.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(vcf)
    val genome_version

    output:
    tuple val(meta), path("*.vcf")      , emit: vcf
    tuple val(meta), path("*.csv")      , emit: csv
    tuple val(meta), path("*.genes.txt"), emit: txt
    tuple val(meta), path("*.html")     , emit: html
    path "*.version.txt"                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def caller   = ("$vcf".contains("_ivar")) ? "ivar" :  ("$vcf".contains("lofreq")) ? "lofreq" : ''
    def avail_mem = 4
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    snpEff \\
        -Xmx${avail_mem}g \\
        ann \\
        $options.args \\
        $genome_version \\
        ${vcf} \\
        -csvStats ${prefix}_${caller}.snpeff.csv \\
        > ${prefix}_${caller}_ann.vcf
    mv snpEff_summary.html ${prefix}_${caller}.snpeff.summary.html
    echo \$(snpEff -version 2>&1) | sed 's/^.*snpEff //; s/Using.*\$//' > ${software}.version.txt
    """
}
