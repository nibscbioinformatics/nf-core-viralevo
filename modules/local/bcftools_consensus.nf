// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BCFTOOLS_CONSENSUS {
    tag "$vcf"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::bcftools=1.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bcftools:1.11--h7c999a4_0"
    } else {
        container "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
    }

    input:
    tuple path(vcf), path(csi)
    path fasta

    output:
    path "*.fasta", emit: fasta
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def filename = "$vcf".tokenize('_')[0]
    def caller   = ("$vcf".contains("_ivar")) ? "ivar" :  ("$vcf".contains("lofreq")) ? "lofreq" : ''
    """
    cat $fasta | bcftools consensus \\
        $options.args \\
        $vcf > ${filename}_${caller}.fasta

    echo \$(bcftools --version 2>&1) | sed 's/^.*bcftools //; s/ .*\$//' > ${software}.version.txt
    """
}
