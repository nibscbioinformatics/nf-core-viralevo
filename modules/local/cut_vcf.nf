// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUT_VCF {
    tag "$vcf"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    path vcf

    output:
    path("*_cutup.vcf"), emit: vcf

    script:
    def software = getSoftwareName(task.process)
    def filename = "$vcf".tokenize('_')[0]
    def caller   = ("$vcf".contains("_ivar")) ? "ivar" :  ("$vcf".contains("lofreq")) ? "lofreq" : ''
    """
    cut \\
        $options.args \\
        $vcf \\
        > ${filename}_${caller}_cutup.vcf
    """
}
