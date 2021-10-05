// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHANGE_FASTA_NAME {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::perl-bioperl=1.7.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/perl-bioperl:1.7.2--pl526_11"
    } else {
        container "quay.io/biocontainers/perl-bioperl:1.7.2--pl526_11"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: new_fa

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    //def caller   = ("$fasta".contains("_ivar")) ? "ivar" :  ("$fasta".contains("lofreq")) ? "lofreq" : ''
    """
    perl $baseDir/bin/change_fasta_name.pl \
    -fasta ${prefix}.consensus.fasta \
    -name ${prefix}L \
    -out ${prefix}_consensus.fa
    """
}
