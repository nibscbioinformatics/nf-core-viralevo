// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TSV2VCF {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'tsv2vcf', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::perl=5.26.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/perl:5.26.2"
    } else {
        container "quay.io/biocontainers/perl:5.26.2"
    }

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf

    script: // This script is bundled with the pipeline, in nf-core-viralevo/bin/
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    perl $projectDir/bin/ivar2vcf.pl --ivar $tsv --vcf ${prefix}_ivar.vcf
    """
}
