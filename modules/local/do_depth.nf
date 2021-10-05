// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DO_DEPTH {
    tag "Calculating per site coverage"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0"
    } else {
        container "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }

    input:
    path(bam)
    path(bai)

    output:
    path("*_samtools.depth"), optional:true, emit: coverage

    script:
    def software = getSoftwareName(task.process)
    """
    mkdir -p bamfiles/
    cd bamfiles
    ln -s $bam .
    mkdir -p ../baifiles/
    cd ../baifiles
    ln -s $bai .
    cd ../
    bamfiles=`ls bamfiles`
    samtools depth -aa -m 0 \$bamfiles > raw_samtools.depth # output all positions, even where no variants and also no max value on depth
    echo Region Position `echo \$bamfiles | sed 's/_indelqual.bam//g'` > header.txt
    cat header.txt raw_samtools.depth > merged_samtools.depth
    """
}
