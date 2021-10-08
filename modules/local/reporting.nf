// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process REPORTING {
    tag 'Generating report..'
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'tsv2vcf', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "$baseDir/containers/environment.yml"  : null)
    if (workflow.containerEngine == 'singularity' ) {
        container '/usr/share/sequencing/singularity_images/viralevo-reporting-1.0.4.img'
    }

    input:
    val vcfData //ch_filtered vars; changed from val
    file(rmodel)  // as is
    val bamData // samtools index output, ch_indelqual_bai = LOFREQ_INDEX_FLAGSTAT.out.bai
    file(trimsummary) // ch_trimlog
    file(alignmentsummary) //ch_alignlog
    file(samdepth) // do depth out ch_samdepth
    file(varianttable)  // ch_vartable

    output:
    path("analysis_report.html"), emit: html
    path("analysis_report.RData"), emit: rdata

    script:

    // handling here the VCF files and metadata
    def sampleNamesList = []
    def callersList = []
    def vcfList = []
    vcfData.each() { sample,caller,vcf -> //each to iterate on each element of list, map doesnt work
    sampleNamesList.add(sample)
    callersList.add(caller)
    vcfList.add(vcf)
    }
    sampleNames = sampleNamesList.join(",")
    callerLabels = callersList.join(",")
    vcfFiles = vcfList.join(",")

    // handling the BAM files and metadata
    // expects tuple: ID, BAM, BAI in one channel
    def bamSampleList = []
    def bamList = []
    def baiList = []
    bamData.each() { sample,bam,bai ->
    bamSampleList.add(sample)
    bamList.add(bam)
    baiList.add(bai)
    }
    bamSamples = bamSampleList.join(",")
    bamFiles = bamList.join(",")
    baiFiles = baiList.join(",")

    """
    ln -s $baseDir/docs/nibsc_report.css .

    Rscript -e "workdir<-getwd()
    rmarkdown::render('$baseDir/docs/analysis_report.Rmd',
    params = list(
        vcf = \\\"$vcfFiles\\\",
        callers = \\\"$callerLabels\\\",
        samples = \\\"$sampleNames\\\",
        genome = \\\"${params.genome}\\\",
        genemodel = \\\"$rmodel\\\",
        baseDir = \\\"$baseDir\\\",
        bamSamples = \\\"$bamSamples\\\",
        bamFiles = \\\"$bamFiles\\\",
        trimsummarytable = \\\"$trimsummary\\\",
        alignmentsummarytable = \\\"$alignmentsummary\\\",
        samdepthtable = \\\"$samdepth\\\",
        varianttable = \\\"$varianttable\\\",
        noannotation = \\\"${params.noannotation}\\\"
        ),
    knit_root_dir=workdir,
    output_dir=workdir)"
    """

}
