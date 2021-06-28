/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowViralevo.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.gff, params.multiqc_config, params.primer_fasta, params.primer_bed ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check if genome exists in the config file
if (params.virus_reference && params.genome && !params.virus_reference.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

//Creating value channels for reference files
params.gff = params.genome ? params.virus_reference[params.genome].gff ?: null : null
if (params.gff) { ch_annotation = Channel.value(file(params.gff, checkIfExists: true)) }

params.fasta = params.genome ? params.virus_reference[params.genome].fasta ?: null : null
if (params.fasta) { ch_fasta = Channel.value(file(params.fasta, checkIfExists: true)) }

params.phyref = params.genome ? params.virus_reference[params.genome].pyloref ?: null : null
if (params.phyref) { ch_phyloref = Channel.value(file(params.phyref, checkIfExists: true)) }

params.genome_rmodel = params.genome ? params.virus_reference[params.genome].rmodel ?: null : null
if (params.genome_rmodel) { ch_genome_rmodel = Channel.value(file(params.genome_rmodel, checkIfExists: true)) }

//Creating channels for other files
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)

ch_primer_bed = params.primer_bed ? Channel.value(file(params.primer_bed)) : "null"

ch_primer_fasta = params.primer_fasta ? Channel.value(file(params.primer_fasta)) : "null"

ch_ivar_variants_header_mqc = file("$projectDir/assets/headers/ivar_variants_header_mqc.txt", checkIfExists: true)


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config_illumina.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['csv':'']]  )
include { LOFREQ_INDELQUAL      } from '../modules/local/lofreq_indelqual'      addParams( options: modules['lofreq_indelqual']   )
include { LOFREQ_CALLPARALLEL   } from '../modules/local/lofreq_callparallel'   addParams( options: modules['lofreq_callparallel'])
include { SNPEFF_ANN            } from '../modules/local/snpeff_ann'            addParams( options: modules['snpeff_ann' ]        )
include { MAKEVARTABLE          } from '../modules/local/makevartable'          addParams( options: modules['makevartable' ]      )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK      } from '../subworkflows/local/input_check'      addParams( options: [:] )
include { PRIMER_TRIM_IVAR } from '../subworkflows/local/primer_trim_ivar' addParams( ivar_trim_options: modules['ivar_trim'], samtools_options: modules['ivar_trim_sort_bam'] )
include { VARIANTS_IVAR    } from '../subworkflows/local/variants_ivar'    addParams( ivar_variants_options: modules['ivar_variants'], ivar_variants_to_vcf_options: modules['ivar_variants_to_vcf'], tabix_bgzip_options: modules['ivar_tabix_bgzip'], tabix_tabix_options: modules['ivar_tabix_tabix'], bcftools_stats_options: modules['ivar_bcftools_stats']               )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//def multiqc_options = modules['multiqc']
//multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW                       } from '../modules/nf-core/software/fastqc/main'         addParams( options: modules['fastqc_raw']               )
include { CUTADAPT                                   } from '../modules/nf-core/software/cutadapt/main'       addParams( options: modules['cutadapt']                 )
include { FASTQC as FASTQC_TRIMMED                   } from '../modules/nf-core/software/fastqc/main'         addParams( options: modules['fastqc_trimmed']           )
include { BWAMEM2_INDEX                              } from '../modules/nf-core/software/bwamem2/index/main'  addParams( options: modules['bwamem2_index']            )
include { BWAMEM2_MEM                                } from '../modules/nf-core/software/bwamem2/mem/main'    addParams( options: modules['bwamem2_mem']              )
include { SAMTOOLS_INDEX                             } from '../modules/nf-core/software/samtools/index/main' addParams( options: modules['samtools_index']           ) 
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_INDELQUAL } from '../modules/nf-core/software/samtools/index/main' addParams( options: modules['samtools_index_indelqual'] )
include { SAMTOOLS_FAIDX                             } from '../modules/nf-core/software/samtools/faidx/main' addParams( options: modules['samtools_faidx']           )
include { MULTIQC                                    } from '../modules/nf-core/software/multiqc/main'        addParams( options: [:]                                 )

//
// MODULE: Installed directly from nf-core/modules
//
//include { BAM_STATS_SAMTOOLS     } from '../subworkflows/nf-core/bam_stats_samtools'          addParams( options: [:] )
include { BAM_SORT_SAMTOOLS      } from '../subworkflows/nf-core/bam_sort_samtools'           addParams( sort_options: modules['samtools_sort'], index_options: modules['samtools_index'], stats_options: modules['samtools_stats'] )
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard'      addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: modules['picard_markduplicates_samtools'], samtools_stats_options: modules['picard_markduplicates_samtools'] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow VIRALEVO {

    ch_software_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ] }
    .set { ch_fastq }
    
    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW (
        ch_fastq
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_RAW.out.version.first().ifEmpty(null))

    //
    // MODULE: CUTADAPT for adapter and quality trimming
    //
    CUTADAPT (
        ch_fastq, ch_primer_fasta
    )
    ch_trimmed_reads = CUTADAPT.out.reads
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.first().ifEmpty(null)) 
   
    //
    // MODULE: Run FASTQC on trimmed reads
    //
    FASTQC_TRIMMED (
        ch_trimmed_reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMMED.out.version.first().ifEmpty(null))

    //
    // MODULE: Create reference genome index using BWA-MEM
    //
    BWAMEM2_INDEX (
        ch_fasta
    )
    ch_index = BWAMEM2_INDEX.out.index
    ch_software_versions = ch_software_versions.mix(BWAMEM2_INDEX.out.version.first().ifEmpty(null))

    //
    // MODULE: Alignment using BWA-MEM
    //
    BWAMEM2_MEM (
        ch_trimmed_reads, ch_index
    )
    ch_bam = BWAMEM2_MEM.out.bam

    //
    // SUBWORKFLOW: Sort, index and stats on bam files using SAMTOOLS
    //
    BAM_SORT_SAMTOOLS (
        ch_bam
    )
    bam = BAM_SORT_SAMTOOLS.out.bam
    ch_software_versions = ch_software_versions.mix(BAM_SORT_SAMTOOLS.out.version.first().ifEmpty(null))    

    //
    // SUBWORKFLOW: Mark duplicate reads using PICARD and stats
    //
    //MARK_DUPLICATES_PICARD (
    //    bam
    //)
    //ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
    //ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
    //ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
    //ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))
 
    //
    // MODULE: Insert indel quality using LOFREQ
    //
    LOFREQ_INDELQUAL (
        bam, ch_fasta
    )
    ch_indelqual_bam          = LOFREQ_INDELQUAL.out.bam
    ch_software_versions      = ch_software_versions.mix(LOFREQ_INDELQUAL.out.version.first().ifEmpty(null)) 

    //
    // MODULE: Run SAMTOOLS to index indelqual bam
    //
    SAMTOOLS_INDEX_INDELQUAL (
        ch_indelqual_bam
    )
    ch_indelqual_bam_bai =  SAMTOOLS_INDEX_INDELQUAL.out.bai

    //
    // MODULE: Run SAMTOOLS to index reference fasta
    //
    SAMTOOLS_FAIDX (
        ch_fasta
    )
    ch_fasta_fai = SAMTOOLS_FAIDX.out.fai

    //
    // MODULE: Run LOFREQ on indelqual_bam for variant calling
    //
    LOFREQ_CALLPARALLEL (
        ch_indelqual_bam, ch_indelqual_bam_bai, ch_fasta, ch_fasta_fai 
    )
    ch_lofreq_variants = LOFREQ_CALLPARALLEL.out.vcf
    //ch_lofreq_variants.view()

    //
    // SUBWORKFLOW: Run IVAR for primer trimming and SAMTOOLS to sort, index and stats on indelqual_bam files 
    //
    ch_indelqual_bam_and_bai = ch_indelqual_bam.join(ch_indelqual_bam_bai, by: [0])

    PRIMER_TRIM_IVAR (
        ch_indelqual_bam_and_bai, ch_primer_bed
    )
    ch_primer_trimmed_sorted_bam = PRIMER_TRIM_IVAR.out.bam        
    
    //
    // MODULE: Mark duplicate reads using PICARD
    //
    //MARK_DUPLICATES_PICARD_IVAR (
    //    ch_primer_trimmed_sorted_bam
    //)
    //ch_samtools_stats         = MARK_DUPLICATES_PICARD.out.stats
    //ch_samtools_flagstat      = MARK_DUPLICATES_PICARD.out.flagstat
    //ch_samtools_idxstats      = MARK_DUPLICATES_PICARD.out.idxstats
    //ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.picard_version.first().ifEmpty(null))



    //
    // SUBWORKFLOW: Call variants with IVAR, then convertion to vcf, gzip and stats using bcftools
    //
    VARIANTS_IVAR (
        ch_primer_trimmed_sorted_bam, ch_fasta, ch_annotation, ch_ivar_variants_header_mqc
    )
    ch_ivar_variants     = VARIANTS_IVAR.out.vcf_orig
    //ch_ivar_variants.view()
    ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.ivar_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.tabix_version.first().ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(VARIANTS_IVAR.out.bcftools_version.first().ifEmpty(null))

    //
    // Merge variant calls from LOFREQ and IVAR
    //
    merged_ch = ch_lofreq_variants.mix(ch_ivar_variants)
    //merged_ch.view()

    //
    // MODULE: Run snpEff to annotate merged variants
    //
    SNPEFF_ANN (
        merged_ch, params.genome_version
    )
    ch_annotatedfortable = SNPEFF_ANN.out.vcf
    ch_software_versions = ch_software_versions.mix(SNPEFF_ANN.out.version.first().ifEmpty(null))

    //
    // MODULE: Take output from annotated vcf files and generate tables
    //
    //vcf = Channel.fromPath('/Data/Users/rbhuller/tmp/new/results/variants/snpeff')

    //MAKEVARTABLE (
    //    vcf, params.alt_depth_threshold, params.vaf_threshold
    //)    

    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowViralevo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report       = MULTIQC.out.report.toList()
    ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
