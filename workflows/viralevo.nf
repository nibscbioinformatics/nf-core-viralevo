/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowViralevo.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.primer_fasta, params.primer_bed, params.anno, params.phyref, params.genome_rmodel, params.ivar_calling_af_threshold, params.ivar_calling_dp_threshold, params.vaf_threshold, params.alt_depth_threshold, params.noannotation, params.tools ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check if genome exists in the config file
if (params.virus_reference && params.genome && !params.virus_reference.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

//Creating value channels for reference files
params.anno = params.genome ? params.virus_reference[params.genome].gff ?: null : null
if (params.anno) { ch_annotation = Channel.value(file(params.anno, checkIfExists: true)) }

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
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['csv':'']] )
include { LOFREQ_INDELQUAL      } from '../modules/local/lofreq_indelqual'      addParams( options: [:]                          )
include { LOFREQ_CALLPARALLEL   } from '../modules/local/lofreq_callparallel'   addParams( options: [:]                          )
include { SNPEFF_ANN            } from '../modules/local/snpeff_ann'            addParams( options: [ 'args': '-ud 1' ] )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK      } from '../subworkflows/local/input_check'      addParams( options: [:] )
include { PRIMER_TRIM_IVAR } from '../subworkflows/local/primer_trim_ivar' addParams( options: [:] )
include { VARIANTS_IVAR    } from '../subworkflows/local/variants_ivar'    addParams( options: [:] )


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
include { FASTQC as FASTQC_RAW     } from '../modules/nf-core/software/fastqc/main'             addParams( options: modules['fastqc'] )
include { CUTADAPT                 } from '../modules/nf-core/software/cutadapt/main'           addParams( options: [ 'args': '-q 30 --times 40 --minimum-length 50 --error-rate 0.1 --max-n 0' ] )
include { FASTQC as FASTQC_TRIMMED } from '../modules/nf-core/software/fastqc/main'             addParams( options: [:] )
include { BWAMEM2_INDEX            } from '../modules/nf-core/software/bwamem2/index/main'      addParams( options: [:] )
include { BWAMEM2_MEM              } from '../modules/nf-core/software/bwamem2/mem/main'        addParams( options: [:] )
include { SAMTOOLS_INDEX           } from '../modules/nf-core/software/samtools/index/main'     addParams( options: [:] )
include { SAMTOOLS_FAIDX           } from '../modules/nf-core/software/samtools/faidx/main'     addParams( options: [:] )
include { MULTIQC                  } from '../modules/nf-core/software/multiqc/main'            addParams( options: [:] )

//
// MODULE: Installed directly from nf-core/modules
//
include { BAM_STATS_SAMTOOLS     } from '../subworkflows/nf-core/bam_stats_samtools'          addParams( options: [:] )
include { BAM_SORT_SAMTOOLS      } from '../subworkflows/nf-core/bam_sort_samtools'           addParams( options: [:] )
include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard'      addParams( options: [:] )

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

    //
    // MODULE: Run FASTQC on trimmed reads
    //

    //
    // MODULE: Create reference genome index using BWA-MEM
    //

    //
    // MODULE: Alignment using BWA-MEM
    //

    //
    // SUBWORKFLOW: Sort, index and stats on bam files using SAMTOOLS
    //

    //
    // SUBWORKFLOW: Mark duplicate reads using PICARD and stats
    //

    //
    // MODULE: Insert indel quality using LOFREQ
    //

    //
    // MODULE: Run SAMTOOLS to index indelqual bam
    //

    //
    // MODULE: Run SAMTOOLS to index reference fasta
    //

    //
    // MODULE: Run LOFREQ for variant calling
    //

    //
    // MODULE: Run IVAR for primer trimming
    //

    //
    // MODULE: Mark duplicate reads using PICARD
    //

    //
    // MODULE: Run IVAR for variant calling
    //

    //
    // Merge variant calls from LOFREQ and IVAR
    //

    //
    // MODULE: Run snpEff to annotate merged variants
    //


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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

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
