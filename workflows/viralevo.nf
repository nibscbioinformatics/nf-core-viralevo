/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowViralevo.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.fasta, params.gff, params.multiqc_config, params.adapter_fasta, params.primer_bed ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check if genome exists in the config file
if (params.virus_reference && params.genome && !params.virus_reference.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// TODO remove if rest works
// Stage dummy file to be used as an optional input where required
//ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()


// TODO create save reference param (see sarek) to save index files etc
// TODO check if possible to provide cache/path to pre-installed snpeff dataset

// TODO add readgroups to bwa mem using more robust method

// Initialize value channels based on params,
//Creating value channels for reference files
//params.snpeff_db = params.genome ? params.virus_reference[params.genome].snpeff_db ?: null : null
//snpeff_db         = params.snpeff_db         ?: ''
//snpeff_cache      = params.snpeff_cache      ? Channel.fromPath(params.snpeff_cache).collect()      : ch_dummy_file



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

ch_primer_bed = params.primer_bed ? Channel.value(file(params.primer_bed)) : "null"

ch_adapter_fasta = params.adapter_fasta ? Channel.value(file(params.adapter_fasta)) : "null"

ch_genome_version = params.genome_version ? Channel.value(params.genome_version) : "null"


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

// Dont want to rewrite module code, so add readgroup here
// TODO check this works, compare to viralevo
//def bwamem2_options   = modules['bwamem2_mem']
// bwamem2_options.args += "-R '@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:Illumina'"

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['csv':'']]  )
include { TRIMLOG               } from '../modules/local/trimlog'               addParams( options: modules['trimlog'] )
include { LOFREQ_INDELQUAL      } from '../modules/local/lofreq_indelqual'      addParams( options: modules['lofreq_indelqual']   )
include { ALIGNMENTLOG          } from '../modules/local/alignmentlog'          addParams( options: modules['alignmentlog'] )
include { LOFREQ_CALLPARALLEL   } from '../modules/local/lofreq_callparallel'   addParams( options: modules['lofreq_callparallel'])
include { DO_DEPTH              } from '../modules/local/do_depth'              addParams( options: modules['do_depth']          )
include { SNPEFF_BUILD          } from '../modules/local/snpeff_build'          addParams( options: modules['snpeff_build' ]    )
include { SNPEFF_ANN            } from '../modules/local/snpeff_ann'            addParams( options: modules['snpeff_ann' ]        )
include { TSV2VCF               } from '../modules/local/tsv2vcf'               addParams( options: modules['tsv2vcf' ]           )
include { MAKEVARTABLE          } from '../modules/local/makevartable'          addParams( options: modules['makevartable' ]      )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK           } from '../subworkflows/local/input_check'      addParams( options: [:] )
include { LOFREQ_INDEX_FLAGSTAT    } from '../subworkflows/local/lofreq_index_flagstat'      addParams( lofreq_indelqual_options: modules['lofreq_indelqual'], index_options: modules['samtools_index'], flagstat_options: modules['samtools_flagstat'] )


include { PRIMER_TRIM_IVAR } from '../subworkflows/local/primer_trim_ivar' addParams( ivar_trim_options: modules['ivar_trim'], samtools_options: modules['ivar_trim_sort_bam'] )
include { CONSENSUS_FASTA  } from '../subworkflows/local/consensus_fasta'  addParams( cut_vcf_options: modules['cut_vcf'], bcftools_norm_options: modules['bcftools_norm'], bcftools_view_options: modules['bcftools_view'], bcftools_index_options: modules['bcftools_index'], bcftools_consensus_options: modules['bcftools_consensus'] )

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
include { FASTQC as FASTQC_RAW                       } from '../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_raw']               )
include { CUTADAPT                                   } from '../modules/nf-core/modules/cutadapt/main'       addParams( options: modules['cutadapt']                 )
include { FASTQC as FASTQC_TRIMMED                   } from '../modules/nf-core/modules/fastqc/main'         addParams( options: modules['fastqc_trimmed']           )
include { BWAMEM2_INDEX                              } from '../modules/nf-core/modules/bwamem2/index/main'  addParams( options: modules['bwamem2_index']            )
include { BWAMEM2_MEM                                } from '../modules/nf-core/modules/bwamem2/mem/main'    addParams( options: modules['bwamem2_mem']              )  // altered module code
// seperated samtools workflow to sort bwa-mem bam but run stats on indelqual out (following viralevo DSL1)
include { SAMTOOLS_SORT                              } from '../modules/nf-core/modules/samtools/sort/main' addParams( options: modules['samtools_sort']          )

include { SAMTOOLS_INDEX                             } from '../modules/nf-core/modules/samtools/index/main' addParams( options: modules['samtools_index']           )
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_INDELQUAL } from '../modules/nf-core/modules/samtools/index/main' addParams( options: modules['samtools_index_indelqual'] )
include { SAMTOOLS_FAIDX                             } from '../modules/nf-core/modules/samtools/faidx/main' addParams( options: modules['samtools_faidx']           )
include { IVAR_VARIANTS                              } from '../modules/nf-core/modules/ivar/variants/main'  addParams( options: modules['ivar_variants']            )
include { SNPEFF                                     } from '../modules/nf-core/modules/snpeff/main'         addParams( options: modules['snpeff']                   )
include { MULTIQC                                    } from '../modules/nf-core/modules/multiqc/main'        addParams( options: [:]                                 )

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { BAM_STATS_SAMTOOLS     } from '../subworkflows/nf-core/bam_stats_samtools'          addParams( options: [:] )
include { BAM_SORT_SAMTOOLS      } from '../subworkflows/nf-core/bam_sort_samtools'           addParams( sort_options: modules['samtools_sort'], index_options: modules['samtools_index'], stats_options: modules['samtools_stats'] )
//include { MARK_DUPLICATES_PICARD } from '../subworkflows/nf-core/mark_duplicates_picard'      addParams( markduplicates_options: modules['picard_markduplicates'], samtools_index_options: modules['picard_markduplicates_samtools'], samtools_stats_options: modules['picard_markduplicates_samtools'] )

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
    //ch_fastq.view()

///////////////////
// QC & Trimming
//////////////////

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
        ch_fastq,
        ch_adapter_fasta
    )
    ch_trimmed_reads = CUTADAPT.out.reads
    ch_trim_logs     = CUTADAPT.out.log.collect{it[1]}
    ch_software_versions = ch_software_versions.mix(CUTADAPT.out.version.first().ifEmpty(null))

    //
    // Trimming Log: Output CSV table of trimming stats for reading in R
    //
    // TODO there is a dot after the name in output file.. correct to make sure it doesnt impact downstream
    TRIMLOG (
        ch_trim_logs
    )

    //
    // MODULE: Run FASTQC on trimmed reads
    //
    FASTQC_TRIMMED (
        ch_trimmed_reads
    )
    ch_software_versions = ch_software_versions.mix(FASTQC_TRIMMED.out.version.first().ifEmpty(null))


///////////////////
// Index References
//////////////////

    //
    // MODULE: Create reference genome index using BWA-MEM
    //
    BWAMEM2_INDEX (
        ch_fasta
    )
    ch_index = BWAMEM2_INDEX.out.index
    ch_software_versions = ch_software_versions.mix(BWAMEM2_INDEX.out.version.first().ifEmpty(null))


    //
    // MODULE: Run SAMTOOLS to index reference fasta
    //
    SAMTOOLS_FAIDX (
        ch_fasta
    )
    ch_fasta_fai = SAMTOOLS_FAIDX.out.fai


//////////////////
// Alignment & Bam Stats
//////////////////


    //
    // MODULE: Alignment using BWA-MEM
    //
    BWAMEM2_MEM (
        ch_trimmed_reads,
        ch_index
    )
    ch_bam = BWAMEM2_MEM.out.bam

    //
    // MODULE: Sort bam file using samtools sort
    //
    SAMTOOLS_SORT (
        ch_bam
    )
    ch_sorted_bam = SAMTOOLS_SORT.out.bam

    // Incorporated module into subworkflow below
    //
    // MODULE: Insert indel quality using LOFREQ
    //
    //LOFREQ_INDELQUAL (
    //    ch_sorted_bam, ch_fasta
    //)
    //ch_indelqual_bam          = LOFREQ_INDELQUAL.out.bam
    // ch_software_versions      = ch_software_versions.mix(LOFREQ_INDELQUAL.out.version.first().ifEmpty(null))

    //
    // SUBWORKFLOW: LoFreq indelqual, index and stats on indelqual bam files using SAMTOOLS
    //
    LOFREQ_INDEX_FLAGSTAT (
        ch_sorted_bam,
        ch_fasta
    )
    ch_indelqual_bam = LOFREQ_INDEX_FLAGSTAT.out.bam
    ch_indelqual_bai = LOFREQ_INDEX_FLAGSTAT.out.bai
    ch_alignment_stats = LOFREQ_INDEX_FLAGSTAT.out.flagstat.collect{it[1]}
    ch_software_versions = ch_software_versions.mix(LOFREQ_INDEX_FLAGSTAT.out.version.first().ifEmpty(null))

    //
    // MODULE: Alignment Log: Output CSV table of samtools alignment stats for reading in R
    //
    ALIGNMENTLOG (
        ch_alignment_stats
    )

///////////////
// LoFreq Variant Calling
//////////////

    //
    // MODULE: Run LOFREQ on indelqual_bam for variant calling
    //
    // TODO why does this module take in 4 inputs but use only 2?
    LOFREQ_CALLPARALLEL (
        ch_indelqual_bam,
        ch_indelqual_bai,
        ch_fasta,
        ch_fasta_fai
    )
    ch_lofreq_variants = LOFREQ_CALLPARALLEL.out.vcf

    //
    // MODULE: Generate depth for LoFreq variant calls (Only for small genomes)
    //
    DO_DEPTH (
        ch_indelqual_bam.collect{it[1]},
        ch_indelqual_bai.collect{it[1]}
    )

///////////////////
// iVAR Variant Calling
// https://andersen-lab.github.io/ivar/html/manualpage.html
//////////////////

    //
    // SUBWORKFLOW: Run IVAR for primer trimming and SAMTOOLS to sort, index and stats on indelqual_bam files
    //
    ch_indelqual_bam_and_bai = ch_indelqual_bam.join(ch_indelqual_bai, by: [0])

    PRIMER_TRIM_IVAR (
        ch_indelqual_bam_and_bai,
        ch_primer_bed
    )
    ch_primer_trimmed_sorted_bam = PRIMER_TRIM_IVAR.out.bam

    //
    // MODULE: Call variants with IVAR
    //
    IVAR_VARIANTS (
        ch_primer_trimmed_sorted_bam, ch_fasta, ch_annotation
    )
    ch_ivar_variants     = IVAR_VARIANTS.out.tsv
    ch_software_versions = ch_software_versions.mix(IVAR_VARIANTS.out.version.first().ifEmpty(null))

    //
    // MODULE: Convert TSV to VCF
    //
    TSV2VCF (
        ch_ivar_variants
    )
    ch_ivar2vcf = TSV2VCF.out.vcf

    //
    // Merge variant calls from LOFREQ and IVAR
    //
    merged_ch = ch_lofreq_variants.mix(ch_ivar2vcf)

/////////////////////////
// SNPEff Annotation
////////////////////////

    //
    // MODULE: Build SNP Annotation DB from gff and fasta
    //

    if (!params.noannotation) {


    SNPEFF_BUILD (
        ch_fasta,
        ch_annotation
    )
    ch_snpeff_db = SNPEFF_BUILD.out.db
    ch_snpeff_config = SNPEFF_BUILD.out.config

    //
    // MODULE: Run snpEff to annotate merged variants
    //
    SNPEFF_ANN (
        merged_ch,
        ch_snpeff_db,
        ch_snpeff_config,
        ch_fasta
    )
    ch_vcffortable = SNPEFF_ANN.out.vcf.collect{it[1]}
    ch_software_versions = ch_software_versions.mix(SNPEFF_ANN.out.version.first().ifEmpty(null))
    ch_annotatedvcf = ch_vcffortable.flatten().filter( ~/^.*vcf/ ) // flattens into one elemenbt, only vcf files... not really sure what this is for yet
    } else {
        ch_vcffortable = merged_ch.collect{it[1]} // if no annotation process simply use merged table for input in report processes
    }

/////////////////////////
// Reporting
////////////////////////


    //
    // MODULE: Take output from annotated vcf files, generate table and write out a filtered VCF file for each input VCF file
    //
    //vcf = Channel.fromPath('/Data/Users/rbhuller/tmp/new/results/variants/snpeff/vcf')
    //vcf = Channel.fromPath( './results/variants/snpeff/vcf', type: 'dir' )



    MAKEVARTABLE (
        ch_vcffortable,
        params.alt_depth_threshold,
        params.vaf_threshold
    )
    ch_filtered_vcfs = MAKEVARTABLE.out.filteredvars.flatten()
    ch_filtered_vcfs.view()




    //
    // SUBWORKFLOW: Index and stats on indelqual bam files using SAMTOOLS
    //
    //BAM_SORT_SAMTOOLS (
    //    ch_indelqual_bam
    //)
    //bam = BAM_SORT_SAMTOOLS.out.bam
    //ch_software_versions = ch_software_versions.mix(BAM_SORT_SAMTOOLS.out.version.first().ifEmpty(null))

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
    // MODULE: Run SAMTOOLS to index indelqual bam
    //
   // SAMTOOLS_INDEX_INDELQUAL (
    //    ch_indelqual_bam
    //)
    ///ch_indelqual_bam_bai =  SAMTOOLS_INDEX_INDELQUAL.out.bai


    //ch_lofreq_variants.view()

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
    // MODULE: Take output from annotated vcf files, generate table and write out a filtered VCF file for each input VCF file
    //
    //vcf = Channel.fromPath('/Data/Users/rbhuller/tmp/new/results/variants/snpeff/vcf')
    //vcf = Channel.fromPath( './results/variants/snpeff/vcf', type: 'dir' )

    //MAKEVARTABLE (
    //    ch_annotatedfortable, params.alt_depth_threshold, params.vaf_threshold
    //)
    //ch_filtered_vcfs = MAKEVARTABLE.out.filteredvars.flatten()
    //ch_filtered_vcfs.view()

    //
    // Subworkflow: Build a consensus using bcftools from the filtered vcfs
    //
    //CONSENSUS_FASTA (
    //    ch_filtered_vcfs, ch_fasta
    //)
    //CONSENSUS_FASTA.out.view_out.view()

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
