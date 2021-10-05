//
// Variant calling and downstream processing for IVar
//

params.ivar_variants_options        = [:]
params.ivar_variants_to_vcf_options = [:]
params.tabix_bgzip_options          = [:]
params.tabix_tabix_options          = [:]
params.bcftools_stats_options       = [:]

include { IVAR_VARIANTS         } from '../../modules/nf-core/modules/ivar/variants/main'  addParams( options: params.ivar_variants_options        )
include { IVAR_VARIANTS_TO_VCF  } from '../../modules/local/ivar_variants_to_vcf'           addParams( options: params.ivar_variants_to_vcf_options )
include { VCF_BGZIP_TABIX_STATS } from '../nf-core/vcf_bgzip_tabix_stats'                   addParams( bgzip_options: params.tabix_bgzip_options, tabix_options: params.tabix_tabix_options, stats_options: params.bcftools_stats_options )

workflow VARIANTS_IVAR {
    take:
    bam                 // channel: [ val(meta), [ bam ] ]
    fasta               // channel: /path/to/genome.fasta
    gff                 // channel: /path/to/genome.gff
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants

    main:
    IVAR_VARIANTS ( bam, fasta, gff )
    IVAR_VARIANTS_TO_VCF ( IVAR_VARIANTS.out.tsv, ivar_multiqc_header ) // Convert original iVar output to VCF, zip and index
    VCF_BGZIP_TABIX_STATS ( IVAR_VARIANTS_TO_VCF.out.vcf )


    emit:
    tsv                 = IVAR_VARIANTS.out.tsv                      // channel: [ val(meta), [ tsv ] ]
    ivar_version        = IVAR_VARIANTS.out.version                  //    path: *.version.txt

    vcf_orig            = IVAR_VARIANTS_TO_VCF.out.vcf               // channel: [ val(meta), [ vcf ] ]
    log_out             = IVAR_VARIANTS_TO_VCF.out.log               // channel: [ val(meta), [ log ] ]
    multiqc_tsv         = IVAR_VARIANTS_TO_VCF.out.tsv               // channel: [ val(meta), [ tsv ] ]

    vcf                 = VCF_BGZIP_TABIX_STATS.out.vcf              // channel: [ val(meta), [ vcf ] ]
    tbi                 = VCF_BGZIP_TABIX_STATS.out.tbi              // channel: [ val(meta), [ tbi ] ]
    stats               = VCF_BGZIP_TABIX_STATS.out.stats            // channel: [ val(meta), [ txt ] ]
    tabix_version       = VCF_BGZIP_TABIX_STATS.out.tabix_version    //    path: *.version.txt
    bcftools_version    = VCF_BGZIP_TABIX_STATS.out.bcftools_version //    path: *.version.txt

}
