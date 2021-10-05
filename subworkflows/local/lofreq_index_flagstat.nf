/*
 * Index BAM file and run samtools stats, flagstat and idxstats
 */

params.index_options = [:]
params.flagstat_options = [:]

include { LOFREQ_INDELQUAL   } from '../../modules/local/lofreq_indelqual'            addParams( options: params.lofreq_indelqual_options )
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/modules/samtools/index/main'    addParams( options: params.index_options )
include { SAMTOOLS_FLAGSTAT  } from '../../modules/nf-core/modules/samtools/flagstat/main' addParams( options: params.flagstat_options )

workflow LOFREQ_INDEX_FLAGSTAT {
    take:
    ch_sorted_bam
    ch_fasta

    main:
    LOFREQ_INDELQUAL    ( ch_sorted_bam, ch_fasta )
    SAMTOOLS_INDEX      ( LOFREQ_INDELQUAL.out.bam )

    // create new ch of bam and bai files for flagstat input
    LOFREQ_INDELQUAL.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .set { ch_bam_bai }

    SAMTOOLS_FLAGSTAT   ( ch_bam_bai )

    emit:
    bam      = LOFREQ_INDELQUAL.out.bam    // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]

    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    version  = SAMTOOLS_INDEX.out.version       //    path: *.version.txt

}
