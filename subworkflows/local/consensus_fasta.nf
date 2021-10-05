//
// bcftools norm (remove duplicates), view (convert output to gz) and index vcf file and run bcftools consensus to create fasta sequence file
//

params.cut_vcf_options             = [:]
params.bcftools_norm_options       = [:]
params.bcftools_view_options       = [:]
params.bcftools_index_options      = [:]
params.bcftools_consensus_options  = [:]

include { CUT_VCF            } from '../../modules/local/cut_vcf'            addParams( options: params.cut_vcf_options            )
include { BCFTOOLS_NORM      } from '../../modules/local/bcftools_norm'      addParams( options: params.bcftools_norm_options      )
include { BCFTOOLS_VIEW      } from '../../modules/local/bcftools_view'      addParams( options: params.bcftools_view_options      )
include { BCFTOOLS_INDEX     } from '../../modules/local/bcftools_index'     addParams( options: params.bcftools_index_options     )
include { BCFTOOLS_CONSENSUS } from '../../modules/local/bcftools_consensus' addParams( options: params.bcftools_consensus_options )

workflow CONSENSUS_FASTA {
    take:
    vcf // channel: [ vcf ]
    fasta // path : fasta

    main:

    //
    // Cut vcf file
    //
    CUT_VCF ( vcf )

    //
    // Remove duplicates using BCFTOOLS
    //
    BCFTOOLS_NORM ( CUT_VCF.out.vcf )

    //
    // Convert vcf to vcf.gz using BCFTOOLS
    //
    BCFTOOLS_VIEW ( BCFTOOLS_NORM.out.vcf )
    BCFTOOLS_VIEW.out.vcf.map { file ->
        def meta = [:]
        meta.id =  file.baseName[0..-5]
        return tuple(meta, file) } // changed meta.id to metals wo
        .groupTuple()
        .set { ch_vcf }

    //
    // Index vcf file using BCFTOOLS
    //
    BCFTOOLS_INDEX ( BCFTOOLS_VIEW.out.vcf )
    //BCFTOOLS_INDEX.out.index.view()
    BCFTOOLS_INDEX.out.index.map { file ->
        def meta = [:]
        meta.id =  file.baseName[0..-8]
        return tuple(meta, file) }
        .groupTuple()
        .set { ch_index }

    vcf_csi = ch_vcf.join(ch_index, by: [0])  //.toList()

    //ch_final = vcf_csi.mix(fasta)
    //rm_id = vcf_csi.remove([0])
    //ch_final = vcf_csi.minus(rm_id)
    //ch_final.view()
    //BCFTOOLS_INDEX.out.index
    //    .join(BCFTOOLS_VIEW.out.vcf, by: [0], remainder: true)
    //    .set { ch_bam_bai }

    //
    // Build consensus fasta using BCFTOOLS
    //
    BCFTOOLS_CONSENSUS ( vcf_csi, fasta )

    emit:
    cut_vcf          = CUT_VCF.out.vcf              // channel: [ vcf   ]
    norm_out         = BCFTOOLS_NORM.out.vcf        // channel: [ vcf   ]
    view_out         = BCFTOOLS_VIEW.out.vcf        // channel: [ vcf   ]
    index_out        = BCFTOOLS_INDEX.out.index     // channel: [ vcf   ]
    consensus_out    = BCFTOOLS_CONSENSUS.out.fasta // channel: [ vcf,fasta ]
    bcftools_version = BCFTOOLS_NORM.out.version   // path: *.version.txt
}
