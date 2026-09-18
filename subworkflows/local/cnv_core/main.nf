include { SUBSET_PAIR_VCF     } from '../../../modules/local/subset_pair_vcf'
include { PREPROCESS_INTERVALS } from '../../../modules/local/preprocess_intervals'
include { COLLECT_READ_COUNTS as COLLECT_READ_COUNTS_TUMOR  } from '../../../modules/local/collect_read_counts'
include { COLLECT_READ_COUNTS as COLLECT_READ_COUNTS_NORMAL } from '../../../modules/local/collect_read_counts'
include { FILTER_VARIANTS     } from '../../../modules/local/filter_variants'
include { ITER_ISEC           } from '../../../modules/local/iter_isec'
include { HETVAR_FILTER       } from '../../../modules/local/hetvar_filter'
include { EAGLE_PHASING       } from '../../../modules/local/eagle_phasing'
include { INDEX_PHASED        } from '../../../modules/local/index_phased'
include { QUERY_HET_AD        } from '../../../modules/local/query_het_ad'
include { AC_CALC             } from '../../../modules/local/ac_calc'
include { TANGENT_XY          } from '../../../modules/local/tangent_xy'
include { SEGMENT             } from '../../../modules/local/segment'

workflow CNV_CORE {
    take:
    joint_vcf_ch     // [ meta, joint_vcf, joint_vcf_index ] -- meta is the single uniform
                     // per-participant map used by every process below: participant_id,
                     // normal_sample_id, tumor_sample_id, sample_sex, agg_bins,
                     // desired_binsize, and every tangent_xy/segment tuning field (Global
                     // Constraints / Task 5's note). It must never be reshaped downstream
                     // (only ever read from), so every .join() on `meta` below matches.
    tumor_cram_ch    // [ meta, cram, crai ] -- same meta as joint_vcf_ch
    normal_cram_ch   // [ meta, cram, crai ] -- same meta as joint_vcf_ch
    ref_fasta_ch     // [ fasta, fai, dict ]
    ref_dir_1000g
    genetic_map
    pon
    ref_clinical
    cohort_clinical

    main:
    SUBSET_PAIR_VCF(joint_vcf_ch)

    PREPROCESS_INTERVALS(ref_fasta_ch)

    COLLECT_READ_COUNTS_TUMOR(
        tumor_cram_ch.combine(PREPROCESS_INTERVALS.out).combine(ref_fasta_ch)
            .map { meta, cram, crai, intervals, fasta, fai, dict -> [ meta, meta.tumor_sample_id, cram, crai, intervals, fasta, fai, dict ] }
    )
    COLLECT_READ_COUNTS_NORMAL(
        normal_cram_ch.combine(PREPROCESS_INTERVALS.out).combine(ref_fasta_ch)
            .map { meta, cram, crai, intervals, fasta, fai, dict -> [ meta, meta.normal_sample_id, cram, crai, intervals, fasta, fai, dict ] }
    )

    FILTER_VARIANTS(SUBSET_PAIR_VCF.out)
    ITER_ISEC(FILTER_VARIANTS.out.combine(ref_dir_1000g))
    HETVAR_FILTER(ITER_ISEC.out)
    EAGLE_PHASING(HETVAR_FILTER.out.map { meta, het, _idx -> [ meta, het ] }.combine(ref_dir_1000g).combine(genetic_map))
    INDEX_PHASED(EAGLE_PHASING.out)

    query_input = HETVAR_FILTER.out
        .join(EAGLE_PHASING.out)
        .join(INDEX_PHASED.out)
        .map { meta, het, het_idx, phased, phased_idx ->
            [ meta, het.find { it.name.contains('_chr1_') }, het, het_idx, phased, phased_idx ]
        }
    QUERY_HET_AD(query_input)

    AC_CALC(
        COLLECT_READ_COUNTS_TUMOR.out
            .join(QUERY_HET_AD.out)
            .map { meta, tumor_counts, hetinfo, ad, phased -> [ meta, tumor_counts, phased, ad, hetinfo ] }
    )

    TANGENT_XY(
        COLLECT_READ_COUNTS_TUMOR.out
            .combine(pon).combine(ref_clinical).combine(cohort_clinical)
    )

    SEGMENT(
        AC_CALC.out.map { meta, _rd, ac, _plots -> [ meta, ac ] }
            .join(TANGENT_XY.out)
    )

    emit:
    read_depth     = AC_CALC.out.map { meta, rd, _ac, _plots -> [ meta, rd ] }
    seg_file       = SEGMENT.out.map { meta, seg, _cts, _plots, _scts -> [ meta, seg ] }
    ac_cts         = SEGMENT.out.map { meta, _seg, cts, _plots, _scts -> [ meta, cts ] }
    segment_plots  = SEGMENT.out.map { meta, _seg, _cts, plots, _scts -> [ meta, plots ] }
    segment_cts    = SEGMENT.out.map { meta, _seg, _cts, _plots, scts -> [ meta, scts ] }
}
