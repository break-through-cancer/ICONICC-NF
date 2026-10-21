include { SUBSET_PAIR_VCF     } from '../../../modules/local/subset_pair_vcf'
include { FILTER_VARIANTS     } from '../../../modules/local/filter_variants'
include { ITER_ISEC           } from '../../../modules/local/iter_isec'
include { HETVAR_FILTER       } from '../../../modules/local/hetvar_filter'
include { EAGLE_PHASING       } from '../../../modules/local/eagle_phasing'
include { INDEX_PHASED        } from '../../../modules/local/index_phased'
include { QUERY_HET_AD        } from '../../../modules/local/query_het_ad'
include { AC_CALC             } from '../../../modules/local/ac_calc'
include { TANGENT_XY          } from '../../../modules/local/tangent_xy'
include { SEGMENT             } from '../../../modules/local/segment'

// Partitions chromosomes into ITER_ISEC/EAGLE_PHASING's shards -- round-robin for
// 'groups' (already size-ordered, so this balances shards), one-per-chromosome otherwise.
def buildChrGroups(mode, count, chromosomes) {
    if (!(mode in ['groups', 'per_chromosome'])) {
        error "CNV_CORE: cnv_scatter_mode must be 'groups' or 'per_chromosome', got '${mode}'"
    }
    if (mode == 'per_chromosome') {
        return chromosomes.collect { [it] }
    }
    def n = count as int
    def groups = (0..<n).collect { [] }
    chromosomes.eachWithIndex { chr, idx -> groups[idx % n] << chr }
    return groups.findAll { !it.isEmpty() }
}

// Extracts the chromosome token from the shared "_chr<token>_" filename convention.
def chrOf(path) {
    def m = (path.name =~ /_chr([^_]+)_/)
    if (!m) {
        error "CNV_CORE: could not extract a chromosome token from filename '${path.name}' (expected '_chr<token>_')"
    }
    return m[0][1]
}

// Fresh channel per call (queue channels split across consumers, don't broadcast); each
// group wrapped as a 1-tuple so combine() doesn't spread its elements into the tuple.
def chrGroupChannel() {
    Channel.fromList(buildChrGroups(params.cnv_scatter_mode, params.cnv_scatter_count, params.scatter_chromosomes))
        .map { [it] }
}

workflow CNV_CORE {
    take:
    joint_vcf_ch     // [ meta, joint_vcf, joint_vcf_index ] -- meta is the single uniform
                     // per-participant map used by every process below: participant_id,
                     // normal_sample_id, tumor_sample_id, sample_sex, agg_bins,
                     // desired_binsize, and every tangent_xy/segment tuning field. It must
                     // never be reshaped downstream (only ever read from), so every
                     // .join() on `meta` below matches.
    tumor_counts_ch  // [ meta, counts_file ] -- same meta as joint_vcf_ch. Pre-computed by
                     // SAMPLE_READ_COUNTS (once per unique sample, before pairing) and
                     // joined into this pair's meta by main.nf before calling CNV_CORE --
                     // this subworkflow no longer computes read counts itself.
    ref_dir_1000g
    genetic_map
    pon
    ref_clinical
    cohort_clinical

    main:
    SUBSET_PAIR_VCF(joint_vcf_ch)

    FILTER_VARIANTS(SUBSET_PAIR_VCF.out)

    // Scatter-gather for ITER_ISEC/EAGLE_PHASING (SPEC.md Sec14.1); HETVAR_FILTER/INDEX_PHASED stay single-task.
    ITER_ISEC(
        FILTER_VARIANTS.out.combine(ref_dir_1000g).combine(chrGroupChannel())
    )

    // groupTuple()+flatten() reassembles all shards' files into the flat list HETVAR_FILTER expects.
    HETVAR_FILTER(
        ITER_ISEC.out.groupTuple()
            .map { meta, vcf_lists, idx_lists -> [ meta, vcf_lists.flatten(), idx_lists.flatten() ] }
    )

    // EAGLE_PHASING needs no module change -- just route it a per-group subset by filename.
    eagle_input = HETVAR_FILTER.out.combine(ref_dir_1000g).combine(genetic_map)
        .combine(chrGroupChannel())
        .map { meta, het_vcfs, het_vcf_indices, ref_dir, gmap, chr_group ->
            def wanted = chr_group as Set
            [
                meta,
                het_vcfs.findAll { wanted.contains(chrOf(it)) },
                het_vcf_indices.findAll { wanted.contains(chrOf(it)) },
                ref_dir,
                gmap
            ]
        }
    EAGLE_PHASING(eagle_input)

    // Gathered once, reused by both INDEX_PHASED and QUERY_HET_AD below.
    eagle_phased_ch = EAGLE_PHASING.out.groupTuple()
        .map { meta, phased_lists -> [ meta, phased_lists.flatten() ] }

    INDEX_PHASED(eagle_phased_ch)

    // QUERY_HET_AD derives its sample list from whichever het VCF it's handed first.
    query_input = HETVAR_FILTER.out
        .join(eagle_phased_ch)
        .join(INDEX_PHASED.out)
    QUERY_HET_AD(query_input)

    AC_CALC(
        tumor_counts_ch
            .join(QUERY_HET_AD.out)
            .map { meta, tumor_counts, hetinfo, ad, phased -> [ meta, tumor_counts, phased, ad, hetinfo ] }
    )

    TANGENT_XY(
        tumor_counts_ch
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
