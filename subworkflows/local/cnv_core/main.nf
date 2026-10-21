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

// Partitions `chromosomes` (params.scatter_chromosomes) into the shards ITER_ISEC and
// EAGLE_PHASING scatter over. 'per_chromosome': one shard per chromosome -- maximum
// parallelism. 'groups': round-robin into `count` shards -- chromosomes are already size-
// ordered (autosomes descending, then X), so interleaving by index balances large and
// small chromosomes across shards without needing GATK's own interval-size-aware
// splitting (there's no genomic interval here to split, just a list of whole-chromosome
// names). No process needed since there's no file or reference dictionary to read.
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

// Extracts the chromosome token from this pipeline's shared "_chr<token>_" filename
// convention (e.g. TEST_PT_chr10_het.vcf.gz -> '10') -- the same convention every
// module's own script already parses via basename/sed.
def chrOf(path) {
    def m = (path.name =~ /_chr([^_]+)_/)
    if (!m) {
        error "CNV_CORE: could not extract a chromosome token from filename '${path.name}' (expected '_chr<token>_')"
    }
    return m[0][1]
}

// Fresh chromosome-groups channel, called once per consumer (ITER_ISEC, then again for
// EAGLE_PHASING) rather than shared -- a multi-item queue channel's emissions split
// across consumers instead of broadcasting to each, so reusing one Channel.fromList()
// would silently hand each process only part of the chromosome groups. Must be a
// top-level def, not a closure assigned inside the workflow body: Nextflow's workflow-
// block parser does not resolve a bare (non-def) closure variable as a callable
// function within the same workflow (confirmed empirically -- "`x` is not defined").
//
// Each group is wrapped as the single element of a 1-tuple (.map { [it] }) before
// combine() sees it, then unwrapped by the destructuring at each call site --
// combine() spreads a list-valued channel item's own elements into the combined tuple
// (confirmed empirically: combining a 5-chromosome group produced 5 separate
// positional values, not one val(chr_group) list, breaking the consuming process's
// input arity), so passing a whole chromosome list through combine() needs this guard.
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

    // Scatter-gather for ITER_ISEC and EAGLE_PHASING (SPEC.md Sec14.1): both already loop
    // per-chromosome sequentially inside one task each, with no cross-chromosome state
    // (every bcftools/eagle invocation only reads that chromosome's own slice against
    // that chromosome's own reference-panel file), so scattering into parallel
    // per-group tasks and gathering the file list before the next stage is exactly
    // equivalent to today's sequential-loop result. HETVAR_FILTER/INDEX_PHASED stay
    // single-task deliberately -- real timing on a production pair showed ~4min/~1min
    // for all 23 chromosomes combined, where per-task scatter overhead would likely cost
    // more than it saves. Same chromosome groups feed both (chrGroupChannel(), above),
    // since both scatter over the same unit -- a whole chromosome.
    ITER_ISEC(
        FILTER_VARIANTS.out.combine(ref_dir_1000g).combine(chrGroupChannel())
    )

    // Every ITER_ISEC shard emits its own List<Path> (one entry per chromosome in its
    // group); groupTuple() collects all shards for a participant, then .flatten() turns
    // the resulting list-of-lists back into one flat list -- the shape HETVAR_FILTER
    // already expects, unchanged since before scattering.
    HETVAR_FILTER(
        ITER_ISEC.out.groupTuple()
            .map { meta, vcf_lists, idx_lists -> [ meta, vcf_lists.flatten(), idx_lists.flatten() ] }
    )

    // EAGLE_PHASING needs no module change to scatter: it already loops over whatever
    // het VCFs it's handed and derives each chromosome from the filename itself, so
    // routing it a per-group subset of HETVAR_FILTER's (already-gathered) output is a
    // pure channel-wiring change. chrOf() extracts the same "_chr<token>_" naming
    // convention the module's own script already parses via basename/sed.
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

    // Gathered once and reused below for both INDEX_PHASED and QUERY_HET_AD, same
    // flatten-after-groupTuple shape as HETVAR_FILTER's gather above.
    eagle_phased_ch = EAGLE_PHASING.out.groupTuple()
        .map { meta, phased_lists -> [ meta, phased_lists.flatten() ] }

    INDEX_PHASED(eagle_phased_ch)

    // QUERY_HET_AD derives its sample list from whichever het VCF it's handed first --
    // every chromosome carries the same sample list, so there's no need to single out
    // chr1 specifically (doing so previously required staging that one file twice, an
    // exact filename collision caught against a real run).
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
