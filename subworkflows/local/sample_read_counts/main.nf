include { PREPROCESS_INTERVALS } from '../../../modules/local/preprocess_intervals'
include { SPLIT_INTERVALS      } from '../../../modules/local/split_intervals'
include { COLLECT_READ_COUNTS  } from '../../../modules/local/collect_read_counts'
include { GATHER_READ_COUNTS   } from '../../../modules/local/gather_read_counts'

workflow SAMPLE_READ_COUNTS {
    take:
    samples_ch    // [ meta(patient,sample,sex,status), cram, crai ]
    ref_fasta_ch  // [ fasta, fai, dict ]

    main:
    // .first() hardens this as a value channel (same reasoning as
    // joint_genotyping/main.nf's own ref_fasta_ch/dbsnp_ch hardening): main.nf builds it
    // via Channel.of([fasta, fai, dict]), a queue channel, and it now has three
    // consumers below instead of the original two. A no-op if it's already a value
    // channel.
    ref_fasta_ch = ref_fasta_ch.first()

    // COLLECT_READ_COUNTS runs over EVERY sample this subworkflow is handed, tumor and
    // normal alike, and that is deliberate even though nothing downstream reads the
    // normals' counts today: main.nf filters `sample_counts` to status == 1 before the
    // join into CNV_CORE, and CNV_CORE's take: interface has no normal_counts input at
    // all (SPEC.md Sec11.3). So for a 1-normal/N-tumor patient the normal's read-counting
    // is computed exactly once and then discarded.
    //
    // This is a known cost, not an oversight -- please do not "optimize" it away without
    // reading this first. Two reasons it stays: (1) this subworkflow is intentionally
    // role-agnostic -- it read-counts samples, and knows nothing about pairing, so
    // special-casing status here would push pair-level knowledge down into a sample-level
    // component purely to skip work; (2) the WDL pipeline this is a port of also collects
    // read counts for the normal without feeding them into the allelic-CN step, so
    // computing them uniformly preserves WDL parity and keeps the door open for a future
    // consumer (e.g. a normal-counts-based QC or tangent step) without re-plumbing.
    PREPROCESS_INTERVALS(ref_fasta_ch)

    // Scatter-gather (SPEC.md real-world finding, 2026-10-20): a real, high-coverage
    // CRAM pushed COLLECT_READ_COUNTS's wall-clock into the tens-of-hours range running
    // the whole genome in one task. SPLIT_INTERVALS shards the preprocessed interval
    // list -- either into read_counts_scatter_count size-balanced groups (the default;
    // GATK's own SplitIntervals already does this sensibly, alt/decoy contigs included,
    // no special-casing needed) or one shard per params.scatter_chromosomes entry
    // (read_counts_scatter_mode = 'per_chromosome', for maximum single-sample
    // parallelism when running one sample at a time rather than a large cohort, where
    // job-count at scale is the bigger risk -- see SPEC.md). GATHER_READ_COUNTS below
    // reassembles the per-shard outputs into the one merged file downstream consumers
    // (TANGENT_XY, AC_CALC) already expect.
    SPLIT_INTERVALS(
        ref_fasta_ch,
        PREPROCESS_INTERVALS.out,
        params.read_counts_scatter_mode,
        params.read_counts_scatter_count,
        params.scatter_chromosomes
    )

    // SPLIT_INTERVALS emits one List<Path> (all N shards for the single, shared
    // PREPROCESS_INTERVALS call); .flatten() turns that into N channel items so
    // .combine() below fans every sample out across every shard (N_samples x N_shards
    // total COLLECT_READ_COUNTS calls), not N_samples calls each holding a list of N.
    shards_ch = SPLIT_INTERVALS.out.flatten()

    COLLECT_READ_COUNTS(
        samples_ch.combine(shards_ch).combine(ref_fasta_ch)
            .map { meta, cram, crai, shard, fasta, fai, dict -> [ meta, meta.sample, cram, crai, shard, fasta, fai, dict ] }
    )

    // Groups each sample's N shard outputs back together by meta (identical object per
    // sample across the fan-out above, so grouping is exact) before merging them back
    // into one file per sample.
    GATHER_READ_COUNTS(COLLECT_READ_COUNTS.out.groupTuple())

    emit:
    sample_counts = GATHER_READ_COUNTS.out // [ meta, counts_file ]
}
