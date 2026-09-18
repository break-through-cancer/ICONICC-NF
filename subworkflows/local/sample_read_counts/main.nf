include { PREPROCESS_INTERVALS } from '../../../modules/local/preprocess_intervals'
include { COLLECT_READ_COUNTS  } from '../../../modules/local/collect_read_counts'

workflow SAMPLE_READ_COUNTS {
    take:
    samples_ch    // [ meta(patient,sample,sex,status), cram, crai ]
    ref_fasta_ch  // [ fasta, fai, dict ]

    main:
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

    COLLECT_READ_COUNTS(
        samples_ch.combine(PREPROCESS_INTERVALS.out).combine(ref_fasta_ch)
            .map { meta, cram, crai, intervals, fasta, fai, dict -> [ meta, meta.sample, cram, crai, intervals, fasta, fai, dict ] }
    )

    emit:
    sample_counts = COLLECT_READ_COUNTS.out // [ meta, counts_file ]
}
