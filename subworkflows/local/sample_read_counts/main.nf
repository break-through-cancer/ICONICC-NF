include { PREPROCESS_INTERVALS } from '../../../modules/local/preprocess_intervals'
include { COLLECT_READ_COUNTS  } from '../../../modules/local/collect_read_counts'

workflow SAMPLE_READ_COUNTS {
    take:
    samples_ch    // [ meta(patient,sample,sex,status), cram, crai ]
    ref_fasta_ch  // [ fasta, fai, dict ]

    main:
    PREPROCESS_INTERVALS(ref_fasta_ch)

    COLLECT_READ_COUNTS(
        samples_ch.combine(PREPROCESS_INTERVALS.out).combine(ref_fasta_ch)
            .map { meta, cram, crai, intervals, fasta, fai, dict -> [ meta, meta.sample, cram, crai, intervals, fasta, fai, dict ] }
    )

    emit:
    sample_counts = COLLECT_READ_COUNTS.out // [ meta, counts_file ]
}
