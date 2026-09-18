include { PARABRICKS_HAPLOTYPECALLER } from '../../../modules/nf-core/parabricks/haplotypecaller/main'
include { INDEX_GVCF                 } from '../../../modules/local/index_gvcf'

workflow SAMPLE_GERMLINE_CALLING {
    take:
    samples_ch           // [ meta(patient,sample,sex,status), cram, crai ]
    ref_fasta_ch         // [ fasta ] -- single file, Parabricks needs no fai/dict
    calling_intervals_ch // [ intervals ] -- value channel, one Picard-style .intervals file

    main:
    // Dedup is free here, the same way it is in SAMPLE_READ_COUNTS: samples_ch already
    // has exactly one row per unique (patient, sample) -- a normal shared across several
    // tumor rows was already collapsed to one row by SAMPLESHEET_TO_PAIRS before this
    // subworkflow ever sees it, so it is called exactly once regardless of how many
    // tumors reference it (SPEC.md Sec12.1).
    calling_input = samples_ch
        .combine(calling_intervals_ch)
        .map { meta, cram, crai, intervals -> [ meta + [id: meta.sample], cram, crai, intervals ] }

    PARABRICKS_HAPLOTYPECALLER(
        calling_input,
        ref_fasta_ch.map { fasta -> [ [id: 'reference'], fasta ] }
    )

    INDEX_GVCF(PARABRICKS_HAPLOTYPECALLER.out.gvcf)

    emit:
    sample_gvcfs = INDEX_GVCF.out // [ meta, gvcf, gvcf_tbi ]
}
