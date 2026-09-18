// main.nf
nextflow.enable.dsl = 2

include { SAMPLESHEET_TO_PAIRS } from './subworkflows/local/samplesheet_to_pairs/main'
include { SAMPLE_READ_COUNTS   } from './subworkflows/local/sample_read_counts/main'
include { CNV_CORE             } from './subworkflows/local/cnv_core/main'

workflow {
    if (!params.input) {
        error "Missing required parameter: input (path to samplesheet CSV). See README.md and assets/samplesheet.csv."
    }

    SAMPLESHEET_TO_PAIRS(params.input)

    ref_fasta_ch = Channel.of([
        file(params.ref_fasta, checkIfExists: true),
        file(params.ref_fasta_fai, checkIfExists: true),
        file(params.ref_fasta_dict, checkIfExists: true)
    ])

    SAMPLE_READ_COUNTS(SAMPLESHEET_TO_PAIRS.out.samples, ref_fasta_ch)

    // Join each pair to its tumor's already-computed (and, when a normal is shared
    // across several pairs, already-deduplicated) read counts. Keyed on a small
    // (patient, sample) map rather than the full meta, since pair_meta and sample_meta
    // are different shapes -- constructing the same-shaped key on both sides is what
    // makes the .join() match, the same principle as every meta-keyed join already in
    // this codebase (see subworkflows/local/cnv_core/main.nf's take: comment).
    tumor_counts_for_pairs = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta, jv, jvi ->
            [ [patient: pair_meta.patient, sample: pair_meta.tumor_sample_id], pair_meta, jv, jvi ]
        }
        .join(
            SAMPLE_READ_COUNTS.out.sample_counts
                .map { sample_meta, counts -> [ [patient: sample_meta.patient, sample: sample_meta.sample], counts ] },
            by: 0, failOnDuplicate: true
        )
        .map { _key, pair_meta, jv, jvi, counts -> [ pair_meta, jv, jvi, counts ] }

    joint_vcf_ch    = tumor_counts_for_pairs.map { pair_meta, jv, jvi, _counts -> [ pair_meta, jv, jvi ] }
    tumor_counts_ch = tumor_counts_for_pairs.map { pair_meta, _jv, _jvi, counts -> [ pair_meta, counts ] }

    CNV_CORE(
        joint_vcf_ch,
        tumor_counts_ch,
        Channel.value(file(params.ref_dir_1000g, checkIfExists: true)),
        Channel.value(file(params.genetic_map, checkIfExists: true)),
        Channel.value(file(params.pon, checkIfExists: true)),
        Channel.value(file(params.ref_clinical, checkIfExists: true)),
        Channel.value(file(params.cohort_clinical, checkIfExists: true))
    )
}
