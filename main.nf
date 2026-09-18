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

    // The VCF side goes straight into CNV_CORE, unjoined. CNV_CORE's phasing chain
    // (SUBSET_PAIR_VCF -> FILTER_VARIANTS -> ITER_ISEC -> HETVAR_FILTER -> EAGLE_PHASING)
    // consumes no read-counts data at all, so routing it through the counts join below
    // would gate hours of phasing behind that pair's tumor COLLECT_READ_COUNTS (a GATK
    // job over a full WGS CRAM) for no reason. Feeding it from `pairs` directly lets
    // phasing start as soon as the samplesheet is parsed, while read-counting runs in
    // parallel; the two sides rejoin inside CNV_CORE. That rejoin still matches because
    // Nextflow's .join() keys on map equality, not object identity -- the pair_meta here
    // and the pair_meta carried through the join below are the same value.
    joint_vcf_ch = SAMPLESHEET_TO_PAIRS.out.pairs // [ pair_meta, joint_vcf, joint_vcf_index ]

    // Join each pair to its tumor's already-computed (and, when a normal is shared
    // across several pairs, already-deduplicated) read counts. Keyed on a small
    // (patient, sample) map rather than the full meta, since pair_meta and sample_meta
    // are different shapes -- constructing the same-shaped key on both sides is what
    // makes the .join() match, the same principle as every meta-keyed join already in
    // this codebase (see subworkflows/local/cnv_core/main.nf's take: comment).
    // sample_counts holds both tumor and normal rows; filtered to status == 1 (tumor)
    // here since only the tumor's counts feed CNV_CORE, and because leaving normals in
    // would let a normal's [patient, sample] key coincidentally collide with a
    // differently-keyed lookup. failOnMismatch: true turns any pair whose tumor has no
    // matching counts entry into a hard, named error instead of a silent drop -- the
    // exact failure mode a .join() defaults to (remainder: false) on a non-match.
    tumor_counts_ch = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta, _jv, _jvi ->
            [ [patient: pair_meta.patient, sample: pair_meta.tumor_sample_id], pair_meta ]
        }
        .join(
            SAMPLE_READ_COUNTS.out.sample_counts
                .filter { sample_meta, _counts -> sample_meta.status == 1 }
                .map { sample_meta, counts -> [ [patient: sample_meta.patient, sample: sample_meta.sample], counts ] },
            by: 0, failOnDuplicate: true, failOnMismatch: true
        )
        .map { _key, pair_meta, counts -> [ pair_meta, counts ] }

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
