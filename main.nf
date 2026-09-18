nextflow.enable.dsl = 2

include { CNV_CORE } from './subworkflows/local/cnv_core/main'

workflow {
    def required_params = [
        'participant_id', 'normal_sample_id', 'tumor_sample_id', 'sample_sex',
        'joint_vcf', 'joint_vcf_index', 'tumor_cram', 'tumor_crai', 'normal_cram', 'normal_crai'
    ]
    def missing = required_params.findAll { !params[it] }
    if (missing) {
        error "Missing required parameter(s): ${missing.join(', ')}. Pass them via --<param_name>."
    }

    // One uniform per-participant meta, reused unchanged everywhere CNV_CORE needs it --
    // see the note in subworkflows/local/cnv_core/main.nf's `take:` block for why this
    // has to stay a single identical map across every channel for .join() to work.
    meta = [
        participant_id: params.participant_id,
        normal_sample_id: params.normal_sample_id,
        tumor_sample_id: params.tumor_sample_id,
        sample_sex: params.sample_sex,
        scale_method: 'median', fix_outliers: false, trailingN: 5, n_dim: 10,
        smooth_region: 10, outlier_SD: 4, smooth_SD: 2, alpha: 0.01, min_width: 2,
        undo_split: 'none', undo_prune: 0.05, undo_SD: 2, recovery_thresh: 20,
        use_het: true, miss_thresh: 1.0, agg_bins: false, desired_binsize: null
    ]

    // checkIfExists on every input: a wrong or unstaged path fails here, at
    // startup, naming the offending file -- instead of being staged as a
    // dangling reference and surfacing much later as an opaque bcftools/htslib
    // error from inside a task container.
    joint_vcf_ch   = Channel.of([ meta,
        file(params.joint_vcf, checkIfExists: true),
        file(params.joint_vcf_index, checkIfExists: true) ])
    tumor_cram_ch  = Channel.of([ meta,
        file(params.tumor_cram, checkIfExists: true),
        file(params.tumor_crai, checkIfExists: true) ])
    normal_cram_ch = Channel.of([ meta,
        file(params.normal_cram, checkIfExists: true),
        file(params.normal_crai, checkIfExists: true) ])
    ref_fasta_ch   = Channel.of([
        file(params.ref_fasta, checkIfExists: true),
        file(params.ref_fasta_fai, checkIfExists: true),
        file(params.ref_fasta_dict, checkIfExists: true) ])

    CNV_CORE(
        joint_vcf_ch, tumor_cram_ch, normal_cram_ch, ref_fasta_ch,
        Channel.value(file(params.ref_dir_1000g, checkIfExists: true)),
        Channel.value(file(params.genetic_map, checkIfExists: true)),
        Channel.value(file(params.pon, checkIfExists: true)),
        Channel.value(file(params.ref_clinical, checkIfExists: true)),
        Channel.value(file(params.cohort_clinical, checkIfExists: true))
    )
}
