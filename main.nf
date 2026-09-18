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

    joint_vcf_ch   = Channel.of([ meta, file(params.joint_vcf), file(params.joint_vcf_index) ])
    tumor_cram_ch  = Channel.of([ meta, file(params.tumor_cram), file(params.tumor_crai) ])
    normal_cram_ch = Channel.of([ meta, file(params.normal_cram), file(params.normal_crai) ])
    ref_fasta_ch   = Channel.of([ file(params.ref_fasta), file(params.ref_fasta_fai), file(params.ref_fasta_dict) ])

    CNV_CORE(
        joint_vcf_ch, tumor_cram_ch, normal_cram_ch, ref_fasta_ch,
        Channel.value(file(params.ref_dir_1000g)),
        Channel.value(file(params.genetic_map)),
        Channel.value(file(params.pon)),
        Channel.value(file(params.ref_clinical)),
        Channel.value(file(params.cohort_clinical))
    )
}
