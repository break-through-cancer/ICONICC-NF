process AC_CALC {
    tag "${meta.participant_id}"
    container params.cnv_r_container
    memory '4 GB'

    input:
    tuple val(meta), path(tumor_counts), path(phased_files, stageAs: 'phased/*'), path(ad_files, stageAs: 'ad/*'), path(hetinfo_files, stageAs: 'hetinfo/*')

    output:
    tuple val(meta), path("*_read_depth.txt"), path("*_allelic_counts.txt"), path("*_plots.pdf")

    script:
    def agg_bins = meta.agg_bins ? 'TRUE' : 'FALSE'
    def binsize_arg = meta.desired_binsize ? "${meta.desired_binsize}" : ''
    """
    set -euo pipefail
    ac_calc.R \\
        ${meta.participant_id} ${tumor_counts} ${agg_bins} ${binsize_arg} \\
        ${meta.normal_sample_id} ${meta.tumor_sample_id} \\
        ${phased_files} ${ad_files} ${hetinfo_files}
    """
}
