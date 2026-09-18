process SEGMENT {
    tag "${meta.participant_id}"
    container params.cnv_r_container
    memory '10 GB'

    input:
    tuple val(meta), path(allelic_counts), path(tangent_normalized)

    output:
    tuple val(meta), path("*.seg.txt"), path("*_processed_counts.txt"), path("*_segmentPlots.pdf"), path("*_segment_cts.txt")

    script:
    """
    set -euo pipefail
    segment.R \\
        --al_cts ${allelic_counts} --tangent_cts ${tangent_normalized} \\
        --sample_sex ${meta.sample_sex} --participant_id ${meta.participant_id} \\
        --smth_reg ${meta.smooth_region} --outlier_SD ${meta.outlier_SD} --smooth_SD ${meta.smooth_SD} \\
        --alpha ${meta.alpha} --min_width ${meta.min_width} --undo_split ${meta.undo_split} \\
        --undo_prune ${meta.undo_prune} --undo_SD ${meta.undo_SD} \\
        --recovery_thresh ${meta.recovery_thresh} --use_hetinfo ${meta.use_het} --miss_thresh ${meta.miss_thresh}
    """
}
