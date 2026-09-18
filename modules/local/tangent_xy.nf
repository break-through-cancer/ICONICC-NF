process TANGENT_XY {
    container params.cnv_r_container
    memory '10 GB'

    input:
    tuple val(meta), path(tumor_counts), path(pon), path(ref_clinical), path(cohort_clinical)

    output:
    tuple val(meta), path("posttangent_log2.rds")

    script:
    """
    set -euo pipefail
    Rscript ${projectDir}/bin/tangent.R \\
        --pon ${pon} --tumor_cts ${tumor_counts} --ref_cli ${ref_clinical} --clin ${cohort_clinical} \\
        --scale-mthd ${meta.scale_method} --fix_outliers ${meta.fix_outliers} \\
        --trailingN ${meta.trailingN} --sample_sex ${meta.sample_sex} --n_dim ${meta.n_dim}
    """
}
