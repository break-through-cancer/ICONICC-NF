process PREPROCESS_INTERVALS {
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 1
    memory '2 GB'

    input:
    tuple path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    output:
    path("wgs.preprocessed.interval_list")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx1500m" PreprocessIntervals \\
        --reference ${ref_fasta} \\
        --padding 0 \\
        --bin-length 25000 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output wgs.preprocessed.interval_list
    """

    stub:
    """
    touch wgs.preprocessed.interval_list
    """
}
