process MAKE_SITES_ONLY_VCF {
    tag "${meta.participant_id}"
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 1
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("${meta.participant_id}.sites_only.vcf.gz"), path("${meta.participant_id}.sites_only.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xms3000m -Xmx3250m" \\
        MakeSitesOnlyVcf \\
        -I ${vcf} \\
        -O ${meta.participant_id}.sites_only.vcf.gz
    """
}
