process INDEX_GVCF {
    tag "${meta.sample}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 1
    memory '2 GB'

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path(gvcf), path("${gvcf}.tbi")

    script:
    """
    set -euo pipefail
    bcftools index --tbi ${gvcf}
    """

    stub:
    """
    touch ${gvcf}.tbi
    """
}
