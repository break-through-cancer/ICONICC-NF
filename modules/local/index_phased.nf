process INDEX_PHASED {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(phased_vcfs)

    output:
    tuple val(meta), path("*_het_eagle2phased.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    for path in ${phased_vcfs}; do
      bcftools index --tbi "\${path}"
    done
    """
}
