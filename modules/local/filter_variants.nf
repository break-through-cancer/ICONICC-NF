process FILTER_VARIANTS {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_HC_var.vcf.gz"), path("*_HC_var.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    bcftools filter -i 'SUM(FORMAT/AD[*:1])>=5' ${vcf} -O z -o ${meta.participant_id}_HC_var.vcf.gz
    bcftools index --tbi ${meta.participant_id}_HC_var.vcf.gz
    """
}
