process SUBSET_PAIR_VCF {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.pair_subset.vcf.gz"), path("*.pair_subset.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    bcftools view -s ${meta.normal_sample_id},${meta.tumor_sample_id} ${vcf} \\
        -O z -o ${meta.participant_id}.pair_subset.vcf.gz
    bcftools index --tbi ${meta.participant_id}.pair_subset.vcf.gz
    """
}
