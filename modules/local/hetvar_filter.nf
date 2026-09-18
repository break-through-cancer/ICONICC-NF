process HETVAR_FILTER {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(chr_vcfs), path(chr_vcf_indices)

    output:
    tuple val(meta), path("*_het.vcf.gz"), path("*_het.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    for chr in ${chr_vcfs}; do
      filename=\$(basename "\${chr}" "_var.vcf.gz")"_het.vcf.gz"
      bcftools filter -i 'GT[1]="het" & FORMAT/AD[1:0]>=5 & FORMAT/AD[1:1]>=5' \\
        "\${chr}" -O z -o "\${filename}"
      bcftools index --tbi "\${filename}"
    done
    """
}
