process HETVAR_FILTER {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(chr_vcfs), path(chr_vcf_indices), val(normal_sm)

    output:
    tuple val(meta), path("*_het.vcf.gz"), path("*_het.vcf.gz.tbi")

    script:
    // bcftools expressions only accept a numeric sample index, not a name -- resolve the
    // normal sample's column position from each VCF's own header instead of assuming a
    // fixed index, since nothing guarantees sample column order stays put upstream.
    """
    set -euo pipefail
    for chr in ${chr_vcfs}; do
      filename=\$(basename "\${chr}" "_var.vcf.gz")"_het.vcf.gz"
      normal_idx=\$(bcftools query -l "\${chr}" | grep -nxF '${normal_sm}' | cut -d: -f1)
      normal_idx=\$((normal_idx - 1))
      bcftools filter -i "GT[\${normal_idx}]=\\"het\\" & FORMAT/AD[\${normal_idx}:0]>=5 & FORMAT/AD[\${normal_idx}:1]>=5" \\
        "\${chr}" -O z -o "\${filename}"
      bcftools index --tbi "\${filename}"
    done
    """
}
