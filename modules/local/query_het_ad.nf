process QUERY_HET_AD {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    memory '4 GB'

    input:
    tuple val(meta), path(chr1_het), path(all_het), path(all_het_indices), path(all_eagle_phased), path(all_eagle_indices)

    output:
    tuple val(meta), path("*_HetInfo.txt"), path("*.AD.chr*.txt"), path("*.phasedGT.chr*.1.txt")

    script:
    """
    set -euo pipefail
    samples=\$(bcftools query -l ${chr1_het})

    for path in ${all_het}; do
      hetinfofile=\$(basename "\${path}" "_het.vcf.gz")"_HetInfo.txt"
      chr=\$(basename "\${path}" "_het.vcf.gz" | sed 's/^.*chr//')
      bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%INFO/DP\\n' "\${path}" > "\${hetinfofile}"
      for sample in \${samples}; do
        bcftools query -s \${sample} -f '[%AD{0}\\t%AD{1}\\n]' "\${path}" > "\${sample}.AD.chr\${chr}.txt"
      done
    done

    for phased in ${all_eagle_phased}; do
      chr=\$(basename "\${phased}" "_het_eagle2phased.vcf.gz" | sed 's/^.*chr//')
      for sample in \${samples}; do
        bcftools query -s \${sample} -f '[%GT\\n]' "\${phased}" > "\${sample}.phasedGT.chr\${chr}.1.txt"
      done
    done
    """
}
