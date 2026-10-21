process ITER_ISEC {
    tag "${meta.participant_id}:${chr_group.join(',')}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'
    disk '2000 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(ref_dir_1000g), val(chr_group)

    output:
    tuple val(meta), path("chr_outputs/*_var.vcf.gz"), path("chr_outputs/*_var.vcf.gz.tbi")

    script:
    def chrs = chr_group.join(' ')
    """
    set -euo pipefail
    mkdir -p tmp/chr_outputs
    tar -xf ${ref_dir_1000g}

    for i in ${chrs}; do
      bcftools isec ${vcf} hg38/1000G/chr\${i}.1000G.genotypes.bcf -w 1 -O z -p tmp
      mv tmp/0002.vcf.gz     tmp/chr_outputs/${meta.participant_id}_chr\${i}_var.vcf.gz
      mv tmp/0002.vcf.gz.tbi tmp/chr_outputs/${meta.participant_id}_chr\${i}_var.vcf.gz.tbi
    done
    mv tmp/chr_outputs chr_outputs
    """

    stub:
    """
    mkdir -p chr_outputs
    for i in ${chr_group.join(' ')}; do
      touch "chr_outputs/${meta.participant_id}_chr\${i}_var.vcf.gz"
      touch "chr_outputs/${meta.participant_id}_chr\${i}_var.vcf.gz.tbi"
    done
    """
}
