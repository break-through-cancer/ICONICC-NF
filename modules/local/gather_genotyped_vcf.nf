process GATHER_GENOTYPED_VCF {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '8 GB'

    input:
    tuple val(meta), path(shard_vcfs), path(shard_tbis)

    output:
    tuple val(meta), path("${meta.participant_id}.genotyped.vcf.gz"), path("${meta.participant_id}.genotyped.vcf.gz.tbi")

    script:
    // SPLIT_INTERVALS's 'groups' mode balances total interval length across shards by
    // bin-packing possibly non-adjacent intervals together (modules/local/split_intervals.nf),
    // so the shards GenotypeGVCFs produces are not in genomic order relative to each other
    // even though each shard's own VCF is internally position-sorted -- the same chromosome
    // can appear in more than one shard, in non-adjacent blocks once concatenated in shard
    // order. bcftools concat refuses that by default ("chromosome block ... is not
    // contiguous"); -a/--allow-overlaps turns that check off. bcftools sort re-establishes
    // true genome-wide order afterward regardless of concat's own output order (same
    // reasoning as GATHER_READ_COUNTS's own dictionary-order sort).
    """
    set -euo pipefail
    bcftools concat --allow-overlaps --output-type u ${shard_vcfs} \\
        | bcftools sort --output-type z --temp-dir . --output-file ${meta.participant_id}.genotyped.vcf.gz
    bcftools index --tbi ${meta.participant_id}.genotyped.vcf.gz
    """

    stub:
    """
    touch ${meta.participant_id}.genotyped.vcf.gz
    touch ${meta.participant_id}.genotyped.vcf.gz.tbi
    """
}
