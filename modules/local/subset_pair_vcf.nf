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
    // parabricks_bqsr -- this pipeline's sole alignment-entry source (README "Cirro
    // deployment") -- embeds each CRAM's read group as SM:${patient}_${sample}, not the
    // bare sample id (subworkflows/local/samplesheet_to_fastqs/main.nf in that repo).
    // GenomicsDBImport/GenotypeGVCFs carry that embedded name straight into the joint
    // VCF's own sample columns, so the subset-by-name here must match the same
    // convention, confirmed against a real run's #CHROM line.
    def normal_sm = "${meta.patient}_${meta.normal_sample_id}"
    def tumor_sm  = "${meta.patient}_${meta.tumor_sample_id}"
    """
    set -euo pipefail
    bcftools view -s ${normal_sm},${tumor_sm} ${vcf} \\
        -O z -o ${meta.participant_id}.pair_subset.vcf.gz
    bcftools index --tbi ${meta.participant_id}.pair_subset.vcf.gz
    """
}
