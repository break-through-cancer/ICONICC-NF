process SUBSET_PAIR_VCF {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.pair_subset.vcf.gz"), path("*.pair_subset.vcf.gz.tbi"), env('NORMAL_SM'), env('TUMOR_SM')

    script:
    // The joint VCF's #CHROM sample columns do not reliably equal the bare
    // normal_sample_id/tumor_sample_id: parabricks_bqsr embeds each CRAM's read group as
    // SM:${patient}_${sample}, but the patient value used at alignment time is not
    // guaranteed to match whatever Cirro's own metadata later records as this sample's
    // patient -- confirmed against two real runs, each embedding a different prefix for
    // the same bare sample id. Resolve the actual two column names from the VCF itself
    // (exact match, or a "_<id>" suffix) instead of reconstructing them, and hard-error
    // on anything but exactly one match per side.
    """
    set -euo pipefail
    match_sample() {
        bcftools query -l ${vcf} | awk -v id="\$1" '
            \$0 == id { print; next }
            { suffix = "_" id; n = length(\$0); m = length(suffix)
              if (n > m && substr(\$0, n-m+1) == suffix) print }
        '
    }
    normal_sm=\$(match_sample '${meta.normal_sample_id}')
    tumor_sm=\$(match_sample '${meta.tumor_sample_id}')
    normal_n=\$(printf '%s\\n' "\$normal_sm" | grep -c .)
    tumor_n=\$(printf '%s\\n' "\$tumor_sm" | grep -c .)
    if [ "\$normal_n" -ne 1 ]; then
        echo "Expected exactly one VCF sample matching normal id '${meta.normal_sample_id}', found: [\$normal_sm]" >&2
        exit 1
    fi
    if [ "\$tumor_n" -ne 1 ]; then
        echo "Expected exactly one VCF sample matching tumor id '${meta.tumor_sample_id}', found: [\$tumor_sm]" >&2
        exit 1
    fi
    bcftools view -s "\$normal_sm,\$tumor_sm" ${vcf} \\
        -O z -o ${meta.participant_id}.pair_subset.vcf.gz
    bcftools index --tbi ${meta.participant_id}.pair_subset.vcf.gz

    # Downstream file-naming (QUERY_HET_AD etc.) derives from these same resolved VCF
    # sample names, not meta.normal_sample_id/tumor_sample_id -- export them for AC_CALC.
    export NORMAL_SM="\$normal_sm"
    export TUMOR_SM="\$tumor_sm"
    """
}
