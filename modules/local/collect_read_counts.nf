process COLLECT_READ_COUNTS {
    tag "${sample_id}:${intervals.baseName}"
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 2
    memory '10 GB'

    input:
    tuple val(meta), val(sample_id), path(cram), path(crai), path(intervals), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    output:
    // Filename includes the shard's own interval_list basename (e.g. chr21, or GATK
    // SplitIntervals' 0000/0001/...) so GATHER_READ_COUNTS can stage every sample's N
    // shards in one task without a filename collision -- every shard otherwise shares
    // the same sample_id and would otherwise all be named identically. meta.sample
    // already equals sample_id, so grouping downstream keys on meta, not a duplicate field.
    tuple val(meta), path("${sample_id}_${intervals.baseName}_counts.txt")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx9G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
        CollectReadCounts \\
        -R ${ref_fasta} \\
        -I ${cram} \\
        --read-index ${crai} \\
        -L ${intervals} \\
        -O ${sample_id}_${intervals.baseName}_counts.txt \\
        --format TSV \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter \\
        --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter \\
        --read-filter PairedReadFilter --read-filter FragmentLengthReadFilter --max-fragment-length 1000 \\
        --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateDifferentStrandReadFilter \\
        --read-filter NotSupplementaryAlignmentReadFilter --read-filter NotSecondaryAlignmentReadFilter \\
        --read-filter FirstOfPairReadFilter --read-filter MappingQualityReadFilter --minimum-mapping-quality 30 \\
        --read-filter OverclippedReadFilter --filter-too-short 50 --seconds-between-progress-updates 100
    """

    stub:
    """
    touch ${sample_id}_${intervals.baseName}_counts.txt
    """
}
