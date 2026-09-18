process COLLECT_READ_COUNTS {
    tag "${sample_id}"
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 2
    memory '10 GB'

    input:
    tuple val(meta), val(sample_id), path(cram), path(crai), path(intervals), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    output:
    tuple val(meta), path("${sample_id}_counts.txt")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx9G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
        CollectReadCounts \\
        -R ${ref_fasta} \\
        -I ${cram} \\
        --read-index ${crai} \\
        -L ${intervals} \\
        -O ${sample_id}_counts.txt \\
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
}
