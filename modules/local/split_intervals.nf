process SPLIT_INTERVALS {
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 1
    memory '4 GB'

    input:
    tuple path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)
    path intervals
    val scatter_mode  // 'groups' or 'per_chromosome'
    val scatter_count // used only when scatter_mode == 'groups'
    val chromosomes   // used only when scatter_mode == 'per_chromosome' -- List<String>, no 'chr' prefix

    output:
    path("shards/*.interval_list")

    script:
    if (!(scatter_mode in ['groups', 'per_chromosome'])) {
        error "SPLIT_INTERVALS: scatter_mode must be 'groups' or 'per_chromosome', got '${scatter_mode}'"
    }
    // 'per_chromosome' only ever shards the chromosomes this pipeline's phasing panel
    // actually covers (params.scatter_chromosomes) -- not every contig in the reference.
    // Filtering the interval_list by contig directly (not GATK) because the input is
    // already a validated Picard interval_list; no engine/reference validation is needed
    // to subset it by the first column.
    if (scatter_mode == 'per_chromosome') {
        def chrList = chromosomes.collect { "chr${it}" }.join(' ')
        """
        set -euo pipefail
        mkdir -p shards
        grep '^@' ${intervals} > header.txt
        for chr in ${chrList}; do
          cat header.txt > "shards/\${chr}.interval_list"
          awk -F'\\t' -v c="\${chr}" '\$1 == c' ${intervals} >> "shards/\${chr}.interval_list"
        done
        """
    } else {
        """
        set -euo pipefail
        mkdir -p shards
        # Both flags are required to stop GATK merging adjacent bins before scattering
        # (confirmed empirically -- --interval-merging-rule alone still merged touching
        # intervals; GATK's own SplitIntervals source carries this exact combination for
        # this exact scenario: "apportion a set of small adjacent intervals ... among X
        # new interval lists").
        gatk SplitIntervals \\
            -R ${ref_fasta} \\
            -L ${intervals} \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \\
            --scatter-count ${scatter_count} \\
            -O shards
        for f in shards/*-scattered.interval_list; do
          mv "\$f" "\${f/-scattered/}"
        done
        """
    }

    stub:
    def n = scatter_mode == 'per_chromosome' ? chromosomes.size() : (scatter_count as int)
    """
    mkdir -p shards
    for i in \$(seq 1 ${n}); do touch "shards/shard_\${i}.interval_list"; done
    """
}
