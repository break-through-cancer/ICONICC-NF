process GATHER_READ_COUNTS {
    tag "${meta.sample}"
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 1
    memory '2 GB'

    input:
    tuple val(meta), path(shard_counts)

    output:
    tuple val(meta), path("${meta.sample}_counts.txt")

    script:
    // Nextflow task completion order is not genomic order, so this cannot assume
    // ${shard_counts} arrives contig-sorted -- every shard shares the identical
    // @HD/@SQ dictionary header (same -R reference for every shard), so the first
    // shard's own @SQ block gives a reliable contig rank (reference order, not
    // alphabetical: chr10 sorts before chr2 alphabetically, which would silently
    // corrupt downstream per-bin CNV math if trusted). Sort every shard's data rows by
    // (contig rank, START) before writing the merged file, rather than trusting file
    // or row arrival order.
    """
    set -euo pipefail
    shards="${shard_counts}"
    first=\$(echo \$shards | awk '{print \$1}')

    grep '^@' "\$first" > ${meta.sample}_counts.txt
    awk '!/^@/{print; exit}' "\$first" >> ${meta.sample}_counts.txt

    grep '^@SQ' "\$first" | awk '{for (i=1;i<=NF;i++) if (\$i ~ /^SN:/) {sub("SN:","",\$i); print \$i}}' > contig_order.txt

    for f in \$shards; do
      awk -v OFS='\t' 'NR==FNR {rank[\$1]=FNR; next} !/^@/ && \$1!="CONTIG" {print rank[\$1], \$0}' contig_order.txt "\$f"
    done | sort -k1,1n -k3,3n | cut -f2- >> ${meta.sample}_counts.txt
    """

    stub:
    """
    touch ${meta.sample}_counts.txt
    """
}
