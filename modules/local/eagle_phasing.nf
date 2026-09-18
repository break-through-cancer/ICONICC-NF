process EAGLE_PHASING {
    tag "${meta.participant_id}"
    container 'docker.io/cbao/eagle@sha256:e837cf05d2821e2bd40d6e6f7dac9ef0a5c1575e744a0e7026d28bfac0cf3d66'
    cpus 16
    memory '128 GB'
    disk '2000 GB'

    input:
    tuple val(meta), path(het_vcfs), path(ref_dir_1000g), path(genetic_map)

    output:
    tuple val(meta), path("*_het_eagle2phased.vcf.gz")

    script:
    """
    set -euo pipefail
    tar -xf ${ref_dir_1000g}

    for path in ${het_vcfs}; do
      filename=\$(basename "\${path}" "_het.vcf.gz")"_het_eagle2phased"
      chr=\$(basename "\${path}" "_het.vcf.gz" | sed 's/^.*chr//')
      eagle \\
        --vcfTarget="\${path}" \\
        --vcfRef=hg38/1000G/chr\${chr}.1000G.genotypes.bcf \\
        --geneticMapFile=${genetic_map} \\
        --chrom=chr\${chr} \\
        --outPrefix="\${filename}" \\
        --numThreads=${task.cpus}
    done
    """
}
