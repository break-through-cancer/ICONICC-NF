include { SPLIT_INTERVALS                                        } from '../../../modules/local/split_intervals'
include { GATK4_GENOMICSDBIMPORT                                 } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATHER_GENOTYPED_VCF                                   } from '../../../modules/local/gather_genotyped_vcf'
include { VARIANT_FILTRATION                                     } from '../../../modules/local/variant_filtration'
include { MAKE_SITES_ONLY_VCF                                    } from '../../../modules/local/make_sites_only_vcf'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_SNP   } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as VARIANTRECALIBRATOR_INDEL } from '../../../modules/nf-core/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as APPLYVQSR_INDEL                     } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as APPLYVQSR_SNP                       } from '../../../modules/nf-core/gatk4/applyvqsr/main'

workflow JOINT_GENOTYPING {
    take:
    gvcf_pairs_ch        // [ pair_meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi ]
    ref_fasta_ch         // [ fasta, fai, dict ]
    calling_intervals_ch // [ intervals ] -- value channel
    dbsnp_ch             // [ dbsnp_vcf, dbsnp_vcf_index ]
    snp_resources_ch     // [ [vcf,vcf,vcf,vcf], [tbi,tbi,tbi,tbi], [label,label,label,label] ]
    indel_resources_ch   // [ [vcf,vcf,vcf], [tbi,tbi,tbi], [label,label,label] ]

    main:
    // .first() hardens these as value channels: each is documented as a single-value
    // channel in the take: comments above, but nothing previously enforced that. If a
    // caller passes a queue channel holding one item instead (e.g. main.nf's ref_fasta_ch
    // built via Channel.of([fasta, fai, dict])), every step after GenomicsDBImport would
    // otherwise run exactly once total regardless of how many pairs are in flight --
    // silently dropping every pair after the first, with no error. .first() is a no-op
    // when the channel is already a value channel.
    ref_fasta_ch       = ref_fasta_ch.first()
    dbsnp_ch           = dbsnp_ch.first()
    snp_resources_ch   = snp_resources_ch.first()
    indel_resources_ch = indel_resources_ch.first()

    // Interval scatter for GenomicsDBImport/GenotypeGVCFs only (SPEC.md Sec12.6, reversed
    // 2026-10-22): a real Cirro run showed these two steps running serially for 246 minutes
    // -- ~40% of the full pipeline's critical path -- for a single tumor/normal pair, against
    // the original no-scatter decision's assumption that scatter "buys little" at that scale.
    // Reuses SAMPLE_READ_COUNTS's own SPLIT_INTERVALS/gather pattern rather than inventing new
    // scatter infrastructure. VARIANT_FILTRATION onward stays a single whole-genome call per
    // pair -- VQSR's Gaussian mixture model needs genome-wide site statistics, so it cannot be
    // scattered the same way.
    SPLIT_INTERVALS(
        ref_fasta_ch,
        calling_intervals_ch,
        params.joint_genotyping_scatter_mode,
        params.joint_genotyping_scatter_count,
        params.scatter_chromosomes
    )
    shards_ch = SPLIT_INTERVALS.out.flatten()

    // Each shard gets its own id (participant + shard name) so GenomicsDBImport/GenotypeGVCFs
    // outputs don't collide across shards, and so the shard's own interval file can be
    // rejoined onto GenomicsDBImport's output below by that id instead of a second combine()
    // (which would re-multiply pairs x shards a second time).
    pair_shard_ch = gvcf_pairs_ch
        .combine(shards_ch)
        .map { meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi, shard ->
            [ meta + [id: "${meta.participant_id}.shard_${shard.baseName}"], normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi, shard ]
        }

    gendb_input = pair_shard_ch
        .map { shard_meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi, shard ->
            [ shard_meta, [normal_gvcf, tumor_gvcf], [normal_tbi, tumor_tbi], shard, [], [] ]
        }

    GATK4_GENOMICSDBIMPORT(gendb_input, false, false, false)

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb
        .join(pair_shard_ch.map { shard_meta, _ng, _nt, _tg, _tt, shard -> [ shard_meta, shard ] }, failOnDuplicate: true, failOnMismatch: true)
        .map { shard_meta, workspace, shard -> [ shard_meta, workspace, [], shard, [] ] }

    GATK4_GENOTYPEGVCFS(
        genotype_input,
        ref_fasta_ch.map { fasta, _fai, _dict -> [ [id: 'reference'], fasta ] },
        ref_fasta_ch.map { _fasta, fai, _dict -> [ [id: 'reference'], fai ] },
        ref_fasta_ch.map { _fasta, _fai, dict -> [ [id: 'reference'], dict ] },
        dbsnp_ch.map { vcf, _tbi -> [ [id: 'reference'], vcf ] },
        dbsnp_ch.map { _vcf, tbi -> [ [id: 'reference'], tbi ] }
    )

    // Regroup each pair's N per-shard VCFs back together (shard_meta - 'id' reconstructs the
    // original pair_meta, the same convention this subworkflow's own emit: block below uses)
    // and reassemble them into one genome-wide VCF before VariantFiltration/VQSR.
    genotyped_shards_ch = GATK4_GENOTYPEGVCFS.out.vcf
        .join(GATK4_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { shard_meta, vcf, tbi -> [ shard_meta - shard_meta.subMap('id'), vcf, tbi ] }
        .groupTuple()

    GATHER_GENOTYPED_VCF(genotyped_shards_ch)

    // ApplyRecalibration's actual input -- and therefore the final output VCF -- is the
    // ExcessHet-filtered VCF with VQSR applied on top, not the raw GenotypeGVCFs output
    // (verified against ICONICC_paired.3.wdl directly, SPEC.md Sec12.3). MakeSitesOnlyVcf's
    // output is training input only and is never the final VCF.
    VARIANT_FILTRATION(GATHER_GENOTYPED_VCF.out)
    MAKE_SITES_ONLY_VCF(VARIANT_FILTRATION.out)

    fasta_bare = ref_fasta_ch.map { fasta, _fai, _dict -> fasta }
    fai_bare   = ref_fasta_ch.map { _fasta, fai, _dict -> fai }
    dict_bare  = ref_fasta_ch.map { _fasta, _fai, dict -> dict }

    VARIANTRECALIBRATOR_SNP(
        MAKE_SITES_ONLY_VCF.out,
        snp_resources_ch.map { vcfs, _tbis, _labels -> vcfs },
        snp_resources_ch.map { _vcfs, tbis, _labels -> tbis },
        snp_resources_ch.map { _vcfs, _tbis, labels -> labels },
        fasta_bare, fai_bare, dict_bare
    )

    VARIANTRECALIBRATOR_INDEL(
        MAKE_SITES_ONLY_VCF.out,
        indel_resources_ch.map { vcfs, _tbis, _labels -> vcfs },
        indel_resources_ch.map { _vcfs, tbis, _labels -> tbis },
        indel_resources_ch.map { _vcfs, _tbis, labels -> labels },
        fasta_bare, fai_bare, dict_bare
    )

    // Chain order verified against ICONICC_paired.3.wdl's ApplyRecalibration task
    // (JointGeno_Tasks.1.wdl:586-634): INDEL runs first, SNP second on INDEL's output --
    // the reverse of sarek's own subworkflow (SPEC.md Sec12.3). Do not swap this order.
    applyvqsr_indel_input = VARIANT_FILTRATION.out
        .join(VARIANTRECALIBRATOR_INDEL.out.recal, failOnDuplicate: true, failOnMismatch: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.idx, failOnDuplicate: true, failOnMismatch: true)
        .join(VARIANTRECALIBRATOR_INDEL.out.tranches, failOnDuplicate: true, failOnMismatch: true)

    APPLYVQSR_INDEL(applyvqsr_indel_input, fasta_bare, fai_bare, dict_bare)

    applyvqsr_snp_input = APPLYVQSR_INDEL.out.vcf
        .join(APPLYVQSR_INDEL.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .join(VARIANTRECALIBRATOR_SNP.out.recal, failOnDuplicate: true, failOnMismatch: true)
        .join(VARIANTRECALIBRATOR_SNP.out.idx, failOnDuplicate: true, failOnMismatch: true)
        .join(VARIANTRECALIBRATOR_SNP.out.tranches, failOnDuplicate: true, failOnMismatch: true)

    APPLYVQSR_SNP(applyvqsr_snp_input, fasta_bare, fai_bare, dict_bare)

    emit:
    joint_vcf = APPLYVQSR_SNP.out.vcf.join(APPLYVQSR_SNP.out.tbi, failOnDuplicate: true, failOnMismatch: true)
        .map { meta, vcf, tbi -> [ meta - meta.subMap('id'), vcf, tbi ] } // [ pair_meta, vcf, tbi ]
}
