include { GATK4_GENOMICSDBIMPORT                                 } from '../../../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                                    } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
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
    // No interval scatter (SPEC.md Sec12.6): every step below runs once per pair over the
    // whole calling_intervals region. This is why there is no BCFTOOLS_SORT/MergeVcfs step
    // -- GenotypeGVCFs already emits one complete, sorted VCF, nothing to reassemble.
    gendb_input = gvcf_pairs_ch
        .combine(calling_intervals_ch)
        .map { meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi, intervals ->
            [ meta + [id: meta.participant_id], [normal_gvcf, tumor_gvcf], [normal_tbi, tumor_tbi], intervals, [], [] ]
        }

    GATK4_GENOMICSDBIMPORT(gendb_input, false, false, false)

    genotype_input = GATK4_GENOMICSDBIMPORT.out.genomicsdb
        .combine(calling_intervals_ch)
        .map { meta, workspace, intervals -> [ meta, workspace, [], intervals, [] ] }

    GATK4_GENOTYPEGVCFS(
        genotype_input,
        ref_fasta_ch.map { fasta, _fai, _dict -> [ [id: 'reference'], fasta ] },
        ref_fasta_ch.map { _fasta, fai, _dict -> [ [id: 'reference'], fai ] },
        ref_fasta_ch.map { _fasta, _fai, dict -> [ [id: 'reference'], dict ] },
        dbsnp_ch.map { vcf, _tbi -> [ [id: 'reference'], vcf ] },
        dbsnp_ch.map { _vcf, tbi -> [ [id: 'reference'], tbi ] }
    )

    // ApplyRecalibration's actual input -- and therefore the final output VCF -- is the
    // ExcessHet-filtered VCF with VQSR applied on top, not the raw GenotypeGVCFs output
    // (verified against ICONICC_paired.3.wdl directly, SPEC.md Sec12.3). MakeSitesOnlyVcf's
    // output is training input only and is never the final VCF.
    VARIANT_FILTRATION(GATK4_GENOTYPEGVCFS.out.vcf.join(GATK4_GENOTYPEGVCFS.out.tbi, failOnDuplicate: true, failOnMismatch: true))
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
