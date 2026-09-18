// main.nf
nextflow.enable.dsl = 2

include { SAMPLESHEET_TO_PAIRS   } from './subworkflows/local/samplesheet_to_pairs/main'
include { SAMPLE_READ_COUNTS     } from './subworkflows/local/sample_read_counts/main'
include { SAMPLE_GERMLINE_CALLING } from './subworkflows/local/sample_germline_calling/main'
include { JOINT_GENOTYPING       } from './subworkflows/local/joint_genotyping/main'
include { CNV_CORE               } from './subworkflows/local/cnv_core/main'

workflow {
    if (!params.input) {
        error "Missing required parameter: input (path to samplesheet CSV). See README.md and assets/samplesheet.csv."
    }

    SAMPLESHEET_TO_PAIRS(params.input)

    ref_fasta_ch = Channel.of([
        file(params.ref_fasta, checkIfExists: true),
        file(params.ref_fasta_fai, checkIfExists: true),
        file(params.ref_fasta_dict, checkIfExists: true)
    ])

    SAMPLE_READ_COUNTS(SAMPLESHEET_TO_PAIRS.out.samples, ref_fasta_ch)

    calling_intervals_ch = Channel.value(file(params.calling_intervals, checkIfExists: true))

    SAMPLE_GERMLINE_CALLING(
        SAMPLESHEET_TO_PAIRS.out.samples,
        Channel.value(file(params.ref_fasta, checkIfExists: true)),
        calling_intervals_ch
    )

    // Join each pair to its normal's and tumor's GVCFs by (patient, sample) -- the same
    // pattern as tumor_counts_ch below, done for both roles (SPEC.md Sec12.1).
    sample_gvcfs_by_key = SAMPLE_GERMLINE_CALLING.out.sample_gvcfs
        .map { meta, gvcf, tbi -> [ [patient: meta.patient, sample: meta.sample], gvcf, tbi ] }

    gvcf_pairs_ch = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta, _jv, _jvi -> [ [patient: pair_meta.patient, sample: pair_meta.normal_sample_id], pair_meta ] }
        .join(sample_gvcfs_by_key, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _key, pair_meta, normal_gvcf, normal_tbi -> [ [patient: pair_meta.patient, sample: pair_meta.tumor_sample_id], pair_meta, normal_gvcf, normal_tbi ] }
        .join(sample_gvcfs_by_key, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _key, pair_meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi -> [ pair_meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi ] }

    dbsnp_ch = Channel.value([
        file(params.dbsnp_vcf, checkIfExists: true),
        file(params.dbsnp_vcf_index, checkIfExists: true)
    ])

    def hapmap_vcf   = file(params.hapmap_resource_vcf, checkIfExists: true)
    def omni_vcf     = file(params.omni_resource_vcf, checkIfExists: true)
    def kg_vcf       = file(params.one_thousand_genomes_resource_vcf, checkIfExists: true)
    def dbsnp_vcf    = file(params.dbsnp_vcf, checkIfExists: true)
    def mills_vcf    = file(params.mills_resource_vcf, checkIfExists: true)
    def axiompoly_vcf = file(params.axiomPoly_resource_vcf, checkIfExists: true)

    snp_resources_ch = Channel.value([
        [ hapmap_vcf, omni_vcf, kg_vcf, dbsnp_vcf ],
        [ file(params.hapmap_resource_vcf_index, checkIfExists: true), file(params.omni_resource_vcf_index, checkIfExists: true), file(params.one_thousand_genomes_resource_vcf_index, checkIfExists: true), file(params.dbsnp_vcf_index, checkIfExists: true) ],
        [
            "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap_vcf.name}",
            "--resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni_vcf.name}",
            "--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${kg_vcf.name}",
            "--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${dbsnp_vcf.name}"
        ]
    ])

    indel_resources_ch = Channel.value([
        [ mills_vcf, axiompoly_vcf, dbsnp_vcf ],
        [ file(params.mills_resource_vcf_index, checkIfExists: true), file(params.axiomPoly_resource_vcf_index, checkIfExists: true), file(params.dbsnp_vcf_index, checkIfExists: true) ],
        [
            "--resource:mills,known=false,training=true,truth=true,prior=12.0 ${mills_vcf.name}",
            "--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 ${axiompoly_vcf.name}",
            "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_vcf.name}"
        ]
    ])

    JOINT_GENOTYPING(
        gvcf_pairs_ch,
        ref_fasta_ch,
        calling_intervals_ch,
        dbsnp_ch,
        snp_resources_ch,
        indel_resources_ch
    )

    joint_vcf_ch = JOINT_GENOTYPING.out.joint_vcf // [ pair_meta, vcf, tbi ]

    tumor_counts_ch = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta, _jv, _jvi ->
            [ [patient: pair_meta.patient, sample: pair_meta.tumor_sample_id], pair_meta ]
        }
        .join(
            SAMPLE_READ_COUNTS.out.sample_counts
                .filter { sample_meta, _counts -> sample_meta.status == 1 }
                .map { sample_meta, counts -> [ [patient: sample_meta.patient, sample: sample_meta.sample], counts ] },
            by: 0, failOnDuplicate: true, failOnMismatch: true
        )
        .map { _key, pair_meta, counts -> [ pair_meta, counts ] }

    CNV_CORE(
        joint_vcf_ch,
        tumor_counts_ch,
        Channel.value(file(params.ref_dir_1000g, checkIfExists: true)),
        Channel.value(file(params.genetic_map, checkIfExists: true)),
        Channel.value(file(params.pon, checkIfExists: true)),
        Channel.value(file(params.ref_clinical, checkIfExists: true)),
        Channel.value(file(params.cohort_clinical, checkIfExists: true))
    )
}
