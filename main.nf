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

    // Join each pair to its normal's and tumor's GVCFs by (patient, sample). The two
    // roles need different join semantics, so this is deliberately NOT the same
    // .join()-on-both-sides pattern as tumor_counts_ch below:
    //  - The tumor side really is 1:1 (SAMPLESHEET_TO_PAIRS.out.pairs has exactly one
    //    row per pair, and a given tumor sample belongs to exactly one pair), so
    //    .join(failOnDuplicate: true, failOnMismatch: true) is correct there, same as
    //    tumor_counts_ch.
    //  - The normal side is 1:N -- SAMPLESHEET_TO_PAIRS fans one normal out to every
    //    one of that patient's tumors (see its own test "one normal + two tumors
    //    produces two pairs sharing the same normal"), so the pairs-side channel can
    //    legitimately repeat the same normal key more than once. A plain .join() would
    //    treat that repetition as a duplicate-key error; .combine(by: 0) is the
    //    operator that allows one right-hand row to fan out to N left-hand rows with
    //    the same key. Unlike .join(..., failOnMismatch: true), .combine() has no way
    //    to fail loudly on an unmatched key -- it just silently drops that row -- so
    //    the explicit key-set guard below (same idiom as SAMPLESHEET_TO_PAIRS's own
    //    orphan-tumor guard) exists specifically to catch that case instead of quietly
    //    losing a pair (SPEC.md Sec12.1: fail loudly, never skip silently).
    normal_gvcfs_by_key = SAMPLE_GERMLINE_CALLING.out.sample_gvcfs
        .filter { meta, _gvcf, _tbi -> meta.status == 0 }
        .map { meta, gvcf, tbi -> [ [patient: meta.patient, sample: meta.sample], gvcf, tbi ] }

    tumor_gvcfs_by_key = SAMPLE_GERMLINE_CALLING.out.sample_gvcfs
        .filter { meta, _gvcf, _tbi -> meta.status == 1 }
        .map { meta, gvcf, tbi -> [ [patient: meta.patient, sample: meta.sample], gvcf, tbi ] }

    normal_gvcf_keys_ch = normal_gvcfs_by_key
        .map { key, _gvcf, _tbi -> key }
        .toList()
        .map { [keys: it] } // wrapped in a Map so .combine() below passes it through as
                             // one opaque element instead of flattening the list

    pairs_by_normal_key = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta -> [ [patient: pair_meta.patient, sample: pair_meta.normal_sample_id], pair_meta ] }
        .combine(normal_gvcf_keys_ch)
        .map { key, pair_meta, wrapper ->
            if (!(key in wrapper.keys)) {
                error "Pair '${pair_meta.participant_id}' references normal sample '${pair_meta.normal_sample_id}' with no matching GVCF from SAMPLE_GERMLINE_CALLING"
            }
            [ key, pair_meta ]
        }

    gvcf_pairs_ch = pairs_by_normal_key
        .combine(normal_gvcfs_by_key, by: 0)
        .map { _key, pair_meta, normal_gvcf, normal_tbi -> [ [patient: pair_meta.patient, sample: pair_meta.tumor_sample_id], pair_meta, normal_gvcf, normal_tbi ] }
        .join(tumor_gvcfs_by_key, by: 0, failOnDuplicate: true, failOnMismatch: true)
        .map { _key, pair_meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi -> [ pair_meta, normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi ] }

    // Resolved once each: dbsnp_vcf/dbsnp_vcf_index feed dbsnp_ch AND appear a second
    // time as the last resource in both snp_resources_ch and indel_resources_ch, so
    // hoisting avoids re-resolving the same two params four times over.
    def hapmap_vcf          = file(params.hapmap_resource_vcf, checkIfExists: true)
    def hapmap_vcf_index    = file(params.hapmap_resource_vcf_index, checkIfExists: true)
    def omni_vcf             = file(params.omni_resource_vcf, checkIfExists: true)
    def omni_vcf_index       = file(params.omni_resource_vcf_index, checkIfExists: true)
    def kg_vcf               = file(params.one_thousand_genomes_resource_vcf, checkIfExists: true)
    def kg_vcf_index         = file(params.one_thousand_genomes_resource_vcf_index, checkIfExists: true)
    def dbsnp_vcf            = file(params.dbsnp_vcf, checkIfExists: true)
    def dbsnp_vcf_index      = file(params.dbsnp_vcf_index, checkIfExists: true)
    def mills_vcf            = file(params.mills_resource_vcf, checkIfExists: true)
    def mills_vcf_index      = file(params.mills_resource_vcf_index, checkIfExists: true)
    def axiompoly_vcf        = file(params.axiomPoly_resource_vcf, checkIfExists: true)
    def axiompoly_vcf_index  = file(params.axiomPoly_resource_vcf_index, checkIfExists: true)

    dbsnp_ch = Channel.value([ dbsnp_vcf, dbsnp_vcf_index ])

    // Each resources_ch below is [ vcfs, tbis, labels ] with all three lists in the
    // same order -- position N of every list must refer to the same resource. One
    // variable per line, in matching order across all three lists, so a transposition
    // is visible at a glance instead of hiding in one long line of file() calls.
    snp_resources_ch = Channel.value([
        [ hapmap_vcf, omni_vcf, kg_vcf, dbsnp_vcf ],
        [ hapmap_vcf_index, omni_vcf_index, kg_vcf_index, dbsnp_vcf_index ],
        [
            "--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap_vcf.name}",
            "--resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni_vcf.name}",
            "--resource:1000G,known=false,training=true,truth=false,prior=10.0 ${kg_vcf.name}",
            "--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ${dbsnp_vcf.name}"
        ]
    ])

    indel_resources_ch = Channel.value([
        [ mills_vcf, axiompoly_vcf, dbsnp_vcf ],
        [ mills_vcf_index, axiompoly_vcf_index, dbsnp_vcf_index ],
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
        .map { pair_meta ->
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
