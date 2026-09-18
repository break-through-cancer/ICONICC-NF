workflow SAMPLESHEET_TO_PAIRS {
    take:
    samplesheet // path to samplesheet CSV (patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index)

    main:
    rows = Channel.fromPath(resolveRelativeToProjectDir(samplesheet), checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> validateSamplesheetRow(row) }

    samples = rows.map { meta, cram, crai, _jv, _jvi -> [ meta, cram, crai ] }

    normal_ch = rows
        .filter { meta, _cram, _crai, _jv, _jvi -> meta.status == 0 }
        .map { meta, _cram, _crai, _jv, _jvi -> [ meta.patient, meta ] }
        .groupTuple()
        .map { patient, metas ->
            if (metas.size() > 1) {
                error "Patient '${patient}' has more than one normal (status=0) sample: ${metas.collect { it.sample }.join(', ')}"
            }
            [ patient, metas[0] ]
        }

    // Broadcast the full set of patients that DO have a normal to every tumor row, so a
    // tumor whose patient has zero normal rows (typo'd patient column, or a genuinely
    // missing normal) is a hard error rather than silently vanishing from `pairs` --
    // `.cross()` below only emits matches, so an unmatched tumor would otherwise produce
    // no pair and no diagnostic at all. This mirrors the existing hard-error-on->1-normal
    // check, per SPEC.md Sec11.1: fail loudly on samplesheet problems, never skip silently.
    // The list is wrapped in a Map before `.combine()`: combine() flattens a bare List
    // value from the other channel into separate trailing tuple positions (confirmed
    // empirically -- with zero normals this produces zero trailing elements instead of
    // one empty list, breaking the downstream closure's arity), but a Map is passed
    // through as a single opaque element.
    normal_patients_ch = normal_ch.map { patient, _meta -> patient }.toList().map { [normals: it] }

    tumor_ch = rows
        .filter { meta, _cram, _crai, _jv, _jvi -> meta.status == 1 }
        .map { meta, _cram, _crai, jv, jvi -> [ meta.patient, meta, jv, jvi ] }
        .combine(normal_patients_ch)
        .map { patient, meta, jv, jvi, wrapper ->
            if (!(patient in wrapper.normals)) {
                error "Tumor sample '${meta.sample}' (patient '${patient}') has no matching normal (status=0) sample in the samplesheet"
            }
            [ patient, meta, jv, jvi ]
        }

    // Sarek's own pattern (workflows/sarek.nf in the vendored reference copy): cross a
    // patient's single normal against every one of that patient's tumors. One normal +
    // N tumors -> N pairs, each referencing the same normal without recomputing anything
    // about it here (the actual compute dedup is in Task 2 -- this step is pure metadata).
    pairs = normal_ch.cross(tumor_ch)
        .map { normal, tumor ->
            def normal_meta = normal[1]
            def tumor_meta  = tumor[1]
            def tumor_jv    = tumor[2]
            def tumor_jvi   = tumor[3]

            if (normal_meta.sex != tumor_meta.sex) {
                error "Sex mismatch for patient '${normal_meta.patient}': normal '${normal_meta.sample}' is '${normal_meta.sex}', tumor '${tumor_meta.sample}' is '${tumor_meta.sex}'"
            }

            def pair_meta = [
                participant_id: "${tumor_meta.sample}_vs_${normal_meta.sample}".toString(),
                patient: normal_meta.patient,
                normal_sample_id: normal_meta.sample,
                tumor_sample_id: tumor_meta.sample,
                sample_sex: normal_meta.sex,
                scale_method: 'median', fix_outliers: false, trailingN: 5, n_dim: 10,
                smooth_region: 10, outlier_SD: 4, smooth_SD: 2, alpha: 0.01, min_width: 2,
                undo_split: 'none', undo_prune: 0.05, undo_SD: 2, recovery_thresh: 20,
                use_het: true, miss_thresh: 1.0, agg_bins: false, desired_binsize: null
            ]
            [ pair_meta, file(resolveRelativeToProjectDir(tumor_jv), checkIfExists: true), file(resolveRelativeToProjectDir(tumor_jvi), checkIfExists: true) ]
        }

    emit:
    samples // [ meta(patient,sample,sex,status), cram, crai ]
    pairs   // [ pair_meta, joint_vcf, joint_vcf_index ]
}

def validateSamplesheetRow(row) {
    if (!row.patient) error "Samplesheet row missing 'patient': ${row}"
    if (!row.sample)  error "Samplesheet row missing 'sample': ${row}"
    if (!(row.sex in ['male', 'female'])) {
        error "Samplesheet row '${row.sample}': sex must be 'male' or 'female', got '${row.sex}'"
    }
    if (!(row.status in ['0', '1'])) {
        error "Samplesheet row '${row.sample}': status must be '0' or '1', got '${row.status}'"
    }
    if (!row.cram) error "Samplesheet row '${row.sample}' missing 'cram'"
    if (!row.crai) error "Samplesheet row '${row.sample}' missing 'crai'"

    def status = row.status.toInteger()
    if (status == 1 && (!row.joint_vcf || !row.joint_vcf_index)) {
        error "Tumor row '${row.sample}' (patient '${row.patient}') missing 'joint_vcf'/'joint_vcf_index' -- required on every tumor row until Phase 2 automates joint-genotyping (SPEC.md §11.2)"
    }

    def meta = [ patient: row.patient, sample: row.sample, sex: row.sex, status: status ]
    def cram = file(resolveRelativeToProjectDir(row.cram), checkIfExists: true)
    def crai = file(resolveRelativeToProjectDir(row.crai), checkIfExists: true)
    [ meta, cram, crai, row.joint_vcf, row.joint_vcf_index ]
}

// Samplesheet paths are authored relative to the pipeline's own directory (this is the
// nf-core/sarek convention: samplesheet-referenced paths are resolved relative to where
// the pipeline lives, not relative to whatever directory a test runner happens to `cd`
// into -- nf-test in particular executes each test from an isolated per-test directory
// under `.nf-test/tests/<hash>/`, not from the repo root, so resolving against Nextflow's
// ambient launch/working directory would break every fixture-driven test). `projectDir`
// is Nextflow's own implicit variable for the directory containing the pipeline's entry
// script, which is stable across both real invocations and nf-test's isolated execution.
// Absolute paths -- including URI-scheme paths like `s3://...`, `gs://...`,
// `https://...` (SPEC.md Sec11.2's canonical samplesheet uses s3:// for every file
// column) -- are passed through unchanged. Note: we must check the ORIGINAL string
// for a leading '/' or a URI scheme, not `file(p).isAbsolute()` -- Nextflow's `file()`
// always returns an absolute Path (it eagerly resolves relative strings against the
// ambient launch/working directory before we ever get to inspect them), so checking
// absoluteness on the already-constructed Path is always true and silently defeats
// this whole function.
def resolveRelativeToProjectDir(p) {
    def s = p.toString()
    def isAbsolute = s.startsWith('/') || s ==~ /^[a-zA-Z][a-zA-Z0-9+.\-]*:\/\/.*/
    isAbsolute ? file(s) : file("${projectDir}/${s}")
}
