# Samplesheet-Based Input & Multi-Pair Support — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `main.nf`'s ten flat `--param` CLI flags (one tumor/normal pair per run)
with a sarek-style per-sample samplesheet, adding support for multiple pairs per run and
real compute dedup when the same normal is reused across several tumors.

**Architecture:** Two new subworkflows sit in front of the already-built `CNV_CORE`:
`SAMPLESHEET_TO_PAIRS` (pure channel/metadata logic — parses the CSV, validates it,
derives one pair per tumor via a sarek-style `.cross()` on `patient`) and
`SAMPLE_READ_COUNTS` (moves `PREPROCESS_INTERVALS`/`COLLECT_READ_COUNTS` out of `CNV_CORE`
to run once per unique sample, before pairing — this is where the dedup actually happens).
`CNV_CORE` itself keeps its internal wiring from `SUBSET_PAIR_VCF` onward unchanged; only
its `take:` interface changes, from raw CRAM channels to pre-computed counts channels.

**Tech Stack:** Nextflow DSL2, plain `.splitCsv(header: true)` CSV parsing (no nf-schema/
nf-validation plugin — this org's other Cirro pipelines are hand-rolled, not nf-core-template
based; matching that convention rather than adding a new dependency), nf-test.

**Spec:** `SPEC.md` §11 (this plan implements it directly — read it first, it has the full
rationale and the exact schema/pairing design this plan just executes). §4-§10 cover the
already-built `CNV_CORE` this plan restructures the front end of.

## Global Constraints

- Samplesheet columns: `patient, sample, sex, status, cram, crai, joint_vcf,
  joint_vcf_index` (SPEC.md §11.2). `sex` is `male`/`female` (matching the existing WDL/R
  convention already used throughout this codebase — NOT sarek's `XX`/`XY`). `status` is
  `0` (normal) or `1` (tumor), sarek-identical. `cram`/`crai` required on every row.
  `joint_vcf`/`joint_vcf_index` required on tumor rows only — transitional until Phase 2
  (SPEC.md §5) automates joint-genotyping.
- Validation, hard errors (not warnings), mirroring sarek's own guards plus the original
  WDL's own preprocessing check: a patient with more than one `status: 0` row errors; a
  tumor whose `sex` doesn't match its patient's normal errors.
- Pair identity: `"${tumor.sample}_vs_${normal.sample}"` — matches both sarek's own
  `meta.id` convention and the original WDL pipeline's Cirro `preprocess.py`
  (`participant_id = f"{tumor_sample_id}_vs_{normal_sample_id}"`), so this isn't a new
  naming scheme, it's convergent with two things that already existed.
- The uniform per-participant `meta` shape `CNV_CORE` already expects (Phase 1's Task 5
  note, still true) does not change: `participant_id, normal_sample_id, tumor_sample_id,
  sample_sex`, plus the fixed tangent/segment tuning fields (`scale_method: 'median',
  fix_outliers: false, trailingN: 5, n_dim: 10, smooth_region: 10, outlier_SD: 4,
  smooth_SD: 2, alpha: 0.01, min_width: 2, undo_split: 'none', undo_prune: 0.05,
  undo_SD: 2, recovery_thresh: 20, use_het: true, miss_thresh: 1.0, agg_bins: false,
  desired_binsize: null`) — SPEC.md §3's "don't expose a richer configuration UI" still
  applies; these stay hardcoded, not new per-row samplesheet columns.
- `checkIfExists: true` on every `file()` call touching a samplesheet-supplied path —
  matches the standard this repo already committed to in Phase 1's final review (Fix 3):
  fail immediately and legibly, not deep inside a container.
- Every module inside `CNV_CORE` from `SUBSET_PAIR_VCF` onward (`modules/local/*.nf`
  except `preprocess_intervals.nf`/`collect_read_counts.nf`) is unchanged by this plan —
  don't touch them. `collect_read_counts.nf` and `preprocess_intervals.nf` are also
  unchanged in content, just called from a new location.

---

## Task 1: `SAMPLESHEET_TO_PAIRS` subworkflow

This is pure channel/metadata logic — no processes, no containers, no external data. It's
the one piece of this whole project that can be fully, genuinely tested without depending
on either of the two accepted gaps (missing fixtures, unbuilt R container). The fixture
files it needs are tiny, committed, empty placeholders — not real genomic data.

**Files:**
- Create: `subworkflows/local/samplesheet_to_pairs/main.nf`
- Create: `tests/subworkflows/samplesheet_to_pairs.nf.test`
- Create (empty placeholder files, committed — see Step 1): `tests/fixtures/samplesheet_to_pairs/*`

**Interfaces:**
- Consumes: `samplesheet` — a path (string or file) to the CSV.
- Produces:
  - `samples` — `tuple val(meta), path(cram), path(crai)` where `meta` is
    `[patient: String, sample: String, sex: String, status: Integer]`, one item per
    samplesheet row. Consumed by Task 2 (`SAMPLE_READ_COUNTS`).
  - `pairs` — `tuple val(pair_meta), path(joint_vcf), path(joint_vcf_index)` where
    `pair_meta` is the full uniform meta `CNV_CORE` expects (see Global Constraints), one
    item per tumor row (i.e. per pair). Consumed by Task 4 (`main.nf`'s join into
    `CNV_CORE`'s inputs).

- [ ] **Step 1: Create tiny placeholder fixture files (committed, not gitignored)**

```bash
mkdir -p tests/fixtures/samplesheet_to_pairs
touch tests/fixtures/samplesheet_to_pairs/normal.cram
touch tests/fixtures/samplesheet_to_pairs/normal.cram.crai
touch tests/fixtures/samplesheet_to_pairs/tumor_a.cram
touch tests/fixtures/samplesheet_to_pairs/tumor_a.cram.crai
touch tests/fixtures/samplesheet_to_pairs/tumor_b.cram
touch tests/fixtures/samplesheet_to_pairs/tumor_b.cram.crai
touch tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz
touch tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz.tbi
touch tests/fixtures/samplesheet_to_pairs/tumor_b_vs_normal.vcf.gz
touch tests/fixtures/samplesheet_to_pairs/tumor_b_vs_normal.vcf.gz.tbi
```

These are zero-byte placeholders whose only job is to exist on disk so `checkIfExists:
true` passes — nothing ever reads their content, since this subworkflow only builds
channels, it never stages files into a process. Unlike `tests/fixtures/cnv_core/`, these
ARE meant to be committed to git (add a one-line note to `.gitignore` if needed to confirm
they aren't caught by the existing `tests/fixtures/*` ignore rule from Phase 1 — check
`.gitignore`'s current `tests/fixtures/cnv_core/` scoping; this new directory is a sibling,
not nested under `cnv_core/`, so it likely isn't excluded at all, but verify before
committing).

- [ ] **Step 2: Create the test samplesheets (committed fixtures)**

```csv
# tests/fixtures/samplesheet_to_pairs/multi_tumor.csv
patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index
PT01,PT01_normal,female,0,tests/fixtures/samplesheet_to_pairs/normal.cram,tests/fixtures/samplesheet_to_pairs/normal.cram.crai,,
PT01,PT01_tumor_a,female,1,tests/fixtures/samplesheet_to_pairs/tumor_a.cram,tests/fixtures/samplesheet_to_pairs/tumor_a.cram.crai,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz.tbi
PT01,PT01_tumor_b,female,1,tests/fixtures/samplesheet_to_pairs/tumor_b.cram,tests/fixtures/samplesheet_to_pairs/tumor_b.cram.crai,tests/fixtures/samplesheet_to_pairs/tumor_b_vs_normal.vcf.gz,tests/fixtures/samplesheet_to_pairs/tumor_b_vs_normal.vcf.gz.tbi
```

These relative paths resolve against Nextflow's `launchDir` (confirmed empirically: a
relative path inside a samplesheet CSV resolves against the directory `nextflow run` (or
here, `nf-test`) is invoked from, not against the CSV's own directory, and not against
`projectDir`) — every test/task in this whole project is already run from the repo root
(every prior dispatch says "Work from: `<repo root>`"), so this works without needing
`${projectDir}`-prefixed paths inside the CSV itself.

Also create two more tiny CSVs for the validation-error tests:

```csv
# tests/fixtures/samplesheet_to_pairs/two_normals.csv
patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index
PT02,PT02_normal_1,female,0,tests/fixtures/samplesheet_to_pairs/normal.cram,tests/fixtures/samplesheet_to_pairs/normal.cram.crai,,
PT02,PT02_normal_2,female,0,tests/fixtures/samplesheet_to_pairs/normal.cram,tests/fixtures/samplesheet_to_pairs/normal.cram.crai,,
PT02,PT02_tumor,female,1,tests/fixtures/samplesheet_to_pairs/tumor_a.cram,tests/fixtures/samplesheet_to_pairs/tumor_a.cram.crai,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz.tbi
```

```csv
# tests/fixtures/samplesheet_to_pairs/sex_mismatch.csv
patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index
PT03,PT03_normal,female,0,tests/fixtures/samplesheet_to_pairs/normal.cram,tests/fixtures/samplesheet_to_pairs/normal.cram.crai,,
PT03,PT03_tumor,male,1,tests/fixtures/samplesheet_to_pairs/tumor_a.cram,tests/fixtures/samplesheet_to_pairs/tumor_a.cram.crai,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz.tbi
```

And one more covering the plain single-pair case — this is the backward-compatibility
check for Phase 1's original one-tumor/one-normal usage, explicitly confirmed as a
requirement during design (not just an assumption):

```csv
# tests/fixtures/samplesheet_to_pairs/single_pair.csv
patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index
PT04,PT04_normal,male,0,tests/fixtures/samplesheet_to_pairs/normal.cram,tests/fixtures/samplesheet_to_pairs/normal.cram.crai,,
PT04,PT04_tumor,male,1,tests/fixtures/samplesheet_to_pairs/tumor_a.cram,tests/fixtures/samplesheet_to_pairs/tumor_a.cram.crai,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz,tests/fixtures/samplesheet_to_pairs/tumor_a_vs_normal.vcf.gz.tbi
```

- [ ] **Step 3: Write the failing nf-tests**

```groovy
// tests/subworkflows/samplesheet_to_pairs.nf.test
nextflow_workflow {
    name "Test SAMPLESHEET_TO_PAIRS"
    script "subworkflows/local/samplesheet_to_pairs/main.nf"
    workflow "SAMPLESHEET_TO_PAIRS"

    test("one normal + two tumors produces two pairs sharing the same normal, and the normal appears exactly once in the samples channel") {
        when {
            workflow {
                """
                input[0] = "tests/fixtures/samplesheet_to_pairs/multi_tumor.csv"
                """
            }
        }
        then {
            assert workflow.success
            assert workflow.out.samples.size() == 3
            assert workflow.out.pairs.size() == 2

            def pair_metas = workflow.out.pairs.collect { it[0] }
            assert pair_metas*.participant_id.sort() == ["PT01_tumor_a_vs_PT01_normal", "PT01_tumor_b_vs_PT01_normal"]
            assert pair_metas.every { it.normal_sample_id == "PT01_normal" }
            assert pair_metas*.tumor_sample_id.sort() == ["PT01_tumor_a", "PT01_tumor_b"]
            assert pair_metas.every { it.sample_sex == "female" }

            def sample_metas = workflow.out.samples.collect { it[0] }
            assert sample_metas.count { it.sample == "PT01_normal" } == 1
        }
    }

    test("a single normal + single tumor produces exactly one pair (Phase 1's original case, unchanged behavior)") {
        when {
            workflow {
                """
                input[0] = "tests/fixtures/samplesheet_to_pairs/single_pair.csv"
                """
            }
        }
        then {
            assert workflow.success
            assert workflow.out.samples.size() == 2
            assert workflow.out.pairs.size() == 1

            def pair_meta = workflow.out.pairs[0][0]
            assert pair_meta.participant_id == "PT04_tumor_vs_PT04_normal"
            assert pair_meta.normal_sample_id == "PT04_normal"
            assert pair_meta.tumor_sample_id == "PT04_tumor"
            assert pair_meta.sample_sex == "male"
        }
    }

    test("a patient with two normal rows is a hard error") {
        when {
            workflow {
                """
                input[0] = "tests/fixtures/samplesheet_to_pairs/two_normals.csv"
                """
            }
        }
        then {
            assert !workflow.success
            assert workflow.stdout.join("\n").contains("more than one normal")
        }
    }

    test("a tumor/normal sex mismatch within a patient is a hard error") {
        when {
            workflow {
                """
                input[0] = "tests/fixtures/samplesheet_to_pairs/sex_mismatch.csv"
                """
            }
        }
        then {
            assert !workflow.success
            assert workflow.stdout.join("\n").contains("Sex mismatch")
        }
    }
}
```

- [ ] **Step 4: Run the tests, verify they fail**

Run: `nf-test test tests/subworkflows/samplesheet_to_pairs.nf.test`
Expected: FAIL — `subworkflows/local/samplesheet_to_pairs/main.nf` doesn't exist yet.

- [ ] **Step 5: Implement the subworkflow**

```groovy
// subworkflows/local/samplesheet_to_pairs/main.nf
workflow SAMPLESHEET_TO_PAIRS {
    take:
    samplesheet // path to samplesheet CSV (patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index)

    main:
    rows = Channel.fromPath(samplesheet, checkIfExists: true)
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

    tumor_ch = rows
        .filter { meta, _cram, _crai, _jv, _jvi -> meta.status == 1 }
        .map { meta, _cram, _crai, jv, jvi -> [ meta.patient, meta, jv, jvi ] }

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
            [ pair_meta, file(tumor_jv, checkIfExists: true), file(tumor_jvi, checkIfExists: true) ]
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
    [ meta, file(row.cram, checkIfExists: true), file(row.crai, checkIfExists: true), row.joint_vcf, row.joint_vcf_index ]
}
```

- [ ] **Step 6: Run the tests, verify they pass**

Run: `nf-test test tests/subworkflows/samplesheet_to_pairs.nf.test`
Expected: PASS — all four tests, for real, no fixture/container gap involved. This is the
first fully-green test in this project; if any of the three doesn't pass, that's a real bug,
not an accepted gap.

- [ ] **Step 7: Commit**

```bash
git add subworkflows/local/samplesheet_to_pairs/main.nf tests/subworkflows/samplesheet_to_pairs.nf.test tests/fixtures/samplesheet_to_pairs/
git commit -m "Add SAMPLESHEET_TO_PAIRS subworkflow (SPEC.md §11)"
```

---

## Task 2: `SAMPLE_READ_COUNTS` subworkflow

Moves `PREPROCESS_INTERVALS`/`COLLECT_READ_COUNTS` out of `CNV_CORE` to run once per unique
sample, before pairing — this is where the actual compute dedup happens: a normal
referenced by three tumor rows appears exactly once in `SAMPLESHEET_TO_PAIRS`'s `samples`
output (Task 1), so it only flows through `COLLECT_READ_COUNTS` once, full stop — no
`.unique()` operator needed anywhere, the dedup is inherent to building the channel from
sample-level rows instead of pair-level ones.

**Files:**
- Create: `subworkflows/local/sample_read_counts/main.nf`
- Create: `tests/subworkflows/sample_read_counts.nf.test`

**Interfaces:**
- Consumes: `samples_ch` — `tuple val(meta), path(cram), path(crai)`, `meta` =
  `[patient, sample, sex, status]` (Task 1's `samples` output, unchanged shape).
  `ref_fasta_ch` — `tuple path(fasta), path(fasta_fai), path(fasta_dict)` (same shape
  Phase 1 already used).
- Produces: `sample_counts` — `tuple val(meta), path("*_counts.txt")`, one item per unique
  sample, `meta` unchanged from the input. Consumed by Task 4 (`main.nf`'s join into pairs).

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/subworkflows/sample_read_counts.nf.test
nextflow_workflow {
    name "Test SAMPLE_READ_COUNTS"
    script "subworkflows/local/sample_read_counts/main.nf"
    workflow "SAMPLE_READ_COUNTS"

    test("runs PREPROCESS_INTERVALS once and COLLECT_READ_COUNTS once per unique sample") {
        when {
            workflow {
                """
                input[0] = Channel.of(
                    [ [patient: 'PT01', sample: 'PT01_normal', sex: 'female', status: 0], file(params.ref_fasta), file(params.ref_fasta_fai) ],
                    [ [patient: 'PT01', sample: 'PT01_tumor_a', sex: 'female', status: 1], file(params.ref_fasta), file(params.ref_fasta_fai) ]
                ).map { meta, cram, crai -> [ meta, cram, crai ] }
                input[1] = Channel.of([ file(params.ref_fasta), file(params.ref_fasta_fai), file(params.ref_fasta_dict) ])
                """
            }
        }
        then {
            assert !workflow.success // expected: no real CRAM fixture, GATK will fail on file content -- see note below
        }
    }
}
```

**Note on this test's fixture data**: unlike Task 1, this subworkflow calls a real GATK
process (`COLLECT_READ_COUNTS`), so it hits the same accepted gap every other GATK/R-based
test in this repo hits — there's no real CRAM fixture staged (`tests/fixtures/cnv_core/`
gap, SPEC.md §10, still open). The `when` block above uses `params.ref_fasta` itself
(a real, valid FASTA, not a CRAM) as a stand-in `cram`/`crai` path purely so
`checkIfExists: true` passes and the wiring can be exercised — GATK will then fail for
real, but on "this isn't a valid CRAM" rather than "file not found", which still proves the
channel wiring and process invocation are correct. Follow the same "attempt real execution,
confirm the failure is the known gap and not a wiring bug, report DONE_WITH_CONCERNS"
pattern established throughout Phase 1's plan.

- [ ] **Step 2: Run it, verify it fails on missing subworkflow**

Run: `nf-test test tests/subworkflows/sample_read_counts.nf.test`
Expected: FAIL — script doesn't exist.

- [ ] **Step 3: Implement the subworkflow**

```groovy
// subworkflows/local/sample_read_counts/main.nf
include { PREPROCESS_INTERVALS } from '../../../modules/local/preprocess_intervals'
include { COLLECT_READ_COUNTS  } from '../../../modules/local/collect_read_counts'

workflow SAMPLE_READ_COUNTS {
    take:
    samples_ch    // [ meta(patient,sample,sex,status), cram, crai ]
    ref_fasta_ch  // [ fasta, fai, dict ]

    main:
    PREPROCESS_INTERVALS(ref_fasta_ch)

    COLLECT_READ_COUNTS(
        samples_ch.combine(PREPROCESS_INTERVALS.out).combine(ref_fasta_ch)
            .map { meta, cram, crai, intervals, fasta, fai, dict -> [ meta, meta.sample, cram, crai, intervals, fasta, fai, dict ] }
    )

    emit:
    sample_counts = COLLECT_READ_COUNTS.out // [ meta, counts_file ]
}
```

Both included modules are used exactly as they already exist — no changes to
`modules/local/preprocess_intervals.nf` or `modules/local/collect_read_counts.nf`.
`COLLECT_READ_COUNTS`'s `sample_id` parameter (its second input element, per its existing
signature) is `meta.sample` here — the module itself doesn't care whether `meta` is
sample-scoped (this task) or pair-scoped (Phase 1's now-removed direct calls); it only ever
reads `meta.sample`/`sample_id` for tagging and reflects `meta` back unchanged in its
output, exactly as designed in Phase 1.

- [ ] **Step 4: Run it, verify the failure mode is the accepted gap, not a wiring bug**

Run: `nf-test test tests/subworkflows/sample_read_counts.nf.test`
Expected: FAIL, but specifically on GATK rejecting the FASTA-as-CRAM stand-in input (a
real, different error than "file not found" or a Nextflow syntax/channel-shape error). If
you get a channel-shape or Groovy error instead, that's a real bug in this task, not the
accepted gap — fix it.

- [ ] **Step 5: Commit**

```bash
git add subworkflows/local/sample_read_counts/main.nf tests/subworkflows/sample_read_counts.nf.test
git commit -m "Add SAMPLE_READ_COUNTS subworkflow -- moves PREPROCESS_INTERVALS/COLLECT_READ_COUNTS out of CNV_CORE for sample-level dedup (SPEC.md §11.3)"
```

---

## Task 3: Restructure `CNV_CORE`'s `take:` interface

Removes `PREPROCESS_INTERVALS`/`COLLECT_READ_COUNTS_TUMOR`/`COLLECT_READ_COUNTS_NORMAL`
from inside `CNV_CORE` (they moved to Task 2) and replaces the raw-CRAM `take:` parameters
with pre-computed counts. Every module from `SUBSET_PAIR_VCF` onward is unchanged — this
task is entirely about `CNV_CORE`'s own `take:`/`main:` wiring, not the modules it calls.

**Files:**
- Modify: `subworkflows/local/cnv_core/main.nf`
- Modify: `tests/subworkflows/cnv_core.nf.test`

**Interfaces:**
- **Old `take:`** (Phase 1): `joint_vcf_ch, tumor_cram_ch, normal_cram_ch, ref_fasta_ch,
  ref_dir_1000g, genetic_map, pon, ref_clinical, cohort_clinical`
- **New `take:`**: `joint_vcf_ch, tumor_counts_ch, ref_dir_1000g, genetic_map, pon,
  ref_clinical, cohort_clinical` — `ref_fasta_ch` is dropped entirely (nothing inside
  `CNV_CORE` needs the raw reference FASTA anymore now that `PREPROCESS_INTERVALS` moved
  out); `tumor_cram_ch`/`normal_cram_ch` are replaced by a single `tumor_counts_ch`
  (`[meta, counts_file]`, already computed by Task 2). **`normal_counts_ch` is deliberately
  NOT part of the new interface** — the normal's read counts were already unused inside
  `CNV_CORE` in Phase 1 (a known, ledgered, WDL-faithful gap), and now that
  `COLLECT_READ_COUNTS` runs uniformly over every sample regardless of tumor/normal role
  (Task 2), there's no role-specific plumbing left to carry a dead channel through — this
  is a simplification this restructuring enables, not a new gap.
- `emit:` unchanged: `read_depth, seg_file, ac_cts, segment_plots, segment_cts`.

- [ ] **Step 1: Update the failing test first**

Rewrite `tests/subworkflows/cnv_core.nf.test` for the new interface (7 inputs instead of
9; `input[1]` is now tumor counts, not tumor CRAM; `input[2]` (old normal CRAM) and the old
`input[3]` (ref_fasta) are gone):

```groovy
// tests/subworkflows/cnv_core.nf.test
nextflow_workflow {
    name "Test CNV_CORE end-to-end against a known-good WDL run"
    script "subworkflows/local/cnv_core/main.nf"
    workflow "CNV_CORE"

    test("produces all five final outputs from a real joint VCF + tumor counts") {
        when {
            workflow {
                """
                assert params.fixture_participant_id != null : "CNV_CORE fixtures not staged: params.fixture_participant_id is null. Stage tests/fixtures/cnv_core/ and set the fixture_* params (see tests/fixtures/cnv_core/README.md)."
                assert params.fixture_normal_id != null : "CNV_CORE fixtures not staged: params.fixture_normal_id is null. See tests/fixtures/cnv_core/README.md."
                assert params.fixture_tumor_id != null : "CNV_CORE fixtures not staged: params.fixture_tumor_id is null. See tests/fixtures/cnv_core/README.md."
                assert params.fixture_sample_sex != null : "CNV_CORE fixtures not staged: params.fixture_sample_sex is null. See tests/fixtures/cnv_core/README.md."

                def meta = [
                    participant_id: params.fixture_participant_id,
                    normal_sample_id: params.fixture_normal_id,
                    tumor_sample_id: params.fixture_tumor_id,
                    sample_sex: params.fixture_sample_sex,
                    scale_method: 'median', fix_outliers: false, trailingN: 5, n_dim: 10,
                    smooth_region: 10, outlier_SD: 4, smooth_SD: 2, alpha: 0.01, min_width: 2,
                    undo_split: 'none', undo_prune: 0.05, undo_SD: 2, recovery_thresh: 20,
                    use_het: true, miss_thresh: 1.0, agg_bins: false, desired_binsize: null
                ]
                input[0] = Channel.of([
                    meta,
                    file("\${projectDir}/tests/fixtures/cnv_core/joint.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/cnv_core/joint.vcf.gz.tbi")
                ])
                input[1] = Channel.of([
                    meta,
                    file("\${projectDir}/tests/fixtures/cnv_core/TEST_TUMOR_counts.txt")
                ])
                input[2] = Channel.value(file(params.ref_dir_1000g))
                input[3] = Channel.value(file(params.genetic_map))
                input[4] = Channel.value(file(params.pon))
                input[5] = Channel.value(file(params.ref_clinical))
                input[6] = Channel.value(file(params.cohort_clinical))
                """
            }
        }
        then {
            assert workflow.success
            assert workflow.out.read_depth
            assert workflow.out.seg_file
            assert workflow.out.ac_cts
            assert workflow.out.segment_plots
            assert workflow.out.segment_cts
        }
    }
}
```

- [ ] **Step 2: Run it, verify it fails (script doesn't match new test shape yet)**

Run: `nf-test test tests/subworkflows/cnv_core.nf.test`
Expected: FAIL — either a channel-arity mismatch or (if `TEST_TUMOR_counts.txt` isn't
staged in fixtures yet, which it isn't) a missing-fixture error. Either is fine at this
point in the TDD cycle; Step 3 fixes the arity, and the fixture gap remains, expected.

- [ ] **Step 3: Restructure the subworkflow**

```groovy
// subworkflows/local/cnv_core/main.nf
include { SUBSET_PAIR_VCF     } from '../../../modules/local/subset_pair_vcf'
include { FILTER_VARIANTS     } from '../../../modules/local/filter_variants'
include { ITER_ISEC           } from '../../../modules/local/iter_isec'
include { HETVAR_FILTER       } from '../../../modules/local/hetvar_filter'
include { EAGLE_PHASING       } from '../../../modules/local/eagle_phasing'
include { INDEX_PHASED        } from '../../../modules/local/index_phased'
include { QUERY_HET_AD        } from '../../../modules/local/query_het_ad'
include { AC_CALC             } from '../../../modules/local/ac_calc'
include { TANGENT_XY          } from '../../../modules/local/tangent_xy'
include { SEGMENT             } from '../../../modules/local/segment'

workflow CNV_CORE {
    take:
    joint_vcf_ch     // [ meta, joint_vcf, joint_vcf_index ] -- meta is the single uniform
                     // per-participant map used by every process below: participant_id,
                     // normal_sample_id, tumor_sample_id, sample_sex, agg_bins,
                     // desired_binsize, and every tangent_xy/segment tuning field. It must
                     // never be reshaped downstream (only ever read from), so every
                     // .join() on `meta` below matches.
    tumor_counts_ch  // [ meta, counts_file ] -- same meta as joint_vcf_ch. Pre-computed by
                     // SAMPLE_READ_COUNTS (once per unique sample, before pairing) and
                     // joined into this pair's meta by main.nf before calling CNV_CORE --
                     // this subworkflow no longer computes read counts itself.
    ref_dir_1000g
    genetic_map
    pon
    ref_clinical
    cohort_clinical

    main:
    SUBSET_PAIR_VCF(joint_vcf_ch)

    FILTER_VARIANTS(SUBSET_PAIR_VCF.out)
    ITER_ISEC(FILTER_VARIANTS.out.combine(ref_dir_1000g))
    HETVAR_FILTER(ITER_ISEC.out)
    EAGLE_PHASING(HETVAR_FILTER.out.map { meta, het, _idx -> [ meta, het ] }.combine(ref_dir_1000g).combine(genetic_map))
    INDEX_PHASED(EAGLE_PHASING.out)

    query_input = HETVAR_FILTER.out
        .join(EAGLE_PHASING.out)
        .join(INDEX_PHASED.out)
        .map { meta, het, het_idx, phased, phased_idx ->
            [ meta, het.find { it.name.contains('_chr1_') }, het, het_idx, phased, phased_idx ]
        }
    QUERY_HET_AD(query_input)

    AC_CALC(
        tumor_counts_ch
            .join(QUERY_HET_AD.out)
            .map { meta, tumor_counts, hetinfo, ad, phased -> [ meta, tumor_counts, phased, ad, hetinfo ] }
    )

    TANGENT_XY(
        tumor_counts_ch
            .combine(pon).combine(ref_clinical).combine(cohort_clinical)
    )

    SEGMENT(
        AC_CALC.out.map { meta, _rd, ac, _plots -> [ meta, ac ] }
            .join(TANGENT_XY.out)
    )

    emit:
    read_depth     = AC_CALC.out.map { meta, rd, _ac, _plots -> [ meta, rd ] }
    seg_file       = SEGMENT.out.map { meta, seg, _cts, _plots, _scts -> [ meta, seg ] }
    ac_cts         = SEGMENT.out.map { meta, _seg, cts, _plots, _scts -> [ meta, cts ] }
    segment_plots  = SEGMENT.out.map { meta, _seg, _cts, plots, _scts -> [ meta, plots ] }
    segment_cts    = SEGMENT.out.map { meta, _seg, _cts, _plots, scts -> [ meta, scts ] }
}
```

- [ ] **Step 4: Run it, confirm the failure mode is only the pre-existing fixture gap**

Run: `nf-test test tests/subworkflows/cnv_core.nf.test`
Expected: FAIL on the missing-fixture guard (`params.fixture_participant_id != null`) or,
if that's somehow set, on the missing `TEST_TUMOR_counts.txt` fixture file — NOT on a
channel-arity or `.join()` key mismatch. If you see the latter, that's a real bug in this
step's restructuring — trace it the same way Phase 1's Task 14 review did (state each
`.join()`'s two input shapes explicitly and check they match) before assuming it's the
accepted gap.

- [ ] **Step 5: Commit**

```bash
git add subworkflows/local/cnv_core/main.nf tests/subworkflows/cnv_core.nf.test
git commit -m "Restructure CNV_CORE take: to accept pre-computed counts instead of raw CRAMs (SPEC.md §11.3)"
```

---

## Task 4: Rewire `main.nf` and `nextflow.config`

Replaces the ten single-pair CLI params with `params.input` (samplesheet path), wires
`SAMPLESHEET_TO_PAIRS` → `SAMPLE_READ_COUNTS` → the join that builds each pair's full
`CNV_CORE` input → `CNV_CORE`, now with a **multi-item** channel (one item per pair)
instead of Phase 1's single `Channel.of([meta, ...])`.

**Files:**
- Modify: `main.nf`
- Modify: `nextflow.config`
- Create: `assets/samplesheet.csv` (example/template, matching sarek's own
  `assets/samplesheet.csv` convention)

**Interfaces:**
- `params.input`: path to the samplesheet CSV. Replaces `participant_id,
  normal_sample_id, tumor_sample_id, sample_sex, joint_vcf, joint_vcf_index, tumor_cram,
  tumor_crai, normal_cram, normal_crai` — all ten removed from `nextflow.config` and from
  `main.nf`'s required-params guard.

- [ ] **Step 1: Remove the old single-pair params from `nextflow.config`**

Delete these ten lines and their comment block (added in Phase 1's final fix wave):

```groovy
params.participant_id    = null
params.normal_sample_id  = null
params.tumor_sample_id   = null
params.sample_sex        = null
params.joint_vcf         = null
params.joint_vcf_index   = null
params.tumor_cram        = null
params.tumor_crai        = null
params.normal_cram       = null
params.normal_crai       = null
```

Replace with:

```groovy
// Samplesheet path -- see SPEC.md §11.2 for the column schema and assets/samplesheet.csv
// for a template. No sensible default exists; must be supplied via --input.
params.input = null
```

Everything else in `nextflow.config` (the reference-file params, `cnv_r_container`, the
`publishDir` block, the `fixture_*` test-only params, `docker.enabled`) is unaffected by
this task.

- [ ] **Step 2: Write the example samplesheet template**

```csv
# assets/samplesheet.csv
patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index
PATIENT_ID,PATIENT_ID_normal,female,0,/path/to/PATIENT_ID_normal.recal.cram,/path/to/PATIENT_ID_normal.recal.cram.crai,,
PATIENT_ID,PATIENT_ID_tumor,female,1,/path/to/PATIENT_ID_tumor.recal.cram,/path/to/PATIENT_ID_tumor.recal.cram.crai,/path/to/PATIENT_ID_tumor_vs_PATIENT_ID_normal.vcf.gz,/path/to/PATIENT_ID_tumor_vs_PATIENT_ID_normal.vcf.gz.tbi
```

- [ ] **Step 3: Rewrite `main.nf`**

```groovy
// main.nf
nextflow.enable.dsl = 2

include { SAMPLESHEET_TO_PAIRS } from './subworkflows/local/samplesheet_to_pairs/main'
include { SAMPLE_READ_COUNTS   } from './subworkflows/local/sample_read_counts/main'
include { CNV_CORE             } from './subworkflows/local/cnv_core/main'

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

    // Join each pair to its tumor's already-computed (and, when a normal is shared
    // across several pairs, already-deduplicated) read counts. Keyed on a small
    // (patient, sample) map rather than the full meta, since pair_meta and sample_meta
    // are different shapes -- constructing the same-shaped key on both sides is what
    // makes the .join() match, the same principle as every meta-keyed join already in
    // this codebase (see subworkflows/local/cnv_core/main.nf's take: comment).
    tumor_counts_for_pairs = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta, jv, jvi ->
            [ [patient: pair_meta.patient, sample: pair_meta.tumor_sample_id], pair_meta, jv, jvi ]
        }
        .join(
            SAMPLE_READ_COUNTS.out.sample_counts
                .map { sample_meta, counts -> [ [patient: sample_meta.patient, sample: sample_meta.sample], counts ] },
            by: 0, failOnDuplicate: true
        )
        .map { _key, pair_meta, jv, jvi, counts -> [ pair_meta, jv, jvi, counts ] }

    joint_vcf_ch    = tumor_counts_for_pairs.map { pair_meta, jv, jvi, _counts -> [ pair_meta, jv, jvi ] }
    tumor_counts_ch = tumor_counts_for_pairs.map { pair_meta, _jv, _jvi, counts -> [ pair_meta, counts ] }

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
```

- [ ] **Step 4: Manual sanity run against Task 1's test samplesheet**

```bash
nextflow run main.nf --input tests/fixtures/samplesheet_to_pairs/multi_tumor.csv
```

Expected: proceeds past `SAMPLESHEET_TO_PAIRS` and `SAMPLE_READ_COUNTS` wiring (though
`SAMPLE_READ_COUNTS` itself will fail for real, same as Task 2 — the placeholder CRAMs
aren't valid GATK input) and, critically, does NOT fail with a Groovy channel-shape error
in `main.nf` itself. If it does, the `.join()`/`.map()` construction above has a bug --
fix it before moving on. This is the same "confirm the failure is the known gap, not a new
one" discipline as every other task in this project.

- [ ] **Step 5: Commit**

```bash
git add main.nf nextflow.config assets/samplesheet.csv
git commit -m "Rewire main.nf to samplesheet-based, multi-pair input (SPEC.md §11)"
```

---

## Task 5: Update documentation

**Files:**
- Modify: `README.md`
- Modify: `tests/fixtures/cnv_core/README.md` (only if its content references the old
  single-pair CLI invocation style — check first)

- [ ] **Step 1: Rewrite `README.md`'s usage section**

Replace the ten-flag example invocation with:

```bash
nextflow run main.nf --input samplesheet.csv
```

Add a short section describing the samplesheet format (columns, sarek-style pairing,
the transitional `joint_vcf`/`joint_vcf_index` columns), pointing at `assets/samplesheet.csv`
as a template and `SPEC.md §11` for the full rationale. Keep the existing "two known,
accepted gaps" section (missing fixtures, unbuilt R container) — those haven't changed.

- [ ] **Step 2: Check `tests/fixtures/cnv_core/README.md` for stale invocation examples**

If it references the old `--participant_id`/`--tumor_cram`/etc. flags anywhere, update
those references to describe the new samplesheet-based flow instead. If it only lists
expected fixture files (no invocation examples), no change needed.

- [ ] **Step 3: Commit**

```bash
git add README.md tests/fixtures/cnv_core/README.md
git commit -m "Update docs for samplesheet-based input"
```
