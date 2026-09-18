# Parabricks Calling + Joint Genotyping (Phase 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute the joint-genotyped, VQSR-filtered VCF the pipeline currently takes as
samplesheet input (`joint_vcf`/`joint_vcf_index`) instead, via Parabricks GPU calling per
sample and a GATK4 joint-genotyping chain per pair.

**Architecture:** Two new subworkflows mirroring the existing sample-scoped/pair-scoped
split (`SAMPLE_READ_COUNTS`/`CNV_CORE`): `SAMPLE_GERMLINE_CALLING` (Parabricks
HaplotypeCaller, once per unique sample) and `JOINT_GENOTYPING` (GenomicsDBImport →
GenotypeGVCFs → ExcessHet filter → VQSR, once per pair, no interval scatter). Four GATK4
tool wrappers are vendored verbatim from the sarek reference checkout in this repo
(`reference/sarek/modules/nf-core/gatk4/*`); the `ExcessHet` filter and sites-only step
have no sarek equivalent and are hand-written local modules, matching this repo's own
Phase 1 convention for GATK4 tasks sarek doesn't wrap.

**Tech Stack:** Nextflow DSL2, GATK4 4.2.6.1 (pinned container digest already used
elsewhere in this repo), Parabricks 4.7.1-1 (NGC container, GPU-only — cannot execute on
this dev machine, confirmed no NVIDIA GPU), bcftools, nf-test.

**Spec:** `SPEC.md` §12 (commits `4f484ad`, `cf59a0b` — the second commit fixed the VQSR
chain order and resolved the tranche-values gap; read the current file, not an earlier
version).

## Global Constraints

- **No interval scatter.** Every GATK4 call in this plan processes the whole
  `params.calling_intervals` region in one shot — no `SplitIntervals`, no
  `meta.num_intervals`, no per-shard merge. (SPEC §12.6)
- **GATK version: 4.2.6.1**, pinned via the exact container digest already used in this
  repo for `PREPROCESS_INTERVALS`/`COLLECT_READ_COUNTS`:
  `broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe`
  (SPEC §12.6). Every vendored GATK4 module's container gets overridden to this digest via
  `nextflow.config`; every hand-written local GATK4 module declares it directly in the
  process definition, matching the existing `modules/local/preprocess_intervals.nf` style.
- **Parabricks version: `nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1`** (sarek's current
  default, vendored module's own container — no override). (SPEC §12.6)
- **bcftools container** (already used 4× in this repo, reuse verbatim):
  `quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a`
- **`ExcessHet` threshold: `54.69`**, exact string `ExcessHet > 54.69`, filter name
  `ExcessHet` — a soft `VariantFiltration` FILTER tag, not a record-dropping hard filter.
  (SPEC §5.1, §12.3)
- **VQSR chain order: INDEL first, then SNP chained on INDEL's output** — verified
  against `ICONICC_paired.3.wdl`'s `ApplyRecalibration` task, the reverse of what sarek's
  own subworkflow does. Do not reorder this. (SPEC §12.3)
- **VQSR filter levels**: INDEL `99.0`, SNP `99.7`. **Tranches, annotations, resource
  labels**: exact literal lists in SPEC §12.3 — copy verbatim, do not substitute GATK's
  tool defaults or sarek's own igenomes.config values. `use_allele_specific_annotations`
  stays `false` — never pass `--use-allele-specific-annotations`.
- **Samplesheet schema**: `joint_vcf`/`joint_vcf_index` columns are removed entirely
  (SPEC §12.4) — this is a breaking change to the CSV format from Phase 1.5.
- **Dedup shape**: calling is deduplicated per unique `(patient, sample)` (matches
  `SAMPLE_READ_COUNTS`'s existing mechanism — one row per unique sample in
  `SAMPLESHEET_TO_PAIRS.out.samples`); joint genotyping stays pairwise per patient, never
  cohort-batched, even when a normal is shared across several tumors. (SPEC §12.1)
- **Testing**: no real fixture data, write code only — same posture as Phase 1/1.5.
  `SAMPLE_GERMLINE_CALLING` tests run via nf-test's `-stub` flag (confirmed against this
  environment: `nextflow help run` lists `-stub-run`/`-stub` as the flag; sarek's own
  vendored tests at `reference/sarek/tests/lib/UTILS.groovy:147` confirm nf-test's
  `options "-stub"` directive is the correct mechanism) and are genuinely green today,
  since stub mode doesn't need valid file content. `JOINT_GENOTYPING` and its component
  GATK4 modules have no real VQSR-scale fixture (a handful of variants cannot fit VQSR's
  Gaussian mixture model regardless of container/network access — this is not a "real
  fixture not staged yet" gap, it is a hard requirement of the algorithm), so those tests
  assert the desired end state (`workflow.success`/`process.success`) with placeholder
  fixtures, documented as a known, accepted gap — exactly the convention already
  established by `tests/subworkflows/cnv_core.nf.test` and
  `tests/subworkflows/sample_read_counts.nf.test`. Never flip these to `!workflow.success`.

---

### Task 1: Vendor `PARABRICKS_HAPLOTYPECALLER` + hand-write `INDEX_GVCF`

**Files:**
- Create: `modules/nf-core/parabricks/haplotypecaller/main.nf`
- Create: `modules/local/index_gvcf.nf`
- Test: `tests/nf-core/parabricks/haplotypecaller.nf.test`
- Test: `tests/local/index_gvcf.nf.test`
- Create: `tests/fixtures/joint_genotyping/placeholder.g.vcf.gz` (0-byte placeholder, same
  convention as `tests/fixtures/samplesheet_to_pairs/*.vcf.gz`)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `PARABRICKS_HAPLOTYPECALLER` process, `INDEX_GVCF` process — `INDEX_GVCF`
  takes `tuple val(meta), path(gvcf)` and emits `tuple val(meta), path(gvcf),
  path("${gvcf}.tbi")`. Task 2 chains these two directly.

- [ ] **Step 1: Vendor the Parabricks module verbatim**

Copy the file byte-for-byet from the vendored sarek reference checkout already in this
repo — do not hand-retype it:

```bash
mkdir -p modules/nf-core/parabricks/haplotypecaller
cp reference/sarek/modules/nf-core/parabricks/haplotypecaller/main.nf \
   modules/nf-core/parabricks/haplotypecaller/main.nf
```

Do not edit this file. Its `conda "${moduleDir}/environment.yml"` line references a file
this repo does not have — that's fine, Nextflow only evaluates the `conda` directive under
`-profile conda`, which this repo never uses (`nextflow.config` only enables `docker`).

- [ ] **Step 2: Write `INDEX_GVCF`**

```groovy
// modules/local/index_gvcf.nf
process INDEX_GVCF {
    tag "${meta.sample}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 1
    memory '2 GB'

    input:
    tuple val(meta), path(gvcf)

    output:
    tuple val(meta), path(gvcf), path("${gvcf}.tbi")

    script:
    """
    set -euo pipefail
    bcftools index --tbi ${gvcf}
    """
}
```

- [ ] **Step 3: Create the placeholder GVCF fixture**

```bash
mkdir -p tests/fixtures/joint_genotyping
touch tests/fixtures/joint_genotyping/placeholder.g.vcf.gz
```

This is a 0-byte stand-in, not real genomic data — the same convention already used for
`tests/fixtures/samplesheet_to_pairs/*.cram`/`*.vcf.gz`.

- [ ] **Step 4: Write the Parabricks module test (stub-run, green today)**

```groovy
// tests/nf-core/parabricks/haplotypecaller.nf.test
nextflow_process {
    name "Test PARABRICKS_HAPLOTYPECALLER (stub)"
    script "modules/nf-core/parabricks/haplotypecaller/main.nf"
    process "PARABRICKS_HAPLOTYPECALLER"
    options "-stub"

    test("emits a GVCF via the module's stub block -- this dev environment has no NVIDIA GPU, confirmed, so real pbrun execution is never attempted here") {
        when {
            params {
                modules {
                    'PARABRICKS_HAPLOTYPECALLER' {
                        ext.args = '--gvcf'
                    }
                }
            }
            process {
                """
                input[0] = [
                    [ id: 'TEST_SAMPLE' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    [],
                    []
                ]
                input[1] = [ [ id: 'reference' ], file(params.ref_fasta) ]
                """
            }
        }
        then {
            assert process.success
            assert process.out.gvcf.size() == 1
        }
    }
}
```

If `params { modules { ... } }` is not how this nf-test version injects `ext.args` for a
bare process test (nf-test's config-injection syntax varies by version and this repo has
no existing process-level `ext.args` test to copy from), set it more simply with a process
`ext` block instead:

```groovy
        when {
            process {
                """
                input[0] = [
                    [ id: 'TEST_SAMPLE' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    [],
                    []
                ]
                input[1] = [ [ id: 'reference' ], file(params.ref_fasta) ]
                """
            }
        }
```

and instead set the arg globally for this test file only via nf-test's
`config "path/to/a/scoped/nextflow.config"` block if the inline approach doesn't parse.
Verify empirically which form this installed nf-test version (`0.9.5`, confirmed via
`nf-test version`) accepts by running the test and reading the actual error before
guessing further — don't assume either form compiles untested.

- [ ] **Step 5: Write the `INDEX_GVCF` module test**

```groovy
// tests/local/index_gvcf.nf.test
nextflow_process {
    name "Test index_gvcf"
    script "modules/local/index_gvcf.nf"
    process "INDEX_GVCF"

    test("indexes a GVCF with bcftools") {
        when {
            process {
                """
                input[0] = [
                    [ patient: 'PT01', sample: 'PT01_normal' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz")
                ]
                """
            }
        }
        then {
            // Known, accepted gap (SPEC.md §12.7, Global Constraints above): the fixture
            // is a 0-byte placeholder, not a real bgzipped GVCF, so bcftools fails with a
            // named "not a BGZF file"-class error before real content ever matters. Real
            // fixture data closes this the same way it closes tests/subworkflows/cnv_core.nf.test.
            // Do not flip this to !process.success.
            assert process.success
        }
    }
}
```

- [ ] **Step 6: Run both tests, confirm expected outcomes**

Run: `nf-test test tests/nf-core/parabricks/haplotypecaller.nf.test tests/local/index_gvcf.nf.test`
Expected: the Parabricks stub test passes. The `INDEX_GVCF` test fails today on the
documented placeholder-fixture gap — confirm the failure message matches what Step 5's
comment describes (a bcftools/BGZF content error), not a wiring or syntax error. If the
failure reason is something else (missing container, syntax error, wrong input arity),
that is a real bug to fix, not the accepted gap.

- [ ] **Step 7: Commit**

```bash
git add modules/nf-core/parabricks/haplotypecaller/main.nf modules/local/index_gvcf.nf \
        tests/nf-core/parabricks/haplotypecaller.nf.test tests/local/index_gvcf.nf.test \
        tests/fixtures/joint_genotyping/placeholder.g.vcf.gz
git commit -m "Vendor PARABRICKS_HAPLOTYPECALLER and add INDEX_GVCF"
```

---

### Task 2: `SAMPLE_GERMLINE_CALLING` subworkflow

**Files:**
- Create: `subworkflows/local/sample_germline_calling/main.nf`
- Test: `tests/subworkflows/sample_germline_calling.nf.test`

**Interfaces:**
- Consumes: `PARABRICKS_HAPLOTYPECALLER`, `INDEX_GVCF` (Task 1).
- Produces: `SAMPLE_GERMLINE_CALLING` workflow, `take: samples_ch` (shape
  `[meta(patient,sample,sex,status), cram, crai]` — identical to `SAMPLE_READ_COUNTS`'s
  first `take:`, so `main.nf` feeds both subworkflows from the same channel), `take:
  ref_fasta_ch` (`[fasta]`, single file — Parabricks needs only the FASTA, no fai/dict,
  confirmed by reading `PARABRICKS_HAPLOTYPECALLER`'s own input declaration), `take:
  calling_intervals_ch` (`[intervals]`, value channel). `emit: sample_gvcfs` — `[meta, gvcf,
  gvcf_tbi]`, one row per unique sample. Task 6 joins this by `[patient, sample]` the same
  way `main.nf` already joins `tumor_counts_ch`.

- [ ] **Step 1: Write the subworkflow**

```groovy
// subworkflows/local/sample_germline_calling/main.nf
include { PARABRICKS_HAPLOTYPECALLER } from '../../../modules/nf-core/parabricks/haplotypecaller/main'
include { INDEX_GVCF                 } from '../../../modules/local/index_gvcf'

workflow SAMPLE_GERMLINE_CALLING {
    take:
    samples_ch           // [ meta(patient,sample,sex,status), cram, crai ]
    ref_fasta_ch         // [ fasta ] -- single file, Parabricks needs no fai/dict
    calling_intervals_ch // [ intervals ] -- value channel, one Picard-style .intervals file

    main:
    // Dedup is free here, the same way it is in SAMPLE_READ_COUNTS: samples_ch already
    // has exactly one row per unique (patient, sample) -- a normal shared across several
    // tumor rows was already collapsed to one row by SAMPLESHEET_TO_PAIRS before this
    // subworkflow ever sees it, so it is called exactly once regardless of how many
    // tumors reference it (SPEC.md Sec12.1).
    calling_input = samples_ch
        .combine(calling_intervals_ch)
        .map { meta, cram, crai, intervals -> [ meta + [id: meta.sample], cram, crai, intervals ] }

    PARABRICKS_HAPLOTYPECALLER(
        calling_input,
        ref_fasta_ch.map { fasta -> [ [id: 'reference'], fasta ] }
    )

    INDEX_GVCF(PARABRICKS_HAPLOTYPECALLER.out.gvcf)

    emit:
    sample_gvcfs = INDEX_GVCF.out // [ meta, gvcf, gvcf_tbi ]
}
```

`meta + [id: meta.sample]` adds the `id` field `PARABRICKS_HAPLOTYPECALLER`'s `tag`/
`prefix` logic requires (`tag "${meta.id}"`, `prefix = task.ext.prefix ?: "${meta.id}"`)
without overwriting any of `SAMPLESHEET_TO_PAIRS`'s existing fields — the merged map keeps
`patient`/`sample`/`sex`/`status` and adds `id`, matching the pattern `main.nf` already
uses when building `tumor_counts_ch`'s join key.

`params.calling_intervals` needs `--gvcf` set as `PARABRICKS_HAPLOTYPECALLER`'s `ext.args`
for this to emit a GVCF instead of a plain VCF — that override lives in `nextflow.config`,
added in Task 6, not here; this subworkflow only wires channels.

- [ ] **Step 2: Write the subworkflow test (stub-run, green today)**

```groovy
// tests/subworkflows/sample_germline_calling.nf.test
nextflow_workflow {
    name "Test SAMPLE_GERMLINE_CALLING"
    script "subworkflows/local/sample_germline_calling/main.nf"
    workflow "SAMPLE_GERMLINE_CALLING"
    options "-stub"

    test("calls every unique sample exactly once, including a normal shared across pairs") {
        when {
            workflow {
                """
                input[0] = Channel.of(
                    [ [patient: 'PT01', sample: 'PT01_normal', sex: 'female', status: 0], file(params.ref_fasta), file(params.ref_fasta_fai) ],
                    [ [patient: 'PT01', sample: 'PT01_tumor_a', sex: 'female', status: 1], file(params.ref_fasta), file(params.ref_fasta_fai) ]
                )
                input[1] = Channel.value(file(params.ref_fasta))
                input[2] = Channel.value(file(params.ref_fasta))
                """
            }
        }
        then {
            assert workflow.success
            assert workflow.out.sample_gvcfs.size() == 2
        }
    }
}
```

`input[2]` (the calling-intervals value) reuses `params.ref_fasta` purely as a
stand-in `path` value to exercise the channel wiring under `-stub` — stub mode never
reads file content, so this is safe and matches `sample_read_counts.nf.test`'s own
precedent of reusing the reference FASTA as a placeholder for an unrelated input type.

- [ ] **Step 3: Run the test**

Run: `nf-test test tests/subworkflows/sample_germline_calling.nf.test`
Expected: PASS — `-stub` mode touches placeholder output files regardless of input
content, so this test is genuinely green, unlike the joint-genotyping-side tests.

- [ ] **Step 4: Commit**

```bash
git add subworkflows/local/sample_germline_calling/main.nf tests/subworkflows/sample_germline_calling.nf.test
git commit -m "Add SAMPLE_GERMLINE_CALLING subworkflow"
```

---

### Task 3: Vendor the four joint-genotyping GATK4 modules

**Files:**
- Create: `modules/nf-core/gatk4/genomicsdbimport/main.nf`
- Create: `modules/nf-core/gatk4/genotypegvcfs/main.nf`
- Create: `modules/nf-core/gatk4/variantrecalibrator/main.nf`
- Create: `modules/nf-core/gatk4/applyvqsr/main.nf`
- Modify: `nextflow.config` (container overrides only — new `params.*` come in Task 6)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: the four GATK4 processes, containers pinned to 4.2.6.1. Task 5 composes them.

- [ ] **Step 1: Vendor all four modules verbatim**

```bash
mkdir -p modules/nf-core/gatk4/genomicsdbimport modules/nf-core/gatk4/genotypegvcfs \
         modules/nf-core/gatk4/variantrecalibrator modules/nf-core/gatk4/applyvqsr

cp reference/sarek/modules/nf-core/gatk4/genomicsdbimport/main.nf    modules/nf-core/gatk4/genomicsdbimport/main.nf
cp reference/sarek/modules/nf-core/gatk4/genotypegvcfs/main.nf       modules/nf-core/gatk4/genotypegvcfs/main.nf
cp reference/sarek/modules/nf-core/gatk4/variantrecalibrator/main.nf modules/nf-core/gatk4/variantrecalibrator/main.nf
cp reference/sarek/modules/nf-core/gatk4/applyvqsr/main.nf           modules/nf-core/gatk4/applyvqsr/main.nf
```

Do not edit any of these four files. Each declares `conda "${moduleDir}/environment.yml"`
and a `topic: versions` output emit — both are harmless here: `conda` is never evaluated
(this repo only runs `docker`), and Nextflow 26 (confirmed installed:
`nextflow help run` reports `26.04.6`) fully supports `topic:` process outputs.

- [ ] **Step 2: Override all four containers to GATK 4.2.6.1**

Add to `nextflow.config`, inside the existing `process { ... }` block (after the existing
`withName:` overrides):

```groovy
    withName: 'GATK4_GENOMICSDBIMPORT' {
        container = 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
        ext.args  = '--genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader'
    }
    withName: 'GATK4_GENOTYPEGVCFS' {
        container = 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    }
    withName: 'VARIANTRECALIBRATOR_SNP' {
        container = 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
        ext.args  = '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP --max-gaussians 6 --trust-all-polymorphic'
    }
    withName: 'VARIANTRECALIBRATOR_INDEL' {
        container = 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
        ext.args  = '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL --max-gaussians 4 --trust-all-polymorphic'
    }
    withName: 'APPLYVQSR_INDEL' {
        container = 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
        ext.args  = '--truth-sensitivity-filter-level 99.0 -mode INDEL'
    }
    withName: 'APPLYVQSR_SNP' {
        container = 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
        ext.args  = '--truth-sensitivity-filter-level 99.7 -mode SNP'
    }
    withName: 'PARABRICKS_HAPLOTYPECALLER' {
        ext.args = '--gvcf'
    }
```

`VARIANTRECALIBRATOR_SNP`/`_INDEL` and `APPLYVQSR_SNP`/`_INDEL` are the process aliases
Task 5's subworkflow will use (`include { GATK4_VARIANTRECALIBRATOR as
VARIANTRECALIBRATOR_SNP }`, etc.) — `withName:` matches the alias, not the underlying
module name, so these selectors only take effect once Task 5 wires the aliases. Do not
add a bare `withName: 'GATK4_VARIANTRECALIBRATOR'` selector — it would apply the SNP
`ext.args` to both aliases indiscriminately.

- [ ] **Step 3: Write a module-level smoke test per vendored module**

These follow the same accepted-gap pattern as Task 1's `INDEX_GVCF` test — a placeholder
fixture produces a named, documented failure today. One example (the pattern repeats for
the other three):

```groovy
// tests/nf-core/gatk4/genomicsdbimport.nf.test
nextflow_process {
    name "Test GATK4_GENOMICSDBIMPORT"
    script "modules/nf-core/gatk4/genomicsdbimport/main.nf"
    process "GATK4_GENOMICSDBIMPORT"

    test("imports a normal+tumor GVCF pair into a GenomicsDB workspace") {
        when {
            process {
                """
                input[0] = [
                    [ id: 'TEST_PAIR' ],
                    [ file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"), file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz") ],
                    [ file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz.tbi"), file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz.tbi") ],
                    file(params.calling_intervals),
                    [],
                    []
                ]
                input[1] = false
                input[2] = false
                input[3] = false
                """
            }
        }
        then {
            // Known, accepted gap (Global Constraints above, SPEC.md Sec12.7): the GVCF
            // fixtures are 0-byte placeholders, so GenomicsDBImport fails on invalid VCF
            // content before any real genomic data would matter. Do not flip to !process.success.
            assert process.success
        }
    }
}
```

`params.calling_intervals` does not exist yet (it's added in Task 6) — add a placeholder
`params.calling_intervals = null` to `nextflow.config` in this task (Task 6 replaces it
with the real S3 path) so the test below can reference the param name now; leaving it
unset would make the test fail on a missing-param error instead of the documented
GVCF-content gap.

```groovy
// tests/nf-core/gatk4/genotypegvcfs.nf.test
nextflow_process {
    name "Test GATK4_GENOTYPEGVCFS"
    script "modules/nf-core/gatk4/genotypegvcfs/main.nf"
    process "GATK4_GENOTYPEGVCFS"

    test("genotypes a GenomicsDB workspace against calling_intervals") {
        when {
            process {
                """
                input[0] = [
                    [ id: 'TEST_PAIR' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    [],
                    file(params.calling_intervals),
                    []
                ]
                input[1] = [ [ id: 'reference' ], file(params.ref_fasta) ]
                input[2] = [ [ id: 'reference' ], file(params.ref_fasta_fai) ]
                input[3] = [ [ id: 'reference' ], file(params.ref_fasta_dict) ]
                input[4] = [ [ id: 'reference' ], file(params.dbsnp_vcf) ]
                input[5] = [ [ id: 'reference' ], file(params.dbsnp_vcf_index) ]
                """
            }
        }
        then {
            // Known, accepted gap (Global Constraints above): placeholder is not a real
            // GenomicsDB workspace. Do not flip to !process.success.
            assert process.success
        }
    }
}
```

`params.dbsnp_vcf`/`params.dbsnp_vcf_index`/`params.calling_intervals` don't exist until
Task 6 either — add all three as `null` placeholders in this task's `nextflow.config`
edit, alongside the `calling_intervals` placeholder already noted above.

```groovy
// tests/nf-core/gatk4/variantrecalibrator.nf.test
nextflow_process {
    name "Test GATK4_VARIANTRECALIBRATOR"
    script "modules/nf-core/gatk4/variantrecalibrator/main.nf"
    process "GATK4_VARIANTRECALIBRATOR"

    test("fits a recalibration model against SNP-shaped resources") {
        when {
            process {
                """
                input[0] = [
                    [ id: 'TEST_PAIR' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz")
                ]
                input[1] = [ file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta) ]
                input[2] = [ file(params.ref_fasta_fai), file(params.ref_fasta_fai), file(params.ref_fasta_fai), file(params.ref_fasta_fai) ]
                input[3] = [
                    '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 ref_fasta',
                    '--resource:omni,known=false,training=true,truth=true,prior=12.0 ref_fasta',
                    '--resource:1000G,known=false,training=true,truth=false,prior=10.0 ref_fasta',
                    '--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 ref_fasta'
                ]
                input[4] = file(params.ref_fasta)
                input[5] = file(params.ref_fasta_fai)
                input[6] = file(params.ref_fasta_dict)
                """
            }
        }
        then {
            // Known, accepted gap (Global Constraints above): VQSR's Gaussian mixture
            // model cannot fit on a placeholder input regardless of container/network
            // access -- this is a hard algorithmic requirement, not a missing fixture
            // that staging real data alone would fix. Do not flip to !process.success.
            assert process.success
        }
    }
}
```

```groovy
// tests/nf-core/gatk4/applyvqsr.nf.test
nextflow_process {
    name "Test GATK4_APPLYVQSR"
    script "modules/nf-core/gatk4/applyvqsr/main.nf"
    process "GATK4_APPLYVQSR"

    test("applies a recalibration model to a VCF") {
        when {
            process {
                """
                input[0] = [
                    [ id: 'TEST_PAIR' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz")
                ]
                input[1] = file(params.ref_fasta)
                input[2] = file(params.ref_fasta_fai)
                input[3] = file(params.ref_fasta_dict)
                """
            }
        }
        then {
            // Known, accepted gap (Global Constraints above): placeholder recal/tranches
            // files, not a real VariantRecalibrator model. Do not flip to !process.success.
            assert process.success
        }
    }
}
```

- [ ] **Step 4: Run all four module tests**

Run: `nf-test test tests/nf-core/gatk4/`
Expected: all four fail today on the documented placeholder-fixture gap. Read each
failure message and confirm it names an invalid-VCF-content or similar GATK error, not a
missing container, syntax error, or wrong input arity — those would be real bugs.

- [ ] **Step 5: Commit**

```bash
git add modules/nf-core/gatk4/ nextflow.config tests/nf-core/gatk4/
git commit -m "Vendor the joint-genotyping GATK4 modules, pinned to 4.2.6.1"
```

---

### Task 4: Hand-write `VARIANT_FILTRATION` and `MAKE_SITES_ONLY_VCF`

**Files:**
- Create: `modules/local/variant_filtration.nf`
- Create: `modules/local/make_sites_only_vcf.nf`
- Test: `tests/local/variant_filtration.nf.test`
- Test: `tests/local/make_sites_only_vcf.nf.test`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `VARIANT_FILTRATION` (`tuple val(meta), path(vcf), path(tbi)` →
  `tuple val(meta), path("*.variant_filtered.vcf.gz"), path("*.variant_filtered.vcf.gz.tbi")`),
  `MAKE_SITES_ONLY_VCF` (same input shape → `tuple val(meta), path("*.sites_only.vcf.gz"),
  path("*.sites_only.vcf.gz.tbi")`). Task 5 chains `GenotypeGVCFs → VARIANT_FILTRATION →
  MAKE_SITES_ONLY_VCF → VariantRecalibrator`, and separately feeds `VARIANT_FILTRATION`'s
  full (not sites-only) output into `ApplyVQSR` — per SPEC §12.3, the sites-only VCF is
  training input only, never the final output.

Neither of these has a sarek module to vendor (confirmed: no `variantfiltration` or
`makesitesonlyvcf`/`sitesonly` directory anywhere under `reference/sarek/modules/`) —
hand-write them following this repo's own Phase 1 convention for GATK4 tasks sarek
doesn't wrap (`modules/local/preprocess_intervals.nf`), reading the exact CLI invocation
directly from `HardFilterAndMakeSitesOnlyVcf` (`JointGeno_Tasks.1.wdl:263-306`).

- [ ] **Step 1: Write `VARIANT_FILTRATION`**

```groovy
// modules/local/variant_filtration.nf
process VARIANT_FILTRATION {
    tag "${meta.participant_id}"
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("${meta.participant_id}.variant_filtered.vcf.gz"), path("${meta.participant_id}.variant_filtered.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xms3000m -Xmx3250m" \\
        VariantFiltration \\
        --filter-expression "ExcessHet > 54.69" \\
        --filter-name ExcessHet \\
        -O ${meta.participant_id}.variant_filtered.vcf.gz \\
        -V ${vcf}
    """
}
```

GATK4's VCF-writing tools default `--create-output-variant-index` to `true`, so the
`.tbi` sidecar is produced automatically — no separate indexing step needed, matching
`PreprocessIntervals`/`CollectReadCounts`'s existing precedent of never hand-indexing a
GATK4 tool's own VCF output.

- [ ] **Step 2: Write `MAKE_SITES_ONLY_VCF`**

```groovy
// modules/local/make_sites_only_vcf.nf
process MAKE_SITES_ONLY_VCF {
    tag "${meta.participant_id}"
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 1
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)

    output:
    tuple val(meta), path("${meta.participant_id}.sites_only.vcf.gz"), path("${meta.participant_id}.sites_only.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xms3000m -Xmx3250m" \\
        MakeSitesOnlyVcf \\
        -I ${vcf} \\
        -O ${meta.participant_id}.sites_only.vcf.gz
    """
}
```

- [ ] **Step 3: Write both module tests**

```groovy
// tests/local/variant_filtration.nf.test
nextflow_process {
    name "Test variant_filtration"
    script "modules/local/variant_filtration.nf"
    process "VARIANT_FILTRATION"

    test("tags ExcessHet > 54.69") {
        when {
            process {
                """
                input[0] = [
                    [ participant_id: 'TEST_TUMOR_vs_TEST_NORMAL' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz")
                ]
                """
            }
        }
        then {
            // Known, accepted gap: placeholder fixture, not real VCF content. Do not flip to !process.success.
            assert process.success
        }
    }
}
```

```groovy
// tests/local/make_sites_only_vcf.nf.test
nextflow_process {
    name "Test make_sites_only_vcf"
    script "modules/local/make_sites_only_vcf.nf"
    process "MAKE_SITES_ONLY_VCF"

    test("produces a sites-only VCF for VQSR training") {
        when {
            process {
                """
                input[0] = [
                    [ participant_id: 'TEST_TUMOR_vs_TEST_NORMAL' ],
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz")
                ]
                """
            }
        }
        then {
            // Known, accepted gap: placeholder fixture, not real VCF content. Do not flip to !process.success.
            assert process.success
        }
    }
}
```

- [ ] **Step 4: Run both tests**

Run: `nf-test test tests/local/variant_filtration.nf.test tests/local/make_sites_only_vcf.nf.test`
Expected: both fail today on the documented placeholder-fixture gap (a GATK "not a valid
VCF" or BGZF-class error). Confirm the failure reason matches before moving on.

- [ ] **Step 5: Commit**

```bash
git add modules/local/variant_filtration.nf modules/local/make_sites_only_vcf.nf \
        tests/local/variant_filtration.nf.test tests/local/make_sites_only_vcf.nf.test
git commit -m "Add VARIANT_FILTRATION and MAKE_SITES_ONLY_VCF local modules"
```

---

### Task 5: `JOINT_GENOTYPING` subworkflow

**Files:**
- Create: `subworkflows/local/joint_genotyping/main.nf`
- Test: `tests/subworkflows/joint_genotyping.nf.test`

**Interfaces:**
- Consumes: `GATK4_GENOMICSDBIMPORT`, `GATK4_GENOTYPEGVCFS`, `VARIANTRECALIBRATOR_SNP`/
  `_INDEL`, `APPLYVQSR_INDEL`/`_SNP` (Task 3), `VARIANT_FILTRATION`, `MAKE_SITES_ONLY_VCF`
  (Task 4).
- Produces: `JOINT_GENOTYPING` workflow. `take: gvcf_pairs_ch` — `[pair_meta, normal_gvcf,
  normal_tbi, tumor_gvcf, tumor_tbi]`, one row per pair (Task 6 builds this by joining
  `SAMPLE_GERMLINE_CALLING.out.sample_gvcfs` against each pair's `normal_sample_id`/
  `tumor_sample_id`). `take: ref_fasta_ch` — `[fasta, fai, dict]`. `take:
  calling_intervals_ch` — `[intervals]`, value channel. `take: dbsnp_ch` — `[dbsnp_vcf,
  dbsnp_vcf_index]`. `take: snp_resources_ch` — `[[hapmap_vcf, omni_vcf, kg_vcf,
  dbsnp_vcf], [hapmap_tbi, omni_tbi, kg_tbi, dbsnp_tbi], [label1, label2, label3,
  label4]]` (a 3-element list: file list, index list, label list, in that fixed order —
  see Step 1 for why). `take: indel_resources_ch` — same 3-element shape, 3 resources
  instead of 4. `emit: joint_vcf` — `[pair_meta, vcf, tbi]`, replacing what `main.nf`
  currently reads straight from the samplesheet.

- [ ] **Step 1: Write the subworkflow**

```groovy
// subworkflows/local/joint_genotyping/main.nf
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
```

`meta - meta.subMap('id')` strips the synthetic `id` field added for `tag`/`prefix`
purposes back off before emitting, so the channel `main.nf` receives carries exactly the
same `pair_meta` shape `CNV_CORE` already expects (participant_id, normal_sample_id,
tumor_sample_id, sample_sex, and every tuning field) — not that map plus a stray `id` key.

- [ ] **Step 2: Write the subworkflow test**

```groovy
// tests/subworkflows/joint_genotyping.nf.test
nextflow_workflow {
    name "Test JOINT_GENOTYPING"
    script "subworkflows/local/joint_genotyping/main.nf"
    workflow "JOINT_GENOTYPING"

    test("produces one joint-genotyped, VQSR-filtered VCF per pair") {
        when {
            workflow {
                """
                def pair_meta = [
                    participant_id: 'TEST_TUMOR_vs_TEST_NORMAL',
                    patient: 'TEST_PATIENT',
                    normal_sample_id: 'TEST_NORMAL',
                    tumor_sample_id: 'TEST_TUMOR',
                    sample_sex: 'female'
                ]
                input[0] = Channel.of([
                    pair_meta,
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/joint_genotyping/placeholder.g.vcf.gz")
                ])
                input[1] = Channel.value([ file(params.ref_fasta), file(params.ref_fasta_fai), file(params.ref_fasta_dict) ])
                input[2] = Channel.value(file(params.ref_fasta))
                input[3] = Channel.value([ file(params.ref_fasta), file(params.ref_fasta) ])
                input[4] = Channel.value([
                    [ file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta) ],
                    [ file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta) ],
                    [ '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 x', '--resource:omni,known=false,training=true,truth=true,prior=12.0 x', '--resource:1000G,known=false,training=true,truth=false,prior=10.0 x', '--resource:dbsnp,known=true,training=false,truth=false,prior=7.0 x' ]
                ])
                input[5] = Channel.value([
                    [ file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta) ],
                    [ file(params.ref_fasta), file(params.ref_fasta), file(params.ref_fasta) ],
                    [ '--resource:mills,known=false,training=true,truth=true,prior=12.0 x', '--resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 x', '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 x' ]
                ])
                """
            }
        }
        then {
            // Known, accepted gap (Global Constraints above): placeholder GVCF fixtures
            // and stand-in resource files, not real WGS-scale joint-called data -- VQSR's
            // Gaussian mixture model cannot fit on this input regardless of container or
            // network access. Real fixture data closes this the same way it closes
            // tests/subworkflows/cnv_core.nf.test. Do not flip to !workflow.success.
            assert workflow.success
        }
    }
}
```

- [ ] **Step 3: Run the test**

Run: `nf-test test tests/subworkflows/joint_genotyping.nf.test`
Expected: FAIL today on the documented placeholder-fixture/VQSR-data gap. Read the
failure and confirm it's a GATK content/statistical error at some point in the chain
(GenomicsDBImport rejecting invalid VCF content, or VariantRecalibrator refusing to fit a
model on too few variants), not a channel-arity or missing-input error — the latter would
mean a wiring mistake in Step 1, not the accepted gap.

- [ ] **Step 4: Commit**

```bash
git add subworkflows/local/joint_genotyping/main.nf tests/subworkflows/joint_genotyping.nf.test
git commit -m "Add JOINT_GENOTYPING subworkflow"
```

---

### Task 6: Wire `main.nf` and add the new reference-input params

**Files:**
- Modify: `main.nf`
- Modify: `nextflow.config`

**Interfaces:**
- Consumes: `SAMPLE_GERMLINE_CALLING` (Task 2), `JOINT_GENOTYPING` (Task 5).
- Produces: `main.nf`'s `joint_vcf_ch` is now computed, not read from the samplesheet.

- [ ] **Step 1: Add the new params to `nextflow.config`**

Replace the placeholder `params.calling_intervals = null` added in Task 3 Step 3, and add
the VQSR resource params, right after the existing `params.pon`/`ref_clinical`/
`cohort_clinical` block:

```groovy
// Phase 2 (SPEC.md §12): calling + joint-genotyping reference inputs.
params.calling_intervals = 's3://gatk-test-data/intervals/hg38.even.handcurated.20k.intervals'

params.dbsnp_vcf              = 's3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf'
params.dbsnp_vcf_index        = 's3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx'
params.hapmap_resource_vcf       = 's3://broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz'
params.hapmap_resource_vcf_index = 's3://broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi'
params.omni_resource_vcf         = 's3://broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz'
params.omni_resource_vcf_index   = 's3://broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi'
params.one_thousand_genomes_resource_vcf       = 's3://broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
params.one_thousand_genomes_resource_vcf_index = 's3://broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi'
params.mills_resource_vcf        = 's3://broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
params.mills_resource_vcf_index  = 's3://broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi'
params.axiomPoly_resource_vcf       = 's3://broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz'
params.axiomPoly_resource_vcf_index = 's3://broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi'
```

- [ ] **Step 2: Wire `main.nf`**

Replace the current `joint_vcf_ch` pass-through and the block that builds
`tumor_counts_ch`, keeping everything else (including `tumor_counts_ch` itself)
unchanged:

```groovy
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
```

`SAMPLESHEET_TO_PAIRS.out.pairs` still emits `[pair_meta, joint_vcf, joint_vcf_index]` at
this point in the plan — Task 7 removes those two trailing fields once the samplesheet
schema itself changes. `gvcf_pairs_ch`'s construction above only reads `pair_meta` from
that tuple (`_jv, _jvi` are discarded), so it does not need to change again once Task 7
lands; do these two tasks in either order without breaking the other.

- [ ] **Step 3: Run the existing single-pair test to confirm nothing else broke**

Run: `nf-test test tests/subworkflows/`
Expected: `samplesheet_to_pairs.nf.test` still passes (untouched by this task).
`sample_read_counts.nf.test` and `cnv_core.nf.test` still fail on their pre-existing,
already-documented gaps — unchanged by this task. `sample_germline_calling.nf.test` and
`joint_genotyping.nf.test` behave as Tasks 2/5 already established. There is no
`main.nf`-level nf-test yet in this repo, so this task has no new integration test of its
own — Task 8 documents this as unchanged from Phase 1.5's posture, not a new gap this
task introduces.

- [ ] **Step 4: Commit**

```bash
git add main.nf nextflow.config
git commit -m "Wire main.nf to compute the joint VCF via SAMPLE_GERMLINE_CALLING + JOINT_GENOTYPING"
```

---

### Task 7: Drop `joint_vcf`/`joint_vcf_index` from the samplesheet schema

**Files:**
- Modify: `subworkflows/local/samplesheet_to_pairs/main.nf`
- Modify: `assets/samplesheet.csv`
- Modify: `tests/fixtures/samplesheet_to_pairs/multi_tumor.csv`
- Modify: `tests/fixtures/samplesheet_to_pairs/orphan_tumor.csv`
- Modify: `tests/fixtures/samplesheet_to_pairs/sex_mismatch.csv`
- Modify: `tests/fixtures/samplesheet_to_pairs/single_pair.csv`
- Modify: `tests/fixtures/samplesheet_to_pairs/two_normals.csv`
- Modify: `tests/subworkflows/samplesheet_to_pairs.nf.test`

**Interfaces:**
- Consumes: nothing from earlier tasks (this task is independent of Tasks 1-6; SPEC §12.4
  calls out this exact ordering-independence).
- Produces: `SAMPLESHEET_TO_PAIRS.out.pairs` now emits a bare `pair_meta` map per item
  (not a list/tuple) — no VCF/index.
  Any caller expecting the old 3-tuple (Task 6's `main.nf`, already written to discard the
  trailing two fields regardless) keeps working either way.

- [ ] **Step 1: Update `validateSamplesheetRow()` and the `pairs` construction**

In `subworkflows/local/samplesheet_to_pairs/main.nf`:

```groovy
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

    def meta = [ patient: row.patient, sample: row.sample, sex: row.sex, status: status ]
    def cram = file(resolveRelativeToProjectDir(row.cram), checkIfExists: true)
    def crai = file(resolveRelativeToProjectDir(row.crai), checkIfExists: true)
    [ meta, cram, crai ]
}
```

And update `rows`/`samples`/`tumor_ch`/`pairs` to match the new 3-element row shape
(dropping every `jv, jvi` field):

```groovy
    rows = Channel.fromPath(resolveRelativeToProjectDir(samplesheet), checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> validateSamplesheetRow(row) }

    samples = rows.map { meta, cram, crai -> [ meta, cram, crai ] }

    normal_ch = rows
        .filter { meta, _cram, _crai -> meta.status == 0 }
        .map { meta, _cram, _crai -> [ meta.patient, meta ] }
        .groupTuple()
        .map { patient, metas ->
            if (metas.size() > 1) {
                error "Patient '${patient}' has more than one normal (status=0) sample: ${metas.collect { it.sample }.join(', ')}"
            }
            [ patient, metas[0] ]
        }

    normal_patients_ch = normal_ch.map { patient, _meta -> patient }.toList().map { [normals: it] }

    tumor_ch = rows
        .filter { meta, _cram, _crai -> meta.status == 1 }
        .map { meta, _cram, _crai -> [ meta.patient, meta ] }
        .combine(normal_patients_ch)
        .map { patient, meta, wrapper ->
            if (!(patient in wrapper.normals)) {
                error "Tumor sample '${meta.sample}' (patient '${patient}') has no matching normal (status=0) sample in the samplesheet"
            }
            [ patient, meta ]
        }

    pairs = normal_ch.cross(tumor_ch)
        .map { normal, tumor ->
            def normal_meta = normal[1]
            def tumor_meta  = tumor[1]

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
            pair_meta
        }
        .ifEmpty { error "Samplesheet produced zero tumor/normal pairs -- check that it has at least one status=1 (tumor) row with a matching status=0 (normal) row for the same patient" }

    emit:
    samples // [ meta(patient,sample,sex,status), cram, crai ]
    pairs   // bare pair_meta map per item, not a list/tuple
```

Also update the `take:` doc comment at the top of the file from
`// path to samplesheet CSV (patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index)`
to `// path to samplesheet CSV (patient,sample,sex,status,cram,crai)`.

- [ ] **Step 2: Update `main.nf`'s reads of `SAMPLESHEET_TO_PAIRS.out.pairs`**

`pairs` now emits a bare `pair_meta`, not a 3-tuple. Update the two places in `main.nf`
(Task 6) that destructure it:

```groovy
    joint_vcf_ch = JOINT_GENOTYPING.out.joint_vcf // unchanged from Task 6

    gvcf_pairs_ch = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta -> [ [patient: pair_meta.patient, sample: pair_meta.normal_sample_id], pair_meta ] }
        // ...rest unchanged from Task 6, this line's .map { pair_meta, _jv, _jvi -> ... } becomes .map { pair_meta -> ... }

    tumor_counts_ch = SAMPLESHEET_TO_PAIRS.out.pairs
        .map { pair_meta ->
            [ [patient: pair_meta.patient, sample: pair_meta.tumor_sample_id], pair_meta ]
        }
        // ...rest unchanged from Task 6
```

If Task 6 was implemented before this task, grep `main.nf` for `_jv, _jvi` and remove
both occurrences now, matching the pattern above exactly.

- [ ] **Step 3: Update the samplesheet template and every test fixture CSV**

`assets/samplesheet.csv`:

```
patient,sample,sex,status,cram,crai
PATIENT_ID,PATIENT_ID_normal,female,0,/path/to/PATIENT_ID_normal.recal.cram,/path/to/PATIENT_ID_normal.recal.cram.crai
PATIENT_ID,PATIENT_ID_tumor,female,1,/path/to/PATIENT_ID_tumor.recal.cram,/path/to/PATIENT_ID_tumor.recal.cram.crai
```

For each of `tests/fixtures/samplesheet_to_pairs/{multi_tumor,orphan_tumor,sex_mismatch,single_pair,two_normals}.csv`:
drop the header's trailing `,joint_vcf,joint_vcf_index` and drop each data row's trailing
two comma-separated fields (whatever they currently hold — empty for normal rows, a file
path for tumor rows). Read each file first, then edit it — do not guess at a row's exact
current content before viewing it.

- [ ] **Step 4: Update `tests/subworkflows/samplesheet_to_pairs.nf.test`**

The existing assertions read `workflow.out.pairs[0][0]` / `pair_metas.collect { it[0] }` —
this indexing assumed `pairs` emitted a 3-tuple with `pair_meta` at position 0. Since
`pairs` now emits a bare `pair_meta`, update every such access to drop the `[0]`:

```groovy
            def pair_metas = workflow.out.pairs.collect { it }
```

and

```groovy
            def pair_meta = workflow.out.pairs[0]
```

Every other assertion in that file (`pair_metas*.participant_id`, `.normal_sample_id`,
etc.) already reads fields off whatever `pair_meta` resolves to, and needs no further
change once the two indexing lines above are fixed.

- [ ] **Step 5: Run the updated test**

Run: `nf-test test tests/subworkflows/samplesheet_to_pairs.nf.test`
Expected: PASS — all 5 existing test cases (multi-tumor dedup, single-pair backward
compat, duplicate-normal error, sex-mismatch error, orphan-tumor error), same as before
this task, now against the trimmed schema.

- [ ] **Step 6: Commit**

```bash
git add subworkflows/local/samplesheet_to_pairs/main.nf assets/samplesheet.csv \
        tests/fixtures/samplesheet_to_pairs/ tests/subworkflows/samplesheet_to_pairs.nf.test \
        main.nf
git commit -m "Drop joint_vcf/joint_vcf_index from the samplesheet schema"
```

---

### Task 8: Update `README.md`

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: the finished state of Tasks 1-7 (this task only documents, no code).

- [ ] **Step 1: Update the "Two logical halves" description**

Replace:

```markdown
1. Per-sample germline calling (Parabricks) + joint genotyping/VQSR (reuses
   nf-core/sarek's `BAM_JOINT_CALLING_GERMLINE_GATK` subworkflow directly).
```

with:

```markdown
1. Per-sample germline calling (`SAMPLE_GERMLINE_CALLING`, Parabricks HaplotypeCaller)
   + joint genotyping/VQSR (`JOINT_GENOTYPING`, GATK4 modules vendored from
   nf-core/sarek and composed with this pipeline's own no-scatter orchestration —
   see `SPEC.md` §12).
```

- [ ] **Step 2: Update the samplesheet column table**

Remove the `joint_vcf` / `joint_vcf_index` row from the column table (Task 7 dropped
them from the schema entirely):

```markdown
| Column | Required | Meaning |
| --- | --- | --- |
| `patient` | always | Groups a patient's normal and tumor rows for pairing. |
| `sample` | always | Sample ID, unique within the samplesheet. |
| `sex` | always | `female` or `male`. |
| `status` | always | `0` = normal, `1` = tumor. |
| `cram` / `crai` | always | The sample's recalibrated CRAM from `sarek_align`. |
```

- [ ] **Step 3: Add the new reference-input defaults to the params table**

Append to the existing "Reference defaults that may be worth overriding" table:

```markdown
| `--calling_intervals` | `s3://gatk-test-data/intervals/hg38.even.handcurated.20k.intervals` |
| `--dbsnp_vcf` / `--dbsnp_vcf_index` | `s3://broad-references/hg38/v0/...` |
| `--hapmap_resource_vcf` / `--omni_resource_vcf` / `--one_thousand_genomes_resource_vcf` / `--mills_resource_vcf` / `--axiomPoly_resource_vcf` (each with a matching `_index`) | `s3://broad-references/hg38/v0/...` |
```

- [ ] **Step 4: Add the third accepted gap to "Known gaps in this phase"**

Append a third numbered item, matching the existing two's style:

```markdown
**3. No real Parabricks/NGC execution has been validated.** This dev environment (and
the dev's own machine, confirmed) has no NVIDIA GPU, so `SAMPLE_GERMLINE_CALLING`'s
wiring is tested via nf-test's `-stub` flag rather than a real `pbrun` run. Real
Parabricks validation can only happen on Cirro/AWS against an actual GPU instance —
where it was always going to run in production anyway. `JOINT_GENOTYPING`'s GATK4 chain
also has no real fixture: VQSR's Gaussian mixture model cannot fit on a handful of
variants regardless of container or network access, so its tests assert the desired end
state against placeholder fixtures, the same convention as gap 1 above.
```

- [ ] **Step 5: Commit**

```bash
git add README.md
git commit -m "Document Phase 2 in README.md"
```

---

## Self-Review Notes (completed during plan authoring, not a task)

- **Spec coverage**: §12.1 (architecture/dedup) → Tasks 2, 5, 6. §12.2
  (`SAMPLE_GERMLINE_CALLING`) → Tasks 1, 2. §12.3 (`JOINT_GENOTYPING` chain, exact
  ext.args) → Tasks 3, 4, 5. §12.4 (schema change) → Task 7. §12.5 (new params) → Task 6.
  §12.6 (decision log: no scatter, vendoring approach, versions) → reflected throughout,
  no dedicated task needed since it's already-made decisions, not new code. §12.7
  (testing posture) → reflected in every task's test steps and Task 8's README update.
- **Placeholder scan**: every `ext.args` string, container digest, and CLI flag above is
  copied from a verified source (the WDL, the production Cirro config, or this repo's own
  existing pinned digest) — none are guessed. The one open question (Task 1 Step 4's
  nf-test config-injection syntax) is flagged for empirical verification during
  implementation rather than asserted as fact, consistent with this project's established
  "verify, don't assume" practice — it is not a logic gap, just an syntax detail this
  plan cannot resolve without running the installed nf-test version's actual parser.
- **Type consistency**: `pair_meta`'s shape (participant_id, patient, normal_sample_id,
  tumor_sample_id, sample_sex, plus 15 tuning fields) is unchanged from Phase 1.5 in
  Task 7's rewrite, and Task 5's `JOINT_GENOTYPING` strips its synthetic `id` field back
  off before emitting so `CNV_CORE`'s existing `take:` interface needs no changes.
  `SAMPLE_GERMLINE_CALLING.out.sample_gvcfs`'s meta (Task 2) keeps `patient`/`sample`
  from `SAMPLESHEET_TO_PAIRS.out.samples` unchanged, so Task 6's join-by-`[patient,
  sample]` key matches on both sides, the same pattern already proven for
  `tumor_counts_ch`.
