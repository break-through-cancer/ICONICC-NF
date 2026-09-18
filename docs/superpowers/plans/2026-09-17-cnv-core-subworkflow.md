# ICONICC → Nextflow Port, Phase 1: CNV Core Subworkflow — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up the repo scaffold and hand-port the 12-task CNV/allelic-copy-number core
(`ICON_Tasks.2.wdl`, §4 of the spec) as a self-contained, independently-testable Nextflow
subworkflow, fed by a known-good joint-genotyped VCF and tumor/normal CRAMs from an existing
completed WDL run — no GPU calling or joint-genotyping code needed yet.

**Architecture:** One process per WDL task (mostly 1:1; `ac_calc_agg`/`ac_calc_noagg` are
consolidated into a single parameterized process — see Task 12). Processes chain into a
`CNV_CORE` subworkflow under `subworkflows/local/`. `main.nf` invokes it directly against
fixture inputs for now; Phase 3 will replace those with the calling+joint-genotyping
subworkflow's real output channel.

**Tech Stack:** Nextflow DSL2, bcftools 1.10.2, GATK 4.2.6.1, Eagle2, R (custom container,
Task 2), nf-test for process/subworkflow testing. Follows this org's existing hand-rolled
Nextflow/Cirro conventions (see `break-through-cancer/btc-checksum` and
`break-through-cancer/ancestry-pipeline`), not the nf-core pipeline template — this repo does
vendor nf-core modules/subworkflows directly, but only where the spec calls for reusing sarek
(Phase 2's joint-genotyping half, not this one).

**Spec:** `SPEC.md` (repo root), §4, §4.1, §4.2, §6, §7. Reference source pulled for this port
lives in `reference/` (gitignored — not part of this repo): the actual WDL
(`reference/btc-wdl-pipelines/ICON_Tasks.2.wdl`, `ICONICC_paired.3.wdl` lines 610-788), the
three R scripts (`reference/r-scripts/*.R`), and the task-mapping TSV
(`reference/iconicc_sarek_task_mapping.tsv`).

## Global Constraints

- No richer Cirro-facing configuration than today: don't add new user-exposed params beyond
  what §3 already allows. Every optional WDL param that Cirro's `process-input.json` never
  set (see below) is pinned as a Nextflow `params.*` default in code, not exposed.
- Production defaults confirmed by reading the current WDL repo's
  `cirro/process-input.json` directly (not just the spec's summary): `desired_binsize` is
  never set → the `ac_calc_noagg` code path is the one that actually runs in production.
  `intervals`/`blacklist_intervals` are never set → `PreprocessIntervals` runs in
  whole-genome mode. `padding`/`readcts_binlength` are never set → task defaults `0`/`25000`
  apply. `include`/`exclude`/`softFilter` are never set → omit these bcftools flags entirely.
- The `hetvarFilter` positional-indexing bug fix (§4.2) is mandatory and must run before any
  other CNV-half processing: subset+reorder the joint VCF to `-s
  ${normal_sample_id},${tumor_sample_id}` immediately after it's handed off.
- GATK-based CNV tasks (`PreprocessIntervals`, `CollectReadCounts`) standardize on
  **GATK 4.2.6.1** (resolved decision, SPEC.md §7): container
  `broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe`
  (digest of `broadinstitute/gatk:4.2.6.1`, verified live — this is the same image as
  `us.gcr.io/broad-gatk/gatk:4.2.6.1`, confirmed identical digest).
- bcftools tasks keep the WDL's exact pinned version for reproducibility: container
  `quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a`
  (digest of `quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2`).
- Eagle's WDL container (`cbao/eagle:latest`) is unpinned, the same class of problem flagged
  for the R image in SPEC.md §7 but not explicitly asked about. Resolve it the same way:
  pin to the digest resolved live during this plan's authoring —
  `docker.io/cbao/eagle@sha256:e837cf05d2821e2bd40d6e6f7dac9ef0a5c1575e744a0e7026d28bfac0cf3d66`.
  Flag to the dev in the PR description that this is a third-party, not-BTC-maintained image
  pinned by digest rather than rebuilt — different treatment than the R image (resolved
  decision: rebuild from scratch) because Eagle itself isn't BTC-authored, unlike the R
  scripts.
- CNV R scripts (`ac_calc.R`, `tangent.R`, `segment.R`) run on a **freshly-built, versioned**
  image per the resolved decision (SPEC.md §7) — not a re-pin of `wchukwu/r-docker:latest`.
  Dockerfile lives in this repo at `containers/iconicc-cnv-r/` (Task 2); building and
  publishing it is a manual step the dev does outside this plan (deliberate scope cut — see
  Task 2 — to avoid an autonomous PR/merge into the shared `btc-containers` registry).
- The three R scripts themselves are vendored into this repo's `bin/` directory (Nextflow
  puts `bin/` on `PATH` automatically) rather than fetched from S3 per-run via a `params.*`
  File input the way the WDL did — same scripts, version-locked with the pipeline code
  instead of a mutable external S3 object. This is a deliberate simplification from the
  WDL's approach, not a spec requirement; call it out in the PR description.

---

## Task 1: Repo scaffold

**Files:**
- Create: `nextflow.config`
- Create: `main.nf`
- Create: `README.md`
- Modify: `.gitignore` (already has `/reference/`, `.nextflow*`, `work/`, `results/`, add
  `tests/fixtures/`)
- Test: none (structural) — verified by running the stub pipeline

**Interfaces:**
- Produces: a `main.nf` that Nextflow can parse and run (`workflow { }` entrypoint), which
  later tasks extend by adding an `include` for `CNV_CORE` and wiring it up.

- [ ] **Step 1: Write `nextflow.config`**

```groovy
nextflow.enable.dsl = 2

params.outdir = './results'

process {
    publishDir = [
        path: { "${params.outdir}/${task.process.tokenize(':')[-1].toLowerCase()}" },
        mode: 'copy',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
}

docker {
    enabled = true
}
```

- [ ] **Step 2: Write a stub `main.nf`**

```groovy
nextflow.enable.dsl = 2

workflow {
    log.info "ICONICC-NF scaffold OK"
}
```

- [ ] **Step 3: Verify it runs**

Run: `nextflow run main.nf`
Expected: completes successfully, logs `ICONICC-NF scaffold OK`.

- [ ] **Step 4: Write `README.md`**

```markdown
# ICONICC-NF

Nextflow port of ICONICC ("HapCNA"), a tumor/normal paired WGS pipeline doing germline
joint-genotyping followed by a haplotype-based allelic copy-number (CNV) analysis.

See `SPEC.md` for the full design. Two logical halves:

1. Per-sample germline calling (Parabricks) + joint genotyping/VQSR (reuses
   nf-core/sarek's `BAM_JOINT_CALLING_GERMLINE_GATK` subworkflow directly).
2. The CNV/allelic-analysis core, hand-ported from `ICON_Tasks.2.wdl` — the actual
   science; see `subworkflows/local/cnv_core/`.

## Local usage

```bash
nextflow run main.nf --help
```

## Testing

```bash
nf-test test
```
```

- [ ] **Step 5: Update `.gitignore`**

Add `tests/fixtures/` (real genomic data, not committed) to the existing `.gitignore`.

- [ ] **Step 6: Commit**

```bash
git add nextflow.config main.nf README.md .gitignore
git commit -m "Scaffold ICONICC-NF Nextflow pipeline"
```

---

## Task 2: Write the CNV R container Dockerfile (in this repo)

This produces the Dockerfile every R-based process in this plan depends on (`ac_calc`,
`tangent_xy`, `segment` — Tasks 12-14), so it has to land before those, but it's independent
of the bcftools/GATK chain (Tasks 3-11) and can happen in parallel with them.

**Deliberate scope cut (dev's call, not the spec's):** keep this entirely inside the
`ICONICC-NF` repo — no cross-repo PR into the shared `break-through-cancer/btc-containers`
registry. Building and publishing the image (`docker build` + push to wherever the dev wants
— GHCR, ECR, `btc-containers` later) is a manual step the dev does on their own machine when
ready; it's a side effect on shared infrastructure outside this repo, not something to
automate here. To keep that decision from blocking Tasks 12-14, the R-based modules
reference a **configurable** `params.cnv_r_container` (default: a local tag) rather than a
hardcoded registry digest — the dev overrides it at runtime (`--cnv_r_container ...`) once
the image is built and pushed somewhere real.

**Files:**
- Create: `containers/iconicc-cnv-r/Dockerfile`
- Create: `containers/iconicc-cnv-r/README.md` (build/publish instructions for the dev)
- Modify: `nextflow.config` (add `params.cnv_r_container` default)

**Interfaces:**
- Produces: `params.cnv_r_container` (default `'iconicc-cnv-r:latest'`), referenced by
  every `container` directive in Tasks 12-14 instead of a literal image string.

**Package inventory** (from grepping `library()`/`require()` calls across all three R
scripts — see `reference/r-scripts/*.R`):
- Base R (no install needed): `tools`, `parallel`
- CRAN: `data.table`, `reshape2`, `scales`, `latex2exp`, `ggplot2`, `useful`, `tidyverse`,
  `pracma`, `optparse`, `dplyr`, `gridExtra`, `Hmisc`
- Bioconductor: `rhdf5`, `rCGH`, `DNAcopy`, `GenomicRanges`
- GitHub (not on CRAN/Bioconductor): `gUtils` (`mskilab-org/gUtils`)

- [ ] **Step 1: Write the Dockerfile**

Pin the whole R environment to the Bioconductor 3.20 / R 4.4 release train (current at
spec-writing time) via its base image, and CRAN packages to a frozen Posit Package Manager
snapshot dated to the same window, so `install.packages()` calls are themselves
reproducible instead of always resolving "latest CRAN":

```dockerfile
# containers/iconicc-cnv-r/Dockerfile
FROM bioconductor/bioconductor_docker:RELEASE_3_20

RUN echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2025-01-15"))' \
    >> /usr/local/lib/R/etc/Rprofile.site

RUN Rscript -e 'install.packages(c( \
      "data.table", "reshape2", "scales", "latex2exp", "ggplot2", "useful", \
      "tidyverse", "pracma", "optparse", "dplyr", "gridExtra", "Hmisc", "remotes" \
    ))'

RUN Rscript -e 'BiocManager::install(c("rhdf5", "rCGH", "DNAcopy", "GenomicRanges"), \
      update = FALSE, ask = FALSE)'

RUN Rscript -e 'remotes::install_github( \
      "mskilab-org/gUtils@095377ab04195aecf6ebb8397072b24864c81599")'
```

- [ ] **Step 2: Write a short README for the dev**

```markdown
# containers/iconicc-cnv-r

R environment for the CNV/allelic-analysis core (`ac_calc.R`, `tangent.R`, `segment.R`).
See `SPEC.md` §7 for the pinning rationale.

## Build and publish (manual — do this yourself when ready)

```bash
docker build -t iconicc-cnv-r:latest containers/iconicc-cnv-r
# then push wherever you want this to live, e.g.:
#   docker tag iconicc-cnv-r:latest ghcr.io/<you>/iconicc-cnv-r:latest
#   docker push ghcr.io/<you>/iconicc-cnv-r:latest
```

Then point the pipeline at it:

```bash
nextflow run main.nf --cnv_r_container ghcr.io/<you>/iconicc-cnv-r:latest ...
```

Until you've built and published it, `params.cnv_r_container` defaults to the local tag
`iconicc-cnv-r:latest` — fine for local Docker testing, not resolvable in a remote executor.
```

- [ ] **Step 3: Wire the default param**

```groovy
// append to nextflow.config
params.cnv_r_container = 'iconicc-cnv-r:latest'
```

- [ ] **Step 4: Commit**

```bash
git add containers/iconicc-cnv-r/Dockerfile containers/iconicc-cnv-r/README.md nextflow.config
git commit -m "Add CNV R container Dockerfile (build/publish is a manual step for the dev)"
```

---

## Task 3: Stage validation fixtures + `subset_pair_vcf` module (§4.2 fix)

**This task has an external dependency the executor cannot self-serve**: a completed,
known-good ICONICC WDL run to validate against. The Cirro project's data bucket
(`btc-gbm`, project `BTC-GBM-Production`) stores datasets under opaque hashed keys with no
human-browsable dataset index reachable from this environment's tools — asking the dev
directly is faster and more reliable than guessing.

**Files:**
- Create: `modules/local/subset_pair_vcf.nf`
- Create: `tests/local/subset_pair_vcf.nf.test`
- Create (fixture staging, not code): `tests/fixtures/cnv_core/*`

**Interfaces:**
- Consumes: `tuple val(meta), path(vcf), path(vcf_index)` where `meta` is
  `[participant_id: String, normal_sample_id: String, tumor_sample_id: String]`
- Produces: `tuple val(meta), path("*.pair_subset.vcf.gz"), path("*.pair_subset.vcf.gz.tbi")`
  — every downstream CNV-core process consumes this exact tuple shape for `meta`.

- [ ] **Step 1: Ask the dev for fixture data**

Ask: "For validating the CNV-core port in isolation, I need a completed ICONICC WDL run's
`output_vcf_file`/`output_vcf_index_file` (the joint-genotyped, VQSR-filtered VCF) plus
both samples' `recal.cram`/`.crai` from `sarek_align`. Cirro dataset ID and project, or
direct S3 paths, either works." Record the real `participant_id`/`normal_sample_id`/
`tumor_sample_id` for that run — later steps need the real IDs, not placeholders.

- [ ] **Step 2: Stage the fixtures under fixed, predictable names**

```bash
mkdir -p tests/fixtures/cnv_core
aws s3 cp <output_vcf_file_path>       tests/fixtures/cnv_core/joint.vcf.gz
aws s3 cp <output_vcf_index_file_path> tests/fixtures/cnv_core/joint.vcf.gz.tbi
aws s3 cp <normal_recal_cram_path>     tests/fixtures/cnv_core/normal.recal.cram
aws s3 cp <normal_recal_crai_path>     tests/fixtures/cnv_core/normal.recal.cram.crai
aws s3 cp <tumor_recal_cram_path>      tests/fixtures/cnv_core/tumor.recal.cram
aws s3 cp <tumor_recal_crai_path>      tests/fixtures/cnv_core/tumor.recal.cram.crai
cat > tests/fixtures/cnv_core/sample_info.env <<EOF
PARTICIPANT_ID=<real participant_id>
NORMAL_SAMPLE_ID=<real normal_sample_id>
TUMOR_SAMPLE_ID=<real tumor_sample_id>
SAMPLE_SEX=<real sample_sex>
EOF
```

Every later nf-test in this plan references these fixed paths — this is the one and only
place real dataset identifiers get plugged in.

- [ ] **Step 3: Write the failing nf-test**

```groovy
// tests/local/subset_pair_vcf.nf.test
nextflow_process {
    name "Test subset_pair_vcf"
    script "modules/local/subset_pair_vcf.nf"
    process "SUBSET_PAIR_VCF"

    test("subsets and reorders to normal,tumor") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT', normal_sample_id: 'TEST_NORMAL', tumor_sample_id: 'TEST_TUMOR'],
                    file("\${projectDir}/tests/fixtures/cnv_core/joint.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/cnv_core/joint.vcf.gz.tbi")
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).exists()
            assert path(process.out[0].get(0)[2]).exists()
        }
    }
}
```

- [ ] **Step 4: Run it, verify it fails**

Run: `nf-test test tests/local/subset_pair_vcf.nf.test`
Expected: FAIL — `modules/local/subset_pair_vcf.nf` doesn't exist yet.

- [ ] **Step 5: Implement the module**

```groovy
// modules/local/subset_pair_vcf.nf
process SUBSET_PAIR_VCF {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*.pair_subset.vcf.gz"), path("*.pair_subset.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    bcftools view -s ${meta.normal_sample_id},${meta.tumor_sample_id} ${vcf} \\
        -O z -o ${meta.participant_id}.pair_subset.vcf.gz
    bcftools index --tbi ${meta.participant_id}.pair_subset.vcf.gz
    """
}
```

Note the fixed `-s normal,tumor` order — this is the §4.2 fix. It pins column 0 to normal
and column 1 to tumor for every downstream bcftools positional-index expression
(`hetvarFilter`'s `GT[1]`/`AD[1:0]`/`AD[1:1]` in Task 7), regardless of what order the
joint-genotyping subworkflow happened to emit samples in.

- [ ] **Step 6: Run it, verify it passes**

Run: `nf-test test tests/local/subset_pair_vcf.nf.test`
Expected: PASS

- [ ] **Step 7: Commit**

```bash
git add modules/local/subset_pair_vcf.nf tests/local/subset_pair_vcf.nf.test
git commit -m "Add subset_pair_vcf module, fixing hetvarFilter's positional-index bug (SPEC.md §4.2)"
```

(Fixtures under `tests/fixtures/` stay untracked, per `.gitignore` — real genomic data
doesn't belong in git. Anyone re-running these tests re-stages fixtures via Task 3 Step 2.)

---

## Task 4: `preprocess_intervals` module

**Files:**
- Create: `modules/local/preprocess_intervals.nf`
- Create: `tests/local/preprocess_intervals.nf.test`

**Interfaces:**
- Consumes: `tuple path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)`
- Produces: `path("*.preprocessed.interval_list")` — consumed by Task 5
  (`collect_read_counts`) as its `intervals` input.

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/preprocess_intervals.nf.test
nextflow_process {
    name "Test preprocess_intervals"
    script "modules/local/preprocess_intervals.nf"
    process "PREPROCESS_INTERVALS"

    test("produces a WGS-wide preprocessed interval list") {
        when {
            process {
                """
                input[0] = [
                    file(params.ref_fasta),
                    file(params.ref_fasta_fai),
                    file(params.ref_fasta_dict)
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)).exists()
            assert path(process.out[0].get(0)).text.contains("@SQ")
        }
    }
}
```

Reference resources (`ref_fasta`, `.fai`, `.dict`) are the same
`s3://broad-references/hg38/v0/Homo_sapiens_assembly38.*` files SPEC.md §5 already lists —
wire them as `params.ref_fasta` etc. with those S3 paths as defaults in `nextflow.config`
(added in this step since `preprocess_intervals` is the first process to need them):

```groovy
// append to nextflow.config
params.ref_fasta      = 's3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta'
params.ref_fasta_fai  = 's3://broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai'
params.ref_fasta_dict = 's3://broad-references/hg38/v0/Homo_sapiens_assembly38.dict'
```

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/preprocess_intervals.nf.test`
Expected: FAIL — module doesn't exist.

- [ ] **Step 3: Implement the module**

WGS mode only (production never sets `intervals`/`blacklist_intervals` — see Global
Constraints):

```groovy
// modules/local/preprocess_intervals.nf
process PREPROCESS_INTERVALS {
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 1
    memory '2 GB'

    input:
    tuple path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    output:
    path("wgs.preprocessed.interval_list")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx1500m" PreprocessIntervals \\
        --reference ${ref_fasta} \\
        --padding 0 \\
        --bin-length 25000 \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --output wgs.preprocessed.interval_list
    """
}
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/preprocess_intervals.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/preprocess_intervals.nf tests/local/preprocess_intervals.nf.test nextflow.config
git commit -m "Add preprocess_intervals module (WGS mode)"
```

---

## Task 5: `collect_read_counts` module

Called twice (tumor, normal) — a single parameterized process, not two near-duplicates.

**Files:**
- Create: `modules/local/collect_read_counts.nf`
- Create: `tests/local/collect_read_counts.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), val(sample_id), path(cram), path(crai), path(intervals),
  path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)`. **`meta` here is the same
  uniform per-participant map used by every other process in this subworkflow**
  (`participant_id`, `normal_sample_id`, `tumor_sample_id`, `sample_sex`, plus the R-script
  tuning fields listed in Task 14) — it is passed through unchanged, never reshaped, so
  every downstream `.join()` on `meta` matches by structural equality. `sample_id` is a
  *separate* plain string value (`meta.tumor_sample_id` or `meta.normal_sample_id`,
  supplied by the caller), used only for output naming — it does not go into `meta` itself,
  precisely so this process's two aliased call sites (`COLLECT_READ_COUNTS_TUMOR`,
  `COLLECT_READ_COUNTS_NORMAL` — Task 14) don't have to hand it two structurally different
  metas.
- Produces: `tuple val(meta), path("*_counts.txt")` — Task 12 (`ac_calc`) and Task 13
  (`tangent_xy`) consume the tumor call's output via a plain `meta`-keyed `.join()`; the
  normal call's output is only consumed inside this same subworkflow's wiring (Task 14),
  not by any R script per §4.1 (tangent normalizes tumor against an external PoN, not the
  matched normal).

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/collect_read_counts.nf.test
nextflow_process {
    name "Test collect_read_counts"
    script "modules/local/collect_read_counts.nf"
    process "COLLECT_READ_COUNTS"

    test("produces a per-bin TSV read-counts file for the tumor CRAM") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT', normal_sample_id: 'TEST_NORMAL', tumor_sample_id: 'TEST_TUMOR'],
                    'TEST_TUMOR',
                    file("\${projectDir}/tests/fixtures/cnv_core/tumor.recal.cram"),
                    file("\${projectDir}/tests/fixtures/cnv_core/tumor.recal.cram.crai"),
                    file("\${projectDir}/tests/fixtures/cnv_core/wgs.preprocessed.interval_list"),
                    file(params.ref_fasta),
                    file(params.ref_fasta_fai),
                    file(params.ref_fasta_dict)
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).exists()
        }
    }
}
```

(This test needs `wgs.preprocessed.interval_list` staged in fixtures too — copy the output
of a manual `nextflow run` of Task 4's module into `tests/fixtures/cnv_core/` once, or run
Task 4's process standalone and save its output there. Note this in the test file as an
inline comment since it's a real setup dependency, not a placeholder.)

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/collect_read_counts.nf.test`
Expected: FAIL — module doesn't exist.

- [ ] **Step 3: Implement the module**

Full read-filter set copied verbatim from `ICON_Tasks.2.wdl`'s `CollectReadCounts` task —
this is exact reproduction, not simplification. CRAM goes in directly (GATK4 tools accept
CRAM+reference natively — no `CramToBamTask` equivalent needed here either, extending the
same simplification SPEC.md §5 already applies to Parabricks):

```groovy
// modules/local/collect_read_counts.nf
process COLLECT_READ_COUNTS {
    tag "${sample_id}"
    container 'broadinstitute/gatk@sha256:21c3cb43b7d11891ed4b63cc7274f20187f00387cfaa0433b3e7991b5be34dbe'
    cpus 2
    memory '10 GB'

    input:
    tuple val(meta), val(sample_id), path(cram), path(crai), path(intervals), path(ref_fasta), path(ref_fasta_fai), path(ref_fasta_dict)

    output:
    tuple val(meta), path("${sample_id}_counts.txt")

    script:
    """
    set -euo pipefail
    gatk --java-options "-Xmx9G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \\
        CollectReadCounts \\
        -R ${ref_fasta} \\
        -I ${cram} \\
        --read-index ${crai} \\
        -L ${intervals} \\
        -O ${sample_id}_counts.txt \\
        --format TSV \\
        --interval-merging-rule OVERLAPPING_ONLY \\
        --read-filter PassesVendorQualityCheckReadFilter --read-filter HasReadGroupReadFilter \\
        --read-filter NotDuplicateReadFilter --read-filter MappingQualityAvailableReadFilter \\
        --read-filter PairedReadFilter --read-filter FragmentLengthReadFilter --max-fragment-length 1000 \\
        --read-filter MateOnSameContigOrNoMappedMateReadFilter --read-filter MateDifferentStrandReadFilter \\
        --read-filter NotSupplementaryAlignmentReadFilter --read-filter NotSecondaryAlignmentReadFilter \\
        --read-filter FirstOfPairReadFilter --read-filter MappingQualityReadFilter --minimum-mapping-quality 30 \\
        --read-filter OverclippedReadFilter --filter-too-short 50 --seconds-between-progress-updates 100
    """
}
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/collect_read_counts.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/collect_read_counts.nf tests/local/collect_read_counts.nf.test
git commit -m "Add collect_read_counts module"
```

---

## Task 6: `filter_variants` module

**Files:**
- Create: `modules/local/filter_variants.nf`
- Create: `tests/local/filter_variants.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(vcf), path(vcf_index)` (the `subset_pair_vcf` output)
- Produces: `tuple val(meta), path("*_HC_var.vcf.gz"), path("*_HC_var.vcf.gz.tbi")` —
  consumed by Task 7 (`iter_isec`)

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/filter_variants.nf.test
nextflow_process {
    name "Test filter_variants"
    script "modules/local/filter_variants.nf"
    process "FILTER_VARIANTS"

    test("keeps only variants with combined ALT depth >= 5") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT'],
                    file("\${projectDir}/tests/fixtures/cnv_core/pair_subset.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/cnv_core/pair_subset.vcf.gz.tbi")
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).exists()
        }
    }
}
```

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/filter_variants.nf.test`
Expected: FAIL

- [ ] **Step 3: Implement the module**

```groovy
// modules/local/filter_variants.nf
process FILTER_VARIANTS {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_index)

    output:
    tuple val(meta), path("*_HC_var.vcf.gz"), path("*_HC_var.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    bcftools filter -i 'SUM(FORMAT/AD[*:1])>=5' ${vcf} -O z -o ${meta.participant_id}_HC_var.vcf.gz
    bcftools index --tbi ${meta.participant_id}_HC_var.vcf.gz
    """
}
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/filter_variants.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/filter_variants.nf tests/local/filter_variants.nf.test
git commit -m "Add filter_variants module"
```

---

## Task 7: `iter_isec` module

**Files:**
- Create: `modules/local/iter_isec.nf`
- Create: `tests/local/iter_isec.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(vcf), path(vcf_index), path(ref_dir_1000g)`
- Produces: `tuple val(meta), path("chr_outputs/*_var.vcf.gz"),
  path("chr_outputs/*_var.vcf.gz.tbi")` — per-chromosome private-variant VCFs, consumed by
  Task 8 (`hetvar_filter`) and (chr1 only) Task 10 (`query_het_ad`, for its sample-name
  lookup).

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/iter_isec.nf.test
nextflow_process {
    name "Test iter_isec"
    script "modules/local/iter_isec.nf"
    process "ITER_ISEC"

    test("produces one private-variant VCF per autosome+X") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT'],
                    file("\${projectDir}/tests/fixtures/cnv_core/HC_var.vcf.gz"),
                    file("\${projectDir}/tests/fixtures/cnv_core/HC_var.vcf.gz.tbi"),
                    file(params.ref_dir_1000g)
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).size() == 23  // chr1-22 + chrX
        }
    }
}
```

Add the reference tarball param to `nextflow.config` in this step:

```groovy
// append to nextflow.config
params.ref_dir_1000g = 's3://btc-refdata/iconicc/reference-genome-1000G.tar'
```

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/iter_isec.nf.test`
Expected: FAIL

- [ ] **Step 3: Implement the module**

Reproduces the WDL task exactly, including its known quirk of overwriting a single
`_isec.vcf.gz` scalar file every loop iteration (harmless — nothing downstream consumes
that particular output; only the per-chromosome `chr_outputs/` files matter, and the WDL
itself never flagged this as a bug to fix):

```groovy
// modules/local/iter_isec.nf
process ITER_ISEC {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'
    disk '2000 GB'

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(ref_dir_1000g)

    output:
    tuple val(meta), path("chr_outputs/*_var.vcf.gz"), path("chr_outputs/*_var.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    mkdir -p tmp/chr_outputs
    tar -xf ${ref_dir_1000g}

    for i in {1..22} X; do
      bcftools isec ${vcf} hg38/1000G/chr\${i}.1000G.genotypes.bcf -w 1 -O z -p tmp
      mv tmp/0002.vcf.gz     tmp/chr_outputs/${meta.participant_id}_chr\${i}_var.vcf.gz
      mv tmp/0002.vcf.gz.tbi tmp/chr_outputs/${meta.participant_id}_chr\${i}_var.vcf.gz.tbi
    done
    mv tmp/chr_outputs chr_outputs
    """
}
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/iter_isec.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/iter_isec.nf tests/local/iter_isec.nf.test nextflow.config
git commit -m "Add iter_isec module"
```

---

## Task 8: `hetvar_filter` module

**Files:**
- Create: `modules/local/hetvar_filter.nf`
- Create: `tests/local/hetvar_filter.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(chr_vcfs), path(chr_vcf_indices)` (flat lists from
  `iter_isec`)
- Produces: `tuple val(meta), path("*_het.vcf.gz"), path("*_het.vcf.gz.tbi")` — consumed by
  Task 9 (`eagle_phasing`) and Task 10 (`query_het_ad`)

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/hetvar_filter.nf.test
nextflow_process {
    name "Test hetvar_filter"
    script "modules/local/hetvar_filter.nf"
    process "HETVAR_FILTER"

    test("keeps only tumor-het sites with adequate depth") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT'],
                    files("\${projectDir}/tests/fixtures/cnv_core/chr_outputs/*_var.vcf.gz"),
                    files("\${projectDir}/tests/fixtures/cnv_core/chr_outputs/*_var.vcf.gz.tbi")
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).size() == 23
        }
    }
}
```

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/hetvar_filter.nf.test`
Expected: FAIL

- [ ] **Step 3: Implement the module**

`GT[1]`/`AD[1:*]` (index 1 = second sample = tumor) is now safe because
`subset_pair_vcf` (Task 3) already pinned column order to `[normal, tumor]`:

```groovy
// modules/local/hetvar_filter.nf
process HETVAR_FILTER {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(chr_vcfs), path(chr_vcf_indices)

    output:
    tuple val(meta), path("*_het.vcf.gz"), path("*_het.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    for chr in ${chr_vcfs}; do
      filename=\$(basename "\${chr}" "_var.vcf.gz")"_het.vcf.gz"
      bcftools filter -i 'GT[1]="het" & FORMAT/AD[1:0]>=5 & FORMAT/AD[1:1]>=5' \\
        "\${chr}" -O z -o "\${filename}"
      bcftools index --tbi "\${filename}"
    done
    """
}
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/hetvar_filter.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/hetvar_filter.nf tests/local/hetvar_filter.nf.test
git commit -m "Add hetvar_filter module"
```

---

## Task 9: `eagle_phasing` module

**Files:**
- Create: `modules/local/eagle_phasing.nf`
- Create: `tests/local/eagle_phasing.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(het_vcfs), path(ref_dir_1000g), path(genetic_map)`
- Produces: `tuple val(meta), path("*_het_eagle2phased.vcf.gz")` — consumed by Task 10
  (`index_phased`, Task 11 `query_het_ad`)

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/eagle_phasing.nf.test
nextflow_process {
    name "Test eagle_phasing"
    script "modules/local/eagle_phasing.nf"
    process "EAGLE_PHASING"

    test("phases het sites against the 1000G reference panel") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT'],
                    files("\${projectDir}/tests/fixtures/cnv_core/het/*_het.vcf.gz"),
                    file(params.ref_dir_1000g),
                    file(params.genetic_map)
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).size() == 23
        }
    }
}
```

Add the genetic map param to `nextflow.config`:

```groovy
// append to nextflow.config
params.genetic_map = 's3://btc-refdata/iconicc/reference-genome-genetic_map_hg38_withX.txt.gz'
```

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/eagle_phasing.nf.test`
Expected: FAIL

- [ ] **Step 3: Implement the module**

```groovy
// modules/local/eagle_phasing.nf
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
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/eagle_phasing.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/eagle_phasing.nf tests/local/eagle_phasing.nf.test nextflow.config
git commit -m "Add eagle_phasing module"
```

---

## Task 10: `index_phased` module

**Files:**
- Create: `modules/local/index_phased.nf`
- Create: `tests/local/index_phased.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(phased_vcfs)`
- Produces: `tuple val(meta), path("*_het_eagle2phased.vcf.gz.tbi")` — consumed by Task 11
  (`query_het_ad`)

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/index_phased.nf.test
nextflow_process {
    name "Test index_phased"
    script "modules/local/index_phased.nf"
    process "INDEX_PHASED"

    test("tabix-indexes each phased VCF") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT'],
                    files("\${projectDir}/tests/fixtures/cnv_core/phased/*_het_eagle2phased.vcf.gz")
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).size() == 23
        }
    }
}
```

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/index_phased.nf.test`
Expected: FAIL

- [ ] **Step 3: Implement the module**

```groovy
// modules/local/index_phased.nf
process INDEX_PHASED {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    cpus 2
    memory '4 GB'

    input:
    tuple val(meta), path(phased_vcfs)

    output:
    tuple val(meta), path("*_het_eagle2phased.vcf.gz.tbi")

    script:
    """
    set -euo pipefail
    for path in ${phased_vcfs}; do
      bcftools index --tbi "\${path}"
    done
    """
}
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/index_phased.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/index_phased.nf tests/local/index_phased.nf.test
git commit -m "Add index_phased module"
```

---

## Task 11: `query_het_ad` module

**Files:**
- Create: `modules/local/query_het_ad.nf`
- Create: `tests/local/query_het_ad.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(chr1_het), path(all_het), path(all_het_indices),
  path(all_eagle_phased), path(all_eagle_indices)`
- Produces: `tuple val(meta), path("*_HetInfo.txt"), path("*.AD.chr*.txt"),
  path("*.phasedGT.chr*.1.txt")` — all three globs feed Task 12 (`ac_calc`) directly.

- [ ] **Step 1: Write the failing nf-test**

```groovy
// tests/local/query_het_ad.nf.test
nextflow_process {
    name "Test query_het_ad"
    script "modules/local/query_het_ad.nf"
    process "QUERY_HET_AD"

    test("extracts per-chr HetInfo, allelic depths, and phased genotypes") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT'],
                    file("\${projectDir}/tests/fixtures/cnv_core/het/TEST_PT_chr1_het.vcf.gz"),
                    files("\${projectDir}/tests/fixtures/cnv_core/het/*_het.vcf.gz"),
                    files("\${projectDir}/tests/fixtures/cnv_core/het/*_het.vcf.gz.tbi"),
                    files("\${projectDir}/tests/fixtures/cnv_core/phased/*_het_eagle2phased.vcf.gz"),
                    files("\${projectDir}/tests/fixtures/cnv_core/phased/*_het_eagle2phased.vcf.gz.tbi")
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).exists()
        }
    }
}
```

- [ ] **Step 2: Run it, verify it fails**

Run: `nf-test test tests/local/query_het_ad.nf.test`
Expected: FAIL

- [ ] **Step 3: Implement the module**

```groovy
// modules/local/query_het_ad.nf
process QUERY_HET_AD {
    tag "${meta.participant_id}"
    container 'quay.io/biocontainers/bcftools@sha256:3130439789ced19fb4ba6e7513d0d12f9bdf3de0618f5b1ce2c76b1835ae347a'
    memory '4 GB'

    input:
    tuple val(meta), path(chr1_het), path(all_het), path(all_het_indices), path(all_eagle_phased), path(all_eagle_indices)

    output:
    tuple val(meta), path("*_HetInfo.txt"), path("*.AD.chr*.txt"), path("*.phasedGT.chr*.1.txt")

    script:
    """
    set -euo pipefail
    samples=\$(bcftools query -l ${chr1_het})

    for path in ${all_het}; do
      hetinfofile=\$(basename "\${path}" "_het.vcf.gz")"_HetInfo.txt"
      chr=\$(basename "\${path}" "_het.vcf.gz" | sed 's/^.*chr//')
      bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT{0}\\t%INFO/DP\\n' "\${path}" > "\${hetinfofile}"
      for sample in \${samples}; do
        bcftools query -s \${sample} -f '[%AD{0}\\t%AD{1}\\n]' "\${path}" > "\${sample}.AD.chr\${chr}.txt"
      done
    done

    for phased in ${all_eagle_phased}; do
      chr=\$(basename "\${phased}" "_het_eagle2phased.vcf.gz" | sed 's/^.*chr//')
      for sample in \${samples}; do
        bcftools query -s \${sample} -f '[%GT\\n]' "\${phased}" > "\${sample}.phasedGT.chr\${chr}.1.txt"
      done
    done
    """
}
```

- [ ] **Step 4: Run it, verify it passes**

Run: `nf-test test tests/local/query_het_ad.nf.test`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add modules/local/query_het_ad.nf tests/local/query_het_ad.nf.test
git commit -m "Add query_het_ad module"
```

---

## Task 12: `ac_calc` module

Consolidates the WDL's `ac_calc_agg`/`ac_calc_noagg` (identical scripts, different arg
lists) into one process branching on whether `desired_binsize` is set — cleaner in
Nextflow than duplicating two near-identical processes, and production only ever exercises
the `noagg` path anyway (Global Constraints).

**Files:**
- Create: `bin/ac_calc.R` (copy of `reference/r-scripts/reference-genome-AC_calc_20231228.R`
  — vendor, don't reference by path)
- Create: `modules/local/ac_calc.nf`
- Create: `tests/local/ac_calc.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(tumor_counts), path(phased_files), path(ad_files),
  path(hetinfo_files)` — `meta` is the same uniform per-participant map threaded through
  the whole subworkflow (see Task 5's note); it must carry `participant_id`,
  `normal_sample_id`, `tumor_sample_id`, `agg_bins`, and `desired_binsize` (`null` when
  unset — the production default per Global Constraints)
- Produces: `tuple val(meta), path("*_read_depth.txt"), path("*_allelic_counts.txt"),
  path("*_plots.pdf")` — `read_depth` is a final pipeline output (SPEC.md §6);
  `allelic_counts` feeds Task 14 (`segment`)

- [ ] **Step 1: Vendor the R script**

```bash
mkdir -p bin
cp reference/r-scripts/reference-genome-AC_calc_20231228.R bin/ac_calc.R
```

- [ ] **Step 2: Write the failing nf-test**

```groovy
// tests/local/ac_calc.nf.test
nextflow_process {
    name "Test ac_calc"
    script "modules/local/ac_calc.nf"
    process "AC_CALC"

    test("computes read depth and allelic counts (noagg path)") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT', normal_sample_id: 'TEST_NORMAL', tumor_sample_id: 'TEST_TUMOR', agg_bins: false, desired_binsize: null],
                    file("\${projectDir}/tests/fixtures/cnv_core/TEST_TUMOR_counts.txt"),
                    files("\${projectDir}/tests/fixtures/cnv_core/query/*.phasedGT.chr*.1.txt"),
                    files("\${projectDir}/tests/fixtures/cnv_core/query/*.AD.chr*.txt"),
                    files("\${projectDir}/tests/fixtures/cnv_core/query/*_HetInfo.txt")
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).exists()
            assert path(process.out[0].get(0)[2]).exists()
        }
    }
}
```

- [ ] **Step 3: Run it, verify it fails**

Run: `nf-test test tests/local/ac_calc.nf.test`
Expected: FAIL

- [ ] **Step 4: Implement the module**

`container` uses `params.cnv_r_container` (Task 2), not a hardcoded image string:

```groovy
// modules/local/ac_calc.nf
process AC_CALC {
    tag "${meta.participant_id}"
    container params.cnv_r_container
    memory '4 GB'

    input:
    tuple val(meta), path(tumor_counts), path(phased_files), path(ad_files), path(hetinfo_files)

    output:
    tuple val(meta), path("*_read_depth.txt"), path("*_allelic_counts.txt"), path("*_plots.pdf")

    script:
    def agg_bins = meta.agg_bins ? 'TRUE' : 'FALSE'
    def binsize_arg = meta.desired_binsize ? "${meta.desired_binsize}" : ''
    """
    set -euo pipefail
    Rscript ${projectDir}/bin/ac_calc.R \\
        ${meta.participant_id} ${tumor_counts} ${agg_bins} ${binsize_arg} \\
        ${meta.normal_sample_id} ${meta.tumor_sample_id} \\
        ${phased_files} ${ad_files} ${hetinfo_files}
    """
}
```

- [ ] **Step 5: Run it, verify it passes**

Run: `nf-test test tests/local/ac_calc.nf.test`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add bin/ac_calc.R modules/local/ac_calc.nf tests/local/ac_calc.nf.test
git commit -m "Add ac_calc module (consolidates WDL's ac_calc_agg/ac_calc_noagg)"
```

---

## Task 13: `tangent_xy` module

**Files:**
- Create: `bin/tangent.R` (copy of `reference/r-scripts/20250409_wdlTangent.R`)
- Create: `modules/local/tangent_xy.nf`
- Create: `tests/local/tangent_xy.nf.test`

**Interfaces:**
- Consumes: `tuple val(meta), path(tumor_counts), path(pon), path(ref_clinical),
  path(cohort_clinical)` — same uniform per-participant `meta` as everywhere else in this
  subworkflow; must carry `sample_sex`, `scale_method`, `fix_outliers`, `trailingN`, `n_dim`
- Produces: `tuple val(meta), path("posttangent_log2.rds")` — consumed by Task 14
  (`segment`)

- [ ] **Step 1: Vendor the R script**

```bash
cp reference/r-scripts/20250409_wdlTangent.R bin/tangent.R
```

- [ ] **Step 2: Write the failing nf-test**

```groovy
// tests/local/tangent_xy.nf.test
nextflow_process {
    name "Test tangent_xy"
    script "modules/local/tangent_xy.nf"
    process "TANGENT_XY"

    test("normalizes tumor depth against the external PoN") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT', normal_sample_id: 'TEST_NORMAL', tumor_sample_id: 'TEST_TUMOR', sample_sex: 'female', scale_method: 'median', fix_outliers: false, trailingN: 5, n_dim: 10],
                    file("\${projectDir}/tests/fixtures/cnv_core/TEST_TUMOR_counts.txt"),
                    file(params.pon),
                    file(params.ref_clinical),
                    file(params.cohort_clinical)
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).exists()
        }
    }
}
```

Add the tangent-specific reference params to `nextflow.config`:

```groovy
// append to nextflow.config
params.pon             = 's3://btc-refdata/iconicc/BTC_pon.rds'
params.ref_clinical     = 's3://btc-refdata/iconicc/reference_clinical.txt'
params.cohort_clinical  = 's3://btc-refdata/iconicc/btc_clinical.txt'
```

- [ ] **Step 3: Run it, verify it fails**

Run: `nf-test test tests/local/tangent_xy.nf.test`
Expected: FAIL

- [ ] **Step 4: Implement the module**

```groovy
// modules/local/tangent_xy.nf
process TANGENT_XY {
    container params.cnv_r_container
    memory '10 GB'

    input:
    tuple val(meta), path(tumor_counts), path(pon), path(ref_clinical), path(cohort_clinical)

    output:
    tuple val(meta), path("posttangent_log2.rds")

    script:
    """
    set -euo pipefail
    Rscript ${projectDir}/bin/tangent.R \\
        --pon ${pon} --tumor_cts ${tumor_counts} --ref_cli ${ref_clinical} --clin ${cohort_clinical} \\
        --scale-mthd ${meta.scale_method} --fix_outliers ${meta.fix_outliers} \\
        --trailingN ${meta.trailingN} --sample_sex ${meta.sample_sex} --n_dim ${meta.n_dim}
    """
}
```

- [ ] **Step 5: Run it, verify it passes**

Run: `nf-test test tests/local/tangent_xy.nf.test`
Expected: PASS

- [ ] **Step 6: Commit**

```bash
git add bin/tangent.R modules/local/tangent_xy.nf tests/local/tangent_xy.nf.test nextflow.config
git commit -m "Add tangent_xy module"
```

---

## Task 14: `segment` module + `CNV_CORE` subworkflow wiring

**Files:**
- Create: `bin/segment.R` (copy of `reference/r-scripts/20250318_wdlSegment.R`)
- Create: `modules/local/segment.nf`
- Create: `subworkflows/local/cnv_core/main.nf`
- Create: `tests/local/segment.nf.test`
- Create: `tests/subworkflows/cnv_core.nf.test`
- Modify: `main.nf` (call `CNV_CORE` against the staged fixtures)

**Interfaces:**
- `SEGMENT` consumes: `tuple val(meta), path(allelic_counts), path(tangent_normalized)`
- `SEGMENT` produces: `tuple val(meta), path("*.seg.txt"), path("*_processed_counts.txt"),
  path("*_segmentPlots.pdf"), path("*_segment_cts.txt")` — these four plus `ac_calc`'s
  `read_depth` are the pipeline's final outputs (SPEC.md §6): `read_depth`, `seg_file`,
  `ac_cts`, `segment_plots`, `segment_cts`.
- `CNV_CORE` subworkflow — **take:** `tuple val(meta), path(joint_vcf),
  path(joint_vcf_index)` plus the per-sample CRAM/CRAI tuples and the static reference
  params listed throughout Tasks 3-13; **emit:** `read_depth`, `seg_file`, `ac_cts`,
  `segment_plots`, `segment_cts` — this is the exact channel set Phase 3's `main.nf` wires
  up once the calling+joint-genotyping subworkflow (Phase 2) exists.

- [ ] **Step 1: Vendor the R script**

```bash
cp reference/r-scripts/20250318_wdlSegment.R bin/segment.R
```

- [ ] **Step 2: Write the failing nf-test for `SEGMENT`**

```groovy
// tests/local/segment.nf.test
nextflow_process {
    name "Test segment"
    script "modules/local/segment.nf"
    process "SEGMENT"

    test("segments allelic counts against tangent-normalized depth") {
        when {
            process {
                """
                input[0] = [
                    [participant_id: 'TEST_PT', sample_sex: 'female', smooth_region: 10, outlier_SD: 4, smooth_SD: 2, alpha: 0.01, min_width: 2, undo_split: 'none', undo_prune: 0.05, undo_SD: 2, recovery_thresh: 20, use_het: true, miss_thresh: 1.0],
                    file("\${projectDir}/tests/fixtures/cnv_core/TEST_PT_allelic_counts.txt"),
                    file("\${projectDir}/tests/fixtures/cnv_core/posttangent_log2.rds")
                ]
                """
            }
        }
        then {
            assert process.success
            assert path(process.out[0].get(0)[1]).exists()
        }
    }
}
```

- [ ] **Step 3: Run it, verify it fails**

Run: `nf-test test tests/local/segment.nf.test`
Expected: FAIL

- [ ] **Step 4: Implement `SEGMENT`**

```groovy
// modules/local/segment.nf
process SEGMENT {
    tag "${meta.participant_id}"
    container params.cnv_r_container
    memory '10 GB'

    input:
    tuple val(meta), path(allelic_counts), path(tangent_normalized)

    output:
    tuple val(meta), path("*.seg.txt"), path("*_processed_counts.txt"), path("*_segmentPlots.pdf"), path("*_segment_cts.txt")

    script:
    """
    set -euo pipefail
    Rscript ${projectDir}/bin/segment.R \\
        --al_cts ${allelic_counts} --tangent_cts ${tangent_normalized} \\
        --sample_sex ${meta.sample_sex} --participant_id ${meta.participant_id} \\
        --smth_reg ${meta.smooth_region} --outlier_SD ${meta.outlier_SD} --smooth_SD ${meta.smooth_SD} \\
        --alpha ${meta.alpha} --min_width ${meta.min_width} --undo_split ${meta.undo_split} \\
        --undo_prune ${meta.undo_prune} --undo_SD ${meta.undo_SD} \\
        --recovery_thresh ${meta.recovery_thresh} --use_hetinfo ${meta.use_het} --miss_thresh ${meta.miss_thresh}
    """
}
```

- [ ] **Step 5: Run it, verify it passes**

Run: `nf-test test tests/local/segment.nf.test`
Expected: PASS

- [ ] **Step 6: Wire the `CNV_CORE` subworkflow**

```groovy
// subworkflows/local/cnv_core/main.nf
include { SUBSET_PAIR_VCF     } from '../../../modules/local/subset_pair_vcf'
include { PREPROCESS_INTERVALS } from '../../../modules/local/preprocess_intervals'
include { COLLECT_READ_COUNTS as COLLECT_READ_COUNTS_TUMOR  } from '../../../modules/local/collect_read_counts'
include { COLLECT_READ_COUNTS as COLLECT_READ_COUNTS_NORMAL } from '../../../modules/local/collect_read_counts'
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
                     // desired_binsize, and every tangent_xy/segment tuning field (Global
                     // Constraints / Task 5's note). It must never be reshaped downstream
                     // (only ever read from), so every .join() on `meta` below matches.
    tumor_cram_ch    // [ meta, cram, crai ] -- same meta as joint_vcf_ch
    normal_cram_ch   // [ meta, cram, crai ] -- same meta as joint_vcf_ch
    ref_fasta_ch     // [ fasta, fai, dict ]
    ref_dir_1000g
    genetic_map
    pon
    ref_clinical
    cohort_clinical

    main:
    SUBSET_PAIR_VCF(joint_vcf_ch)

    PREPROCESS_INTERVALS(ref_fasta_ch)

    COLLECT_READ_COUNTS_TUMOR(
        tumor_cram_ch.combine(PREPROCESS_INTERVALS.out).combine(ref_fasta_ch)
            .map { meta, cram, crai, intervals, fasta, fai, dict -> [ meta, meta.tumor_sample_id, cram, crai, intervals, fasta, fai, dict ] }
    )
    COLLECT_READ_COUNTS_NORMAL(
        normal_cram_ch.combine(PREPROCESS_INTERVALS.out).combine(ref_fasta_ch)
            .map { meta, cram, crai, intervals, fasta, fai, dict -> [ meta, meta.normal_sample_id, cram, crai, intervals, fasta, fai, dict ] }
    )

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
        COLLECT_READ_COUNTS_TUMOR.out
            .join(QUERY_HET_AD.out)
            .map { meta, tumor_counts, hetinfo, ad, phased -> [ meta, tumor_counts, phased, ad, hetinfo ] }
    )

    TANGENT_XY(
        COLLECT_READ_COUNTS_TUMOR.out
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

- [ ] **Step 7: Write the subworkflow-level nf-test against the real fixtures**

```groovy
// tests/subworkflows/cnv_core.nf.test
nextflow_workflow {
    name "Test CNV_CORE end-to-end against a known-good WDL run"
    script "subworkflows/local/cnv_core/main.nf"
    workflow "CNV_CORE"

    test("produces all five final outputs from a real joint VCF + CRAM pair") {
        when {
            workflow {
                """
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
                    file("\${projectDir}/tests/fixtures/cnv_core/tumor.recal.cram"),
                    file("\${projectDir}/tests/fixtures/cnv_core/tumor.recal.cram.crai")
                ])
                input[2] = Channel.of([
                    meta,
                    file("\${projectDir}/tests/fixtures/cnv_core/normal.recal.cram"),
                    file("\${projectDir}/tests/fixtures/cnv_core/normal.recal.cram.crai")
                ])
                input[3] = Channel.of([ file(params.ref_fasta), file(params.ref_fasta_fai), file(params.ref_fasta_dict) ])
                input[4] = Channel.value(file(params.ref_dir_1000g))
                input[5] = Channel.value(file(params.genetic_map))
                input[6] = Channel.value(file(params.pon))
                input[7] = Channel.value(file(params.ref_clinical))
                input[8] = Channel.value(file(params.cohort_clinical))
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

Load `params.fixture_participant_id`/`fixture_normal_id`/`fixture_tumor_id`/
`fixture_sample_sex` from `tests/fixtures/cnv_core/sample_info.env` (staged in Task 3) via
an `nf-test.config` `params-file`, or inline them directly once the real values are known
from Task 3 Step 1 — either way, this is the one test in the plan that needs the real
dataset's identifiers, not synthetic `TEST_*` placeholders, since it's checking actual
scientific output, not just wiring.

- [ ] **Step 8: Run it, verify it fails, then passes**

Run: `nf-test test tests/subworkflows/cnv_core.nf.test`
First run expected: FAIL (subworkflow file doesn't exist until Step 6 lands — if executing
strictly TDD-style, write this test before Step 6's implementation). After Step 6:
Expected: PASS.

- [ ] **Step 9: Wire `main.nf` to run `CNV_CORE` standalone against the fixtures**

```groovy
// main.nf
nextflow.enable.dsl = 2

include { CNV_CORE } from './subworkflows/local/cnv_core/main'

workflow {
    // One uniform per-participant meta, reused unchanged everywhere CNV_CORE needs it --
    // see the note in subworkflows/local/cnv_core/main.nf's `take:` block for why this
    // has to stay a single identical map across every channel for .join() to work.
    meta = [
        participant_id: params.participant_id,
        normal_sample_id: params.normal_sample_id,
        tumor_sample_id: params.tumor_sample_id,
        sample_sex: params.sample_sex,
        scale_method: 'median', fix_outliers: false, trailingN: 5, n_dim: 10,
        smooth_region: 10, outlier_SD: 4, smooth_SD: 2, alpha: 0.01, min_width: 2,
        undo_split: 'none', undo_prune: 0.05, undo_SD: 2, recovery_thresh: 20,
        use_het: true, miss_thresh: 1.0, agg_bins: false, desired_binsize: null
    ]

    joint_vcf_ch   = Channel.of([ meta, file(params.joint_vcf), file(params.joint_vcf_index) ])
    tumor_cram_ch  = Channel.of([ meta, file(params.tumor_cram), file(params.tumor_crai) ])
    normal_cram_ch = Channel.of([ meta, file(params.normal_cram), file(params.normal_crai) ])
    ref_fasta_ch   = Channel.of([ file(params.ref_fasta), file(params.ref_fasta_fai), file(params.ref_fasta_dict) ])

    CNV_CORE(
        joint_vcf_ch, tumor_cram_ch, normal_cram_ch, ref_fasta_ch,
        Channel.value(file(params.ref_dir_1000g)),
        Channel.value(file(params.genetic_map)),
        Channel.value(file(params.pon)),
        Channel.value(file(params.ref_clinical)),
        Channel.value(file(params.cohort_clinical))
    )
}
```

(Phase 3 will replace `joint_vcf_ch`/`tumor_cram_ch`/`normal_cram_ch` here with the real
calling+joint-genotyping subworkflow's output — this version exists purely to let the CNV
half run and be validated standalone, per SPEC.md §10 point 1.)

- [ ] **Step 10: Run the real pipeline end-to-end against the fixture, by hand**

```bash
nextflow run main.nf \
    --participant_id      "$(grep PARTICIPANT_ID tests/fixtures/cnv_core/sample_info.env | cut -d= -f2)" \
    --normal_sample_id    "$(grep NORMAL_SAMPLE_ID tests/fixtures/cnv_core/sample_info.env | cut -d= -f2)" \
    --tumor_sample_id     "$(grep TUMOR_SAMPLE_ID tests/fixtures/cnv_core/sample_info.env | cut -d= -f2)" \
    --sample_sex          "$(grep SAMPLE_SEX tests/fixtures/cnv_core/sample_info.env | cut -d= -f2)" \
    --joint_vcf tests/fixtures/cnv_core/joint.vcf.gz --joint_vcf_index tests/fixtures/cnv_core/joint.vcf.gz.tbi \
    --tumor_cram tests/fixtures/cnv_core/tumor.recal.cram --tumor_crai tests/fixtures/cnv_core/tumor.recal.cram.crai \
    --normal_cram tests/fixtures/cnv_core/normal.recal.cram --normal_crai tests/fixtures/cnv_core/normal.recal.cram.crai
```

Compare `results/segment/*.seg.txt` and `results/ac_calc/*_read_depth.txt` against the
original WDL run's equivalent outputs for the same patient (SPEC.md §10 point 1's actual
purpose) — general concordance, not necessarily bit-identical given Task 2's fresh R
container may resolve slightly different transitive package versions than
`wchukwu/r-docker:latest` did. Flag any material divergence to the dev before proceeding to
Phase 2.

- [ ] **Step 11: Commit**

```bash
git add bin/segment.R modules/local/segment.nf subworkflows/local/cnv_core/main.nf \
        tests/local/segment.nf.test tests/subworkflows/cnv_core.nf.test main.nf
git commit -m "Wire CNV_CORE subworkflow end-to-end; validate against a known-good WDL run"
```

---

## What Phase 1 does NOT cover (by design)

- Any GPU/Parabricks calling, joint-genotyping, or VQSR — Phase 2.
- Cirro wiring (`.cirro/process-form.json` etc.) — Phase 3, after both halves exist.
- The final end-to-end run comparing this port's full output against the original WDL
  pipeline for the same patient (SPEC.md §10 point 4) — needs Phase 2 wired in first.

Once this plan's tasks are done and Step 10's manual concordance check looks reasonable,
return to brainstorming for Phase 2 (calling + joint-genotyping subworkflow, SPEC.md §5) as
its own plan.
