# ICONICC → Nextflow Port: Technical Spec

Status: design/planning notes, no code written yet. Written to hand off into a fresh
coding session. Everything below was verified against the actual WDL source, the actual
R scripts, and the actual nf-core/sarek source (not just docs) unless explicitly flagged
as unverified.

## 1. Background

ICONICC ("HapCNA" — the actual Cromwell workflow name) is a tumor/normal paired WGS
pipeline doing germline joint-genotyping followed by a bespoke haplotype-based allelic
copy-number (CNV) analysis. Current implementation:

- Repo: `break-through-cancer/btc-wdl-pipelines`, branch `iconicc`
- Main workflow: `ICONICC_paired.3.wdl` (workflow name `HapCNA`)
- Task libraries: `ICON_Tasks.2.wdl` (16 tasks), `JointGeno_Tasks.1.wdl` (19 tasks)
- Cirro wiring: `cirro/{process-form.json, process-input.json, preprocess.py,
  process-compute.config, process-output.json}` on the same branch
- Runs on Cromwell/AWS via Cirro, in project `BTC-GBM-Production`
  (`e4f45b6a-e2b0-4ddb-a1d5-b276c414af05`)

Trigger for this spec: a real run (`ab0cbb5c-9961-4345-882b-77e663808970`, "ICONICC test
on P3 S1.L2 recalibrated crams") failed with
`java.nio.file.ProviderNotFoundException: Provider "s3" not found` inside
`HaplotypeCaller`/`CollectReadCounts` — root cause was Cromwell's AWS backend silently
skipping localization for files that live in the execution's own scratch bucket
(`CramToBamTask`'s output), then GATK trying to read the resulting `s3://...` path
directly, which vanilla GATK 4.2.0.0 can't do (no S3 NIO provider, only GCS). That
investigation is what prompted evaluating a Nextflow port instead of patching the WDL.

## 2. Architecture Decision: Two Pipelines, Not One

**Pipeline A — unchanged, already exists.** `sarek_align`
(`process-nf-core-sarek-align-3-2` in Cirro's `nf-core/sarek_align/3.2/`) runs nf-core/sarek
through alignment/BQSR. Stays exactly as-is; not part of this port. Its
`preprocessing/recalibrated/<sample>/<sample>.recal.cram` (+`.crai`) output is this port's
input. Reason to keep it separate: the recalibrated CRAM has independent reuse value
outside ICONICC (other analyses want it too), so it earns its own Cirro process.

**Pipeline B — new, this port.** Bundles three logical stages into one Nextflow pipeline,
because none of their intermediate artifacts (a per-sample GVCF, a per-pair joint VCF)
have standalone reuse value outside this one analysis — there's no reason to force a user
to manually chain three separate Cirro dataset launches for what is really one workflow:

1. Per-sample germline calling (Parabricks HaplotypeCaller, GPU) — replaces
   `CramToBamTask` + `HaplotypeCaller` + `MergeVCFs` + `GenerateSampleMapFile` from
   `ICON_Tasks.2.wdl`.
2. Joint genotyping + VQSR (reuse sarek's `BAM_JOINT_CALLING_GERMLINE_GATK` subworkflow
   directly, imported as a subworkflow — NOT via sarek's top-level `tools`/`joint_germline`
   CLI params, see §5) — replaces all 19 tasks in `JointGeno_Tasks.1.wdl`.
3. The bespoke CNV/allelic-analysis core — hand-ported from `ICON_Tasks.2.wdl`'s
   remaining 12 tasks, see §4.

Even though these three stages are bundled behind one Cirro launch, **write them as two
internal Nextflow subworkflows** (calling+joint-genotyping, and CNV) rather than one flat
script — this keeps the CNV logic independently testable (feed it a known-good VCF
without re-running GPU calling) and lets `-resume` skip the expensive front end if only
the CNV/R tail needs a fix.

Mixed GPU (stage 1) + CPU (stages 2-3) processes in one Nextflow pipeline is normal and
already has precedent in this org: see `btc-wdl-pipelines`'s sibling repo
`Cirro-pipelines` branch `feature/sarek-parabricks-gpu`, commit `7da9376`, which added
`accelerator` overrides to `sarek_align`/`sarek_call_variants`'s
`process-compute.config` to force GPU requests unconditionally (Cirro can't pass
`-profile gpu` conditionally per run).

## 3. Cirro-facing scope for v1: keep it minimal

Confirmed with the dev: the *current* ICONICC Cirro form only exposes `input_type`,
`participant_id`, `tumor_sample_id`/`normal_sample_id`, the BAM/BAI file pickers, and
`sample_sex`. None of the GATK/VQSR internals below are user-adjustable today, and most
optional WDL params are never touched. **Do not build a richer configuration UI for v1.**
Reproduce the current fixed behavior behind an equally thin form; pin everything else in
code (compute-config `ext.args` overrides), not exposed params. Contamination filtering
(`-contamination`, defaults to `0` in the WDL, never exposed via Cirro, no known issues)
is explicitly out of scope — don't try to replicate it.

## 4. The CNV Core — what actually gets hand-ported

From `ICON_Tasks.2.wdl`, these 12 tasks are ICONICC-unique and have no sarek/nf-core
equivalent. This is the real porting work:

`PreprocessIntervals`, `CollectReadCounts` (×2, tumor+normal), `Filter`, `iter_isec`,
`hetvarFilter`, `eagle_phasing`, `index`, `query`, `ac_calc_agg`/`ac_calc_noagg`,
`tangent_XY`, `segment`.

### 4.1 The actual science (verified against the real R scripts, not just the WDL)

Pulled and read the live scripts from `s3://btc-refdata/iconicc/`:
`reference-genome-AC_calc_20231228.R` (`ac_rscript`), `20250409_wdlTangent.R`
(`tangent_rscript`), `20250318_wdlSegment.R` (`segment_rscript`).

This is a genuine haplotype-based allelic copy-number method (matches the workflow's own
name, "HapCNA"), not boilerplate GATK usage:

- **Genotype/phasing comes from the *normal*; allele depth comes from the *tumor*,
  recombined at the same coordinates.** In `ac_calc.R`:
  ```r
  read_data(phased_GT_directory, normal_prefix, ".phasedGT.%s.1.txt", ...)   # NORMAL
  read_data(allelic_depths_directory, tumor_prefix,  ".AD.%s.txt", ...)      # TUMOR
  ```
  `eagle_phasing` phases germline-heterozygous sites (against the normal) using the
  1000G reference panel; `ac_calc.R` then measures allelic imbalance in the *tumor's*
  reads at those same phased positions.
- **Depth-based copy ratio (`tangent_XY`) normalizes the tumor against an external
  panel of normals (PoN)**, not this patient's own matched normal — `tangent.R` only
  takes `--pon` and `--tumor_cts`; the patient's own normal never enters the depth-ratio
  calculation at all. Two different data sources for two different signals, deliberately
  kept separate.

### 4.2 Known bug to fix during the port, not blindly replicate

`hetvarFilter`'s bcftools expression uses **positional** VCF-column indexing:
`GT[1]="het" & FORMAT/AD[1:0]>=5 & FORMAT/AD[1:1]>=5` (task comment: "assumes the first
sample is germline"). This only works because `GenerateSampleMapFile` happens to write
`[normal, tumor]` in that order today, and GATK's GenomicsDB/joint-genotyping sample
order is not a guaranteed contract. `ac_calc.R` itself is robust (matches samples **by
name** via `tumor_prefix`/`normal_prefix`), so only this one upstream selection step is
fragile.

**Fix in the port**: immediately after the joint-calling stage hands off its VCF, run
`bcftools view -s normal_sample_id,tumor_sample_id` — this both subsets to this patient's
pair (see §5, sarek's joint-calling is cohort-wide, not pairwise) and pins the column
order `hetvarFilter` needs, in one step. Do this before any other CNV-half processing.

## 5. Joint genotyping — reuse sarek's subworkflow directly, not the top-level pipeline

`nf-core/sarek`'s `subworkflows/local/bam_joint_calling_germline_gatk/main.nf`
(`BAM_JOINT_CALLING_GERMLINE_GATK`) is structurally identical to
`JointGeno_Tasks.1.wdl`: `GenomicsDBImport → GenotypeGVCFs → BCFTOOLS_SORT/merge →
VariantRecalibrator(SNP+INDEL) → ApplyVQSR(SNP then INDEL, chained)`. Import it directly
as a subworkflow.

**Do not go through sarek's packaged `sarek_call_variants` Cirro process / top-level
`tools`+`joint_germline` params for this**, for two reasons:
1. Its `tools=parabricks_haplotypecaller` path has **no `joint_germline` branch wired up
   at all** in `BAM_VARIANT_CALLING_GERMLINE_ALL` (verified in source) — selecting
   Parabricks there silently gives you two independent per-sample VCFs, never the joint
   VCF this pipeline needs.
2. Sarek's joint-germline mode genotypes the **whole cohort together**, not
   pairwise-per-patient like ICONICC does today. Feed it just the one patient's two
   samples per pipeline run to match today's execution granularity (don't batch across
   patients).

**Correct wiring for this port**: call `PARABRICKS_HAPLOTYPECALLER` (from
`nf-core/modules`, `ext.args` including `--gvcf`) directly per sample, feeding it
`sarek_align`'s `recal.cram` — Parabricks' `--in-bam` accepts CRAM directly per NVIDIA's
own docs ("takes an aligned BAM/CRAM file"), so no `CramToBamTask` equivalent is needed
at all. Feed both resulting GVCFs into `BAM_JOINT_CALLING_GERMLINE_GATK` imported
directly.

Reference resource files needed for VQSR — already known-good S3 paths, reuse as-is
(from the current `cirro/process-input.json`):
```
ref_fasta / ref_fasta_index / ref_dict     s3://broad-references/hg38/v0/Homo_sapiens_assembly38.*
dbsnp_vcf(+idx)                             s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf(.idx)
hapmap_resource_vcf(+idx)                   s3://broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz(.tbi)
omni_resource_vcf(+idx)                     s3://broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz(.tbi)
one_thousand_genomes_resource_vcf(+idx)     s3://broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz(.tbi)
mills_resource_vcf(+idx)                    s3://broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz(.tbi)
axiomPoly_resource_vcf(+idx)                s3://broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz(.tbi)
```

### 5.1 GATK/VQSR fixed-value overrides needed (per §3, pin these in code, don't expose)

ICONICC's actual defaults (`ICONICC_paired.3.wdl` input block) vs sarek's hardcoded
`conf/modules/joint_germline.config` values:

| Param | ICONICC default | Sarek hardcoded | Action |
|---|---|---|---|
| SNP truth-sensitivity filter level | `99.7` | `99.9` | Override via `ext.args` on `GATK4_APPLYVQSR_SNP` |
| INDEL truth-sensitivity filter level | `99.0` | `99.9` | Override via `ext.args` on `GATK4_APPLYVQSR_INDEL` |
| SNP VQSR annotations | `QD, MQRankSum, ReadPosRankSum, FS, MQ, SOR, DP` (7) | `QD, MQ, MQRankSum, ReadPosRankSum, FS, SOR` (missing `DP`) | Override `ext.args` on `VARIANTRECALIBRATOR_SNP` to add `-an DP` |
| INDEL VQSR annotations | `FS, ReadPosRankSum, MQRankSum, QD, SOR, DP` | same 6, different order | No action — sets match |
| `ExcessHet` hard-filter (54.69) | present (`HardFilterAndMakeSitesOnlyVcf`) | step doesn't exist in sarek | No equivalent at all — decide whether to add a standalone `VariantFiltration` step or drop it (needs a call: is 54.69 a deliberately-tuned value or an inherited Broad-template default nobody revisited?) |
| `use_gnarly_genotyper` | defaults `false` | N/A, not implemented in sarek | Moot at ICONICC's actual default |
| `use_allele_specific_annotations` | defaults `false` | N/A, sarek's VQSR args are fixed non-AS | Moot at ICONICC's actual default |
| Tranche value lists | explicit custom lists | none passed, GATK tool built-in default applies | Not verified identical — diff directly if exact reproducibility matters |

Override mechanism: `process { withName: 'VARIANTRECALIBRATOR_SNP' { ext.args = '...' } }`
in this pipeline's own compute config — same override pattern already proven for the
Parabricks GPU work (commit `7da9376` above).

**Resolved (confirmed with the dev)**: `99.7`/`99.0`/the `DP` annotation/`54.69` were
deliberately tuned for this assay, not inherited defaults. Replicate all of them exactly —
override `ext.args` on `VARIANTRECALIBRATOR_SNP`/`_INDEL` and `GATK4_APPLYVQSR_SNP`/`_INDEL`
as specified above, and add the `ExcessHet` hard filter (see below) rather than dropping it.

**`ExcessHet` hard filter — resolved**: add it as a standalone step. Sarek's
`BAM_JOINT_CALLING_GERMLINE_GATK` has no equivalent, so insert a `GATK4_VARIANTFILTRATION`
process (filter expression `ExcessHet > 54.69`, matching `HardFilterAndMakeSitesOnlyVcf`)
between `GenotypeGVCFs`/merge and `VariantRecalibrator`, in the same place the WDL runs it.

### 5.2 HaplotypeCaller: Parabricks vs GATK4 — verified against the full NVIDIA option reference

Full options pulled directly from
`https://docs.nvidia.com/clara/parabricks/tool-reference/tools/haplotypecaller`
(rendered page, not the doc-search summaries, which were unreliable/inconsistent across
several attempts — always re-verify against the live page if this matters later).

| ICONICC flag | Parabricks equivalent | Verdict |
|---|---|---|
| `-ERC GVCF` (`make_gvcf`) | `--gvcf` | Matches |
| `-GQB 10 -GQB 20 ... -GQB 90` | `-GQB` / `--gvcf-gq-bands` — "Exclusive upper bounds for reference confidence GQ bands... in increasing order" | Matches exactly, same flag name |
| `-G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation` | `-G` / `--annotation-group` — "supported annotation groups are StandardAnnotation, StandardHCAnnotation, and AS_StandardAnnotation" | **Matches exactly**, same three group names |
| `-bamout` (`make_bamout`, defaults `false`) | `--htvc-bam-output` | Matches (different name, same function); moot anyway since default is off |
| `-contamination` | Not present anywhere in the full ~40-flag reference | No equivalent — explicitly out of scope per §3 |

Net: 4 of 5 match directly; the one gap is the one already agreed to not worry about.

**Known accuracy caveat** (NVIDIA's own "Source of Mismatches" docs, quoted in full):
Parabricks is *concordant*, not bit-identical, with GATK4, for four specific reasons —
separate RNG streams (parallelism), Java `Math.log10()` vs C++ `cmath` `log10` giving
different haplotype-penalty scores (can select different haplotypes, not just round
differently), AVX/SIMD floating-point differences vs GPU kernels, and (in Parabricks'
favor) deterministic hash-iteration order vs GATK's Java `HashMap` non-determinism. Worth
knowing before claiming bit-for-bit reproducibility with historical WDL runs.

## 6. Data-flow contract (the actual interface boundary)

Traced directly from `ICONICC_paired.3.wdl` lines 655-660 — this is the *entire*
dependency surface between joint-genotyping and the CNV half:
```
File output_vcf_file       = select_first([FinalGatherVcf.output_vcf, ApplyRecalibration.recalibrated_vcf])
File output_vcf_index_file = select_first([FinalGatherVcf.output_vcf_index, ApplyRecalibration.recalibrated_vcf_index])
```
Consumed purely via `bcftools filter`/`bcftools isec` downstream.

Full input set for the new pipeline (Pipeline B, §2):
```
from sarek_align:   <normal_sample>.recal.cram + .crai
                     <tumor_sample>.recal.cram + .crai
generated in-pipeline: joint-genotyped, VQSR-filtered VCF + .tbi (§5)
static, unchanged:  tumor_sample_id, normal_sample_id, participant_id, sample_sex,
                     ref_dir1000G      s3://btc-refdata/iconicc/reference-genome-1000G.tar
                     geneticMap        s3://btc-refdata/iconicc/reference-genome-genetic_map_hg38_withX.txt.gz
                     haplotype_database s3://btc-refdata/iconicc/Homo_sapiens_assembly38.haplotype_database.txt
                     pon               s3://btc-refdata/iconicc/BTC_pon.rds
                     ref_clinical      s3://btc-refdata/iconicc/reference_clinical.txt
                     cohort_clinical   s3://btc-refdata/iconicc/btc_clinical.txt
                     ac_rscript        s3://btc-refdata/iconicc/reference-genome-AC_calc_20231228.R
                     tangent_rscript   s3://btc-refdata/iconicc/20250409_wdlTangent.R
                     segment_rscript   s3://btc-refdata/iconicc/20250318_wdlSegment.R
```

Final outputs (from `ICONICC_paired.3.wdl`'s own `output{}` block, carry over unchanged):
`read_depth`, `seg_file`, `ac_cts`, `segment_plots`, `segment_cts`.

## 7. Container/version notes — resolved

- Existing WDL mixes GATK **4.2.0.0** (`CramToBamTask` only, via `gitc_docker`) and
  **4.2.6.1** (everywhere else). `CramToBamTask` itself is moot (replaced by Parabricks,
  §5.2). **Resolved**: standardize the CNV-half's own `gatk`-based tasks
  (`PreprocessIntervals`, `CollectReadCounts`) on **GATK 4.2.6.1**, matching the majority
  of existing tasks — one pinned version across the whole CNV half.
- CNV R tasks currently run on `wchukwu/r-docker:latest` — an **unpinned tag**.
  **Resolved**: build a fresh, versioned Dockerfile for the port rather than reverse-pinning
  the existing image. Identify the actual R packages the three scripts (`ac_calc.R`,
  `tangent.R`, `segment.R`) import, pin explicit versions (CRAN/Bioconductor snapshot or
  renv lockfile), and build+tag it as part of this repo rather than depending on an
  external, unversioned image.

## 8. Things explicitly deferred / not needed for v1

- Fingerprinting (`SelectFingerprintSiteVariants`, `CrossCheckFingerprint`,
  `GatherPicardMetrics`, `GetFingerprintingIntervalIndices`) and Picard variant-calling
  metrics (`CollectVariantCallingMetrics`, `GatherVariantCallingMetrics`) have **no sarek
  equivalent at all**. Nf-core/modules has ready-made wrappers for two related pieces if
  ever needed later — `gatk4/selectvariants` and `picard/collectvariantcallingmetrics`
  (both near-drop-in, confirmed against their actual `main.nf` source) — but there's no
  `accumulatevariantcallingmetrics`/gather equivalent; that one's a short custom process
  if ever needed. Not required unless sample-swap detection or these specific QC metrics
  turn out to be load-bearing for validation sign-off.
- `GnarlyGenotyper` and the scattered/model-based SNP VQSR path
  (`SNPsVariantRecalibratorCreateModel`, `GatherTranches`) — moot, ICONICC defaults both
  off (`use_gnarly_genotyper = false`).
- A full task-by-task mapping TSV (all 23 generic tasks from both WDL files → sarek
  module/subworkflow, with notes) already exists at
  `/Users/DLT/Downloads/iconicc_sarek_task_mapping.tsv` — reference it directly rather
  than re-deriving.

## 9. Rough sizing

Full 1:1 hand-port of both WDL files: large — ~34 tasks including several
mutually-exclusive branch-pairs (Gnarly vs GenotypeGVCFs, scattered-model vs classic SNP
VQSR, scattered vs solo fingerprinting) each needing separate implementation and
validation in Nextflow's channel-branching idioms.

This architecture (two pipelines, reuse sarek wherever possible): meaningfully smaller —
the ~19-task `JointGeno_Tasks.1.wdl` and 4 generic `ICON_Tasks.2.wdl` tasks are replaced
by importing existing, already-tested sarek modules/subworkflows. Actual new authoring is
concentrated in: the `PARABRICKS_HAPLOTYPECALLER` wiring (small), the VQSR `ext.args`
overrides from §5.1 (small, config-only), and the 12-task CNV core from §4 (the real
work — genuine science, not boilerplate, needs careful validation against existing WDL
run outputs before trusting it).

## 10. Suggested validation order

1. Stand up the CNV-half subworkflow (§4) first, in isolation, fed with an **existing**
   completed WDL run's `output_vcf_file`/CRAMs — de-risks the CNV port without also
   debugging a new GPU-calling front end at the same time.
2. Stand up calling+joint-genotyping (§5) separately, validate its output VCF against a
   known-good existing WDL joint VCF (sample presence/order via `bcftools query -l`,
   general concordance).
3. Wire the two together, confirm the `bcftools view -s normal,tumor` reorder step (§4.2)
   produces the same column order the CNV half expects.
4. Full end-to-end run against a real patient pair, compare final CNV segment/allelic-count
   outputs against the existing WDL pipeline's outputs for the same patient.

## 11. Samplesheet-based input & multi-pair support (post-Phase-1 revision)

Phase 1 shipped `main.nf` taking one tumor/normal pair via ten flat `--param` flags
(§ commit history, not this doc — see `docs/superpowers/plans/2026-09-17-cnv-core-subworkflow.md`).
This section replaces that with a samplesheet, adds support for multiple pairs per run, and
adds real compute dedup when the same normal is reused across several tumors (e.g. one
patient's PBMC germline vs. several tumor samples) — confirmed with the dev as a real
requirement, not just a metadata convenience.

### 11.1 Why sample-level, not pair-level

Modeled directly on nf-core/sarek's own samplesheet convention
(`assets/schema_input.json`, `assets/samplesheet.csv` in the vendored sarek source under
`reference/sarek/` — not part of this repo, pulled for reference during Phase 1). Sarek's
sheet is **one row per sample** (patient/sample/sex/status + file), not one row per pair —
pairing is derived downstream via a `.cross()` on `patient` between samples with
`status: 0` (normal) and `status: 1` (tumor) (`workflows/sarek.nf:360-375` in the vendored
copy). For a patient with one normal and N tumors, this produces exactly N pairs, each
naturally reusing the same normal — the comment in sarek's own source literally reads
`// Use cross to combine normal with all tumor samples, i.e. multi tumor samples from
recurrences`. This is the same scenario the dev asked about (one PBMC germline vs. several
tumors), so the pattern transfers directly rather than needing to be invented.

Pairing this way also isn't a new convention for this project specifically: the *original*
WDL pipeline's own Cirro `preprocess.py` (`reference/btc-wdl-pipelines/cirro/preprocess.py`,
`yield_single_inputs()`) already derived `participant_id = f"{tumor_sample_id}_vs_
{normal_sample_id}"` when auto-pairing multiple tumor/normal files from one dataset
selection. Sarek's `meta.id = "${tumor.sample}_vs_${normal.sample}"` is the identical
convention. Adopting it here converges on something that already existed in two independent
places, rather than introducing a third naming scheme.

**Validation**: mirror sarek's own guard — a patient with more than one `status: 0` row is
a hard error (`workflows/sarek.nf` calls this out explicitly: "has more than one sample with
normal status"). Supporting multiple normals per patient was not requested and would make
"which normal pairs with which tumor" ambiguous; out of scope. Also carry over the original
WDL Cirro `preprocess.py`'s own check (`yield_single_inputs()`): a tumor whose `sex` doesn't
match its patient's normal is a hard error, not a silent skip as the WDL version did — a
mismatch indicates a sample mix-up, and failing loudly is more appropriate for a
samplesheet a human curated directly than for a preprocessing script auto-pairing files.

### 11.2 Samplesheet schema

CSV, one row per sample:

```
patient,sample,sex,status,cram,crai,joint_vcf,joint_vcf_index
PT01,PT01_normal,XX,0,s3://.../PT01_normal.recal.cram,s3://.../PT01_normal.recal.cram.crai,,
PT01,PT01_tumor_A,XX,1,s3://.../PT01_tumor_A.recal.cram,s3://.../PT01_tumor_A.recal.cram.crai,s3://.../PT01_tumor_A_vs_PT01_normal.vcf.gz,s3://.../PT01_tumor_A_vs_PT01_normal.vcf.gz.tbi
PT01,PT01_tumor_B,XX,1,s3://.../PT01_tumor_B.recal.cram,s3://.../PT01_tumor_B.recal.cram.crai,s3://.../PT01_tumor_B_vs_PT01_normal.vcf.gz,s3://.../PT01_tumor_B_vs_PT01_normal.vcf.gz.tbi
```

- `patient`, `sample`, `sex` (`male`/`female`, matching the WDL's existing convention —
  not sarek's `XX`/`XY`, to avoid a needless translation layer downstream in the R scripts,
  which already expect `male`/`female`), `status` (`0`=normal, `1`=tumor): required on every
  row, sarek-identical semantics.
- `cram`/`crai`: required on every row — the `sarek_align` output this pipeline has always
  consumed (SPEC.md §1, §5).
- `joint_vcf`/`joint_vcf_index`: **required on tumor rows, ignored on normal rows.**
  Transitional — joint-genotyping is pairwise per patient (§5: "feed it just the one
  patient's two samples per pipeline run", not cohort-wide), so a patient with 3 tumors
  needs 3 separate joint-genotyping runs today, one per tumor row, each against the shared
  normal. Once Phase 2 (§5) is built, these two columns are computed by the pipeline instead
  of supplied, and disappear from the schema — the schema is expected to change again then,
  not a stable long-term contract.

Schema validation: a `assets/schema_input.json` following sarek's own JSON-Schema
convention (referenced from `nextflow_schema.json` via nf-schema, matching this org's other
Nextflow/Cirro pipelines where applicable) — enumerated as its own plan task, not written
inline in this spec.

### 11.3 Pairing and dedup architecture

`main.nf` (or a thin subworkflow it delegates to) does, in order:

1. Parse the samplesheet into a per-sample channel: `[meta(patient, sample, sex, status),
   cram, crai]`.
2. Feed **every** sample (both normal and tumor) through a new sample-level subworkflow —
   `subworkflows/local/sample_read_counts/` — running `PREPROCESS_INTERVALS` once globally
   (unchanged from Phase 1) and `COLLECT_READ_COUNTS` once per unique `(patient, sample)`.
   This is the actual dedup: a normal referenced by 3 tumor rows only gets its read counts
   computed once, because it only appears once in the sample-level channel.
3. Branch the sample-level channel by `status`, `.cross()` on `patient` (sarek's exact
   pattern, §11.1) to build one pair-level `meta` per tumor: `{participant_id:
   "${tumor.sample}_vs_${normal.sample}", normal_sample_id, tumor_sample_id, sample_sex,
   ...tuning fields}` — the same uniform per-participant `meta` shape `CNV_CORE` already
   expects (Phase 1 plan, Task 5's note), constructed once per pair instead of once per
   pipeline invocation.
4. Join each pair's `meta` to its normal's and tumor's already-computed read-counts (by
   `(patient, sample)`) and to its `joint_vcf`/`joint_vcf_index` (from the tumor row) — one
   `CNV_CORE`-shaped input tuple per pair.
5. Call `CNV_CORE` with a **multi-item** channel (one item per pair) instead of Phase 1's
   single-item `Channel.of([meta, ...])`. Nothing inside `CNV_CORE` from `SUBSET_PAIR_VCF`
   onward needs to change — every module in it already keys on `meta` and Nextflow
   processes handle N-items-in-flight natively; Phase 1's per-task reviews already confirmed
   this wiring is correct in general, not just for exactly one item.

**`CNV_CORE`'s `take:` interface changes**: `tumor_cram_ch`/`normal_cram_ch` (raw CRAMs) are
replaced by `tumor_counts_ch`/`normal_counts_ch` (pre-computed, already-deduplicated
`COLLECT_READ_COUNTS` outputs, `[meta, counts_file]`). `ref_fasta_ch` is dropped from
`CNV_CORE`'s `take:` entirely — `PREPROCESS_INTERVALS` moves out with `COLLECT_READ_COUNTS`
into the new sample-level subworkflow, since both are sample/reference-scoped, not
pair-scoped. This is a breaking change to `CNV_CORE`'s external interface from Phase 1 —
expected and acceptable, since Phase 1 was explicitly a standalone validation harness
(§10 point 1), not a frozen contract.

**Single-pair case**: a two-row samplesheet (one normal, one tumor for one patient)
produces exactly one item through the cross in step 3 — identical behavior to Phase 1,
just expressed through a samplesheet instead of ten CLI flags.
