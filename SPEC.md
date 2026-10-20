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

**Verified against a real Cirro run**: the joint VCF's actual `#CHROM` sample columns are
not the bare `normal_sample_id`/`tumor_sample_id` — `parabricks_bqsr` (this pipeline's
sole alignment-entry source) embeds each CRAM's read group as `SM:${patient}_${sample}`,
which GenomicsDBImport/GenotypeGVCFs carry straight through. `SUBSET_PAIR_VCF` subsets by
`${patient}_${normal_sample_id}`/`${patient}_${tumor_sample_id}` to match. This is a fixed
property of `parabricks_bqsr`'s own read-group convention, not something either pipeline
is expected to change — if that convention ever changes, `SUBSET_PAIR_VCF` must change
with it.

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
| Tranche value lists | explicit custom lists | none passed, GATK tool built-in default applies | **Now verified** (§12.3): use ICONICC's exact literal lists, not sarek's tool-default behavior |

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
PT01,PT01_normal,female,0,s3://.../PT01_normal.recal.cram,s3://.../PT01_normal.recal.cram.crai,,
PT01,PT01_tumor_A,female,1,s3://.../PT01_tumor_A.recal.cram,s3://.../PT01_tumor_A.recal.cram.crai,s3://.../PT01_tumor_A_vs_PT01_normal.vcf.gz,s3://.../PT01_tumor_A_vs_PT01_normal.vcf.gz.tbi
PT01,PT01_tumor_B,female,1,s3://.../PT01_tumor_B.recal.cram,s3://.../PT01_tumor_B.recal.cram.crai,s3://.../PT01_tumor_B_vs_PT01_normal.vcf.gz,s3://.../PT01_tumor_B_vs_PT01_normal.vcf.gz.tbi
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

Schema validation (**resolved**): hand-rolled Groovy inside the `SAMPLESHEET_TO_PAIRS`
subworkflow (`subworkflows/local/samplesheet_to_pairs/main.nf`), not a JSON-Schema file.
The sheet is read with plain `.splitCsv(header: true)` and each row passes through a
`validateSamplesheetRow()` helper that hard-errors, naming the offending row, on a missing
`patient`/`sample`/`cram`/`crai`, a `sex` outside `male`/`female`, a `status` outside
`0`/`1`, or a tumor row missing `joint_vcf`/`joint_vcf_index`. The pair-level checks
(>1 normal per patient, a tumor with no matching normal, tumor/normal sex mismatch, and a
samplesheet yielding zero pairs at all) live in the same file. Sarek's
`assets/schema_input.json` + `nextflow_schema.json` + nf-schema route was considered and
deliberately not taken: this org's other Nextflow/Cirro pipelines are hand-rolled rather
than nf-core-template based, so matching that convention avoids introducing an nf-schema
dependency and a second place where the schema is defined. Neither
`assets/schema_input.json` nor `nextflow_schema.json` exists in this repo, by design.

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
4. Join each pair's `meta` to its **tumor's** already-computed read-counts (by
   `(patient, sample)`) — one `[pair_meta, tumor_counts]` tuple per pair. Only the tumor's
   counts are joined: `CNV_CORE` consumes no normal counts (see the `take:` note below), so
   the sample-level counts channel is filtered to `status == 1` before this join. The
   normal's counts are still computed in step 2 — once per normal, role-uniformly — and
   then simply go unconsumed; that cost is deliberate and documented in
   `subworkflows/local/sample_read_counts/main.nf`.
   The pair's `joint_vcf`/`joint_vcf_index` (from the tumor row) are *not* routed through
   this join — they go straight from the samplesheet into `CNV_CORE`, because the phasing
   chain (`SUBSET_PAIR_VCF` → … → `EAGLE_PHASING`) needs no counts data and must not be
   gated behind that pair's hours-long `COLLECT_READ_COUNTS`. The two sides rejoin inside
   `CNV_CORE` on `meta`, which matches by map equality regardless of which channel built it.
5. Call `CNV_CORE` with a **multi-item** channel (one item per pair) instead of Phase 1's
   single-item `Channel.of([meta, ...])`. Nothing inside `CNV_CORE` from `SUBSET_PAIR_VCF`
   onward needs to change — every module in it already keys on `meta` and Nextflow
   processes handle N-items-in-flight natively; Phase 1's per-task reviews already confirmed
   this wiring is correct in general, not just for exactly one item.

**`CNV_CORE`'s `take:` interface changes**: `tumor_cram_ch`/`normal_cram_ch` (raw CRAMs) are
replaced by a single `tumor_counts_ch` (pre-computed, already-deduplicated
`COLLECT_READ_COUNTS` output, `[meta, counts_file]`). There is deliberately **no**
`normal_counts_ch`: nothing from `SUBSET_PAIR_VCF` onward ever read the normal's counts —
already true in Phase 1, and WDL-faithful, since the WDL likewise runs `CollectReadCounts`
on the normal without feeding it into the allelic-CN step (§4) — so now that
`COLLECT_READ_COUNTS` runs role-uniformly over every sample there is no role-specific
plumbing left that would justify threading a dead channel through the interface. The final
`take:` is therefore `joint_vcf_ch, tumor_counts_ch, ref_dir_1000g, genetic_map, pon,
ref_clinical, cohort_clinical`. `ref_fasta_ch` is dropped from
`CNV_CORE`'s `take:` entirely — `PREPROCESS_INTERVALS` moves out with `COLLECT_READ_COUNTS`
into the new sample-level subworkflow, since both are sample/reference-scoped, not
pair-scoped. This is a breaking change to `CNV_CORE`'s external interface from Phase 1 —
expected and acceptable, since Phase 1 was explicitly a standalone validation harness
(§10 point 1), not a frozen contract.

**Single-pair case**: a two-row samplesheet (one normal, one tumor for one patient)
produces exactly one item through the cross in step 3 — identical behavior to Phase 1,
just expressed through a samplesheet instead of ten CLI flags.

## 12. Parabricks calling + joint genotyping (Phase 2)

Automates what §11.2 called transitional: the samplesheet's `joint_vcf`/`joint_vcf_index`
columns, supplied by hand until now, get computed by the pipeline instead. This closes
the gap between §5's design (already resolved: VQSR tuning values, the `ExcessHet` hard
filter, GATK 4.2.6.1, the Parabricks/GATK4 HaplotypeCaller flag mapping) and actual code.

### 12.1 Architecture

Two new subworkflows, mirroring the sample-scoped/pair-scoped split §11.3 already
established for `SAMPLE_READ_COUNTS`/`CNV_CORE`:

```
SAMPLESHEET_TO_PAIRS  (unchanged, minus joint_vcf/joint_vcf_index columns — see §12.4)
        │
        ├─→ samples ──→ SAMPLE_READ_COUNTS        (existing, unchanged)
        │
        └─→ samples ──→ SAMPLE_GERMLINE_CALLING   (new, sample-scoped, §12.2)
                              │
                              ▼
                         [meta(patient,sample,...), gvcf, gvcf_tbi]
                              │
        pairs ────────────────┤  join normal's GVCF + tumor's GVCF by (patient, sample)
                              ▼
                       JOINT_GENOTYPING            (new, pair-scoped, §12.3)
                              │
                              ▼
                       [pair_meta, vcf, tbi]  ──→ replaces joint_vcf_ch in main.nf
                                                        │
                                                        ▼
                                                    CNV_CORE (unchanged)
```

**Dedup shape** (confirmed with the dev): calling is the dedup point, joint genotyping is
not. `SAMPLE_GERMLINE_CALLING` gets dedup for free the same way `SAMPLE_READ_COUNTS` does
— `SAMPLESHEET_TO_PAIRS.out.samples` already has one row per unique sample, so a normal
referenced by 3 tumor rows is called once. `JOINT_GENOTYPING` stays pairwise per patient,
matching both the WDL's existing behavior and §5's "feed it just the one patient's two
samples per pipeline run" — a shared normal is joint-genotyped separately against each of
its tumors, never cohort-batched. `JOINT_GENOTYPING`'s input is therefore a join of
`SAMPLE_GERMLINE_CALLING`'s output against each pair's `normal_sample_id` and
`tumor_sample_id`, keyed on `[patient, sample]` — the same join pattern `main.nf` already
uses for `tumor_counts_ch` (§11.3 step 4), done for both roles and combined into one
`[normal_gvcf, normal_tbi, tumor_gvcf, tumor_tbi]` tuple per pair. A missing GVCF on either
side is a hard error (`failOnMismatch: true`), matching every other meta-keyed join in this
codebase.

### 12.2 `SAMPLE_GERMLINE_CALLING` (new, sample-scoped)

Runs `PARABRICKS_HAPLOTYPECALLER` (vendored from `nf-core/modules`, verbatim — container
`nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1`, sarek's current default; no org precedent
to match, see decision log below) once per unique `(patient, sample)`, `ext.args` including
`--gvcf` per §5.2's flag mapping. Feeds it the sample's `recal.cram` directly — Parabricks
accepts CRAM as `--in-bam` per NVIDIA's own docs, so no `CramToBamTask` equivalent is
needed, exactly as §5 already specified. Restricts calling to the WDL's own
`unpadded_intervals_file` (see §12.5) rather than the true unrestricted whole genome, so
assembly-gap regions aren't wastefully processed — but does **not** split that file into
shards (see decision log: no interval scatter for v1). Output: `[meta(patient, sample,
sex, status), gvcf, gvcf_tbi]` — one GVCF per unique sample, `.tbi`-indexed via a
local `INDEX_GVCF` module (bcftools-based, `modules/local/index_gvcf.nf`) alongside
the Parabricks call, since `PARABRICKS_HAPLOTYPECALLER` itself emits `.g.vcf.gz`
without an index.

### 12.3 `JOINT_GENOTYPING` (new, pair-scoped)

Composes vendored nf-core modules directly — `GATK4_GENOMICSDBIMPORT`,
`GATK4_GENOTYPEGVCFS`, `GATK4_VARIANTFILTRATION`, `GATK4_VARIANTRECALIBRATOR`,
`GATK4_APPLYVQSR` — rather than importing sarek's `BAM_JOINT_CALLING_GERMLINE_GATK`
subworkflow as a black box (decision log below explains why). Chain, per pair:

```
GenomicsDBImport(normal_gvcf, tumor_gvcf, calling_intervals)
  → GenotypeGVCFs                                     (one call, whole calling_intervals region)
  → VariantFiltration  --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet
  → MakeSitesOnlyVcf                                   (training input only, not the final VCF)
  → VariantRecalibrator(SNP, on sites-only)  ┐
  → VariantRecalibrator(INDEL, on sites-only)┘
  → ApplyVQSR(INDEL)  → ApplyVQSR(SNP, chained on INDEL's output)
  → [pair_meta, vcf, tbi]
```

Chain order verified directly against `ICONICC_paired.3.wdl`'s `ApplyRecalibration` task
(`JointGeno_Tasks.1.wdl:586-634`): **INDEL runs first, SNP second, chained on INDEL's
output** — the reverse of sarek's own `BAM_JOINT_CALLING_GERMLINE_GATK` comment ("Step 1:
ApplyVQSR to SNPs... Step 2: ...run ApplyVQSR_INDEL"), which this design initially copied
without independently re-checking against the WDL. Caught while pulling exact `ext.args`
values for the implementation plan — corrected here since each pass only tags its own
variant type's FILTER field and passes the other type through untouched, so the final
tags are equivalent either way, but WDL-order match is followed for exact parity with
what has actually been running in production.

Verified directly against `ICONICC_paired.3.wdl` (not inferred from sarek's subworkflow,
which doesn't do this): `ApplyRecalibration`'s actual input — and therefore the *final*
`output_vcf_file`/`output_vcf_index_file` that feeds the CNV half (§6) — is the
`ExcessHet`-filtered VCF with VQSR applied on top, not the raw `GenotypeGVCFs` output. The
`ExcessHet` filter is a soft `VariantFiltration` tag (matching the WDL's own
`HardFilterAndMakeSitesOnlyVcf` task — the FILTER column is annotated, not the record
dropped), and `MakeSitesOnlyVcf`'s output is used only as `VariantRecalibrator`'s training
input, per the WDL. `ext.args`/tranche overrides for `VariantRecalibrator`/`ApplyVQSR`
(SNP/INDEL) are exactly §5.1's already-resolved values — no new decision needed there.

No `BCFTOOLS_SORT`/`MergeVcfs` step: those exist in sarek's subworkflow purely to
reassemble interval-scattered shards, and per the no-scatter decision below,
`GenotypeGVCFs` runs once per pair over the whole `calling_intervals` region and already
emits one complete, sorted VCF.

**Exact `ext.args`/resource values** — read directly from `ICONICC_paired.3.wdl`'s own
input-block defaults (lines 56-59) and confirmed unoverridden in the production Cirro
config (`cirro/process-input.json` has no `*_tranche_values`/`*_annotation_values`/
`*_filter_level` keys at all, so the WDL's embedded defaults are what production actually
runs), resolving §5.1's "Tranche value lists... Not verified identical" row:

```
VARIANTRECALIBRATOR_SNP:
  -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5
  -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0
  -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP
  -mode SNP --max-gaussians 6 --trust-all-polymorphic
  resource labels (prior values are literal WDL constants, GATK expects a trailing ".0"):
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0   hapmap_3.3.hg38.vcf.gz
    --resource:omni,known=false,training=true,truth=true,prior=12.0     1000G_omni2.5.hg38.vcf.gz
    --resource:1000G,known=false,training=true,truth=false,prior=10.0  1000G_phase1.snps.high_confidence.hg38.vcf.gz
    --resource:dbsnp,known=true,training=false,truth=false,prior=7.0   Homo_sapiens_assembly38.dbsnp138.vcf

VARIANTRECALIBRATOR_INDEL:
  -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0
  -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0
  -tranche 91.0 -tranche 90.0
  -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP
  -mode INDEL --max-gaussians 4 --trust-all-polymorphic
  resource labels:
    --resource:mills,known=false,training=true,truth=true,prior=12.0    Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0    Homo_sapiens_assembly38.dbsnp138.vcf

APPLYVQSR_INDEL (first):  --truth-sensitivity-filter-level 99.0 -mode INDEL
APPLYVQSR_SNP   (second): --truth-sensitivity-filter-level 99.7 -mode SNP
```

`use_allele_specific_annotations` stays `false` (§5.1, already resolved) — the
`--use-allele-specific-annotations` flag is never passed, on either `VariantRecalibrator`
call. `max-gaussians` (6 for SNP, 4 for INDEL) and `--trust-all-polymorphic` are WDL
defaults carried over unconditionally — no decision needed, just literal reproduction.
Resource file basenames match §12.5's params, following sarek's own convention
(`conf/igenomes.config`) of hardcoding the label string with the reference file's literal
basename rather than deriving it from a `path` value at runtime.

### 12.4 Samplesheet schema change

`joint_vcf`/`joint_vcf_index` columns are removed entirely from the schema §11.2 defined —
exactly the outcome §11.2 already flagged as expected once this phase landed.
`validateSamplesheetRow()` drops the tumor-row check for those two columns. This is a
breaking schema change (any hand-written samplesheet from Phase 1.5 needs its two trailing
columns removed); acceptable since §11.2 called the schema unstable pending this phase and
there is no external consumer yet.

### 12.5 New reference inputs

New `params.*`, all resolved from the production Cirro config
(`reference/btc-wdl-pipelines/cirro/process-input.json`) rather than guessed:

```
calling_intervals                          s3://gatk-test-data/intervals/hg38.even.handcurated.20k.intervals
dbsnp_vcf(+idx)                             s3://broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf(.idx)
hapmap_resource_vcf(+idx)                   s3://broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz(.tbi)
omni_resource_vcf(+idx)                     s3://broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz(.tbi)
one_thousand_genomes_resource_vcf(+idx)     s3://broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz(.tbi)
mills_resource_vcf(+idx)                    s3://broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz(.tbi)
axiomPoly_resource_vcf(+idx)                s3://broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz(.tbi)
```

`calling_intervals` is the WDL's own `unpadded_intervals_file` — the un-scattered,
gap-excluding base interval list `SplitIntervalList` shards in the original WDL. Used
unsplit for both `PARABRICKS_HAPLOTYPECALLER` and `GATK4_GENOTYPEGVCFS` (§12.2/§12.3),
per the no-scatter decision below. The WDL's `eval_interval_list` and
`scattered_calling_intervals_list` are not needed: the former only feeds
`CollectVariantCallingMetrics`, already deferred (§8); the latter is a pre-split shard
manifest with no purpose once scatter is dropped.

### 12.6 Decision log

**Interval scatter — dropped for v1.** The WDL scatters `HaplotypeCaller` and
`GenotypeGVCFs`/VQSR across ~50 interval shards, a CPU-parallelism trick from the original
Broad WGS-joint-calling template. Parabricks HaplotypeCaller already processes the whole
genome in one GPU job — the reason for that scatter disappears on the calling side — and
this pipeline runs one patient's two samples per invocation, never a cohort, so
joint-genotyping's scatter buys little either. **Resolved (confirmed with the dev)**: run
each stage as a single whole-genome job per pair. This is why §12.3's chain has no
`BCFTOOLS_SORT`/merge step and why `JOINT_GENOTYPING` needs no `num_intervals` bookkeeping.

**Joint-genotyping integration — vendor sarek's modules, write our own orchestration, not
sarek's subworkflow.** Three approaches were weighed: (A) fork sarek's
`BAM_JOINT_CALLING_GERMLINE_GATK` subworkflow file and patch in `ExcessHet` — rejected,
because that subworkflow's wiring exists mainly to handle interval-scatter bookkeeping
(`meta.num_intervals`, `groupTuple(by: intervals)`) this port just decided to drop, so
forking it means immediately deleting half of what was forked. (B) hand-write the whole
GATK4 chain from raw processes, ignoring sarek entirely — rejected, throws away sarek's
already-tested module wiring (correct args, containers, output globs) for no reason, and
contradicts §2's explicit direction to reuse sarek rather than re-porting all 19
`JointGeno_Tasks.1.wdl` tasks by hand. **Resolved (C, confirmed with the dev)**: vendor the
individual nf-core modules unmodified into `modules/nf-core/`, and write a new, simple
local subworkflow (`JOINT_GENOTYPING`) that calls them in WDL order with our own no-scatter
wiring — reusing exactly the part of sarek worth reusing (the GATK tool wrappers) while
writing only the orchestration ourselves, which is now cheap with scatter out of scope.

**Parabricks version — sarek's current default.** No known org precedent to match (asked
directly); `nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1` is what the vendored sarek
modules currently pin, and matching the modules we're vendoring is the safer default over
guessing at an org-wide standard that may not exist yet.

**GATK4 version for the joint-genotyping modules — pinned to 4.2.6.1 for now.** The
vendored nf-core modules default to whatever GATK4 version sarek currently ships, likely
newer than the WDL's 4.2.6.1. **Resolved (confirmed with the dev)**: override each
vendored GATK4 module's container to 4.2.6.1 via `nextflow.config`, matching the WDL
exactly and keeping one GATK version across the whole pipeline (consistent with §7's
CNV-half resolution). The dev flagged this as provisional — a candidate for an org-wide
GATK version bump later, revisited outside this port's scope.

### 12.7 Testing & accepted gaps

Same posture as Phase 1/1.5: no real fixture data, write code only. Adds a third accepted
gap alongside the two already open (§10 point 1's fixture data, the unpushed R container):
no real Parabricks/NGC execution has been validated, since neither the dev's nor this
session's environment has GPU/NGC access. `PARABRICKS_HAPLOTYPECALLER`'s vendored module
already ships a `stub {}` block (verified by reading it directly), so `SAMPLE_GERMLINE_CALLING`'s
wiring gets tested via nf-test's `-stub-run` rather than a real GPU run. The GATK4 modules
in `JOINT_GENOTYPING` use public containers, so those may run for real in nf-test if
Docker/network access allows during implementation — verify empirically rather than assume
(this pipeline's own precedent: the CRLF/shebang bug and the `--clin` flag both turned out
differently than assumed until actually checked against the real thing).

## 13. Cirro wiring (Phase 3)

Deploys this pipeline as an actual Cirro process, matching how `sarek_align`/
`sarek_call_variants` are already deployed (both live in `CirroBio/Cirro-pipelines`,
inspected directly for real precedent rather than designed from scratch). `.cirro/` files
live in this same repo, not a separate one — confirmed with the dev during Phase 2's
wrap-up.

### 13.1 Upstream chain and parent process

Real pipeline chain, now confirmed end-to-end: `sarek_align` → `parabricks_bqsr`
(`break-through-cancer/parabricks-bqsr`, a separate repo, built independently of this
one) → this pipeline. `parabricks_bqsr` applies BQSR genome-wide (no interval
restriction, so unmapped reads still get recalibrated) and optionally quantizes quality
scores to replicate GATK `ApplyBQSR`'s `--static-quantized-quals`; its real, confirmed
output shape is `preprocessing/recalibrated/<sample>/<sample>.recal.cram(+.crai)` — a real
completed run was inspected directly via the Cirro Python SDK (project
`af54fa61-9e6c-4b87-9922-3c2c073aeaaa`, dataset `a89df334-e188-4b55-a50a-8d4ac3d6e08e`,
process `dlt-parabricks-bqsr`, status `COMPLETED`) to confirm this rather than assume it
from the spec alone.

**Resolved**: ICONICC-NF declares `parabricks_bqsr`'s alignment-entry process (the one
that accepts an already-aligned BAM/CRAM and emits the recalibrated+quantized CRAM above
— not its separate fastq-entry process, which handles alignment itself upstream of that)
as its sole Cirro parent process. `sarek_align` is not declared as a parent directly —
every real input to this pipeline has gone through BQSR+optional quantization first, so
there is no supported path that skips `parabricks_bqsr`.

### 13.2 Samplesheet assembly (`preprocess.py`)

Real precedent: `nf-core/sarek_call_variants`'s own `preprocess.py`
(`CirroBio/Cirro-pipelines`, `nf-core/sarek_call_variants/3/`) already solves this exact
problem for sarek's own `patient,sex,status,sample` samplesheet shape — read directly
rather than designed blind. Two Cirro-native building blocks it draws on:

- **`ds.files`**: the selected input dataset's own file listing, passed into
  `preprocess.py` as `ds.files["file"]`. Unlike the `sarek_call_variants` precedent, this
  pipeline does **not** rely on a committed `fileMappingRules` entry to do the indexed
  BAM/CRAM pairing — `fileMappingRules` lives in `process-definition.json`, which §13.5
  already establishes this repo deliberately does not commit. Instead, `preprocess.py`
  does its own file discovery: `alignment_pairs()` matches each path against a regex
  (`RECAL_FILE`, matching `preprocessing/recalibrated/<sample>/<file>.(bam|cram)$`) and
  keeps only samples that also have a same-named `.bai`/`.crai` alongside, dropping any
  sample missing its index.
- **`ds.samplesheet`**: Cirro's own per-sample metadata table, which the user fills in
  directly in Cirro's UI when registering samples, passed into `preprocess.py` as a plain
  DataFrame. `preprocess.py` calls `.set_index("sample")` on it itself, which only makes
  sense if `sample` arrives as an ordinary **column**, not an already-set index — so that
  is the shape this code assumes, not an index Cirro hands over pre-set. Already uses
  `patient`/`sex`/`status` as literal column names in the existing `sarek_call_variants`
  precedent — the exact three columns `SAMPLESHEET_TO_PAIRS` needs beyond `cram`/`crai`.
  **Verified against a real failed Cirro launch**: `sample` does arrive as an ordinary
  column as assumed, but the other column names are not reliable. The project's own
  declared common-metadata schema does define `patient` — but a sample whose metadata was
  carried over from a fastq-ingestion samplesheet instead exposes `patient_id`, and the
  first real launch failed on exactly that (`KeyError: 'patient'`). Separately, `sex` has
  been observed recorded as karyotype (`XX`/`XY`), not literal `male`/`female`.
  `preprocess.py` now accepts both column names (preferring `patient` when both are
  present) and maps `XX`/`XY` to `female`/`male` before validating.

`preprocess.py`'s job: join the two into `SAMPLESHEET_TO_PAIRS`'s exact CSV shape
(`patient,sample,sex,status,cram,crai`), write it out, and point `--input` at it via
`process-input.json`.

**Resolved (validation duplication — confirmed with the dev)**: `preprocess.py` also
re-implements `SAMPLESHEET_TO_PAIRS`'s own pairing hard-errors (duplicate normal, orphan
tumor, sex mismatch) so a bad samplesheet fails before Cirro launches compute, not after.
This is a deliberate two-places-same-rules tradeoff, not an oversight — accepted
explicitly in exchange for faster, cheaper failure. If the Groovy-side rules in
`SAMPLESHEET_TO_PAIRS` ever change, `preprocess.py`'s copy needs the same change; there is
no shared source of truth between the two.

### 13.3 GPU compute-config

Real, current precedent from `parabricks-bqsr`'s own `.cirro/alignment/process-compute.config`
(`break-through-cancer/parabricks-bqsr`, `fastq-entry` branch — read directly, not
inferred from the older `sarek_align`/§2 precedent, which turned out to describe a
different, instance-family-based mechanism this org has since moved away from):

```groovy
withName: 'PARABRICKS_HAPLOTYPECALLER' {
    accelerator = 1
    queue       = env('PW_ONDEMAND_JOB_QUEUE')
}
```

GPU sizing is a **queue selection**, not a compute-config instance-family choice — the
queue name is an environment variable Cirro injects, pointing at an AWS Batch job queue
tied to a GPU-capable Compute Environment the admin provisions separately, outside this
repo. `accelerator = 1` matches a single-GPU request, taken from the alignment-entry
config above (the directly analogous case — a single already-aligned sample through one
GPU-requiring Parabricks step). Worth noting rather than smoothing over: the sibling
repo's *other* config, for its fastq-entry process, sets no `accelerator` directive at all
for its own GPU-requiring processes (`PARABRICKS_FQ2BAM`/`_PART`/`MARKDUP`/`APPLYBQSR`) —
only `queue`. Read directly, not assumed; the two real configs genuinely disagree on
whether `accelerator` needs setting at all, and this spec does not resolve why. Start from
`accelerator = 1` here since it matches the more directly analogous precedent, verify it
actually schedules correctly once real Cirro access exists to test against, and treat the
value (and whether it's needed at all) as revisable — the dev already expects to tweak
this. No instance-family or vCPU-ceiling details belong in this repository's
`nextflow.config`/`process-compute.config` at all —
that is 100% the admin's queue-provisioning responsibility.

### 13.4 Form scope

Thin, per §2/§3's already-resolved "keep it minimal" directive: no GATK/VQSR internals
exposed (those are pinned in code via `nextflow.config` `ext.args`, already true since
Phase 2). The form only needs the input-dataset selector (driving §13.2's file discovery)
and `--outdir`-equivalent Cirro plumbing — no new user-facing parameters beyond what
`usesSampleSheet: true` and Cirro's own dataset-selection UI already provide for free.

### 13.5 Repo layout additions

**Correction during plan-writing**: this section originally listed a committed
`process-definition.json`, modeled on the public nf-core forks (`sarek_align`/
`sarek_call_variants`), which do commit one. Checking `parabricks_bqsr`'s actual repo
tree directly (`break-through-cancer/parabricks-bqsr`, both the `fastq`/`alignment`
`.cirro/` subdirectories, full tree listing, not just the two directories first
inspected) shows zero `process-definition.json` anywhere — this org's real, current
convention for its own custom (non-nf-core-fork) pipelines is to register the process
directly in Cirro's admin UI (id, parent-process declaration, executor version, code
repository pointer), not as a version-controlled file. Only the files below are
committed, matching that real precedent exactly:

```
.cirro/
  process-form.json         -- thin form, §13.4
  process-compute.config    -- GPU queue wiring, §13.3
  process-input.json        -- --input mapped to preprocess.py's assembled CSV
  process-output.json
  preprocess.py             -- samplesheet assembly + validation, §13.2
  test_preprocess.py
```

The parent-process linkage itself (adding ICONICC-NF's registered process id to
`parabricks_bqsr`'s own `childProcessIds`) is consequently a manual, UI-side step in the
*other* repo's Cirro registration, not a file either repo commits — flagged explicitly
as a follow-up, not silently assumed to happen automatically.

### Decision log

- **Parent process**: `parabricks_bqsr`'s alignment-entry process only (not `sarek_align`
  directly, not `parabricks_bqsr`'s fastq-entry process) — every real input has gone
  through BQSR+optional quantization first, so no other parent is a supported path.
- **Validation duplication**: accepted deliberately — `preprocess.py` re-implements
  `SAMPLESHEET_TO_PAIRS`'s pairing checks for faster failure, at the cost of two
  same-rules places to maintain.
- **GPU wiring**: queue-based (`env('PW_ONDEMAND_JOB_QUEUE')` + `accelerator = 1`),
  matching `parabricks_bqsr`'s own real, current pattern — not the older
  instance-family-based approach §2 originally cited, which this org has since moved
  past. Accelerator count may be revisited later; not treated as final here.
