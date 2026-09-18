# CNV-core validation fixtures

This directory is expected to hold real genomic data staged from a completed
ICONICC WDL run, per Task 3 Step 2 of the `cnv-core-subworkflow` plan
(`.superpowers/sdd/2026-09-17-cnv-core-subworkflow/task-3-brief.md`). These
files are gitignored (`tests/fixtures/` in `.gitignore`) and are not, and
must never be, committed to this repository.

**Status: not staged.** The controller asked the dev for a Cirro dataset ID /
S3 paths to a known-good run and was told to defer fixture staging and proceed
with code-only work. Every `nf-test` that reads a fixture will fail until these
files are staged by hand — that is all of `tests/local/` except
`preprocess_intervals.nf.test`, plus `tests/subworkflows/cnv_core.nf.test`.

## Expected files

Two kinds of fixture appear below. **Staged** files are copied out of a
completed WDL run. **Derived** files are produced by running an earlier module
of this pipeline standalone and copying its output here — they let each module
be tested in isolation without re-running its upstream. Derived fixtures are
real setup dependencies, not placeholders.

Filenames containing a sample or participant ID use the `TEST_PT` /
`TEST_TUMOR` / `TEST_NORMAL` placeholders the module tests pass as `meta`.
Where a test globs (`*_het.vcf.gz`), the expected count is 23 files — chr1-22
plus chrX — matching the per-chromosome loop in the corresponding module.

### Run inputs (staged)

| File | What it is |
| --- | --- |
| `joint.vcf.gz` | The joint-genotyped, VQSR-filtered VCF (`output_vcf_file`) from a completed ICONICC WDL run's joint-genotyping subworkflow. |
| `joint.vcf.gz.tbi` | Tabix index for `joint.vcf.gz` (`output_vcf_index_file`). |
| `normal.recal.cram` | The normal sample's recalibrated CRAM from `sarek_align`, for the same run. |
| `normal.recal.cram.crai` | CRAM index for `normal.recal.cram`. |
| `tumor.recal.cram` | The tumor sample's recalibrated CRAM from `sarek_align`, for the same run. |
| `tumor.recal.cram.crai` | CRAM index for `tumor.recal.cram`. |
| `sample_info.env` | Shell-sourceable file recording the real `PARTICIPANT_ID`, `NORMAL_SAMPLE_ID`, `TUMOR_SAMPLE_ID`, and `SAMPLE_SEX` for the staged run, so later tests use real IDs instead of placeholders. Its four values must also be copied into `params.fixture_participant_id` / `fixture_normal_id` / `fixture_tumor_id` / `fixture_sample_sex` in `nextflow.config`, or passed on the nf-test command line — `tests/subworkflows/cnv_core.nf.test` asserts they are non-null before it runs. |

Used by `tests/subworkflows/cnv_core.nf.test` (all of the above),
`tests/local/subset_pair_vcf.nf.test` (`joint.vcf.gz*`), and
`tests/local/collect_read_counts.nf.test` (`tumor.recal.cram*`).

### Interval and variant-subset stage (derived)

| File | Produced by | Needed by |
| --- | --- | --- |
| `wgs.preprocessed.interval_list` | `PREPROCESS_INTERVALS` (GATK `PreprocessIntervals`) against the reference fasta | `collect_read_counts.nf.test` |
| `pair_subset.vcf.gz` | `SUBSET_PAIR_VCF` (normal+tumor sample subset of `joint.vcf.gz`) | `filter_variants.nf.test` |
| `pair_subset.vcf.gz.tbi` | Tabix index for `pair_subset.vcf.gz` | `filter_variants.nf.test` |
| `HC_var.vcf.gz` | `FILTER_VARIANTS` (high-confidence variant filter over `pair_subset.vcf.gz`) | `iter_isec.nf.test` |
| `HC_var.vcf.gz.tbi` | Tabix index for `HC_var.vcf.gz` | `iter_isec.nf.test` |

Note the rename: both producing modules prefix their output with the
participant ID (`SUBSET_PAIR_VCF` emits `TEST_PT.pair_subset.vcf.gz`,
`FILTER_VARIANTS` emits `TEST_PT_HC_var.vcf.gz`), but the tests that consume
these as fixtures reference the bare names above. Drop the prefix when copying.

### Per-chromosome phasing stage (derived, 23 files per glob)

| Path | Produced by | Needed by |
| --- | --- | --- |
| `chr_outputs/TEST_PT_chr{1..22,X}_var.vcf.gz` | `ITER_ISEC` (per-chromosome private variants vs. the 1000G panel) | `hetvar_filter.nf.test` |
| `chr_outputs/TEST_PT_chr{1..22,X}_var.vcf.gz.tbi` | Tabix indices for the above | `hetvar_filter.nf.test` |
| `het/TEST_PT_chr{1..22,X}_het.vcf.gz` | `HETVAR_FILTER` (tumor-het sites with adequate depth) | `eagle_phasing.nf.test`, `query_het_ad.nf.test` |
| `het/TEST_PT_chr{1..22,X}_het.vcf.gz.tbi` | Tabix indices for the above | `query_het_ad.nf.test` |
| `phased/TEST_PT_chr{1..22,X}_het_eagle2phased.vcf.gz` | `EAGLE_PHASING` (Eagle2 phasing against the 1000G panel) | `index_phased.nf.test`, `query_het_ad.nf.test` |
| `phased/TEST_PT_chr{1..22,X}_het_eagle2phased.vcf.gz.tbi` | `INDEX_PHASED` (tabix indices for the phased VCFs) | `query_het_ad.nf.test` |

`query_het_ad.nf.test` additionally reads
`het/TEST_PT_chr1_het.vcf.gz` by exact name, to pull the sample list — it is
the chr1 member of the `het/` set above, not a separate file.

### Allelic-count and segmentation stage (derived)

| Path | Produced by | Needed by |
| --- | --- | --- |
| `query/*_HetInfo.txt` | `QUERY_HET_AD` (per-chromosome het site info) | `ac_calc.nf.test` |
| `query/*.AD.chr*.txt` | `QUERY_HET_AD` (per-sample, per-chromosome allelic depths) | `ac_calc.nf.test` |
| `query/*.phasedGT.chr*.1.txt` | `QUERY_HET_AD` (per-sample, per-chromosome phased genotypes) | `ac_calc.nf.test` |
| `TEST_TUMOR_counts.txt` | `COLLECT_READ_COUNTS` (per-bin tumor read counts) | `ac_calc.nf.test`, `tangent_xy.nf.test` |
| `TEST_PT_allelic_counts.txt` | `AC_CALC` | `segment.nf.test` |
| `posttangent_log2.rds` | `TANGENT_XY` | `segment.nf.test` |

Note on `posttangent_log2.rds`: `bin/tangent.R` hardcodes this name internally,
and `modules/local/tangent_xy.nf` renames the result to
`<participant_id>_posttangent_log2.rds` so published copies do not collide
across participants. The fixture keeps the bare name, which is fine —
`segment.nf.test` passes it by path, and `SEGMENT` does not care what it is
called.

## How to stage them

See Task 3 Step 2 in the brief for the exact `aws s3 cp` commands. In short:
ask the dev for a completed run's `output_vcf_file` / `output_vcf_index_file`
plus both samples' `recal.cram` / `.crai`, note the real
`participant_id` / `normal_sample_id` / `tumor_sample_id` / sample sex, and
copy each file to the fixed name listed above.

Then produce the derived fixtures in pipeline order — intervals and variant
subset first, then the per-chromosome phasing stage, then the allelic-count
stage — by running each module standalone and copying its outputs into the
directory shown above. The three R-based modules (`AC_CALC`, `TANGENT_XY`,
`SEGMENT`) require the `iconicc-cnv-r` image, which is not built yet; see
`containers/iconicc-cnv-r/README.md`.
