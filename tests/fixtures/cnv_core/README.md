# CNV-core validation fixtures

This directory is expected to hold real genomic data staged from a completed
ICONICC WDL run, per Task 3 Step 2 of the `cnv-core-subworkflow` plan
(`.superpowers/sdd/2026-09-17-cnv-core-subworkflow/task-3-brief.md`). These
files are gitignored (`tests/fixtures/` in `.gitignore`) and are not, and
must never be, committed to this repository.

**Status as of this task: not staged.** The controller asked the dev for a
Cirro dataset ID / S3 paths to a known-good run and was told to defer
fixture staging and proceed with code-only work. `nf-test` runs against this
module will fail until these files are staged by hand.

## Expected files

| File | What it is |
| --- | --- |
| `joint.vcf.gz` | The joint-genotyped, VQSR-filtered VCF (`output_vcf_file`) from a completed ICONICC WDL run's joint-genotyping subworkflow. |
| `joint.vcf.gz.tbi` | Tabix index for `joint.vcf.gz` (`output_vcf_index_file`). |
| `normal.recal.cram` | The normal sample's recalibrated CRAM from `sarek_align`, for the same run. |
| `normal.recal.cram.crai` | CRAM index for `normal.recal.cram`. |
| `tumor.recal.cram` | The tumor sample's recalibrated CRAM from `sarek_align`, for the same run. |
| `tumor.recal.cram.crai` | CRAM index for `tumor.recal.cram`. |
| `sample_info.env` | Shell-sourceable file recording the real `PARTICIPANT_ID`, `NORMAL_SAMPLE_ID`, `TUMOR_SAMPLE_ID`, and `SAMPLE_SEX` for the staged run, so later tests use real IDs instead of placeholders. |
| `wgs.preprocessed.interval_list` | Output of Task 4's `PREPROCESS_INTERVALS` module (GATK `PreprocessIntervals`), needed by `tests/local/collect_read_counts.nf.test`. Generate by running that module standalone against the reference fasta and copying its output here — not fetched from an external run. |
| `pair_subset.vcf.gz` | Output of Task 3's `SUBSET_PAIR_VCF` module (normal+tumor sample subset of `joint.vcf.gz`), needed by `tests/local/filter_variants.nf.test`. Generate by running that module standalone and copying its output here — not fetched from an external run. |
| `pair_subset.vcf.gz.tbi` | Tabix index for `pair_subset.vcf.gz`. |

## How to stage them

See Task 3 Step 2 in the brief for the exact `aws s3 cp` commands. In short:
ask the dev for a completed run's `output_vcf_file` / `output_vcf_index_file`
plus both samples' `recal.cram` / `.crai`, note the real
`participant_id` / `normal_sample_id` / `tumor_sample_id` / sample sex, and
copy each file to the fixed name listed above.
