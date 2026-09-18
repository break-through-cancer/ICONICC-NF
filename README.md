# ICONICC-NF

Nextflow port of ICONICC ("HapCNA"), a tumor/normal paired WGS pipeline doing germline
joint-genotyping followed by a haplotype-based allelic copy-number (CNV) analysis.

See `SPEC.md` for the full design. Two logical halves:

1. Per-sample germline calling (Parabricks) + joint genotyping/VQSR (reuses
   nf-core/sarek's `BAM_JOINT_CALLING_GERMLINE_GATK` subworkflow directly).
2. The CNV/allelic-analysis core, hand-ported from `ICON_Tasks.2.wdl` — the actual
   science; see `subworkflows/local/cnv_core/`.

## Local usage

`main.nf` runs the `CNV_CORE` subworkflow standalone. Ten per-run parameters are
required; every reference file has a working default in `nextflow.config` and only
needs an override to point somewhere else.

```bash
nextflow run main.nf \
    --participant_id    PT001 \
    --normal_sample_id  PT001_N \
    --tumor_sample_id   PT001_T \
    --sample_sex        female \
    --joint_vcf         /data/PT001/joint.vcf.gz \
    --joint_vcf_index   /data/PT001/joint.vcf.gz.tbi \
    --tumor_cram        /data/PT001/tumor.recal.cram \
    --tumor_crai        /data/PT001/tumor.recal.cram.crai \
    --normal_cram       /data/PT001/normal.recal.cram \
    --normal_crai       /data/PT001/normal.recal.cram.crai
```

`--sample_sex` takes `female` or `male`. Omitting any of the ten required parameters
fails immediately, naming the ones that are missing; a path that does not exist fails
at startup naming that path, rather than deep inside a container.

Reference defaults that may be worth overriding:

| Parameter | Default |
| --- | --- |
| `--outdir` | `./results` |
| `--cnv_r_container` | `iconicc-cnv-r:latest` (see the gap below) |
| `--ref_fasta` / `--ref_fasta_fai` / `--ref_fasta_dict` | hg38 from `s3://broad-references` |
| `--ref_dir_1000g` | `s3://btc-refdata/iconicc/reference-genome-1000G.tar` |
| `--genetic_map` | `s3://btc-refdata/iconicc/reference-genome-genetic_map_hg38_withX.txt.gz` |
| `--pon` | `s3://btc-refdata/iconicc/BTC_pon.rds` |
| `--ref_clinical` / `--cohort_clinical` | `s3://btc-refdata/iconicc/...` |

Only the five final outputs named in `SPEC.md` §6 are published to `--outdir`:
`read_depth`, `seg_file`, `ac_cts`, `segment_plots`, `segment_cts`. Publishing is
disabled for every upstream process, whose per-chromosome intermediates run to
hundreds of GB per run and are not consumed downstream.

## Testing

```bash
nf-test test
```

## Known gaps in this phase

Two prerequisites are deliberately unresolved. Both are accepted for this phase, and
together they mean most module-level tests and the subworkflow-level test cannot fully
pass yet. The pipeline code itself is complete and reviewed.

**1. Real fixture data is not staged.** `tests/fixtures/cnv_core/` holds only a README;
no patient genomic data has been copied into it, and the directory is gitignored so that
none ever lands in this repository. This blocks every test that reads a fixture: all of
`tests/local/` except `preprocess_intervals.nf.test` (which needs only the reference
fasta), plus `tests/subworkflows/cnv_core.nf.test`. See
**`tests/fixtures/cnv_core/README.md`** for the full list of expected files and how to
stage them.

**2. The R container is not built.** `AC_CALC`, `TANGENT_XY`, and `SEGMENT` all run in
`params.cnv_r_container`, which defaults to the local tag `iconicc-cnv-r:latest`. The
Dockerfile exists but the image has never been built or published, so those three
processes cannot run anywhere yet, locally or remotely. See
**`containers/iconicc-cnv-r/README.md`** for the build and publish commands.
