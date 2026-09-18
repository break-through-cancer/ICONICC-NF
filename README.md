# ICONICC-NF

Nextflow port of ICONICC ("HapCNA"), a tumor/normal paired WGS pipeline doing germline
joint-genotyping followed by a haplotype-based allelic copy-number (CNV) analysis.

See `SPEC.md` for the full design. Two logical halves:

1. Per-sample germline calling (Parabricks) + joint genotyping/VQSR (reuses
   nf-core/sarek's `BAM_JOINT_CALLING_GERMLINE_GATK` subworkflow directly).
2. The CNV/allelic-analysis core, hand-ported from `ICON_Tasks.2.wdl` — the actual
   science; see `subworkflows/local/cnv_core/`.

## Local usage

`main.nf` takes a single samplesheet, driving the germline read-counting and
`CNV_CORE` subworkflows for every tumor/normal pair it describes. Every reference
file has a working default in `nextflow.config` and only needs an override to point
somewhere else.

```bash
nextflow run main.nf --input samplesheet.csv
```

`--input` is the only required parameter; omitting it fails immediately with a
message naming the missing parameter and pointing back at this README and
`assets/samplesheet.csv`. A path inside the samplesheet that does not exist also
fails at startup naming that path, rather than deep inside a container.

### Samplesheet format

CSV, one row per sample (not per pair) — modeled on nf-core/sarek's own
samplesheet convention. See `assets/samplesheet.csv` for a filled-in template;
full rationale in `SPEC.md` §11.

| Column | Required | Meaning |
| --- | --- | --- |
| `patient` | always | Groups a patient's normal and tumor rows for pairing. |
| `sample` | always | Sample ID, unique within the samplesheet. |
| `sex` | always | `female` or `male`. |
| `status` | always | `0` = normal, `1` = tumor. |
| `cram` / `crai` | always | The sample's recalibrated CRAM from `sarek_align`. |
| `joint_vcf` / `joint_vcf_index` | tumor rows only | Transitional: the joint-genotyped, VQSR-filtered VCF for that tumor's pair with its patient's normal. Ignored on normal rows. Once Phase 2 (`SPEC.md` §5) lands, these two columns are computed by the pipeline instead of supplied, and will disappear from the schema. |

Pairing is derived automatically: each patient must have exactly one `status: 0`
row, which is crossed against every one of that patient's `status: 1` rows,
producing one pair per tumor (a normal referenced by several tumors is paired,
and its read counts computed, exactly once). A patient with more than one normal,
a tumor with no matching normal, or a tumor/normal sex mismatch all fail fast with
a named error rather than pairing silently incorrectly.

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
