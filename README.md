# ICONICC-NF

Nextflow port of ICONICC ("HapCNA"), a tumor/normal paired WGS pipeline doing germline
joint-genotyping followed by a haplotype-based allelic copy-number (CNV) analysis.

See `SPEC.md` for the full design. Two logical halves:

1. Per-sample germline calling (`SAMPLE_GERMLINE_CALLING`, Parabricks HaplotypeCaller)
   + joint genotyping/VQSR (`JOINT_GENOTYPING`, GATK4 modules vendored from
   nf-core/sarek and composed with this pipeline's own no-scatter orchestration —
   see `SPEC.md` §12).
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

Pairing is derived automatically: each patient must have exactly one `status: 0`
row, which is crossed against every one of that patient's `status: 1` rows,
producing one pair per tumor (a normal referenced by several tumors is paired,
and its read counts computed, exactly once). A patient with more than one normal,
a tumor with no matching normal, a tumor/normal sex mismatch, or a samplesheet that
yields no pairs at all (for example, no tumor rows) all fail fast with a named error
rather than pairing silently incorrectly or succeeding with no work.

**Path resolution:** relative paths — both the `--input` samplesheet path itself and
any relative path inside the CSV — resolve against the pipeline's own directory
(Nextflow's `projectDir`), not against the directory you launch from; absolute paths,
local or `s3://`, are recommended for real runs to avoid the ambiguity.

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
| `--calling_intervals` | `s3://gatk-test-data/intervals/hg38.even.handcurated.20k.intervals` |
| `--dbsnp_vcf` / `--dbsnp_vcf_index` | `s3://broad-references/hg38/v0/...` |
| `--hapmap_resource_vcf` / `--omni_resource_vcf` / `--one_thousand_genomes_resource_vcf` / `--mills_resource_vcf` / `--axiomPoly_resource_vcf` (each with a matching `_index`) | `s3://broad-references/hg38/v0/...` |

Only the five final outputs named in `SPEC.md` §6 are published to `--outdir`:
`read_depth`, `seg_file`, `ac_cts`, `segment_plots`, `segment_cts`. Publishing is
disabled for every upstream process, whose per-chromosome intermediates run to
hundreds of GB per run and are not consumed downstream.

## Testing

```bash
nf-test test
```

## Known gaps in this phase

Three prerequisites are deliberately unresolved. All three are accepted for this phase,
and together they mean most module-level tests and the subworkflow-level test cannot
fully pass yet. The pipeline code itself is complete and reviewed.

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

**3. No real Parabricks/NGC execution has been validated.** This dev environment (and
the dev's own machine, confirmed) has no NVIDIA GPU, so `SAMPLE_GERMLINE_CALLING`'s
wiring is tested via nf-test's `-stub` flag rather than a real `pbrun` run. Real
Parabricks validation can only happen on Cirro/AWS against an actual GPU instance —
where it was always going to run in production anyway. `JOINT_GENOTYPING`'s GATK4 chain
also has no real fixture: VQSR's Gaussian mixture model cannot fit on a handful of
variants regardless of container or network access, so its tests assert the desired end
state against placeholder fixtures, the same convention as gap 1 above.
