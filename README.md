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

`main.nf` takes a single samplesheet, driving the germline read-counting,
`SAMPLE_GERMLINE_CALLING`, `JOINT_GENOTYPING`, and `CNV_CORE` subworkflows for every
tumor/normal pair it describes. Every reference file has a working default in
`nextflow.config` and only needs an override to point somewhere else.

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
| `--cnv_r_container` | `ghcr.io/break-through-cancer/iconicc-cnv-r:latest` |
| `--ref_fasta` / `--ref_fasta_fai` / `--ref_fasta_dict` | hg38 from `s3://broad-references` |
| `--ref_dir_1000g` | `s3://btc-refdata/iconicc/reference-genome-1000G.tar` |
| `--genetic_map` | `s3://btc-refdata/iconicc/reference-genome-genetic_map_hg38_withX.txt.gz` |
| `--pon` | `s3://btc-refdata/iconicc/BTC_pon.rds` |
| `--ref_clinical` / `--cohort_clinical` | `s3://btc-refdata/iconicc/...` |
| `--calling_intervals` | `s3://gatk-test-data/intervals/hg38.even.handcurated.20k.intervals` |
| `--dbsnp_vcf` / `--dbsnp_vcf_index` | `s3://broad-references/hg38/v0/...` |
| `--hapmap_resource_vcf` / `--omni_resource_vcf` / `--one_thousand_genomes_resource_vcf` / `--mills_resource_vcf` / `--axiomPoly_resource_vcf` (each with a matching `_index`) | `s3://broad-references/hg38/v0/...` |
| `--haplotypecaller_memory_gb` | `72` (then scaled by `task.attempt` on each OOM retry; see below) |
| `--scatter_chromosomes` | `chr1-22 + X` (hg38; matches `--ref_dir_1000g`'s own per-chromosome files) |
| `--read_counts_scatter_mode` | `groups` (`groups` or `per_chromosome`; see below) |
| `--read_counts_scatter_count` | `8` (shards per sample, `groups` mode only) |
| `--cnv_scatter_mode` | `groups` (`groups` or `per_chromosome`; see below) |
| `--cnv_scatter_count` | `4` (shards per pair, `groups` mode only) |

`PARABRICKS_HAPLOTYPECALLER` and the GATK joint-genotyping chain
(`GATK4_GENOMICSDBIMPORT`, `GATK4_GENOTYPEGVCFS`, `VARIANTRECALIBRATOR_SNP`/`_INDEL`,
`APPLYVQSR_SNP`/`_INDEL`) all scale their memory as `base * task.attempt` rather than a
flat value, and `.cirro/process-compute.config` retries `PARABRICKS_HAPLOTYPECALLER` on
exactly the exit codes a container memory-kill produces (including AWS Batch's
no-exit-code case, reported to Nextflow as `2147483647`) — a real OOM on real (non-toy)
coverage data retries automatically with more memory rather than failing identically
every attempt. Override `--haplotypecaller_memory_gb` if even the scaled default isn't
enough for a given sample's coverage.

`SAMPLE_READ_COUNTS` scatters each sample's `COLLECT_READ_COUNTS` call across multiple
parallel shards rather than one whole-genome task (real high-coverage CRAMs pushed that
single call into a tens-of-hours wall-clock — see `SPEC.md` §14). Two mutually exclusive
modes: `groups` (default) splits into `--read_counts_scatter_count` size-balanced shards
via GATK's own `SplitIntervals`, bounding total job count for cohort-scale runs; `per_chromosome`
shards one-per-entry in `--scatter_chromosomes` instead, for maximum parallelism on a
single-sample run where job count isn't the concern. `--scatter_chromosomes` is also the
chromosome universe `--cnv_scatter_mode`/`--cnv_scatter_count` use — override it for
other species/reference builds.

`ITER_ISEC`/`EAGLE_PHASING` (the CNV chain's two slowest steps — real timing on a
production pair: 63.7 min and 65.7 min respectively, for all 23 chromosomes sequentially
in one task) scatter the same way, over the same `--scatter_chromosomes` universe, via
`--cnv_scatter_mode`/`--cnv_scatter_count`. `HETVAR_FILTER`/`INDEX_PHASED` stay
single-task deliberately (3.8 min / 0.8 min for all 23 chromosomes combined — scattering
them would likely add more overhead than it would save; see `SPEC.md` §14.1).

The five final outputs named in `SPEC.md` §6 (`read_depth`, `seg_file`, `ac_cts`,
`segment_plots`, `segment_cts`) are published to `--outdir`, plus each sample's indexed
germline GVCF (from `INDEX_GVCF`, under `haplotypecaller/`) — a standard GATK-format GVCF
(`--gvcf` mode, same GQ confidence bands and annotation groups any GATK joint-genotyping
pipeline expects), reusable beyond this pipeline's own CNV analysis. This is a genuine
addition beyond the original WDL pipeline's own output set, not something it ever
published. Publishing stays disabled for every other upstream process, whose
per-chromosome intermediates run to hundreds of GB per run and are not consumed
downstream.

## Testing

```bash
nf-test test
```

`.cirro/preprocess.py`'s own offline unit tests run separately, without `nf-test` or the
Cirro SDK:

```bash
python3 .cirro/test_preprocess.py -v
```

## Cirro deployment

This pipeline is deployed as a Cirro process launched from a completed
`parabricks_bqsr` dataset (that pipeline's alignment-entry process specifically --
`break-through-cancer/parabricks-bqsr`, not its separate fastq-entry process). See
`SPEC.md` §13 for the full design.

`.cirro/preprocess.py` auto-discovers each selected sample's recalibrated CRAM/CRAI
pair and joins it against Cirro's own per-sample metadata table -- set `patient`,
`sex` (`male`/`female`), and `status` (`Normal`/`Tumor`, case-insensitive) for every
sample directly in Cirro's own sample-metadata UI before launching. A normal sample
referenced by multiple tumor samples for the same patient is supported and deduplicated
automatically, same as this pipeline's own samplesheet convention (§11.1).

**Manual follow-up, not automated by any file in this repo**: this org's convention for
its own custom pipelines does not commit a `process-definition.json` (unlike the public
nf-core forks), so registering this process happens entirely in Cirro's admin UI. That
registration needs: a process id; `usesSampleSheet: true`; the executor/engine version;
a pointer to this code repository; and the `.cirro/process-compute.config` compute
profile attached. Separately, once a process id is assigned, add that id to
`parabricks_bqsr`'s own `childProcessIds` list so a user can launch this pipeline
directly from a completed `parabricks_bqsr` dataset in Cirro's UI -- a one-line edit in
the *other* repo's own Cirro registration, not something either repo's committed files
accomplish on their own.

## Known gaps in this phase

Two prerequisites are deliberately unresolved. Both are accepted for this phase, and
together they mean most module-level tests and the subworkflow-level test cannot fully
pass yet. The pipeline code itself is complete and reviewed. (The R container gap from
earlier in this phase is resolved: `params.cnv_r_container` now defaults to the
published `ghcr.io/break-through-cancer/iconicc-cnv-r:latest`. See
`containers/iconicc-cnv-r/README.md` for rebuild/republish commands.)

**1. Real fixture data is not staged.** `tests/fixtures/cnv_core/` holds only a README;
no patient genomic data has been copied into it, and the directory is gitignored so that
none ever lands in this repository. This blocks every test that reads a fixture: all of
`tests/local/` except `preprocess_intervals.nf.test` (which needs only the reference
fasta), plus `tests/subworkflows/cnv_core.nf.test`. See
**`tests/fixtures/cnv_core/README.md`** for the full list of expected files and how to
stage them.

**2. No real Parabricks/NGC execution has been validated.** This dev environment (and
the dev's own machine, confirmed) has no NVIDIA GPU, so `SAMPLE_GERMLINE_CALLING`'s
wiring is tested via nf-test's `-stub` flag rather than a real `pbrun` run. Real
Parabricks validation can only happen on Cirro/AWS against an actual GPU instance —
where it was always going to run in production anyway. `JOINT_GENOTYPING`'s GATK4 chain
also has no real fixture: VQSR's Gaussian mixture model cannot fit on a handful of
variants regardless of container or network access, so its tests assert the desired end
state against placeholder fixtures, the same convention as gap 1 above.
