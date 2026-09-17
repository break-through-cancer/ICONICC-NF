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
