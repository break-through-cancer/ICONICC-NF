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
