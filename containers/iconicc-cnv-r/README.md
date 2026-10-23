# containers/iconicc-cnv-r

R environment for the CNV/allelic-analysis core (`ac_calc.R`, `tangent.R`, `segment.R`).
See `SPEC.md` §7 for the pinning rationale.

## Published

`params.cnv_r_container` defaults to `ghcr.io/break-through-cancer/iconicc-cnv-r:latest`,
a public image. A commit-pinned tag is also published alongside `latest` for
reproducibility (tagged with the short SHA of the commit that last touched the
Dockerfile).

## Rebuilding and republishing

```bash
docker build -t ghcr.io/break-through-cancer/iconicc-cnv-r:latest containers/iconicc-cnv-r
docker push ghcr.io/break-through-cancer/iconicc-cnv-r:latest
```

To point a run at a different image instead (e.g. a local build for testing):

```bash
nextflow run main.nf --cnv_r_container iconicc-cnv-r:latest ...
```
