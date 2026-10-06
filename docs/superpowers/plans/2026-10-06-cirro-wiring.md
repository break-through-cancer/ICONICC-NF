# Cirro Wiring (Phase 3) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deploy ICONICC-NF as an actual Cirro process, launchable directly from a
completed `parabricks_bqsr` dataset, matching how `sarek_align`/`sarek_call_variants`
are already deployed.

**Architecture:** A new `.cirro/` directory holds the five Cirro-facing config/script
files this org's own custom (non-nf-core-fork) pipelines commit — `preprocess.py` is the
one with real logic (auto-discovers recalibrated CRAM/CRAI pairs from the input dataset,
joins Cirro's own per-sample `patient`/`sex`/`status` metadata, re-validates the same
pairing rules `SAMPLESHEET_TO_PAIRS` already enforces, writes the final samplesheet CSV).
Everything else is thin, declarative config.

**Tech Stack:** Python 3 (first Python code in this repo — `pandas`, stdlib
`unittest`/`importlib`, no new runtime dependency since the real Cirro SDK is stubbed
out for tests), Nextflow/Groovy config (unchanged elsewhere), JSON.

**Spec:** `SPEC.md` §13 (commits `cef0d8e`, `81b0914` — the second commit corrected
§13.5's file list after a local interruption lost that edit before it could be
committed the first time; read the current file, not an earlier version).

## Global Constraints

- **No committed `process-definition.json`.** Verified directly against
  `break-through-cancer/parabricks-bqsr`'s full repo tree (`gh api
  repos/break-through-cancer/parabricks-bqsr/git/trees/fastq-entry?recursive=true`) —
  zero instances anywhere. This org's real convention for its own custom pipelines is
  registering the process directly in Cirro's admin UI, not as a version-controlled
  file. Do not create one.
- **Parent process: `parabricks_bqsr`'s alignment-entry process only** — not
  `sarek_align` directly, not `parabricks_bqsr`'s separate fastq-entry process. Every
  real input to this pipeline has gone through BQSR + optional quantization first.
- **Samplesheet shape is fixed**: `patient,sample,sex,status,cram,crai` — exactly what
  `subworkflows/local/samplesheet_to_pairs/main.nf`'s `validateSamplesheetRow()` already
  expects (`sex` is the literal strings `male`/`female`, not sarek's `XX`/`XY`; `status`
  is `0`=normal/`1`=tumor). Do not invent a different column order or encoding.
- **Validation duplication is deliberate** (confirmed with the dev): `preprocess.py`
  re-implements `SAMPLESHEET_TO_PAIRS`'s own pairing hard-errors (duplicate normal,
  orphan tumor, sex mismatch) for fast pre-launch failure. There is no shared source of
  truth between the Groovy and Python copies — if one changes, the other needs the same
  change by hand.
- **GPU wiring is queue-based, not instance-family-based**: `accelerator = 1` +
  `queue = env('PW_ONDEMAND_JOB_QUEUE')` on `PARABRICKS_HAPLOTYPECALLER`, copied from
  `parabricks_bqsr`'s own `.cirro/alignment/process-compute.config` — the real, current
  pattern this org uses, not the older instance-family approach SPEC.md §2 originally
  described. The accelerator count is explicitly not treated as final (SPEC.md §13.3) —
  implement it as specified, don't try to resolve the sibling repo's own internal
  inconsistency (its fastq-entry config sets no `accelerator` directive at all for a
  heavier workload; this plan doesn't explain why, and isn't blocked on explaining it).
- **Form stays thin**: no GATK/VQSR/reference-file parameters exposed — every one of
  those already has a working default in `nextflow.config` (SPEC.md §3, confirmed
  already true for this pipeline). `process-form.json` needs no properties at all for
  v1.
- **Testing**: `preprocess.py`'s logic (the pairing/validation function) gets a real,
  standalone test suite (`test_preprocess.py`) that stubs the `cirro` package rather than
  requiring real Cirro access — following the exact pattern already proven in
  `CirroBio/Cirro-pipelines`' `nf-core/sarek_call_variants/3/test_preprocess.py` (stub
  `sys.modules['cirro.helpers.preprocess_dataset']` before dynamically loading
  `preprocess.py` via `importlib.util.spec_from_file_location`). Real end-to-end Cirro
  launch validation is an accepted gap for this plan — a real test launch is a manual,
  human step (the dev administers the actual Cirro project), not something any test in
  this repo can exercise.

---

### Task 1: Cirro process scaffolding

**Files:**
- Create: `.cirro/process-form.json`
- Create: `.cirro/process-compute.config`
- Create: `.cirro/process-input.json`
- Create: `.cirro/process-output.json`

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces: the four static/declarative Cirro config files. Task 2's `preprocess.py`
  writes its output CSV to the exact path `process-input.json`'s `"input"` key
  references (`$.dataset.s3|/config/manifest.csv` — a Cirro expression meaning "this
  dataset's own S3 prefix, `/config/manifest.csv` relative to it", not a literal
  repo-root path).

- [ ] **Step 1: Write `process-form.json`**

```json
{
  "form": {
    "type": "object",
    "properties": {}
  },
  "ui": {}
}
```

No properties: every parameter this pipeline needs already has a working default in
`nextflow.config` (SPEC.md §3's "do not build a richer configuration UI" directive,
already true for every phase of this pipeline) — there is nothing left to expose.

- [ ] **Step 2: Write `process-compute.config`**

```groovy
process {
    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 2

    withName: 'PARABRICKS_HAPLOTYPECALLER' {
        accelerator   = 1
        queue         = env('PW_ONDEMAND_JOB_QUEUE')
        errorStrategy = { task.exitStatus in ((130..145) + 104 + (175..177) + 2147483647) ? 'retry' : 'finish' }
        maxRetries    = 2
    }
}
```

Copied directly from `break-through-cancer/parabricks-bqsr`'s own
`.cirro/alignment/process-compute.config` (real, current precedent — SPEC.md §13.3), the
GPU-requiring-process-specific retry predicate included (it is broader than the
top-level default: it also retries on exit codes 175-177 and the sentinel
`2147483647`, both AWS-Batch-specific signals for a reclaimed spot GPU instance, not
just the generic 130-145/104 range every other process gets). This file is *only*
applied when Cirro launches the pipeline with this compute profile attached — it is
never merged into the committed `nextflow.config`, by design (a Phase 2 ruling in this
same repo already established that GPU/queue-specific directives must stay out of the
shared config; this file is the Cirro-native place for them instead).

- [ ] **Step 3: Write `process-input.json`**

```json
{
  "input": "$.dataset.s3|/config/manifest.csv",
  "outdir": "$.dataset.dataPath"
}
```

`preprocess.py` (Task 2) writes its assembled samplesheet to exactly this path. No other
param needs mapping here — every reference file (`ref_fasta`, `ref_dir_1000g`,
`calling_intervals`, the VQSR resource files, etc.) already has a working S3 default in
`nextflow.config`, and Cirro passes through any param it doesn't explicitly override.

- [ ] **Step 4: Write `process-output.json`**

```json
{
  "commands": []
}
```

Matches `parabricks_bqsr`'s own minimal precedent — no special output-viewer wiring
needed for this pipeline's five final outputs (SPEC.md §6).

- [ ] **Step 5: Commit**

```bash
git add .cirro/process-form.json .cirro/process-compute.config .cirro/process-input.json .cirro/process-output.json
git commit -m "Add Cirro process scaffolding (form, compute-config, input/output mapping)"
```

---

### Task 2: `preprocess.py` + `test_preprocess.py`

**Files:**
- Create: `.cirro/preprocess.py`
- Create: `.cirro/test_preprocess.py`

**Interfaces:**
- Consumes: nothing from other tasks (independent of Task 1's static config files).
- Produces: `build_samplesheet(files, samplesheet)` — a pure function, not touching the
  real Cirro SDK at all, so it is directly unit-testable. `files` is an iterable of path
  strings (matching `ds.files["file"]`'s shape — a pandas Series of strings iterates as
  strings). `samplesheet` is a `pandas.DataFrame` with at least `sample`, `patient`,
  `sex`, `status` columns (matching `ds.samplesheet`'s real shape, confirmed directly
  against `CirroBio/Cirro-pipelines`' `sarek_call_variants/3/preprocess.py`, which
  already relies on this same Cirro-native per-sample metadata table). Returns a
  `pandas.DataFrame` with columns `["patient", "sample", "sex", "status", "cram",
  "crai"]`, or raises `ValueError` with a descriptive message on any invalid input.

- [ ] **Step 1: Write the failing tests**

```python
# .cirro/test_preprocess.py
"""Offline tests for ICONICC-NF's Cirro preprocess script.

Runs without the Cirro SDK or boto3: both are stubbed before preprocess is loaded,
matching the pattern established in CirroBio/Cirro-pipelines'
nf-core/sarek_call_variants/3/test_preprocess.py. Plain `python3 -m unittest`.
"""
import sys
import types
import unittest
from pathlib import Path

import pandas as pd


def _load_preprocess():
    helpers = types.ModuleType("cirro.helpers.preprocess_dataset")
    helpers.PreprocessDataset = object
    for name in ["cirro", "cirro.helpers"]:
        sys.modules.setdefault(name, types.ModuleType(name))
    sys.modules["cirro.helpers.preprocess_dataset"] = helpers

    import importlib.util
    path = Path(__file__).with_name("preprocess.py")
    spec = importlib.util.spec_from_file_location("iconicc_nf_preprocess", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


preprocess = _load_preprocess()


def samplesheet_df(rows):
    """rows: list of (sample, patient, sex, status) tuples."""
    return pd.DataFrame(
        [dict(sample=s, patient=p, sex=x, status=st) for s, p, x, st in rows]
    )


def recal_paths(sample):
    return [
        f"preprocessing/recalibrated/{sample}/{sample}.recal.cram",
        f"preprocessing/recalibrated/{sample}/{sample}.recal.cram.crai",
    ]


class BuildSamplesheetTests(unittest.TestCase):

    def test_single_pair(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "female", "Normal"),
            ("PT01_tumor", "PT01", "female", "Tumor"),
        ])
        result = preprocess.build_samplesheet(paths, sheet)
        self.assertEqual(sorted(result["sample"]), ["PT01_normal", "PT01_tumor"])
        self.assertEqual(
            result.loc[result["sample"] == "PT01_normal", "status"].iloc[0], 0
        )
        self.assertEqual(
            result.loc[result["sample"] == "PT01_tumor", "status"].iloc[0], 1
        )
        self.assertEqual(list(result.columns), ["patient", "sample", "sex", "status", "cram", "crai"])

    def test_shared_normal_across_two_tumors(self):
        paths = (
            recal_paths("PT01_normal")
            + recal_paths("PT01_tumor_a")
            + recal_paths("PT01_tumor_b")
        )
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "female", "Normal"),
            ("PT01_tumor_a", "PT01", "female", "Tumor"),
            ("PT01_tumor_b", "PT01", "female", "Tumor"),
        ])
        result = preprocess.build_samplesheet(paths, sheet)
        self.assertEqual(len(result), 3)

    def test_missing_index_drops_sample(self):
        paths = ["preprocessing/recalibrated/PT01_normal/PT01_normal.recal.cram"]
        sheet = samplesheet_df([("PT01_normal", "PT01", "female", "Normal")])
        with self.assertRaises(ValueError):
            preprocess.build_samplesheet(paths, sheet)

    def test_no_recalibrated_files_raises(self):
        paths = ["preprocessing/markduplicates/PT01_normal/PT01_normal.md.cram"]
        sheet = samplesheet_df([("PT01_normal", "PT01", "female", "Normal")])
        with self.assertRaises(ValueError):
            preprocess.build_samplesheet(paths, sheet)

    def test_duplicate_normal_is_hard_error(self):
        paths = (
            recal_paths("PT01_normal_a")
            + recal_paths("PT01_normal_b")
            + recal_paths("PT01_tumor")
        )
        sheet = samplesheet_df([
            ("PT01_normal_a", "PT01", "female", "Normal"),
            ("PT01_normal_b", "PT01", "female", "Normal"),
            ("PT01_tumor", "PT01", "female", "Tumor"),
        ])
        with self.assertRaisesRegex(ValueError, "more than one normal"):
            preprocess.build_samplesheet(paths, sheet)

    def test_orphan_tumor_is_hard_error(self):
        paths = recal_paths("PT01_tumor")
        sheet = samplesheet_df([("PT01_tumor", "PT01", "female", "Tumor")])
        with self.assertRaisesRegex(ValueError, "no matching normal"):
            preprocess.build_samplesheet(paths, sheet)

    def test_sex_mismatch_is_hard_error(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "female", "Normal"),
            ("PT01_tumor", "PT01", "male", "Tumor"),
        ])
        with self.assertRaisesRegex(ValueError, "Sex mismatch"):
            preprocess.build_samplesheet(paths, sheet)

    def test_zero_pairs_is_hard_error(self):
        paths = recal_paths("PT01_normal")
        sheet = samplesheet_df([("PT01_normal", "PT01", "female", "Normal")])
        with self.assertRaisesRegex(ValueError, "zero tumor/normal pairs"):
            preprocess.build_samplesheet(paths, sheet)

    def test_sample_missing_from_cirro_metadata_is_hard_error(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = samplesheet_df([("PT01_normal", "PT01", "female", "Normal")])
        with self.assertRaisesRegex(ValueError, "No Cirro sample metadata"):
            preprocess.build_samplesheet(paths, sheet)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m unittest .cirro.test_preprocess -v` (from the repo root; if `pandas`
isn't installed in whatever Python runs this, create a throwaway venv first:
`python3 -m venv /tmp/iconicc-preprocess-venv && /tmp/iconicc-preprocess-venv/bin/pip
install pandas` — this mirrors the venv approach already used once this session for
unrelated Cirro SDK work — then run with that venv's `python3`).
Expected: FAIL with `ModuleNotFoundError` or `AttributeError: module 'iconicc_nf_preprocess' has no attribute 'build_samplesheet'` — `preprocess.py` does not exist yet.

- [ ] **Step 3: Write `preprocess.py`**

```python
#!/usr/bin/env python3
"""Cirro preprocess script for ICONICC-NF.

Auto-discovers recalibrated CRAM/BAM + index pairs from the selected input dataset
(parabricks_bqsr's own output shape, preprocessing/recalibrated/<sample>/...), joins
Cirro's own per-sample patient/sex/status metadata, and assembles the exact samplesheet
subworkflows/local/samplesheet_to_pairs/main.nf expects. See SPEC.md §13.2.
"""
import re

import pandas as pd

RECAL_FILE = re.compile(
    r"(?:^|/)preprocessing/recalibrated/(?P<sample>[^/]+)/[^/]+\.(?P<ext>bam|cram)$"
)
INDEX_SUFFIX = {"bam": "bai", "cram": "crai"}
COLUMNS = ["patient", "sample", "sex", "status", "cram", "crai"]


def alignment_pairs(paths):
    """Every indexed recalibrated CRAM/BAM in the input dataset, one row per sample.

    Only preprocessing/recalibrated/<sample>/: parabricks_bqsr's own, and only, output
    stage (SPEC.md §13.1 -- that pipeline's alignment-entry process is this pipeline's
    sole declared parent). A sample with no matching index is dropped --
    SAMPLESHEET_TO_PAIRS requires both cram and crai.
    """
    present = set(paths)
    rows = []
    for path in paths:
        match = RECAL_FILE.search(path)
        if not match:
            continue
        ext = match["ext"]
        index = f"{path}.{INDEX_SUFFIX[ext]}"
        if index not in present:
            continue
        rows.append(dict(sample=match["sample"], cram=path, crai=index))
    return pd.DataFrame(rows, columns=["sample", "cram", "crai"])


def build_samplesheet(files, samplesheet: pd.DataFrame) -> pd.DataFrame:
    """Build ICONICC-NF's samplesheet (patient,sample,sex,status,cram,crai).

    Re-implements subworkflows/local/samplesheet_to_pairs/main.nf's own pairing
    hard-errors (duplicate normal, orphan tumor, sex mismatch) so a bad samplesheet
    fails before Cirro launches compute, not after -- SPEC.md §13.2's deliberate
    validation-duplication decision. If those rules change in that Groovy file, update
    this function to match; there is no shared source of truth between the two.
    """
    pairs = alignment_pairs(files)
    if pairs.empty:
        raise ValueError(
            "No indexed recalibrated CRAM/BAM found under "
            "preprocessing/recalibrated/<sample>/ in the input dataset. This pipeline "
            "expects parabricks_bqsr's own output shape."
        )

    meta = samplesheet.set_index("sample")
    missing = sorted(set(pairs["sample"]) - set(meta.index))
    if missing:
        raise ValueError(
            "No Cirro sample metadata (patient/sex/status) found for: "
            + ", ".join(missing)
        )

    manifest = pairs.assign(
        patient=pairs["sample"].map(meta["patient"]),
        sex=pairs["sample"].map(meta["sex"]),
        status=pairs["sample"].map(meta["status"]),
    )

    # Normalize status the same way sarek_call_variants does -- "Normal"/"normal" -> 0,
    # "Tumor"/"tumor" -> 1 -- then validate every row resolved to a real 0/1.
    manifest = manifest.replace(to_replace=dict(
        status=dict(normal=0, Normal=0, tumor=1, Tumor=1)
    ))
    bad_status = manifest.loc[~manifest["status"].isin([0, 1]), "sample"].tolist()
    if bad_status:
        raise ValueError(
            "status must be Normal/0 or Tumor/1 for: " + ", ".join(bad_status)
        )
    manifest = manifest.assign(status=manifest["status"].astype(int))

    # Normalize sex -- SAMPLESHEET_TO_PAIRS accepts only the literal strings
    # 'male'/'female' (not sarek's XX/XY), case-insensitively here for user tolerance.
    manifest = manifest.assign(sex=manifest["sex"].astype(str).str.strip().str.lower())
    bad_sex = manifest.loc[~manifest["sex"].isin(["male", "female"]), "sample"].tolist()
    if bad_sex:
        raise ValueError("sex must be 'male' or 'female' for: " + ", ".join(bad_sex))

    # Mirror SAMPLESHEET_TO_PAIRS's pairing hard-errors exactly: >1 normal per patient,
    # a tumor whose patient has no normal, a sex mismatch within a pair.
    normals = manifest.loc[manifest["status"] == 0]
    tumors = manifest.loc[manifest["status"] == 1]

    dup_normal_patients = normals["patient"][normals["patient"].duplicated()].unique()
    for patient in dup_normal_patients:
        samples = normals.loc[normals["patient"] == patient, "sample"].tolist()
        raise ValueError(
            f"Patient '{patient}' has more than one normal (status=0) sample: "
            + ", ".join(samples)
        )

    normal_by_patient = normals.set_index("patient")
    for _, tumor in tumors.iterrows():
        patient = tumor["patient"]
        if patient not in normal_by_patient.index:
            raise ValueError(
                f"Tumor sample '{tumor['sample']}' (patient '{patient}') has no "
                "matching normal (status=0) sample in the samplesheet"
            )
        normal = normal_by_patient.loc[patient]
        if normal["sex"] != tumor["sex"]:
            raise ValueError(
                f"Sex mismatch for patient '{patient}': normal "
                f"'{normal['sample']}' is '{normal['sex']}', tumor "
                f"'{tumor['sample']}' is '{tumor['sex']}'"
            )

    if tumors.empty:
        raise ValueError(
            "Samplesheet produced zero tumor/normal pairs -- check that the input "
            "dataset has at least one status=1 (tumor) sample with a matching "
            "status=0 (normal) sample for the same patient"
        )

    return manifest.reindex(columns=COLUMNS)


if __name__ == "__main__":
    from cirro.helpers.preprocess_dataset import PreprocessDataset

    ds = PreprocessDataset.from_running()
    manifest = build_samplesheet(ds.files["file"], ds.samplesheet)
    ds.logger.info(manifest.to_csv(index=False))
    manifest.to_csv(ds.params["input"], index=False)
    ds.logger.info(f"Wrote {manifest.shape[0]} row(s) to {ds.params['input']}")
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m unittest .cirro.test_preprocess -v` (same venv as Step 2 if one was
needed)
Expected: `OK`, all 9 tests passing.

- [ ] **Step 5: Commit**

```bash
git add .cirro/preprocess.py .cirro/test_preprocess.py
git commit -m "Add preprocess.py: Cirro samplesheet assembly + pairing validation"
```

---

### Task 3: Documentation

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes: the finished state of Tasks 1-2 (this task only documents, no code).

- [ ] **Step 1: Add a "Cirro deployment" section to `README.md`**

Insert after the existing "Testing" section, before "Known gaps in this phase":

```markdown
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

**Manual follow-up, not automated by any file in this repo**: once this process is
registered in Cirro (process id assigned via Cirro's admin UI -- this org's convention
for its own custom pipelines does not commit a `process-definition.json`, unlike the
public nf-core forks), add that process id to `parabricks_bqsr`'s own `childProcessIds`
list so a user can launch this pipeline directly from a completed `parabricks_bqsr`
dataset in Cirro's UI. This is a one-line edit in the *other* repo's own Cirro
registration, not something either repo's committed files accomplish on their own.
```

- [ ] **Step 2: Commit**

```bash
git add README.md
git commit -m "Document Cirro deployment (Phase 3)"
```

---

## Self-Review Notes (completed during plan authoring, not a task)

- **Spec coverage**: §13.1 (parent process) → documented in Task 3, enforced by Task 2's
  `RECAL_FILE` regex only matching the `recalibrated` stage. §13.2 (preprocess.py) →
  Task 2. §13.3 (GPU compute-config) → Task 1 Step 2. §13.4 (thin form) → Task 1 Step 1.
  §13.5 (repo layout, corrected) → Tasks 1-2's file list matches it exactly, no
  `process-definition.json`.
- **Placeholder scan**: every file's content above is complete, real code or config —
  none are guessed; the `process-input.json`/`process-compute.config` content is copied
  directly from a real, inspected sibling-repo file, and `preprocess.py`'s validation
  logic is a direct line-by-line port of `SAMPLESHEET_TO_PAIRS`'s own Groovy rules,
  re-read from the current file during this plan's authoring, not from memory.
- **Type consistency**: `build_samplesheet`'s return columns
  (`patient,sample,sex,status,cram,crai`) match `SAMPLESHEET_TO_PAIRS`'s own expected
  CSV header exactly (re-verified against the current file's `take:` comment during
  planning). `status` is emitted as Python `int` (0/1); confirmed empirically (not
  assumed) that `pandas.DataFrame.to_csv()` serializes an `int`-typed column as plain
  `"0"`/`"1"` digit strings, no decimal point — matching
  `validateSamplesheetRow()`'s exact-string check `row.status in ['0', '1']` exactly.
  A `float`-typed status column would have serialized as `"0.0"`/`"1.0"` and silently
  failed that check; the explicit `.astype(int)` cast in Step 3's code is what prevents
  this, and is now verified necessary, not just defensive.
