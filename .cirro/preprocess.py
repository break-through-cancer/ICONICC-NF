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
