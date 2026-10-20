"""Offline tests for ICONICC-NF's Cirro preprocess script.

Runs without the Cirro SDK or boto3: both are stubbed before preprocess is loaded,
matching the pattern established in CirroBio/Cirro-pipelines'
nf-core/sarek_call_variants/3/test_preprocess.py. Run with
`python3 .cirro/test_preprocess.py -v` -- `.cirro/` is not a Python package, so
`python3 -m unittest` (bare, or as `python3 -m unittest .cirro.test_preprocess`) does not
discover or run these tests.
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
        # PT01's sample has no .crai alongside it, so it must be dropped by
        # alignment_pairs() -- but PT02's pair is fully indexed and valid, so the
        # function should still succeed and return only PT02's pair, rather than
        # raising for the unrelated "zero pairs at all" reason.
        paths = (
            ["preprocessing/recalibrated/PT01_normal/PT01_normal.recal.cram"]
            + recal_paths("PT02_normal")
            + recal_paths("PT02_tumor")
        )
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "female", "Normal"),
            ("PT02_normal", "PT02", "female", "Normal"),
            ("PT02_tumor", "PT02", "female", "Tumor"),
        ])
        result = preprocess.build_samplesheet(paths, sheet)
        self.assertEqual(sorted(result["sample"]), ["PT02_normal", "PT02_tumor"])
        self.assertNotIn("PT01_normal", result["sample"].tolist())

    def test_status_accepts_string_zero_one_and_uppercase(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "female", "0"),
            ("PT01_tumor", "PT01", "female", "1"),
        ])
        result = preprocess.build_samplesheet(paths, sheet)
        self.assertEqual(
            result.loc[result["sample"] == "PT01_normal", "status"].iloc[0], 0
        )
        self.assertEqual(
            result.loc[result["sample"] == "PT01_tumor", "status"].iloc[0], 1
        )

        paths2 = recal_paths("PT02_normal") + recal_paths("PT02_tumor")
        sheet2 = samplesheet_df([
            ("PT02_normal", "PT02", "female", "NORMAL"),
            ("PT02_tumor", "PT02", "female", "TUMOR"),
        ])
        result2 = preprocess.build_samplesheet(paths2, sheet2)
        self.assertEqual(
            result2.loc[result2["sample"] == "PT02_normal", "status"].iloc[0], 0
        )
        self.assertEqual(
            result2.loc[result2["sample"] == "PT02_tumor", "status"].iloc[0], 1
        )

    def test_bad_status_is_still_hard_error(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "female", "unknown"),
            ("PT01_tumor", "PT01", "female", "Tumor"),
        ])
        with self.assertRaisesRegex(ValueError, "status must be Normal/0 or Tumor/1"):
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

    def test_patient_id_column_is_accepted(self):
        # Real Cirro sample metadata carried over from a fastq ingestion samplesheet
        # uses "patient_id", not the project schema's "patient".
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = pd.DataFrame([
            dict(sample="PT01_normal", patient_id="PT01", sex="female", status="Normal"),
            dict(sample="PT01_tumor", patient_id="PT01", sex="female", status="Tumor"),
        ])
        result = preprocess.build_samplesheet(paths, sheet)
        self.assertEqual(sorted(result["patient"].unique()), ["PT01"])

    def test_patient_column_preferred_over_patient_id(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = pd.DataFrame([
            dict(sample="PT01_normal", patient="PT01", patient_id="wrong", sex="female", status="Normal"),
            dict(sample="PT01_tumor", patient="PT01", patient_id="wrong", sex="female", status="Tumor"),
        ])
        result = preprocess.build_samplesheet(paths, sheet)
        self.assertEqual(sorted(result["patient"].unique()), ["PT01"])

    def test_sex_accepts_xx_xy_karyotype(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "XY", "Normal"),
            ("PT01_tumor", "PT01", "xy", "Tumor"),
        ])
        result = preprocess.build_samplesheet(paths, sheet)
        self.assertEqual(sorted(result["sex"].unique()), ["male"])

    def test_bad_sex_is_still_hard_error(self):
        paths = recal_paths("PT01_normal") + recal_paths("PT01_tumor")
        sheet = samplesheet_df([
            ("PT01_normal", "PT01", "unknown", "Normal"),
            ("PT01_tumor", "PT01", "female", "Tumor"),
        ])
        with self.assertRaisesRegex(ValueError, "sex must be 'male'/'female' or 'XX'/'XY'"):
            preprocess.build_samplesheet(paths, sheet)


class ProcessFormTests(unittest.TestCase):
    """The Cirro launch form/input mapping, not preprocess.py -- same file as
    break-through-cancer/parabricks-fq2bam-bqsr's own .cirro/align/test_align_preprocess.py,
    which mixes the two for the same reason (one process, one Cirro folder)."""

    def test_haplotypecaller_memory_override_is_offered_and_mapped(self):
        import json
        here = Path(__file__).parent
        field = json.loads((here / "process-form.json").read_text())["form"]["properties"]["haplotypecaller_memory_gb"]
        self.assertEqual(field["type"], "integer")
        self.assertNotIn("default", field)
        mapping = json.loads((here / "process-input.json").read_text())
        self.assertEqual(mapping["haplotypecaller_memory_gb"], "$.dataset.params.haplotypecaller_memory_gb")


if __name__ == "__main__":
    unittest.main()
