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
