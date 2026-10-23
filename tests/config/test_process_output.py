"""Offline check that .cirro/process-output.json's dashboard-visualization commands
target the real published paths/columns for segment.R's two main tabular outputs.
Run with `python3 tests/config/test_process_output.py -v`.
"""
import json
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


def load():
    return json.loads((ROOT / ".cirro" / "process-output.json").read_text())


class ProcessOutputTests(unittest.TestCase):

    def test_seg_txt_command(self):
        data = load()
        cmd = next(c for c in data["commands"] if c["params"]["target"] == "segments.parquet")
        self.assertEqual(cmd["command"], "hot.Parquet")
        self.assertEqual(cmd["params"]["source"], "$data_directory/segment/*.seg.txt")
        cols = {c["col"] for c in cmd["params"]["cols"]}
        self.assertEqual(cols, {"ID", "chrom", "loc.start", "loc.end", "num.mark", "seg_nature", "seg.mean"})

    def test_processed_counts_txt_command(self):
        data = load()
        cmd = next(c for c in data["commands"] if c["params"]["target"] == "processed_counts.parquet")
        self.assertEqual(cmd["command"], "hot.Parquet")
        self.assertEqual(cmd["params"]["source"], "$data_directory/segment/*_processed_counts.txt")
        cols = {c["col"] for c in cmd["params"]["cols"]}
        # Confirmed against a real 'segment/*_processed_counts.txt' header (dataset 7f8935b2).
        self.assertEqual(cols, {
            "CONTIG", "START", "END", "loci", "COUNT", "NORMALIZED_COUNT", "AF1", "AF2",
            "snp_count", "C1", "C2", "log2_tangent", "AFMAX", "AFMIN", "AF_ratio",
            "cum_start", "cum_end", "AR_segment", "tangent_segment",
        })


if __name__ == "__main__":
    unittest.main()
