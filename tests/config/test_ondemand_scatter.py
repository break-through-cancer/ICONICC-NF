"""Offline check that the JOINT_GENOTYPING scatter steps run on-demand, not Spot.
Run with `python3 tests/config/test_ondemand_scatter.py -v`.
"""
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]

ONDEMAND_PROCESSES = ("GATK4_GENOMICSDBIMPORT", "GATK4_GENOTYPEGVCFS")


def with_name_block(text, name):
    for selector, body in re.findall(r"withName: '([^']*)' \{(.*?)\n    \}", text, re.S):
        if selector == name:
            return body
    raise AssertionError(f"no withName: '{name}' block found")


class OnDemandScatterTests(unittest.TestCase):

    def test_scatter_steps_use_ondemand_queue(self):
        text = (ROOT / ".cirro" / "process-compute.config").read_text()
        for name in ONDEMAND_PROCESSES:
            block = with_name_block(text, name)
            self.assertIn("PW_ONDEMAND_JOB_QUEUE", block)


if __name__ == "__main__":
    unittest.main()
