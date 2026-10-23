"""Offline check that JOINT_GENOTYPING's pair-level processes name outputs by
participant_id, not meta.id -- the pair-level meta this subworkflow uses from
GATHER_GENOTYPED_VCF onward never carries an 'id' key (only the internal per-shard
meta does, stripped again before this point), so `meta.id` resolves to null and
produces filenames like "null.recalibrated.vcf.gz". Run with
`python3 tests/config/test_joint_genotyping_prefixes.py -v`.
"""
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]

PAIR_LEVEL_PROCESSES = (
    "VARIANTRECALIBRATOR_SNP", "VARIANTRECALIBRATOR_INDEL",
    "APPLYVQSR_INDEL", "APPLYVQSR_SNP",
)


def with_name_block(text, name):
    for selector, body in re.findall(r"withName: '([^']*)' \{(.*?)\n    \}", text, re.S):
        if selector == name:
            return body
    raise AssertionError(f"no withName: '{name}' block found")


class JointGenotypingPrefixTests(unittest.TestCase):

    def test_pair_level_processes_use_participant_id_not_meta_id(self):
        text = (ROOT / "nextflow.config").read_text()
        for name in PAIR_LEVEL_PROCESSES:
            block = with_name_block(text, name)
            prefix = [l for l in block.splitlines() if l.strip().startswith("ext.prefix")]
            self.assertTrue(prefix, f"{name}: no ext.prefix set, falls back to the module's own meta.id default")
            self.assertIn("meta.participant_id", prefix[0], f"{name}: {prefix[0].strip()}")
            self.assertNotIn("${meta.id}", prefix[0], f"{name}: {prefix[0].strip()}")


if __name__ == "__main__":
    unittest.main()
