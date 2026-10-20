"""Offline checks that OOM-killed GATK/Parabricks tasks retry with more memory, not the
same fixed ceiling. Mirrors break-through-cancer/parabricks-fq2bam-bqsr's own
tests/config/test_oom_retry.py -- same convention, same repo family. Run with
`python3 tests/config/test_oom_retry.py -v`.
"""
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]

MEMORY_SCALED_PROCESSES = (
    "GATK4_GENOMICSDBIMPORT", "GATK4_GENOTYPEGVCFS",
    "VARIANTRECALIBRATOR_SNP", "VARIANTRECALIBRATOR_INDEL",
    "APPLYVQSR_INDEL", "APPLYVQSR_SNP",
    "PARABRICKS_HAPLOTYPECALLER",
)


def with_name_block(text, name):
    for selector, body in re.findall(r"withName: '([^']*)' \{(.*?)\n    \}", text, re.S):
        if selector == name:
            return body
    raise AssertionError(f"no withName: '{name}' block found")


class OomRetryTests(unittest.TestCase):

    def test_memory_grows_with_each_attempt(self):
        text = (ROOT / "nextflow.config").read_text()
        for name in MEMORY_SCALED_PROCESSES:
            block = with_name_block(text, name)
            memory = [l for l in block.splitlines() if l.strip().startswith("memory")]
            self.assertTrue(memory, f"{name}: no memory directive found")
            self.assertIn("task.attempt", memory[0], f"{name}: {memory[0].strip()}")

    def test_haplotypecaller_memory_override_param_exists_and_is_wired(self):
        text = (ROOT / "nextflow.config").read_text()
        self.assertIn("params.haplotypecaller_memory_gb", text)
        block = with_name_block(text, "PARABRICKS_HAPLOTYPECALLER")
        memory = [l for l in block.splitlines() if l.strip().startswith("memory")][0]
        self.assertIn("params.haplotypecaller_memory_gb", memory)

    def test_cirro_runs_haplotypecaller_on_demand_and_retries_memory_kills(self):
        # AWS Batch reports a container memory kill without an exit code; Nextflow
        # records Integer.MAX_VALUE for it.
        text = (ROOT / ".cirro" / "process-compute.config").read_text()
        block = with_name_block(text, "PARABRICKS_HAPLOTYPECALLER")
        self.assertIn("PW_ONDEMAND_JOB_QUEUE", block)
        self.assertIn("2147483647", block)


if __name__ == "__main__":
    unittest.main()
