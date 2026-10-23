"""Offline checks for the optional germline-GVCF / joint-genotyped-VCF publishing toggles.
Run with `python3 tests/config/test_vcf_publishing.py -v`.
"""
import json
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]


def with_name_block(text, name):
    for selector, body in re.findall(r"withName: '([^']*)' \{(.*?)\n    \}", text, re.S):
        if selector == name:
            return body
    raise AssertionError(f"no withName: '{name}' block found")


class VcfPublishingTests(unittest.TestCase):

    def test_params_declared_with_expected_defaults(self):
        text = (ROOT / "nextflow.config").read_text()
        self.assertIn("params.publish_germline_gvcf = false", text)
        self.assertIn("params.publish_joint_genotyped_vcf = true", text)

    def test_index_gvcf_publishdir_gated_off_by_default(self):
        text = (ROOT / "nextflow.config").read_text()
        block = with_name_block(text, "INDEX_GVCF")
        self.assertIn("enabled: params.publish_germline_gvcf", block)

    def test_applyvqsr_snp_publishdir_gated_on_by_default(self):
        text = (ROOT / "nextflow.config").read_text()
        block = with_name_block(text, "APPLYVQSR_SNP")
        self.assertIn("enabled: params.publish_joint_genotyped_vcf", block)
        self.assertIn("joint_genotyped_vcf", block)

    def test_cirro_form_exposes_both_checkboxes(self):
        form = json.loads((ROOT / ".cirro" / "process-form.json").read_text())
        props = form["form"]["properties"]
        self.assertEqual(props["publish_germline_gvcf"]["type"], "boolean")
        self.assertEqual(props["publish_germline_gvcf"]["default"], False)
        self.assertEqual(props["publish_joint_genotyped_vcf"]["type"], "boolean")
        self.assertEqual(props["publish_joint_genotyped_vcf"]["default"], True)

    def test_cirro_input_wires_both_params(self):
        process_input = json.loads((ROOT / ".cirro" / "process-input.json").read_text())
        self.assertEqual(process_input["publish_germline_gvcf"], "$.dataset.params.publish_germline_gvcf")
        self.assertEqual(process_input["publish_joint_genotyped_vcf"], "$.dataset.params.publish_joint_genotyped_vcf")


if __name__ == "__main__":
    unittest.main()
