"""Tests for update_state.py: completed rows are recorded with their settings hash; others are reported."""

import json
import tempfile
import unittest
from pathlib import Path

from testlib import run_script

import gwas_pipeline_common as common


class UpdateStateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.state = self.tmp / ".standardise_state.tsv"

    def record(self, name, with_plot=True):
        out = self.tmp / name
        out.mkdir()
        paths = {key: out / f"{key}.out" for key in ["standard_vcf", "standard_index", "metadata_json", "qc_json", "report_html"]}
        for key, path in paths.items():
            if key != "report_html" or with_plot:
                path.write_text("x", encoding="utf-8")
        record = {"row_hash": f"hash_{name}", "study_id": "s", "output_id": name, "GWAS_location": "loc",
                  **{k: str(v) for k, v in paths.items()}, "qc": {"settings_hash": f"settings_{name}"}}
        path = self.tmp / f"{name}.row_record.json"
        path.write_text(json.dumps(record), encoding="utf-8")
        return path

    def update(self, *records, report="true"):
        summary = self.tmp / "summary.json"
        run_script("update_state.py", "--state-file", self.state, "--report", report, "--output", summary,
                   *(["--records", *records] if records else []))
        return json.loads(summary.read_text(encoding="utf-8"))

    def test_complete_records_are_saved_with_their_settings_hash(self):
        summary = self.update(self.record("a"))
        self.assertEqual((summary["updated"], summary["updated_row_hashes"]), (1, ["hash_a"]))
        self.assertEqual(common.read_state(self.state)["hash_a"]["settings_hash"], "settings_a")

    def test_incomplete_and_missing_records_are_ignored_with_reasons(self):
        summary = self.update(self.record("noplot", with_plot=False), self.tmp / "absent.row_record.json")
        self.assertEqual((summary["updated"], summary["ignored"]), (0, 2))
        reasons = {item["record"]: item["reason"] for item in summary["ignored_records"]}
        self.assertEqual(reasons[str(self.tmp / "absent.row_record.json")], "record file not found")
        self.assertEqual(reasons[str(self.tmp / "noplot.row_record.json")], "outputs incomplete")

    def test_existing_rows_are_kept(self):
        self.update(self.record("a"))
        self.update(self.record("b"))
        self.assertEqual(sorted(common.read_state(self.state)), ["hash_a", "hash_b"])

    def test_no_records_still_writes_a_summary(self):
        summary = self.update()
        self.assertEqual((summary["updated"], summary["ignored"]), (0, 0))


if __name__ == "__main__":
    unittest.main()
