"""manifest.html: 1 row per GWAS output folder, status kept across runs, relative links, counts from the QC JSON."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

import testlib  # noqa: F401  (puts pipeline/bin on sys.path)

import registry

STATE_COLUMNS = ["row_hash", "study_id", "output_id", "gwas_location", "standard_vcf", "standard_index", "metadata_json",
                 "qc_json", "report_html", "settings_hash", "updated_at"]
RUN_COLUMNS = ["GWAS_location", "phenotype", "ancestry", "author", "year", "PMID", "sex", "source_type", "study_yaml",
               "study_id", "output_id", "row_hash", "standardise_dir"]


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


class RegistryTests(unittest.TestCase):
    def setUp(self):
        self.work = Path(tempfile.mkdtemp()) / "GWAS"
        self.logs = self.work.parent / "logs"
        self.logs.mkdir(parents=True)
        self.done = self.output("locke_2015_25673413", "GWAScatalog_GCST002783")
        self.failed = self.work / "omara_2018_30093612" / "standardise" / "GWAScatalog_GCST006464"

    def output(self, study: str, output_id: str) -> Path:
        folder = self.work / study / "standardise" / output_id
        folder.mkdir(parents=True)
        (folder / "standardised.gwas.vcf.gz").write_bytes(b"x")
        (folder / "standardised.report.html").write_text("<html></html>", encoding="utf-8")
        (folder / "standardised.metadata.json").write_text(json.dumps(
            {"phenotype": "BMI", "ancestry": "EUR", "author": "Locke", "year": "2015", "PMID": "25673413",
             "GWAS_location": "https://www.ebi.ac.uk/gwas/studies/GCST002783"}), encoding="utf-8")
        (folder / "standardised.qc.json").write_text(json.dumps(
            {"input_rows": 2555510, "output_rows": 2554903, "lambda_gc": 1.234, "completed_at": "2026-11-12T10:00:00+00:00",
             "coordinate_source": "dbsnp", "qc_flags": {"n_p_sig": "more than 1,000 variants at P < 5e-8"}, "warnings": [],
             "sample_header": {"StudyType": "Continuous", "TotalControls": "236781", "AFFilled": "54408"}}), encoding="utf-8")
        return folder

    def state(self, rows: list[dict]) -> None:
        write_tsv(self.work / registry.STATE, STATE_COLUMNS, rows)

    def run_logs(self, rows: list[dict], completed: list[str]) -> None:
        write_tsv(self.logs / "run_manifest.tsv", RUN_COLUMNS, rows)
        (self.logs / "run_manifest_summary.json").write_text(json.dumps({"rows_to_run": len(rows)}), encoding="utf-8")
        (self.logs / "state_update_summary.json").write_text(json.dumps({"updated_row_hashes": completed}), encoding="utf-8")

    def test_page_lists_completed_and_failed_rows_with_relative_links(self):
        vcf = str(self.done / "standardised.gwas.vcf.gz")
        self.state([{"row_hash": "a1", "study_id": "locke_2015_25673413", "output_id": "GWAScatalog_GCST002783",
                     "standard_vcf": vcf, "updated_at": "2026-11-12T10:00:00+00:00"}])
        self.run_logs([{"row_hash": "b2", "study_id": "omara_2018_30093612", "output_id": "GWAScatalog_GCST006464",
                        "phenotype": "endometrial cancer", "ancestry": "EUR", "standardise_dir": str(self.failed)}], completed=[])
        page = registry.write_page(self.work, self.logs).read_text(encoding="utf-8")
        self.assertIn("1 complete, 1 failed", page)
        self.assertIn('href="locke_2015_25673413/standardise/GWAScatalog_GCST002783/standardised.report.html"', page)
        self.assertIn("2,554,903 of 2,555,510 (99.98%)", page)
        self.assertIn("more than 1,000 variants at P &lt; 5e-8", page)
        self.assertIn("rsIDs", page)  # build column for rsID-only inputs
        self.assertIn('href="https://pubmed.ncbi.nlm.nih.gov/25673413/"', page)
        self.assertIn("endometrial cancer", page)
        self.assertNotIn(str(self.work), page.split("Work root:")[1].split("</p>")[1])  # no absolute paths in the table

    def test_failed_rows_are_remembered_and_cleared_by_a_later_completion(self):
        failed_row = {"row_hash": "b2", "study_id": "omara_2018_30093612", "output_id": "GWAScatalog_GCST006464",
                      "standardise_dir": str(self.failed)}
        self.run_logs([failed_row], completed=[])
        registry.write_page(self.work, self.logs)
        (self.logs / "run_manifest.tsv").unlink()  # a later run that did not touch this row
        page = registry.write_page(self.work, self.logs).read_text(encoding="utf-8")
        self.assertIn("0 complete, 1 failed", page)
        self.output("omara_2018_30093612", "GWAScatalog_GCST006464")
        self.run_logs([failed_row], completed=["b2"])
        page = registry.write_page(self.work, self.logs).read_text(encoding="utf-8")
        self.assertIn("1 complete, 0 failed", page)

    def test_1_row_per_output_folder_when_its_manifest_row_changed(self):
        vcf = str(self.done / "standardised.gwas.vcf.gz")
        self.state([{"row_hash": "old", "study_id": "locke_2015_25673413", "output_id": "GWAScatalog_GCST002783",
                     "standard_vcf": vcf, "updated_at": "2026-11-11T10:00:00+00:00"},
                    {"row_hash": "new", "study_id": "locke_2015_25673413", "output_id": "GWAScatalog_GCST002783",
                     "standard_vcf": vcf, "updated_at": "2026-11-12T10:00:00+00:00"}])
        page = registry.write_page(self.work).read_text(encoding="utf-8")
        self.assertEqual(page.count("<tr class="), 1)

    def test_a_completed_row_whose_outputs_were_deleted_is_marked(self):
        vcf = self.done / "standardised.gwas.vcf.gz"
        self.state([{"row_hash": "a1", "study_id": "locke_2015_25673413", "output_id": "GWAScatalog_GCST002783",
                     "standard_vcf": str(vcf), "updated_at": "2026-11-12T10:00:00+00:00"}])
        vcf.unlink()
        page = registry.write_page(self.work).read_text(encoding="utf-8")
        self.assertIn("1 with outputs missing", page)


if __name__ == "__main__":
    unittest.main()
