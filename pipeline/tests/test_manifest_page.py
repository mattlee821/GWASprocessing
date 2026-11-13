"""manifest.html in the work root (pipeline/bin/registry.py): every GWAS this install has processed or tried to, 1 row
per output folder, with its manifest fields, status, counts, QC flags, warnings and relative links.

1 test; each README rule is a subTest. Synthetic data only, no network.
"""

from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from testlib import run_script

import registry

STATE_COLUMNS = ["row_hash", "study_id", "output_id", "gwas_location", "standard_vcf", "standard_index", "metadata_json",
                 "qc_json", "report_html", "settings_hash", "updated_at"]
RUN_COLUMNS = ["GWAS_location", "phenotype", "ancestry", "author", "year", "PMID", "sex", "source_type", "study_yaml",
               "study_id", "output_id", "row_hash", "standardise_dir"]
LOCKE = ("locke_2015_25673413", "GWAScatalog_GCST002783")
OMARA = ("omara_2018_30093612", "GWAScatalog_GCST006464")


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


class Install:
    """A work root and a log root, with helpers to write outputs, the state file and the last run's logs."""

    def __init__(self, base: Path):
        self.work = base / "GWAS"
        self.logs = base / "logs"
        self.logs.mkdir(parents=True)
        self.work.mkdir()

    def folder(self, study: str, output_id: str) -> Path:
        return self.work / study / "standardise" / output_id

    def outputs(self, study: str, output_id: str) -> Path:
        folder = self.folder(study, output_id)
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "standardised.gwas.vcf.gz").write_bytes(b"x")
        (folder / "standardised.report.html").write_text("<html></html>", encoding="utf-8")
        (folder / "standardised.metadata.json").write_text(json.dumps(
            {"phenotype": "BMI", "ancestry": "ALL", "sex": "combined", "author": "Locke", "year": "2015", "PMID": "25673413",
             "GWAS_location": "https://www.ebi.ac.uk/gwas/studies/GCST002783"}), encoding="utf-8")
        (folder / "standardised.qc.json").write_text(json.dumps(
            {"input_rows": 2555510, "output_rows": 2554903, "lambda_gc": 1.234, "completed_at": "2026-11-12T10:00:00+00:00",
             "coordinate_source": "dbsnp", "qc_flags": {"n_p_sig": "more than 1,000 variants at P < 5e-8"},
             "warnings": ["weak frequency correlation"],
             "sample_header": {"StudyType": "Continuous", "TotalControls": "236781", "AFFilled": "54408"}}), encoding="utf-8")
        return folder

    def state(self, rows: list[dict]) -> None:
        write_tsv(self.work / registry.STATE, STATE_COLUMNS, rows)

    def run_logs(self, rows: list[dict], completed: list[str]) -> None:
        write_tsv(self.logs / "run_manifest.tsv", RUN_COLUMNS, rows)
        (self.logs / "run_manifest_summary.json").write_text(json.dumps({"rows_to_run": len(rows)}), encoding="utf-8")
        (self.logs / "state_update_summary.json").write_text(json.dumps({"updated_row_hashes": completed}), encoding="utf-8")

    def page(self, logs: bool = True) -> str:
        return registry.write_page(self.work, self.logs if logs else None).read_text(encoding="utf-8")


def state_row(row_hash: str, folder: Path, updated_at: str = "2026-11-12T10:00:00+00:00") -> dict:
    study, output_id = folder.parts[-3], folder.parts[-1]
    return {"row_hash": row_hash, "study_id": study, "output_id": output_id,
            "standard_vcf": str(folder / "standardised.gwas.vcf.gz"), "updated_at": updated_at}


def run_row(row_hash: str, folder: Path, **fields: str) -> dict:
    return {"row_hash": row_hash, "study_id": folder.parts[-3], "output_id": folder.parts[-1],
            "standardise_dir": str(folder), **fields}


class ManifestPageTest(unittest.TestCase):
    def test_manifest_page(self):
        base = Path(tempfile.mkdtemp(prefix="gwas_test_manifest_page_"))
        counter = iter(range(1000))

        def install() -> Install:
            return Install(base / f"install{next(counter)}")

        with self.subTest("completed and failed rows are listed with links relative to the page"):
            i = install()
            done = i.outputs(*LOCKE)
            i.state([state_row("a1", done)])
            failed = i.folder(*OMARA)
            i.run_logs([run_row("b2", failed, phenotype="endometrial cancer", ancestry="EUR", sex="female", author="OMara",
                                year="2018", PMID="30093612", GWAS_location="https://www.ebi.ac.uk/gwas/studies/GCST006464")],
                       completed=[])
            page = i.page()
            self.assertTrue((i.work / "manifest.html").exists())
            self.assertIn("1 complete, 1 failed", page)
            prefix = "locke_2015_25673413/standardise/GWAScatalog_GCST002783"
            for target in (prefix, f"{prefix}/standardised.report.html", f"{prefix}/standardised.gwas.vcf.gz"):
                self.assertIn(f'href="{target}"', page)
            table = page.split("Work root:")[1].split("</p>")[1]
            self.assertNotIn(str(i.work), table, "no absolute paths in the table")

        with self.subTest("each row shows its manifest fields, status, variants in and out, sample size, QC flags and warnings"):
            self.assertIn('<span class="status complete">complete</span>', page)
            self.assertIn('<span class="status failed">failed</span>', page)
            for text in ("endometrial cancer", "OMara 2018", "female", "https://www.ebi.ac.uk/gwas/studies/GCST006464",
                         "2,554,903 of 2,555,510 (99.98%)", "N 236,781", "54,408", "1.234",
                         "more than 1,000 variants at P &lt; 5e-8", "weak frequency correlation", "rsIDs",
                         'href="https://pubmed.ncbi.nlm.nih.gov/25673413/"', 'href="https://pubmed.ncbi.nlm.nih.gov/30093612/"'):
                self.assertIn(text, page)

        with self.subTest("failed rows stay listed until a later run completes them"):
            i = install()
            failed = i.folder(*OMARA)
            i.run_logs([run_row("b2", failed)], completed=[])
            i.page()
            for name in ("run_manifest.tsv", "run_manifest_summary.json", "state_update_summary.json"):
                (i.logs / name).unlink()  # a later run that did not touch this row
            self.assertIn("0 complete, 1 failed", i.page())
            i.outputs(*OMARA)
            i.run_logs([run_row("b2", failed)], completed=["b2"])
            self.assertIn("1 complete, 0 failed", i.page())

        with self.subTest("1 row per output folder, even when its manifest row changed"):
            i = install()
            done = i.outputs(*LOCKE)
            i.state([state_row("old", done, "2026-11-11T10:00:00+00:00"), state_row("new", done, "2026-11-12T10:00:00+00:00")])
            page = i.page(logs=False)
            self.assertEqual(page.count("<tr class="), 1)
            self.assertIn("1 complete, 0 failed", page)

        with self.subTest("a completed row whose outputs were deleted is marked"):
            i = install()
            done = i.outputs(*LOCKE)
            i.state([state_row("a1", done)])
            (done / "standardised.gwas.vcf.gz").unlink()
            page = i.page(logs=False)
            self.assertIn("0 complete, 0 failed, 1 with outputs missing", page)
            self.assertIn('<span class="status outputs-missing">outputs missing</span>', page)

        with self.subTest("registry.py writes the page from the command line"):
            i = install()
            i.state([state_row("a1", i.outputs(*LOCKE))])
            proc = run_script("registry.py", "--work-root", i.work, "--log-root", i.logs)
            self.assertIn(str(i.work / "manifest.html"), proc.stdout)
            self.assertIn("1 complete", (i.work / "manifest.html").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
