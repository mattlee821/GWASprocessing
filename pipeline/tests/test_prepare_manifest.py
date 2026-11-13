"""Tests for prepare_manifest.py: row derivation, validation, per-row JSON files and settings-aware skipping."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

from testlib import EXAMPLES, example_reference, run_script

MANIFEST_HEADER = "GWAS_location\tphenotype\tancestry\tauthor\tyear\tPMID\n"


class PrepareManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.raw = EXAMPLES / "tables/generic.tsv"
        self.config_root = self.tmp / "manifests"
        self.config_root.mkdir()

    def manifest(self, body, header=MANIFEST_HEADER):
        path = self.config_root / "manifest.tsv"
        path.write_text(header + body, encoding="utf-8")
        return path

    def prepare(self, manifest, *extra, check=True):
        args = ["--manifest", manifest, "--work-root", self.tmp / "work", "--state-file", self.tmp / ".standardise_state.tsv",
                "--output", self.tmp / "run.tsv", "--summary", self.tmp / "summary.json", "--rows-dir", self.tmp / "rows",
                "--config-root", self.config_root, "--reference-root", example_reference(), *extra]
        return run_script("prepare_manifest.py", *args, check=check)

    def rows(self):
        with open(self.tmp / "run.tsv", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def row_jsons(self):
        return [json.loads(p.read_text(encoding="utf-8")) for p in sorted((self.tmp / "rows").glob("*.json"))]

    def summary(self):
        return json.loads((self.tmp / "summary.json").read_text(encoding="utf-8"))

    def complete_outputs(self, row, settings_hash):
        standardise_dir = Path(row["standardise_dir"])
        standardise_dir.mkdir(parents=True, exist_ok=True)
        files = {key: standardise_dir / name for key, name in [
            ("standard_vcf", "standardised.gwas.vcf.gz"), ("standard_index", "standardised.gwas.vcf.gz.tbi"),
            ("metadata_json", "standardised.metadata.json"), ("qc_json", "standardised.qc.json"), ("report_html", "standardised.report.html")]}
        for path in files.values():
            path.write_text("ok", encoding="utf-8")
        columns = ["row_hash", "study_id", "output_id", "gwas_location", "standard_vcf", "standard_index", "metadata_json",
                   "qc_json", "report_html", "settings_hash", "updated_at"]
        values = [row["row_hash"], row["study_id"], row["output_id"], row["GWAS_location"], *map(str, files.values()), settings_hash, "now"]
        (self.tmp / ".standardise_state.tsv").write_text("\t".join(columns) + "\n" + "\t".join(values) + "\n", encoding="utf-8")

    def test_prepare_manifest_derives_ids_and_skip_state(self):
        manifest = self.manifest(f"{self.raw}\tBMI\tall\tLocke\t2015\t25673413\n")
        self.prepare(manifest)
        rows = self.rows()
        self.assertEqual(len(rows), 1)
        self.assertEqual((rows[0]["study_id"], rows[0]["source_type"], rows[0]["sex"]), ("locke_2015_25673413", "local", ""))
        self.assertTrue(rows[0]["output_id"].startswith("local_generic_"))
        self.assertTrue(rows[0]["output_prefix"].endswith(f"standardise/{rows[0]['output_id']}/standardised"))
        settings_hash = self.row_jsons()[0]["settings_hash"]
        self.complete_outputs(rows[0], settings_hash)
        self.prepare(manifest)
        self.assertEqual(self.rows(), [])
        self.assertEqual(self.row_jsons(), [])
        self.assertEqual(self.summary()["rows_skipped_complete"], 1)

    def test_changed_settings_make_a_completed_row_eligible_again(self):
        manifest = self.manifest(f"{self.raw}\tBMI\tall\tLocke\t2015\t25673413\n")
        self.prepare(manifest)
        row, settings_hash = self.rows()[0], self.row_jsons()[0]["settings_hash"]
        self.complete_outputs(row, settings_hash)
        self.prepare(manifest)
        self.assertEqual(len(self.rows()), 0, "unchanged settings must not regenerate a complete row")
        (self.config_root / "locke_2015_25673413.yaml").write_text("input_build: hg19\n", encoding="utf-8")
        self.prepare(manifest)
        self.assertEqual(len(self.rows()), 1, "a new study YAML must regenerate the row")

    def test_rows_are_written_as_json_files_carrying_every_field(self):
        manifest = self.manifest(f'{self.raw}\tBMI "adjusted" $HOME `x`\tEUR\tLocke\t2015\t25673413\n')
        self.prepare(manifest)
        (row,) = self.row_jsons()
        self.assertEqual(row["phenotype"], 'BMI "adjusted" $HOME `x`')
        self.assertEqual(row["GWAS_location"], str(self.raw))
        self.assertEqual(len(row["settings_hash"]), 16)

    def test_manifest_optional_fields_are_preserved_and_hashed(self):
        header = "GWAS_location\tphenotype\tancestry\tsex\tauthor\tyear\tPMID\tsource_type\tinput_build\tpopulation\tdelimiter\tformat\tstudy_yaml\n"
        (self.config_root / "custom.yaml").write_text("description: custom\n", encoding="utf-8")
        manifest = self.manifest(f"{self.raw}\tBMI\tall\tfemale\tLocke\t2015\t25673413\tlocal\thg19\tALL\t\\s+\ttable\tcustom.yaml\n", header)
        self.prepare(manifest)
        row = self.rows()[0]
        self.assertEqual((row["sex"], row["input_build"], row["population"], row["delimiter"], row["format"], row["study_yaml"]),
                         ("female", "hg19", "ALL", "\\s+", "table", "custom.yaml"))
        first_hash = row["row_hash"]
        self.manifest(f"{self.raw}\tBMI\tall\tLocke\t2015\t25673413\n")
        self.prepare(manifest)
        self.assertNotEqual(first_hash, self.rows()[0]["row_hash"])

    def test_prepare_manifest_row_limit_defers_remaining_rows(self):
        manifest = self.manifest(f"{self.raw}\tBMI\tall\tLocke\t2015\t25673413\n{self.raw}\tWHR\tall\tPulit\t2018\t30239722\n")
        self.prepare(manifest, "--row-limit", "1")
        summary = self.summary()
        self.assertEqual((len(self.rows()), len(self.row_jsons())), (1, 1))
        self.assertEqual((summary["rows_eligible"], summary["rows_to_run"], summary["rows_deferred_due_to_limit"]), (2, 1, 1))

    def test_only_study_filters_rows(self):
        manifest = self.manifest(f"{self.raw}\tBMI\tall\tLocke\t2015\t25673413\n{self.raw}\tWHR\tall\tPulit\t2018\t30239722\n")
        self.prepare(manifest, "--only-study", "pulit_2018_30239722")
        self.assertEqual([r["study_id"] for r in self.rows()], ["pulit_2018_30239722"])
        self.assertEqual(self.summary()["rows_filtered"], 1)

    def test_invalid_ancestry_stops_the_run_before_any_row_starts(self):
        manifest = self.manifest(f"{self.raw}\tBMI\tEUR\tLocke\t2015\t25673413\n{self.raw}\tWHR\tunknown\tPulit\t2018\t30239722\n")
        proc = self.prepare(manifest, check=False)
        self.assertEqual(proc.returncode, 78)
        self.assertIn("row 3", proc.stderr)
        self.assertIn("ancestry", proc.stderr)

    def test_sex_and_format_must_be_known_values(self):
        header = "GWAS_location\tphenotype\tancestry\tsex\tauthor\tyear\tPMID\tformat\n"
        manifest = self.manifest(f"{self.raw}\tBMI\tEUR\tFemale\tLocke\t2015\t25673413\tvcf\n"
                                 f"{self.raw}\tWHR\tEUR\twomen\tPulit\t2018\t30239722\tparquet\n", header)
        proc = self.prepare(manifest, check=False)
        self.assertEqual(proc.returncode, 78)
        self.assertIn("row 3: sex 'women' must be combined, male or female", proc.stderr)
        self.assertIn("row 3: format 'parquet' must be table or vcf", proc.stderr)
        self.assertNotIn("row 2", proc.stderr)

    def test_two_rows_writing_the_same_outputs_are_rejected(self):
        manifest = self.manifest("ieu-a-2\tBMI\tEUR\tLocke\t2015\t25673413\nieu-a-2\tBMI (copy)\tEUR\tLocke\t2015\t25673413\n")
        proc = self.prepare(manifest, check=False)
        self.assertEqual(proc.returncode, 78)
        self.assertIn("same output directory", proc.stderr)

    def test_synapse_rows_are_rejected_up_front(self):
        proc = self.prepare(self.manifest("syn26123456\tBMI\tEUR\tLocke\t2015\t25673413\n"), check=False)
        self.assertEqual(proc.returncode, 78)
        self.assertIn("Synapse", proc.stderr)

    def test_relative_local_paths_resolve_against_the_manifest_directory(self):
        (self.config_root / "data").mkdir()
        local = self.config_root / "data" / "gwas.tsv"
        local.write_text(self.raw.read_text(encoding="utf-8"), encoding="utf-8")
        self.prepare(self.manifest("data/gwas.tsv\tBMI\tEUR\tLocke\t2015\t25673413\n"))
        self.assertEqual(self.row_jsons()[0]["GWAS_location"], str(local.resolve()))


if __name__ == "__main__":
    unittest.main()
