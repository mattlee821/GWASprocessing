"""Tests for stage_source.py: each staging route writes the raw file and a stage record from the row JSON."""

import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from testlib import EXAMPLES, LocalHttpServer, run_script, write_stub


class StageSourceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def row(self, location, source_type, **extra):
        row = {"GWAS_location": str(location), "phenotype": 'BMI "adj"', "ancestry": "EUR", "author": "Locke", "year": "2015",
               "PMID": "25673413", "sex": "", "source_type": source_type, "input_build": "", "population": "", "delimiter": "",
               "format": "", "study_yaml": "", "study_id": "locke_2015_25673413", "output_id": "out1", "row_hash": "h1",
               "raw_dir": str(self.tmp / "raw"), "standardise_dir": str(self.tmp / "std"),
               "output_prefix": str(self.tmp / "std" / "standardised"), "settings_hash": "s" * 16}
        row.update(extra)
        path = self.tmp / "row.json"
        path.write_text(json.dumps(row), encoding="utf-8")
        return path

    def stage(self, row_json, *extra, env=None, check=True):
        out = self.tmp / "stage_record.json"
        proc = run_script("stage_source.py", "--row-json", row_json, "--config-root", self.tmp, "--reference-root",
                          EXAMPLES / "reference", "--output", out, *extra, env=env, check=check)
        return proc, (json.loads(out.read_text(encoding="utf-8")) if out.exists() else None)

    def test_local_file_is_copied_and_refreshed_when_the_source_changes(self):
        source = self.tmp / "gwas.tsv"
        source.write_text("SNP\tP\nrs1\t0.1\n", encoding="utf-8")
        _, record = self.stage(self.row(source, "local"))
        self.assertEqual(Path(record["raw_path"]).read_text(encoding="utf-8"), "SNP\tP\nrs1\t0.1\n")
        source.write_text("SNP\tP\nrs1\t0.1\nrs2\t0.2\n", encoding="utf-8")
        _, record = self.stage(self.row(source, "local"))
        self.assertEqual(Path(record["raw_path"]).read_text(encoding="utf-8"), "SNP\tP\nrs1\t0.1\nrs2\t0.2\n")

    def test_direct_url_is_downloaded(self):
        served = self.tmp / "served"
        served.mkdir()
        shutil.copy(EXAMPLES / "tables/generic.tsv", served / "generic.tsv")
        with LocalHttpServer(served) as server:
            _, record = self.stage(self.row(f"{server.url}/generic.tsv", "direct_url"))
        self.assertEqual(Path(record["raw_path"]).read_bytes(), (EXAMPLES / "tables/generic.tsv").read_bytes())

    def test_s3_and_gs_sources_use_the_cloud_cli(self):
        bin_dir = self.tmp / "bin"
        log = self.tmp / "cli.log"
        for tool in ("aws", "gsutil"):
            write_stub(bin_dir, tool, f'echo "{tool} $*" >> "{log}"\ndest="${{@: -1}}"\nprintf "SNP\\tP\\n" > "$dest/gwas.tsv"\n')
        env = {**os.environ, "PATH": f"{bin_dir}{os.pathsep}{os.environ['PATH']}"}
        _, record = self.stage(self.row("s3://bucket/dir/gwas.tsv", "s3"), env=env)
        self.assertTrue(Path(record["raw_path"]).exists())
        _, record = self.stage(self.row("gs://bucket/dir/gwas.tsv", "gs"), env=env)
        self.assertEqual(log.read_text(encoding="utf-8").splitlines(),
                         [f"aws s3 cp s3://bucket/dir/gwas.tsv {self.tmp / 'raw'}", f"gsutil cp gs://bucket/dir/gwas.tsv {self.tmp / 'raw'}"])

    def test_dry_run_records_paths_without_staging(self):
        _, record = self.stage(self.row(self.tmp / "absent.tsv", "local"), "--dry-run", "true")
        self.assertEqual(Path(record["raw_path"]).parent, self.tmp / "raw")
        self.assertFalse(Path(record["raw_path"]).exists())

    def test_record_carries_manifest_fields_and_the_settings_hash(self):
        source = self.tmp / "gwas.tsv"
        source.write_text("SNP\tP\n", encoding="utf-8")
        _, record = self.stage(self.row(source, "local"))
        self.assertEqual((record["phenotype"], record["settings_hash"], record["study_id"]), ('BMI "adj"', "s" * 16, "locke_2015_25673413"))

    def test_missing_local_file_fails_and_unsupported_sources_are_configuration_errors(self):
        proc, _ = self.stage(self.row(self.tmp / "absent.tsv", "local"), check=False)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("not found", proc.stderr)
        proc, _ = self.stage(self.row("syn123", "synapse"), check=False)
        self.assertEqual(proc.returncode, 78)


if __name__ == "__main__":
    unittest.main()
