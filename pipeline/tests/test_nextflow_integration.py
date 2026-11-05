"""Integration test: the real Nextflow workflow driven by src/standardise.sh on a two-row local manifest.

Skipped unless nextflow, java, bgzip and tabix are on PATH. Takes about a minute (three Nextflow runs).
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from testlib import EXAMPLES, SRC, read_vcf_text, write_table

import gwas_pipeline_common as common

TOOLS = ("nextflow", "java", "bgzip", "tabix")


@unittest.skipUnless(all(shutil.which(t) for t in TOOLS), "nextflow, java and htslib are required")
class NextflowIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        manifests = self.tmp / "manifests"
        manifests.mkdir()
        good = shutil.copy(EXAMPLES / "tables/generic.tsv", self.tmp / "good.tsv")
        bad = write_table(self.tmp / "bad.tsv", ["CHR", "POS", "SNP", "EA", "OA", "Estimate", "SE", "P"],
                          [[1, 1000, "rs1", "A", "G", 0.1, 0.02, 1e-6]])
        self.manifest = manifests / "manifest.tsv"
        self.manifest.write_text(
            "GWAS_location\tphenotype\tancestry\tauthor\tyear\tPMID\n"
            f'{good}\tBMI "adjusted" $HOME\tEUR\tLocke\t2015\t1\n'
            f"{bad}\tWHR\tEUR\tBad\t2020\t2\n",
            encoding="utf-8",
        )
        self.env = {
            **os.environ,
            "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}",
            "GWAS_WORK_ROOT": str(self.tmp / "GWAS"), "GWAS_REFERENCE_ROOT": str(EXAMPLES / "reference"),
            "GWAS_CONFIG_ROOT": str(manifests), "GWAS_TEMP_ROOT": str(self.tmp / "temp"),
            "GWAS_LOG_ROOT": str(self.tmp / "temp" / "logs"), "GWAS_NEXTFLOW_WORK": str(self.tmp / "work"),
            "GWAS_NEXTFLOW_LOG": str(self.tmp / "temp" / ".nextflow.log"), "GWAS_QCPLOT": "true", "NXF_ANSI_LOG": "false",
        }

    def run_pipeline(self, *extra):
        cmd = ["bash", str(SRC / "standardise.sh"), "--manifest", str(self.manifest), "--profile", "local", "--row-limit", "0", *extra]
        return subprocess.run(cmd, env=self.env, capture_output=True, text=True, timeout=900)

    def state(self):
        return common.read_state(self.tmp / "GWAS" / ".standardise_state.tsv")

    def test_rows_run_fail_independently_skip_when_complete_and_rerun_when_settings_change(self):
        first = self.run_pipeline()
        self.assertEqual(first.returncode, 1, first.stdout + first.stderr)
        self.assertIn("1 row(s) failed", first.stdout)
        self.assertIn("bad_2020_2/", first.stdout)
        (row,) = self.state().values()
        self.assertEqual(row["study_id"], "locke_2015_1")
        vcf = Path(row["standard_vcf"])
        for key in ("standard_vcf", "standard_index", "metadata_json", "qc_json", "plot_tiff"):
            self.assertTrue(Path(row[key]).stat().st_size > 0, key)
        self.assertIn('##phenotype=BMI "adjusted" $HOME', read_vcf_text(vcf))
        first_hash = row["settings_hash"]

        second = self.run_pipeline()
        self.assertEqual(second.returncode, 1)
        summary = json.loads((self.tmp / "temp" / "logs" / "run_manifest_summary.json").read_text(encoding="utf-8"))
        self.assertEqual((summary["rows_skipped_complete"], summary["rows_to_run"]), (1, 1))

        third = self.run_pipeline("--analysis-qc", "true", "--only-study", "locke_2015_1")
        self.assertEqual(third.returncode, 0, third.stdout + third.stderr)
        (row,) = self.state().values()
        self.assertNotEqual(row["settings_hash"], first_hash)
        self.assertTrue(json.loads(Path(row["qc_json"]).read_text(encoding="utf-8"))["analysis_qc"])


if __name__ == "__main__":
    unittest.main()
