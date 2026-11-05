"""Tests for the GWAS-VCF audit scripts against real pipeline outputs (Python, and R via gwasvcf)."""

import csv
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from testlib import EXAMPLES, PIPELINE, run_standardise

PY_AUDIT = PIPELINE / "tests" / "gwasvcf_guide_audit.py"
R_AUDIT = PIPELINE / "tests" / "gwasvcf_guide_audit.R"


def r_has_gwasvcf() -> bool:
    if not shutil.which("Rscript"):
        return False
    probe = subprocess.run(["Rscript", "-e", "quit(status = !requireNamespace('gwasvcf', quietly = TRUE))"], capture_output=True)
    return probe.returncode == 0


def read_flat(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp())
        result = run_standardise(cls.tmp, EXAMPLES / "tables/generic.tsv", EXAMPLES / "reference")
        assert result.returncode == 0, result.stderr
        cls.vcf = Path(result.record["standard_vcf"])

    def test_python_audit_flattens_a_pipeline_output(self):
        out = self.tmp / "flat.python.tsv"
        proc = subprocess.run([sys.executable, str(PY_AUDIT), str(self.vcf), str(out)], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        rows = read_flat(out)
        self.assertEqual([(r["CHR"], r["POS"], r["SNP"], r["REF"], r["ALT"]) for r in rows],
                         [("1", "1100", "rs1", "G", "A"), ("1", "2100", "rs2", "T", "C")])
        self.assertEqual((rows[0]["BETA"], rows[0]["P"], rows[0]["N"]), ("0.2", "0.001", "1000"))

    def test_python_audit_rejects_records_without_the_required_format_keys(self):
        vcf = self.tmp / "noss.vcf"
        vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                       "1\t100\trs1\tA\tG\t.\tPASS\t.\tES:SE:LP\t0.1:0.01:2\n", encoding="utf-8")
        proc = subprocess.run([sys.executable, str(PY_AUDIT), str(vcf), str(self.tmp / "x.tsv")], capture_output=True, text=True)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("SS", proc.stderr)

    @unittest.skipUnless(r_has_gwasvcf(), "R with the gwasvcf package is required")
    def test_r_gwasvcf_reads_the_output_identically(self):
        py_out, r_out = self.tmp / "flat.python2.tsv", self.tmp / "flat.r.tsv"
        subprocess.run([sys.executable, str(PY_AUDIT), str(self.vcf), str(py_out)], check=True, capture_output=True)
        proc = subprocess.run(["Rscript", str(R_AUDIT), str(self.vcf), str(r_out)], capture_output=True, text=True, timeout=600)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        key = lambda r: (r["CHR"], r["POS"], r["REF"], r["ALT"])  # noqa: E731
        r_rows, py_rows = read_flat(r_out), read_flat(py_out)
        self.assertEqual([key(r) for r in r_rows], [key(r) for r in py_rows])
        self.assertEqual([float(r["BETA"]) for r in r_rows], [float(r["BETA"]) for r in py_rows])


if __name__ == "__main__":
    unittest.main()
