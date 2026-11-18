"""Fixed-size batch manifest splitting, status tracking, raw cleanup, and safe VCF purging (batches.py / src/batches.sh)."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from testlib import run_script

HEADER = ("GWAS_location\tphenotype\tancestry\tsex\tauthor\tyear\tPMID\tsource_type\t"
          "input_build\tpopulation\tdelimiter\tformat\tstudy_yaml\toutput_name\n")


class BatchesTest(unittest.TestCase):
    def test_split_status_clean_and_purge(self):
        with tempfile.TemporaryDirectory(prefix="gwas_batches_") as tmp_str:
            tmp = Path(tmp_str)
            manifest = tmp / "sun.tsv"
            rows = [
                f"syn5100000{i}\tProt {i}\tcombined\tcombined\tSun\t2023\t37794186\tsynapse\t"
                f"hg38\tALL\tspace\ttable\tmanifests/sun_2023_37794186.yaml\tP0000{i}_OID2000{i}"
                for i in range(1, 6)
            ]
            manifest.write_text(HEADER + "\n".join(rows) + "\n", encoding="utf-8")
            batches_dir = tmp / "sun_batches"
            work_root = tmp / "GWAS"
            export_dir = tmp / "exports"
            export_dir.mkdir()

            with self.subTest("split writes fixed-size batch_XX.tsv files and batches_index.tsv"):
                proc = run_script("batches.py", "split", "--manifest", manifest, "--out-dir", batches_dir, "--size", "2")
                self.assertEqual(proc.returncode, 0, proc.stderr)
                files = sorted(batches_dir.glob("batch_*.tsv"))
                self.assertEqual([f.name for f in files], ["batch_01.tsv", "batch_02.tsv", "batch_03.tsv"])
                self.assertEqual(len(files[0].read_text(encoding="utf-8").splitlines()), 3)
                self.assertEqual(len(files[2].read_text(encoding="utf-8").splitlines()), 2)

            with self.subTest("status starts with all batches PENDING"):
                proc = run_script("batches.py", "status", "--batches", batches_dir, "--work-root", work_root)
                self.assertEqual(proc.returncode, 0, proc.stderr)
                self.assertIn("Pending: 3", proc.stdout)
                self.assertIn("batch_01.tsv", proc.stdout)

            # Simulate batch_01 completing standardisation with staged raw files
            for i in (1, 2):
                oid = f"P0000{i}_OID2000{i}"
                raw = work_root / "sun_2023_37794186" / "raw" / oid
                std = work_root / "sun_2023_37794186" / "standardise" / oid
                raw.mkdir(parents=True)
                std.mkdir(parents=True)
                (raw / f"{oid}.combined.txt.gz").write_bytes(b"staged-raw-data")
                (std / f"{oid}.gwas.vcf.gz").write_bytes(b"vcf-bytes")
                (std / f"{oid}.gwas.vcf.gz.tbi").write_bytes(b"tbi-bytes")
                (std / f"{oid}.qc.json").write_text('{"output_rows": 10}\n', encoding="utf-8")
                (std / f"{oid}.metadata.json").write_text('{"phenotype": "Prot"}\n', encoding="utf-8")

            with self.subTest("clean-raw removes raw/<id> while keeping VCFs intact"):
                proc = run_script("batches.py", "status", "--batches", batches_dir, "--work-root", work_root)
                self.assertIn("STANDARDISED (raw staged)", proc.stdout)
                clean = run_script("batches.py", "clean-raw", "--batches", batches_dir, "--batch", "1", "--work-root", work_root)
                self.assertEqual(clean.returncode, 0, clean.stderr)
                self.assertFalse((work_root / "sun_2023_37794186" / "raw" / "P00001_OID20001").exists())
                self.assertTrue((work_root / "sun_2023_37794186" / "standardise" / "P00001_OID20001" / "P00001_OID20001.gwas.vcf.gz").exists())

            with self.subTest("purge refuses when --export-dir is missing protein exports, and refuses unstandardised batches"):
                bad_exp = run_script("batches.py", "purge", "--batches", batches_dir, "--batch", "01",
                                     "--work-root", work_root, "--export-dir", export_dir, check=False)
                self.assertNotEqual(bad_exp.returncode, 0)
                self.assertIn("no exported output", bad_exp.stderr)
                bad_pending = run_script("batches.py", "purge", "--batches", batches_dir, "--batch", "02",
                                         "--work-root", work_root, check=False)
                self.assertNotEqual(bad_pending.returncode, 0)
                self.assertIn("missing .qc.json", bad_pending.stderr)

            with self.subTest("purge deletes VCFs after export verification, keeps .qc.json, and marks batch DONE"):
                for i in (1, 2):
                    (export_dir / f"P0000{i}_OID2000{i}.finemap.tsv").write_text("credset\n", encoding="utf-8")
                ok = run_script("batches.py", "purge", "--batches", batches_dir, "--batch", "01",
                                "--work-root", work_root, "--export-dir", export_dir)
                self.assertEqual(ok.returncode, 0, ok.stderr)
                std1 = work_root / "sun_2023_37794186" / "standardise" / "P00001_OID20001"
                self.assertFalse((std1 / "P00001_OID20001.gwas.vcf.gz").exists())
                self.assertTrue((std1 / "P00001_OID20001.qc.json").exists())
                self.assertTrue((std1 / "P00001_OID20001.purged.json").exists())

                status = run_script("batches.py", "status", "--batches", batches_dir, "--work-root", work_root)
                self.assertIn("DONE (purged)", status.stdout)
                self.assertIn("Done (purged): 1", status.stdout)
                self.assertIn("Next pending batch: " + str(batches_dir / "batch_02.tsv"), status.stdout)


if __name__ == "__main__":
    unittest.main()
