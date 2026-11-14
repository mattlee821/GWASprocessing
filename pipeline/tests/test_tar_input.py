"""A protein delivered as a tar of per-chromosome files (UKB-PPP) is joined into 1 table by stage_source."""

from __future__ import annotations

import gzip
import io
import tarfile
import tempfile
import unittest
from pathlib import Path

import testlib  # noqa: F401  (puts pipeline/bin on the path)
import stage_source
from gwas_pipeline_common import EXIT_DATA

HEADER = "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA\n"


def make_tar(path: Path, members: dict[str, str]) -> None:
    with tarfile.open(path, "w") as tar:
        folder = tarfile.TarInfo("protein/")
        folder.type = tarfile.DIRTYPE
        tar.addfile(folder)
        for name, text in members.items():
            data = gzip.compress(text.encode())
            info = tarfile.TarInfo(f"protein/{name}")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


class TarTableTest(unittest.TestCase):
    def test_members_are_joined_with_one_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            archive = tmp / "A_P1_OID1_v1_Panel.tar"
            make_tar(archive, {
                "combined_chr2_A.gz": HEADER + "2 100 2:1:G:A:imp:v1 G A 0.1 0.9 100 ADD 0.1 0.1 1 0.5 NA\n",
                "combined_chr1_A.gz": HEADER + "1 50 1:1:C:T:imp:v1 C T 0.2 0.9 100 ADD 0.2 0.1 4 1.3 NA\n1 60 1:2:C:T:imp:v1 C T 0.3 0.9 100 ADD 0.3 0.1 9 3.0 NA\n",
            })
            raw = tmp / "raw"
            raw.mkdir()
            combined, info = stage_source.combine_tar_table(archive, raw, "P1_OID1")
            lines = gzip.open(combined, "rt").read().splitlines()
            self.assertEqual(lines[0] + "\n", HEADER)
            self.assertEqual(len(lines), 4)
            self.assertEqual(info["members"], 2)
            self.assertEqual(combined.name, "P1_OID1.combined.txt.gz")
            again, info2 = stage_source.combine_tar_table(archive, raw, "P1_OID1")
            self.assertTrue(info2["reused"])
            self.assertEqual(again, combined)

    def test_a_different_header_stops_the_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            archive = tmp / "B.tar"
            make_tar(archive, {"a.gz": HEADER + "1 1 x A C 0.1 1 1 ADD 0 1 0 0 NA\n", "b.gz": "CHROM GENPOS\n1 2\n"})
            (tmp / "raw").mkdir()
            with self.assertRaises(SystemExit) as raised:
                stage_source.combine_tar_table(archive, tmp / "raw", "B")
            self.assertEqual(raised.exception.code, EXIT_DATA)
            self.assertFalse(list((tmp / "raw").glob("*.part")))


if __name__ == "__main__":
    unittest.main()
