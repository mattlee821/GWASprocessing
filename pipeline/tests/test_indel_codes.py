"""Indels reported as codes (Pietzner: D and I) get their sequences from the dbSNP indel table."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

import testlib  # noqa: F401  (puts pipeline/bin on the path)
sys.path.insert(0, str(testlib.ROOT / "src" / "reference_data"))
import build_dbsnp_indels
import gwas_reference as ref
import gwas_vcf
import standardise_gwas as sg

RECORDS = [  # rs, chromosome code, position, REF, ALT
    (1, 6, 1000, "CT", "C"),        # a deletion
    (2, 6, 2000, "A", "AGG"),       # an insertion
    (4, 6, 4000, "A", "AT"),        # 2 indels under 1 rsID: ambiguous
    (4, 6, 4000, "A", "ATT"),
    (5, 6, 5000, "AG", "T"),        # not a simple indel
    (6, 6, 6000, "GA", "G"),        # reached through a merged rsID (600 -> 6)
    (8, 7, 8000, "TC", "T"),        # on a different chromosome than the study gives
]


def build_index(root: Path) -> ref.DbsnpIndels:
    out = root / "dbsnp" / "GRCh38" / "indels"
    out.mkdir(parents=True)
    bucket = root / "bucket.tsv"
    bucket.write_text("".join(f"{rs}\t{chrom}\t{pos}\t{r}\t{a}\n" for rs, chrom, pos, r, a in RECORDS))
    counts = build_dbsnp_indels.write_index([bucket], out)
    (out / "meta.json").write_text('{"version": 1, "rows": %d}' % counts["rows"])
    return ref.DbsnpIndels.find(root)


def table() -> pd.DataFrame:
    rows = [  # SNP, CHR, EA, OA
        ("rs1", "6", "D", "I"), ("rs2", "6", "I", "D"), ("rs3", "6", "D", "I"), ("rs4", "6", "D", "I"), ("rs5", "6", "D", "I"),
        ("rs600", "6", "d", "i"), ("rs8", "6", "D", "I"), (".", "6", "D", "I"), ("rs9", "6", "A", "G"),
    ]
    return pd.DataFrame(rows, columns=["SNP", "CHR", "EA", "OA"]).assign(POS38=np.nan, IR=False)


class IndelCodesTest(unittest.TestCase):
    def test_recovery_and_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            indels = build_index(Path(tmp))
            dbsnp = SimpleNamespace(has_merges=True, current_rsids=lambda numbers: np.where(np.asarray(numbers) == 600, 6, 0))
            config = {"indel_codes": {"deletion": "D", "insertion": "I"}}
            qc: dict = {}
            out = sg.recover_indel_codes(table(), config, indels, dbsnp, qc)
        got = out.set_index("SNP")
        # the deletion allele is the shorter one; the effect allele is the one the code names
        self.assertEqual((got.loc["rs1", "EA"], got.loc["rs1", "OA"], got.loc["rs1", "POS38"]), ("C", "CT", 1000))
        self.assertEqual((got.loc["rs2", "EA"], got.loc["rs2", "OA"]), ("AGG", "A"))
        self.assertEqual((got.loc["rs600", "EA"], got.loc["rs600", "OA"], got.loc["rs600", "POS38"]), ("G", "GA", 6000))
        self.assertEqual(sorted(got.index[got["IR"]]), ["rs1", "rs2", "rs600"])
        self.assertEqual(got.loc["rs3", "EA"], "D")  # not recovered: left as it was
        self.assertEqual(got.loc["rs9", "EA"], "A")  # not a code: untouched
        self.assertEqual(qc["indel_codes"], {
            "coded_rows": 8, "recovered": 3, "recovered_via_merged_rsid": 1, "no_rsid": 1, "not_in_dbsnp": 1,
            "ambiguous_in_dbsnp": 1, "not_a_simple_indel": 1, "chromosome_conflict": 1})

    def test_off_without_the_yaml_setting(self):
        out = sg.recover_indel_codes(table(), {}, None, None, {})
        self.assertFalse(out["IR"].any())
        self.assertEqual(out["EA"].iloc[0], "D")

    def test_missing_table_is_a_configuration_error(self):
        with self.assertRaises(SystemExit) as raised:
            sg.recover_indel_codes(table(), {"indel_codes": {"deletion": "D", "insertion": "I"}}, None, None, {})
        self.assertEqual(raised.exception.code, 78)

    def test_header_and_flag(self):
        record = {"study_id": "s", "output_id": "o", "PMID": "1", "phenotype": "x", "ancestry": "EUR", "sex": "combined",
                  "author": "a", "year": "2021"}
        header = "\n".join(gwas_vcf.table_header(record, [], ["ES"], {}, "GRCh38", "h", indel_recovered=True))
        self.assertIn("##INFO=<ID=IR,", header)
        self.assertIn("bcftools view -e INFO/IR=1", header)
        plain = "\n".join(gwas_vcf.table_header(record, [], ["ES"], {}, "GRCh38", "h"))
        self.assertNotIn("ID=IR,", plain)


if __name__ == "__main__":
    unittest.main()
