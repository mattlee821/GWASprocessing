"""Tests for gwas_coordinates (variant-ID parsing, rsID numbers) and build detection against the paired
GRCh37/GRCh38 sample (gwas_reference.build_hits, standardise_gwas.detect_build_paired)."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import write_grch38_references

import gwas_coordinates as coords
import gwas_reference
from standardise_gwas import detect_build_paired

OFFSET = 100  # GRCh38 = GRCh37 + 100 in the fixture


class BuildDetectionTests(unittest.TestCase):
    def setUp(self):
        self.positions = [1000 * (i + 1) for i in range(40)]
        root = Path(tempfile.mkdtemp())
        write_grch38_references(root, [], paired=[("1", p, p + OFFSET, "G", "A") for p in self.positions])
        self.paired = gwas_reference.paired_build_sample(root)

    def detect(self, positions, a="A", b="G"):
        n = len(positions)
        rows = pd.DataFrame({"CHR": ["1"] * n, "POS": positions, "EA": [a] * n, "OA": [b] * n})
        return detect_build_paired(rows, self.paired)

    def test_grch37_positions_are_detected(self):
        result = self.detect(self.positions)
        self.assertEqual((result["build"], result["hits_hg19"], result["hits_hg38"]), ("hg19", 40, 0))

    def test_grch38_positions_are_detected(self):
        result = self.detect([p + OFFSET for p in self.positions])
        self.assertEqual((result["build"], result["hits_hg19"], result["hits_hg38"]), ("hg38", 0, 40))

    def test_reverse_strand_alleles_still_match(self):
        self.assertEqual(self.detect(self.positions, a="T", b="C")["build"], "hg19")

    def test_mismatched_alleles_do_not_count(self):
        result = self.detect(self.positions, a="A", b="C")
        self.assertEqual((result["build"], result["hits_hg19"]), ("", 0))

    def test_no_overlap_is_inconclusive(self):
        self.assertEqual(self.detect([123_457 + 7 * i for i in range(40)])["build"], "")

    def test_a_build_needs_ten_times_the_other_builds_hits(self):
        mixed = self.positions[:20] + [p + OFFSET for p in self.positions[20:25]]
        self.assertEqual(self.detect(mixed)["build"], "")


class VariantIdParsingTests(unittest.TestCase):
    def test_common_chr_pos_formats_are_parsed(self):
        parsed = coords.parse_variant_ids(pd.Series(["1:100", "chr2:200:A:G", "X_300_C_T", "rs123", "1:abc", None]))
        self.assertEqual(parsed["CHR"].tolist()[:3], ["1", "2", "X"])
        self.assertTrue(parsed["CHR"].iloc[3:].isna().all())
        self.assertEqual(parsed["POS"].tolist()[:3], [100, 200, 300])
        self.assertTrue(parsed["POS"].iloc[3:].isna().all())

    def test_rs_numbers_require_a_bare_rsid(self):
        np.testing.assert_array_equal(coords.rs_number(pd.Series(["rs12", "RS7", "rs12:100:A:G", "x", None])), [12, 7, -1, -1, -1])


if __name__ == "__main__":
    unittest.main()
