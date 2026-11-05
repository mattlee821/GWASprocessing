"""Tests for gwas_coordinates: liftover outcomes, build detection, variant-ID parsing, rsID coordinate index."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import chain_block, write_chain, write_offset_chains

import gwas_coordinates as coords
from pyliftover import LiftOver


def lifter_from(blocks):
    tmp = Path(tempfile.mkdtemp())
    return LiftOver(str(write_chain(tmp / "test.over.chain.gz", blocks)))


class LiftTests(unittest.TestCase):
    def test_minus_strand_target_is_reported_with_forward_position(self):
        # Block chr1:[7000,8000) -> chr1 minus strand [10000,11000) of a 248,956,422 bp chromosome.
        # 1-based 7500 is offset 499: reverse coordinate 10499, forward 0-based 248956422-1-10499.
        lifter = lifter_from([chain_block("chr1", "chr1", 249250621, 248956422, 7000, 1000, 10000, 1, dst_strand="-")])
        result = coords.lift(pd.Series(["1"]), pd.Series([7500]), lifter)
        self.assertEqual(list(result.status), ["ok"])
        self.assertEqual(list(result.strand), ["-"])
        self.assertEqual(result.pos[0], 248945923)

    def test_landing_on_another_chromosome_is_a_failure(self):
        lifter = lifter_from([chain_block("chr1", "chr2", 249250621, 242193529, 5000, 1000, 90000, 1)])
        result = coords.lift(pd.Series(["1"]), pd.Series([5500]), lifter)
        self.assertEqual(list(result.status), ["other_chrom"])
        self.assertTrue(np.isnan(result.pos[0]))

    def test_unmapped_multiple_and_missing_inputs_are_distinguished(self):
        lifter = lifter_from([
            chain_block("chr1", "chr1", 249250621, 248956422, 0, 1000, 5000, 1),
            chain_block("chr1", "chr1", 249250621, 248956422, 0, 1000, 9000, 2),
            chain_block("chr1", "chr1", 249250621, 248956422, 2000, 1000, 2100, 3),
        ])
        result = coords.lift(pd.Series(["1", "1", "1", None]), pd.Series([500, 2500, 50_000, 10]), lifter)
        self.assertEqual(list(result.status), ["multiple", "ok", "unmapped", "missing"])
        self.assertEqual(result.pos[1], 2600)

    def test_mitochondrial_contig_uses_ucsc_chrM_name(self):
        lifter = lifter_from([chain_block("chrM", "chrM", 16571, 16569, 0, 1000, 0, 1)])
        result = coords.lift(pd.Series(["MT"]), pd.Series([100]), lifter)
        self.assertEqual((result.status[0], result.pos[0]), ("ok", 100))


class BuildDetectionTests(unittest.TestCase):
    def setUp(self):
        positions = [1000 * (i + 1) for i in range(40)]
        self.lookup = pd.DataFrame({"CHR": ["1"] * 40, "POS19": positions, "A1": ["A"] * 40, "A2": ["G"] * 40})
        tmp = Path(tempfile.mkdtemp())
        write_offset_chains(tmp)  # GRCh38 = GRCh37 + 100
        self.lift38to19 = LiftOver(str(tmp / "hg38ToHg19.over.chain.gz"))
        self.positions = positions

    def detect(self, positions, a="A", b="G", lifter="default"):
        n = len(positions)
        return coords.detect_build(pd.Series(["1"] * n), pd.Series(positions), pd.Series([a] * n), pd.Series([b] * n),
                                   self.lookup, self.lift38to19 if lifter == "default" else lifter)

    def test_grch37_positions_are_detected(self):
        result = self.detect(self.positions)
        self.assertEqual((result["build"], result["hits_hg19"], result["hits_hg38"]), ("hg19", 40, 0))

    def test_grch38_positions_are_detected_by_lifting_back(self):
        result = self.detect([p + 100 for p in self.positions])
        self.assertEqual((result["build"], result["hits_hg19"], result["hits_hg38"]), ("hg38", 0, 40))

    def test_reverse_strand_alleles_still_match(self):
        self.assertEqual(self.detect(self.positions, a="T", b="C")["build"], "hg19")

    def test_mismatched_alleles_do_not_count(self):
        result = self.detect(self.positions, a="A", b="C")
        self.assertEqual((result["build"], result["hits_hg19"]), ("", 0))

    def test_no_overlap_is_inconclusive(self):
        self.assertEqual(self.detect([123_457 + 7 * i for i in range(40)])["build"], "")

    def test_without_the_hg38_chain_only_grch37_can_be_confirmed(self):
        self.assertEqual(self.detect(self.positions, lifter=None)["build"], "hg19")
        self.assertEqual(self.detect([p + 100 for p in self.positions], lifter=None)["build"], "")


class VariantIdParsingTests(unittest.TestCase):
    def test_common_chr_pos_formats_are_parsed(self):
        parsed = coords.parse_variant_ids(pd.Series(["1:100", "chr2:200:A:G", "X_300_C_T", "rs123", "1:abc", None]))
        self.assertEqual(parsed["CHR"].tolist()[:3], ["1", "2", "X"])
        self.assertTrue(parsed["CHR"].iloc[3:].isna().all())
        self.assertEqual(parsed["POS"].tolist()[:3], [100, 200, 300])
        self.assertTrue(parsed["POS"].iloc[3:].isna().all())


class CoordinateIndexTests(unittest.TestCase):
    def test_rs_numbers_require_a_bare_rsid(self):
        np.testing.assert_array_equal(coords.rs_number(pd.Series(["rs12", "RS7", "rs12:100:A:G", "x", None])), [12, 7, -1, -1, -1])

    def test_round_trip_keeps_consistent_duplicates_and_drops_conflicts(self):
        tmp = Path(tempfile.mkdtemp())
        coords.write_coordinate_index(tmp, rs=np.array([5, 3, 3, 9, 9]), chrom=np.array(["1", "2", "2", "1", "3"]),
                                      pos=np.array([500, 300, 300, 900, 901]))
        index = coords.CoordinateIndex.load(tmp)
        found = index.lookup(pd.Series(["rs3", "rs5", "rs9", "rs4", "1:100", None]))
        self.assertEqual(found["CHR"].tolist()[:2], ["2", "1"])
        self.assertTrue(found["CHR"].iloc[2:].isna().all())
        self.assertEqual(found["POS19"].tolist()[:2], [300, 500])
        self.assertEqual(index.size, 2)

    def test_missing_index_loads_as_none(self):
        self.assertIsNone(coords.CoordinateIndex.load(Path(tempfile.mkdtemp()) / "absent"))


if __name__ == "__main__":
    unittest.main()
