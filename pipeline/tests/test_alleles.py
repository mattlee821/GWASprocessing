"""Tests for gwas_alleles: 1000 Genomes GRCh38 annotation, EAF fill, palindromic resolution, de-duplication,
allele-frequency concordance, FASTA access and reference-allele alignment."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import write_fasta, write_grch38_references

import gwas_alleles as alleles
import gwas_reference


def frame(**columns):
    """Rows on chromosome 1 at GRCh38 positions 1000, 1001, ... unless POS38 is given."""
    n = len(next(iter(columns.values())))
    base = {"CHR": ["1"] * n, "POS38": [1000 + i for i in range(n)], "SNP": [pd.NA] * n, "EAF": [np.nan] * n,
            "EAF_SOURCE": [""] * n, "BETA": [0.1] * n}
    base.update(columns)
    return pd.DataFrame(base)


def frequencies(rows):
    """1000 Genomes GRCh38 frequencies from rows (chrom, pos38, ALT, REF, ALT frequency)."""
    root = Path(tempfile.mkdtemp())
    write_grch38_references(root, [(c, p, 0, ref, alt, af) for c, p, alt, ref, af in rows])
    return gwas_reference.Frequencies.find(root, "ALL")


def annotate(data, rows):
    return alleles.annotate_with_frequencies(data, frequencies(rows))


class AnnotationTests(unittest.TestCase):
    def test_match_needs_position_and_alleles(self):
        data = frame(POS38=[1000, 1000, 2000], EA=["A", "A", "C"], OA=["G", "C", "T"])
        # Two 1000 Genomes records share position 1000 (A/G and A/C); each input row must get its own allele match.
        annotated = annotate(data, [("1", 1000, "A", "G", 0.2), ("1", 1000, "A", "C", 0.4), ("1", 3000, "C", "T", 0.5)])
        np.testing.assert_allclose(annotated["REF_AF1"].iloc[:2], [0.2, 0.4])
        self.assertEqual((annotated["REF_A1"].iloc[1], annotated["REF_A2"].iloc[1]), ("A", "C"))
        self.assertTrue(np.isnan(annotated["REF_AF1"].iloc[2]))

    def test_reverse_strand_alleles_match(self):
        annotated = annotate(frame(EA=["T"], OA=["C"]), [("1", 1000, "A", "G", 0.2)])
        self.assertAlmostEqual(annotated["REF_AF1"].iloc[0], 0.2, places=6)


class EafFillTests(unittest.TestCase):
    rows = [("1", 1000, "A", "G", 0.8), ("1", 1001, "A", "G", 0.8), ("1", 1002, "A", "G", 0.8), ("1", 1003, "A", "G", 0.8),
            ("1", 1004, "A", "C", 0.8)]

    def test_fill_uses_the_frequency_of_the_effect_allele(self):
        # 1000: EA=ALT. 1001: EA=REF. 1002: GWAS on the reverse strand (T/C), EA=T is ALT's complement.
        data = frame(EA=["A", "G", "T"], OA=["G", "A", "C"])
        filled, count = alleles.fill_eaf_from_reference(annotate(data, self.rows))
        np.testing.assert_allclose(filled["EAF"], [0.8, 0.2, 0.8])
        self.assertEqual((count, filled["EAF_SOURCE"].tolist()), (3, ["reference"] * 3))

    def test_study_frequencies_are_kept_and_unmatched_rows_stay_missing(self):
        data = frame(EA=["A", "A", "A", "A", "A"], OA=["G", "G", "G", "G", "T"], EAF=[0.3, np.nan, np.nan, np.nan, np.nan],
                     EAF_SOURCE=["input", "", "", "", ""], POS38=[1000, 1001, 1002, 9999, 1004])
        filled, count = alleles.fill_eaf_from_reference(annotate(data, self.rows))
        self.assertEqual(filled["EAF"].iloc[0], 0.3)
        self.assertEqual(filled["EAF_SOURCE"].tolist(), ["input", "reference", "reference", "", ""])
        self.assertTrue(np.isnan(filled["EAF"].iloc[3:]).all())
        self.assertEqual(count, 2)


class PalindromicResolutionTests(unittest.TestCase):
    def test_palindromes_are_kept_complemented_or_marked_ambiguous_never_dropped(self):
        rows = [("1", 1000, "A", "T", 0.20), ("1", 1001, "C", "G", 0.15), ("1", 1002, "A", "T", 0.47)]
        data = frame(SNP=["rsP1", "rsP2", "rsP3", "rsP4", "rsN1"], EA=["A", "C", "A", "A", "A"], OA=["T", "G", "T", "T", "G"],
                     EAF=[0.20, 0.85, 0.47, 0.20, 0.30], EAF_SOURCE=["input"] * 5)
        qc = {}
        resolved = alleles.resolve_palindromic_strand(annotate(data, rows), qc)
        self.assertEqual(len(resolved), 5)
        detail = qc["palindromic_resolution"]
        self.assertEqual((detail["kept_same_strand"], detail["flipped"], detail["ambiguous_frequency_near_half"],
                          detail["ambiguous_no_reference"]), (1, 1, 1, 1))
        by_snp = resolved.set_index("SNP")
        self.assertEqual((by_snp.loc["rsP1", "EA"], by_snp.loc["rsP1", "OA"], by_snp.loc["rsP1", "PAL"]), ("A", "T", "forward"))
        self.assertEqual((by_snp.loc["rsP2", "EA"], by_snp.loc["rsP2", "OA"], by_snp.loc["rsP2", "PAL"]), ("G", "C", "complemented"))
        self.assertEqual((by_snp.loc["rsP3", "PAL"], by_snp.loc["rsP4", "PAL"], by_snp.loc["rsN1", "PAL"]), ("ambiguous", "ambiguous", ""))
        self.assertEqual((by_snp.loc["rsP4", "EA"], by_snp.loc["rsP4", "OA"]), ("A", "T"))

    def test_reference_filled_frequency_never_drives_a_flip(self):
        # EAF came from the reference itself, so it carries no strand information: ambiguous, never flipped.
        data = frame(EA=["A"], OA=["T"], EAF=[0.2], EAF_SOURCE=["reference"], BETA=[0.1])
        qc = {}
        resolved = alleles.resolve_palindromic_strand(annotate(data, [("1", 1000, "A", "T", 0.8)]), qc)
        self.assertEqual((len(resolved), resolved["PAL"].iloc[0], resolved["EA"].iloc[0]), (1, "ambiguous", "A"))
        self.assertEqual(qc["palindromic_resolution"]["ambiguous_no_study_eaf"], 1)
        self.assertEqual(qc["palindromic_resolution"]["flipped"], 0)

    def test_without_reference_alleles_every_palindrome_is_ambiguous(self):
        data = frame(EA=["A", "C", "A"], OA=["T", "G", "G"], EAF=[0.2, 0.2, 0.2], EAF_SOURCE=["input"] * 3)
        qc = {}
        resolved = alleles.resolve_palindromic_strand(annotate(data, []), qc)
        self.assertEqual(resolved["PAL"].tolist(), ["ambiguous", "ambiguous", ""])

    def test_tolerance_can_be_widened(self):
        data = frame(EA=["A"], OA=["T"], EAF=[0.38], EAF_SOURCE=["input"])
        rows = [("1", 1000, "A", "T", 0.30)]
        self.assertEqual(alleles.resolve_palindromic_strand(annotate(data, rows), {})["PAL"].iloc[0], "forward")
        self.assertEqual(alleles.resolve_palindromic_strand(annotate(data, rows), {}, 0.15)["PAL"].iloc[0], "ambiguous")


class DeduplicationTests(unittest.TestCase):
    def test_the_same_variant_coded_both_ways_is_kept_once(self):
        data = frame(POS38=[5000, 5000, 5000], EA=["A", "G", "A"], OA=["G", "A", "C"], BETA=[0.1, -0.1, 0.2])
        kept, removed = alleles.deduplicate(data)
        self.assertEqual((removed, kept["BETA"].tolist()), (1, [0.1, 0.2]))


class ConcordanceTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(821)
        self.af = rng.uniform(0.05, 0.95, 2000)
        self.freqs = frequencies([("1", 1000 + i, "A", "G", af) for i, af in enumerate(self.af)])

    def test_concordant_frequencies_correlate_positively(self):
        data = frame(POS38=[1000 + i for i in range(2000)], EA=["A"] * 2000, OA=["G"] * 2000, EAF=self.af, EAF_SOURCE=["input"] * 2000)
        result = alleles.af_concordance(alleles.annotate_with_frequencies(data, self.freqs))
        self.assertEqual(result["n"], 2000)
        self.assertGreater(result["pearson_r"], 0.99)

    def test_swapped_allele_columns_give_a_negative_correlation(self):
        data = frame(POS38=[1000 + i for i in range(2000)], EA=["G"] * 2000, OA=["A"] * 2000, EAF=self.af, EAF_SOURCE=["input"] * 2000)
        self.assertLess(alleles.af_concordance(alleles.annotate_with_frequencies(data, self.freqs))["pearson_r"], -0.99)

    def test_reference_filled_and_palindromic_rows_are_excluded(self):
        data = frame(POS38=[1000, 1001], EA=["A", "A"], OA=["G", "T"], EAF=[0.3, 0.3], EAF_SOURCE=["reference", "input"])
        self.assertEqual(alleles.af_concordance(alleles.annotate_with_frequencies(data, self.freqs))["n"], 0)


class FastaTests(unittest.TestCase):
    def test_appledouble_files_are_not_taken_for_the_fasta(self):
        root = Path(tempfile.mkdtemp())
        (root / "genome").mkdir()
        (root / "genome" / "._hg38.fa").write_bytes(b"\x00\x05\x16\x07")
        self.assertIsNone(alleles.Fasta.find(root))

    def test_bases_are_read_by_position_with_flexible_contig_names(self):
        genome = Path(tempfile.mkdtemp())
        write_fasta(genome, {"chr1": "ACGTACGTAC" * 13, "chrM": "GATTACA"}, width=7)
        fasta = alleles.Fasta(genome / "GRCh38.fa")
        self.assertTrue((genome / "GRCh38.fa.fai").exists())
        self.assertEqual(fasta.bases("1", np.array([1, 2, 8, 9, 130])).tolist(), ["A", "C", "T", "A", "C"])
        self.assertEqual(fasta.bases("1", np.array([3]), length=4).tolist(), ["GTAC"])
        self.assertEqual(fasta.bases("MT", np.array([4])).tolist(), ["T"])
        self.assertEqual(fasta.bases("1", np.array([131, 0])).tolist(), ["", ""])
        self.assertEqual(fasta.bases("2", np.array([1])).tolist(), [""])

    def test_find_locates_the_grch38_fasta_under_the_reference_root(self):
        root = Path(tempfile.mkdtemp())
        self.assertIsNone(alleles.Fasta.find(root))
        write_fasta(root / "genome", {"1": "ACGT"})
        self.assertEqual(alleles.Fasta.find(root).path.name, "GRCh38.fa")


class ReferenceAlignmentTests(unittest.TestCase):
    def test_actions_for_snvs_palindromes_and_indels(self):
        genome = Path(tempfile.mkdtemp())
        # Positions 1..12 on chromosome 1.
        write_fasta(genome, {"1": "GGGGGGGGAAGTCC"})
        fasta = alleles.Fasta(genome / "GRCh38.fa")
        cases = [  # (position, EA, OA, expected action)
            (1, "A", "G", "keep"), (2, "G", "A", "swap"), (3, "T", "C", "complement"), (4, "C", "T", "complement_swap"),
            (5, "A", "T", "drop"), (9, "T", "A", "keep"), (10, "A", "T", "swap"), (11, "C", "G", "keep"),
            (12, "AT", "T", "keep"), (13, "C", "CA", "swap"), (6, "GTT", "A", "indel_unmatched"),
        ]
        result = alleles.reference_alignment(pd.Series(["1"] * len(cases)), pd.Series([c[0] for c in cases]),
                                             pd.Series([c[1] for c in cases]), pd.Series([c[2] for c in cases]), fasta)
        self.assertEqual(list(result["action"]), [c[3] for c in cases])

    def test_indels_are_left_aligned_and_trimmed(self):
        genome = Path(tempfile.mkdtemp())
        write_fasta(genome, {"1": "GCACACATT"})  # a CA repeat at positions 2-7
        fasta = alleles.Fasta(genome / "GRCh38.fa")
        data = frame(POS38=[5, 3, 8], EA=["A", "ACAC", "C"], OA=["ACA", "AC", "T"])
        action = np.array(["keep", "keep", "keep"], dtype=object)
        normalised, moved = alleles.normalise_indels(data, fasta, action)
        # deletion of CA written at the right of the repeat moves to the left anchor G at position 1;
        # an over-long insertion AC>ACAC is trimmed to its shortest left-aligned form; the SNV is untouched
        self.assertEqual(list(zip(normalised["POS38"].astype(int), normalised["OA"], normalised["EA"])),
                         [(1, "GCA", "G"), (1, "G", "GCA"), (8, "T", "C")])
        self.assertEqual(moved, 2)

    def test_alignment_reorients_effect_size_and_frequency(self):
        genome = Path(tempfile.mkdtemp())
        write_fasta(genome, {"1": "AAAAGC"})
        fasta = alleles.Fasta(genome / "GRCh38.fa")
        data = frame(POS38=[5, 6, 5], EA=["G", "A", "T"], OA=["A", "G", "A"], BETA=[0.2, 0.3, 0.4], EAF=[0.1, 0.6, 0.5])
        aligned, summary = alleles.align_table_to_reference(data, fasta)
        self.assertEqual(aligned["EA"].tolist(), ["A", "T"])
        self.assertEqual(aligned["OA"].tolist(), ["G", "C"])
        np.testing.assert_allclose(aligned["BETA"], [-0.2, 0.3])
        np.testing.assert_allclose(aligned["EAF"], [0.9, 0.6])
        self.assertEqual((summary["swapped"], summary["strand_flipped"], summary["dropped_mismatch"]), (1, 1, 1))


if __name__ == "__main__":
    unittest.main()
