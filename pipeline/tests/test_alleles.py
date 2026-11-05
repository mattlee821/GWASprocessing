"""Tests for gwas_alleles: reference lookup, positional annotation, EAF fill, palindromic resolution,
de-duplication, allele-frequency concordance, FASTA access and reference-allele alignment."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import write_fasta, write_lookup

import gwas_alleles as alleles


def frame(**columns):
    n = len(next(iter(columns.values())))
    base = {"CHR": ["1"] * n, "POS19": [1000 + i for i in range(n)], "POS38": [1100 + i for i in range(n)],
            "SNP": [pd.NA] * n, "EAF": [np.nan] * n, "EAF_SOURCE": [""] * n, "BETA": [0.1] * n, "LIFT_STRAND": ["+"] * n}
    base.update(columns)
    return pd.DataFrame(base)


def lookup(rows):
    return pd.DataFrame(rows, columns=["CHR", "SNP", "POS19", "A1", "A2", "AF1"])


class ReferenceLookupTests(unittest.TestCase):
    def test_population_table_is_preferred_and_fallback_is_reported(self):
        root = Path(tempfile.mkdtemp())
        write_lookup(root, [("1", "rs1", 1000, "A", "G", 0.1)], population="ALL")
        self.assertEqual(alleles.read_reference_lookup(root, "EUR")[1], "ALL")
        write_lookup(root, [("1", "rs1", 1000, "A", "G", 0.3)], population="EUR")
        table, used = alleles.read_reference_lookup(root, "EUR")
        self.assertEqual((used, table.loc[0, "AF1"]), ("EUR", 0.3))

    def test_legacy_columns_are_ignored_and_values_normalised(self):
        root = Path(tempfile.mkdtemp())
        write_lookup(root, [("chr1", "RS1", 1000, "a", "g", 0.1, "1:1000", 0.1), ("X", "rs2", 500, "C", "T", 0.3, "X:500", 0.3)],
                     columns=["CHR", "SNP", "POS19", "A1", "A2", "AF1", "SNPID", "EAF_ALL"])
        table, _ = alleles.read_reference_lookup(root, "ALL")
        self.assertNotIn("EAF_ALL", table.columns)
        row = table.iloc[0]
        self.assertEqual((row["CHR"], row["RS"], row["POS19"], row["A1"], row["A2"], row["AF1"]), ("1", 1, 1000, "A", "G", 0.1))
        self.assertEqual(table["KEY"].tolist(), [1 * 10**10 + 1000, 23 * 10**10 + 500])

    def test_missing_tables_give_an_empty_frame(self):
        table, used = alleles.read_reference_lookup(Path(tempfile.mkdtemp()), "EUR")
        self.assertTrue(table.empty)
        self.assertEqual(used, "")


class AnnotationTests(unittest.TestCase):
    def test_match_needs_position_and_alleles_and_supplies_missing_rsids(self):
        data = frame(POS19=[1000, 1000, 2000], EA=["A", "A", "C"], OA=["G", "C", "T"])
        # Two reference variants share position 1000 (A/G and A/C); each input row must get its own allele match.
        ref = lookup([("1", "rs101", 1000, "A", "G", 0.2), ("1", "rs102", 1000, "A", "C", 0.4), ("1", "rs103", 3000, "C", "T", 0.5)])
        annotated = alleles.annotate_with_reference(data, ref)
        self.assertEqual(annotated["REF_SNP"].tolist()[:2], ["rs101", "rs102"])
        self.assertTrue(pd.isna(annotated["REF_SNP"].iloc[2]))
        filled, count = alleles.fill_missing_rsids(annotated)
        self.assertEqual((filled["SNP"].tolist()[:2], count), (["rs101", "rs102"], 2))

    def test_minus_strand_loci_are_matched_in_grch37_orientation(self):
        # After liftover onto the minus strand the GRCh38 alleles are the complement of the GRCh37 ones.
        data = frame(EA=["T"], OA=["C"], LIFT_STRAND=["-"])
        annotated = alleles.annotate_with_reference(data, lookup([("1", "rs1", 1000, "A", "G", 0.2)]))
        self.assertEqual(annotated["REF_SNP"].iloc[0], "rs1")
        ea37, oa37 = alleles.grch37_alleles(annotated)
        self.assertEqual((ea37.iloc[0], oa37.iloc[0]), ("A", "G"))


class EafFillTests(unittest.TestCase):
    def setUp(self):
        self.ref = lookup([("1", "rs1", 1000, "A", "G", 0.8), ("1", "rs2", 1001, "A", "G", 0.8), ("1", "rs3", 1002, "A", "G", 0.8),
                           ("1", "rs4", 1003, "A", "G", 0.8), ("1", "rs5", 1004, "A", "C", 0.8)])

    def test_fill_uses_the_frequency_of_the_effect_allele(self):
        # rs1: EA=A1. rs2: EA=A2. rs3: GWAS on the reverse strand (T/C), EA=T is A1's complement.
        data = frame(EA=["A", "G", "T"], OA=["G", "A", "C"])
        filled, count = alleles.fill_eaf_from_reference(alleles.annotate_with_reference(data, self.ref))
        np.testing.assert_allclose(filled["EAF"], [0.8, 0.2, 0.8])
        self.assertEqual((count, filled["EAF_SOURCE"].tolist()), (3, ["reference"] * 3))

    def test_study_frequencies_are_kept_and_unmatched_rows_stay_missing(self):
        data = frame(EA=["A", "A", "A", "A", "A"], OA=["G", "G", "G", "G", "T"], EAF=[0.3, np.nan, np.nan, np.nan, np.nan],
                     EAF_SOURCE=["input", "", "", "", ""], POS19=[1000, 1001, 1002, 9999, 1004])
        filled, count = alleles.fill_eaf_from_reference(alleles.annotate_with_reference(data, self.ref))
        self.assertEqual(filled["EAF"].iloc[0], 0.3)
        self.assertEqual(filled["EAF_SOURCE"].tolist(), ["input", "reference", "reference", "", ""])
        self.assertTrue(np.isnan(filled["EAF"].iloc[3:]).all())
        self.assertEqual(count, 2)


class PalindromicResolutionTests(unittest.TestCase):
    def test_resolve_palindromic_strand_keep_flip_drop(self):
        ref = lookup([("1", "rsP1", 1000, "A", "T", 0.20), ("1", "rsP2", 1001, "C", "G", 0.15), ("1", "rsP3", 1002, "A", "T", 0.47)])
        data = frame(SNP=["rsP1", "rsP2", "rsP3", "rsP4", "rsN1"], EA=["A", "C", "A", "A", "A"], OA=["T", "G", "T", "T", "G"],
                     EAF=[0.20, 0.85, 0.47, 0.20, 0.30], EAF_SOURCE=["input"] * 5)
        qc = {}
        resolved, dropped = alleles.resolve_palindromic_strand(alleles.annotate_with_reference(data, ref), dict(alleles.ANALYSIS_QC_DEFAULTS), qc)
        self.assertEqual(dropped, 2)
        detail = qc["palindromic_resolution"]
        self.assertEqual((detail["kept_same_strand"], detail["flipped"], detail["dropped_ambiguous"], detail["dropped_no_reference"]), (1, 1, 1, 1))
        by_snp = resolved.set_index("SNP")
        self.assertEqual((by_snp.loc["rsP1", "EA"], by_snp.loc["rsP1", "OA"]), ("A", "T"))
        self.assertEqual((by_snp.loc["rsP2", "EA"], by_snp.loc["rsP2", "OA"]), ("G", "C"))
        self.assertIn("rsN1", by_snp.index)

    def test_reference_filled_frequency_never_drives_a_flip(self):
        # EAF came from the reference itself, so it carries no strand information: drop, never flip.
        data = frame(EA=["A"], OA=["T"], EAF=[0.2], EAF_SOURCE=["reference"], BETA=[0.1])
        qc = {}
        resolved, dropped = alleles.resolve_palindromic_strand(
            alleles.annotate_with_reference(data, lookup([("1", "rs1", 1000, "A", "T", 0.8)])), dict(alleles.ANALYSIS_QC_DEFAULTS), qc)
        self.assertEqual((len(resolved), dropped), (0, 1))
        self.assertEqual(qc["palindromic_resolution"]["dropped_no_study_eaf"], 1)
        self.assertEqual(qc["palindromic_resolution"]["flipped"], 0)

    def test_minus_strand_palindrome_is_judged_in_grch37_orientation(self):
        # GRCh38 alleles T/A at a minus-strand locus are A/T on GRCh37. Study EAF 0.2 for GRCh38 T = GRCh37 A,
        # and the reference frequency of A is 0.2: same strand, so the variant is kept unchanged.
        data = frame(EA=["T"], OA=["A"], EAF=[0.2], EAF_SOURCE=["input"], LIFT_STRAND=["-"])
        qc = {}
        resolved, dropped = alleles.resolve_palindromic_strand(
            alleles.annotate_with_reference(data, lookup([("1", "rs1", 1000, "A", "T", 0.2)])), dict(alleles.ANALYSIS_QC_DEFAULTS), qc)
        self.assertEqual(dropped, 0)
        self.assertEqual((resolved["EA"].iloc[0], resolved["OA"].iloc[0]), ("T", "A"))
        self.assertEqual(qc["palindromic_resolution"]["kept_same_strand"], 1)

    def test_without_reference_alleles_every_palindrome_is_dropped(self):
        data = frame(EA=["A", "C", "A"], OA=["T", "G", "G"], EAF=[0.2, 0.2, 0.2], EAF_SOURCE=["input"] * 3)
        bare = pd.DataFrame(columns=["CHR", "SNP", "POS19", "A1", "A2", "AF1"])
        qc = {}
        resolved, dropped = alleles.resolve_palindromic_strand(alleles.annotate_with_reference(data, bare), dict(alleles.ANALYSIS_QC_DEFAULTS), qc)
        self.assertEqual((dropped, len(resolved)), (2, 1))


class DeduplicationAndQcTests(unittest.TestCase):
    def test_the_same_variant_coded_both_ways_is_kept_once(self):
        data = frame(POS38=[5000, 5000, 5000], EA=["A", "G", "A"], OA=["G", "A", "C"], BETA=[0.1, -0.1, 0.2])
        kept, removed = alleles.deduplicate(data)
        self.assertEqual((removed, kept["BETA"].tolist()), (1, [0.1, 0.2]))

    def test_analysis_qc_filters_run_in_order_with_counts(self):
        ref = lookup([("1", "rs1", 1000, "A", "G", 0.3)])
        data = frame(POS38=[1, 2, 3, 4, 1], POS19=[1000, 1001, 1002, 1003, 1000], EA=["A", "AT", "A", "A", "A"],
                     OA=["G", "G", "T", "G", "G"], EAF=[0.30, 0.30, 0.30, 0.005, 0.30], EAF_SOURCE=["input"] * 5)
        qc = {}
        kept = alleles.apply_analysis_qc(alleles.annotate_with_reference(data, ref), {}, qc)
        self.assertEqual(len(kept), 1)
        self.assertEqual(qc["qc_removed"], {"non_biallelic": 1, "palindromic": 1, "maf_below_0.01": 1, "duplicate": 1})


class ConcordanceTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(821)
        self.af = rng.uniform(0.05, 0.95, 2000)
        self.ref = lookup([("1", f"rs{i}", 1000 + i, "A", "G", af) for i, af in enumerate(self.af)])

    def test_concordant_frequencies_correlate_positively(self):
        data = frame(POS19=[1000 + i for i in range(2000)], EA=["A"] * 2000, OA=["G"] * 2000, EAF=self.af, EAF_SOURCE=["input"] * 2000)
        result = alleles.af_concordance(alleles.annotate_with_reference(data, self.ref))
        self.assertEqual(result["n"], 2000)
        self.assertGreater(result["pearson_r"], 0.99)

    def test_swapped_allele_columns_give_a_negative_correlation(self):
        data = frame(POS19=[1000 + i for i in range(2000)], EA=["G"] * 2000, OA=["A"] * 2000, EAF=self.af, EAF_SOURCE=["input"] * 2000)
        self.assertLess(alleles.af_concordance(alleles.annotate_with_reference(data, self.ref))["pearson_r"], -0.99)

    def test_reference_filled_and_palindromic_rows_are_excluded(self):
        data = frame(POS19=[1000, 1001], EA=["A", "A"], OA=["G", "T"], EAF=[0.3, 0.3], EAF_SOURCE=["reference", "input"])
        self.assertEqual(alleles.af_concordance(alleles.annotate_with_reference(data, self.ref))["n"], 0)


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
            (12, "AT", "T", "keep"), (13, "C", "CA", "swap"), (6, "GTT", "A", "indel_unmatched"), (7, "I", "D", "unalignable"),
        ]
        result = alleles.reference_alignment(pd.Series(["1"] * len(cases)), pd.Series([c[0] for c in cases]),
                                             pd.Series([c[1] for c in cases]), pd.Series([c[2] for c in cases]), fasta)
        self.assertEqual(list(result["action"]), [c[3] for c in cases])

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
