"""README step 8, palindromic SNPs: the strand of each A/T and C/G SNP is inferred from the study frequency against the
1000 Genomes frequency of the same allele (TwoSampleMR harmonise_data, action 2). Same side of 0.5: kept; opposite
sides: alleles complemented with ES and EAF unchanged; either frequency within the tolerance (0.08, palindromic_tolerance)
of 0.5, or either missing: kept as reported and marked ambiguous. INFO/HM records the decision.

The first tables are declared GRCh38 and the reference has no FASTA, so HM shows the palindrome decision alone. The
liftover case uses a minus-strand chain block (seeded genome, 821): GRCh37 chr1 [1000, 2000) maps to the reverse strand of
GRCh38 chr1, so GRCh37 position P is GRCh38 position 3501 - P.
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import (chain_block, make_reference, needs_liftover, run_standardise, write_chain, write_fasta,  # noqa: I001
                     write_grch38_references, write_table)

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
HEADER = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
# 1000 Genomes records (pos38, REF, ALT, ALT frequency)
THOUSAND_GENOMES = [(1000, "A", "T", 0.8), (1100, "C", "G", 0.15), (1200, "A", "T", 0.2), (1300, "A", "T", 0.53),
                    (1500, "A", "T", 0.8), (1600, "A", "T", 0.7), (1700, "G", "A", 0.3)]
TABLE = [
    [1, 1000, "same_side", "A", "T", 0.2, 0.1, 0.02, 1e-6, 5000],  # 1000 Genomes A = 0.2
    [1, 1100, "opposite_sides", "C", "G", 0.15, 0.1, 0.02, 1e-6, 5000],  # 1000 Genomes C = 0.85
    [1, 1200, "study_near_half", "A", "T", 0.45, 0.1, 0.02, 1e-6, 5000],
    [1, 1300, "reference_near_half", "A", "T", 0.2, 0.1, 0.02, 1e-6, 5000],  # 1000 Genomes A = 0.47
    [1, 1400, "no_reference", "A", "T", 0.2, 0.1, 0.02, 1e-6, 5000],
    [1, 1500, "filled", "A", "T", None, 0.1, 0.02, 1e-6, 5000],  # filled from 1000 Genomes: no strand information
    [1, 1600, "outside_tolerance", "A", "T", 0.59, 0.1, 0.02, 1e-6, 5000],  # 0.09 from 0.5; 1000 Genomes A = 0.3
    [1, 1700, "not_palindromic", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
]


def frequencies_reference(root: Path, records: list[tuple]) -> Path:
    return write_grch38_references(root, [("1", p, 0, ref, alt, af) for p, ref, alt, af in records])


def hm(result) -> dict[str, str]:
    return {snp: row["INFO"]["HM"] for snp, row in result.by_snpid().items()}


class Step08PalindromesTest(unittest.TestCase):
    @needs_liftover
    def test_step08_palindromic_snps(self):
        tmp = Path(tempfile.mkdtemp())
        reference = frequencies_reference(tmp / "reference", THOUSAND_GENOMES)
        raw = write_table(tmp / "gwas.tsv", HEADER, TABLE)
        result = run_standardise(tmp, raw, reference, input_build="hg38", ancestry="EUR", name="default")
        self.assertEqual(result.returncode, 0, result.stderr)
        rows = result.by_snpid()
        marks = hm(result)
        detail = result.qc["palindromic_resolution"]

        with self.subTest("palindromes are kept, complemented or marked ambiguous, never dropped"):
            self.assertEqual(sorted(rows), sorted(r[2] for r in TABLE))
            self.assertEqual(detail["palindromic_total"], 7)

        with self.subTest("same side of 0.5: kept as reported (palindrome_forward)"):
            self.assertEqual((rows["same_side"]["REF"], rows["same_side"]["ALT"], marks["same_side"]), ("T", "A", "palindrome_forward"))

        with self.subTest("opposite sides of 0.5: alleles complemented with ES and EAF unchanged (palindrome_complemented)"):
            flipped = rows["opposite_sides"]
            self.assertEqual((flipped["REF"], flipped["ALT"], flipped["FORMAT"]["ES"], flipped["FORMAT"]["AF"]), ("C", "G", "0.1", "0.15"))
            self.assertEqual(marks["opposite_sides"], "palindrome_complemented")

        with self.subTest("strand is inferred per variant (TwoSampleMR action 2), not from 1 file-wide assumption"):
            self.assertEqual((marks["same_side"], marks["outside_tolerance"]), ("palindrome_forward", "palindrome_complemented"))
            self.assertEqual((detail["kept_same_strand"], detail["flipped"]), (1, 2))

        with self.subTest("either frequency within 0.08 of 0.5: kept as reported and marked ambiguous (palindrome_ambiguous)"):
            self.assertEqual((marks["study_near_half"], marks["reference_near_half"]), ("palindrome_ambiguous", "palindrome_ambiguous"))
            self.assertEqual((rows["study_near_half"]["REF"], rows["study_near_half"]["ALT"]), ("T", "A"))
            self.assertEqual(detail["ambiguous_frequency_near_half"], 2)

        with self.subTest("a palindrome without a 1000 Genomes frequency is ambiguous"):
            self.assertEqual(marks["no_reference"], "palindrome_ambiguous")
            self.assertEqual(detail["ambiguous_no_reference"], 1)

        with self.subTest("a reference-filled frequency never drives a flip"):
            filled = rows["filled"]
            self.assertEqual((filled["REF"], filled["ALT"], filled["FORMAT"]["AFS"], marks["filled"]),
                             ("T", "A", "1000G_EUR", "palindrome_ambiguous"))
            self.assertEqual(detail["ambiguous_no_study_eaf"], 1)

        with self.subTest("non-palindromic SNPs are not touched by this step"):
            self.assertFalse(marks["not_palindromic"].startswith("palindrome_"))

        with self.subTest("the tolerance can be widened with palindromic_tolerance"):
            manifests = tmp / "manifests"
            manifests.mkdir()
            (manifests / "test_2026_1.yaml").write_text("palindromic_tolerance: 0.15\n", encoding="utf-8")
            wide = run_standardise(tmp, raw, reference, config_root=manifests, input_build="hg38", ancestry="EUR", name="wide")
            self.assertEqual(wide.returncode, 0, wide.stderr)
            self.assertEqual(hm(wide)["outside_tolerance"], "palindrome_ambiguous")
            self.assertEqual(wide.qc["palindromic_resolution"]["ambiguous_frequency_near_half"], 3)

        with self.subTest("without 1000 Genomes alleles every palindrome is ambiguous"):
            empty = frequencies_reference(tmp / "empty_reference", [])
            bare = run_standardise(tmp, raw, empty, input_build="hg38", ancestry="EUR", name="bare")
            self.assertEqual(bare.returncode, 0, bare.stderr)
            marks_bare = hm(bare)
            self.assertEqual({marks_bare[r[2]] for r in TABLE if r[2] != "not_palindromic"}, {"palindrome_ambiguous"})
            self.assertEqual(len(bare.rows), len(TABLE))
            self.assertEqual(bare.qc["palindromic_resolution"]["ambiguous_no_reference"], 7)

        # --- after a minus-strand liftover ----------------------------------------------------------------------
        rng = np.random.default_rng(821)
        g37 = list("".join(rng.choice(list("ACGT"), size=3000)))
        g37[1501 - 1], g37[1601 - 1] = "A", "T"
        g38 = list("".join(rng.choice(list("ACGT"), size=3000)))
        for o in range(1000):
            g38[3000 - 1 - (500 + o)] = COMPLEMENT[g37[1000 + o]]
        lifted_ref = make_reference(tmp / "liftover_reference", lookup_rows=[], offset_chains=False)
        write_chain(lifted_ref / "genome" / "hg19ToHg38.over.chain.gz",
                    [chain_block("chr1", "chr1", 3000, 3000, 1000, 1000, 500, 1, dst_strand="-")])
        write_fasta(lifted_ref / "genome", {"1": "".join(g37)}, name="GRCh37.fa")
        write_fasta(lifted_ref / "genome", {"1": "".join(g38)}, name="GRCh38.fa")
        # GRCh38 2000 is T (the study's A complemented); 1000 Genomes T = 0.2, the study's A = 0.2: same side
        frequencies_reference(lifted_ref, [(2000, "T", "A", 0.8), (1900, "A", "T", 0.3)])
        raw = write_table(tmp / "b37.tsv", HEADER, [[1, 1501, "resolved", "A", "T", 0.2, 0.1, 0.02, 1e-6, 5000],
                                                    [1, 1601, "ambiguous", "A", "T", 0.5, 0.1, 0.02, 1e-6, 5000]])
        minus = run_standardise(tmp, raw, lifted_ref, input_build="hg19", ancestry="EUR", name="minus")
        self.assertEqual(minus.returncode, 0, minus.stderr)
        lifted = minus.by_snpid()

        with self.subTest("a palindrome resolved after a minus-strand liftover reads complemented (palindrome_complemented)"):
            row = lifted["resolved"]
            self.assertEqual((row["POS"], row["REF"], row["ALT"]), (2000, "T", "A"))
            # the study's A is GRCh38 T, the reference base, so the alleles are then swapped and ES negated
            self.assertEqual((row["INFO"]["HM"], row["FORMAT"]["ES"]), ("palindrome_complemented_swapped", "-0.1"))

        with self.subTest("an ambiguous palindrome complemented by liftover is palindrome_ambiguous_complemented"):
            row = lifted["ambiguous"]
            self.assertEqual((row["POS"], row["REF"], row["ALT"], row["FORMAT"]["ES"]), (1900, "A", "T", "0.1"))
            self.assertEqual(row["INFO"]["HM"], "palindrome_ambiguous_complemented")


if __name__ == "__main__":
    unittest.main()
