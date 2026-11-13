"""README step 7, 1000 Genomes annotation: each variant is matched on GRCh38 position and alleles (either strand) to the
frequencies of the manifest population; FORMAT/AF is filled from 1000 Genomes where the study reports none, with
FORMAT/AFS naming the source; the frequency check stops a row whose effect-allele frequencies anti-correlate with
1000 Genomes over at least 1,000 non-palindromic SNPs, and warns when the agreement is weak.

Tables are declared GRCh38 and the references have no FASTA, so the alleles are written as reported (no swap to the
GRCh38 reference changes AF).
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import run_standardise, write_grch38_references, write_table

HEADER = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
# (pos38, REF, ALT, EUR ALT frequency, EAS ALT frequency); EAS differs everywhere, so the population used is visible
SITES = [(1000, "G", "A", 0.2, 0.6), (1000, "G", "T", 0.05, 0.15), (1100, "G", "A", 0.87449801, 0.3), (1200, "C", "T", 0.9, 0.5),
         (1300, "G", "A", 0.4, 0.7), (1500, "G", "A", 0.25, 0.35)]


def populations_reference(root: Path) -> Path:
    write_grch38_references(root, [("1", p, 0, r, a, eur) for p, r, a, eur, _ in SITES], populations=("EUR", "ALL"))
    write_grch38_references(root, [("1", p, 0, r, a, eas) for p, r, a, _, eas in SITES], populations=("EAS",))
    return root


class Step07ThousandGenomesTest(unittest.TestCase):
    def test_step07_1000genomes_annotation(self):
        tmp = Path(tempfile.mkdtemp())
        reference = populations_reference(tmp / "reference")
        raw = write_table(tmp / "gwas.tsv", HEADER, [
            [1, 1000, "fill_alt", "A", "G", None, 0.1, 0.02, 1e-6, 5000],  # effect allele = 1000 Genomes ALT
            [1, 1000, "second_record", "T", "G", None, 0.1, 0.02, 1e-6, 5000],  # a 2nd 1000 Genomes record at 1000
            [1, 1100, "fill_ref", "G", "A", None, 0.1, 0.02, 1e-6, 5000],  # effect allele = 1000 Genomes REF
            [1, 1200, "study", "T", "C", 0.1234567, 0.1, 0.02, 1e-6, 5000],  # the study's frequency, far from 1000 Genomes
            [1, 1300, "reverse_strand", "T", "C", None, 0.1, 0.02, 1e-6, 5000],  # G/A reported on the reverse strand
            [1, 1500, "other_alleles", "T", "G", None, 0.1, 0.02, 1e-6, 5000],  # position matches, alleles do not
            [1, 1600, "other_position", "A", "G", None, 0.1, 0.02, 1e-6, 5000],
        ])
        eur = run_standardise(tmp, raw, reference, input_build="hg38", ancestry="EUR", name="eur")
        self.assertEqual(eur.returncode, 0, eur.stderr)
        rows = eur.by_snpid()
        fmt = {snp: (row["FORMAT"]["AF"], row["FORMAT"]["AFS"]) for snp, row in rows.items()}

        with self.subTest("frequencies are those of the manifest population"):
            self.assertEqual(fmt["fill_alt"], ("0.2", "1000G_EUR"))
            self.assertEqual((eur.qc["reference_population"], eur.qc["reference_population_used"]), ("EUR", "EUR"))
            eas = run_standardise(tmp, raw, reference, input_build="hg38", ancestry="EAS", name="eas")
            self.assertEqual(eas.returncode, 0, eas.stderr)
            self.assertEqual({s: (r["FORMAT"]["AF"], r["FORMAT"]["AFS"]) for s, r in eas.by_snpid().items()
                              if s in ("fill_alt", "fill_ref")}, {"fill_alt": ("0.6", "1000G_EAS"), "fill_ref": ("0.7", "1000G_EAS")})

        with self.subTest("a manifest population column overrides the population from ancestry"):
            chosen = run_standardise(tmp, raw, reference, input_build="hg38", ancestry="EUR", population="eas", name="population")
            self.assertEqual(chosen.returncode, 0, chosen.stderr)
            self.assertEqual((chosen.by_snpid()["fill_alt"]["FORMAT"]["AF"], chosen.by_snpid()["fill_alt"]["FORMAT"]["AFS"]),
                             ("0.6", "1000G_EAS"))
            self.assertEqual(chosen.qc["reference_population"], "EAS")
            self.assertEqual(chosen.sample["Ancestry"], "EUR")

        with self.subTest("variants are matched on GRCh38 position and alleles, on either strand"):
            self.assertEqual(fmt["reverse_strand"], ("0.4", "1000G_EUR"))
            self.assertEqual(fmt["second_record"], ("0.05", "1000G_EUR"))  # its own allele match, not the G/A record
            self.assertEqual(fmt["other_alleles"], (".", "."))
            self.assertEqual(fmt["other_position"], (".", "."))

        with self.subTest("AF is filled only where the study has none, as the effect allele's frequency"):
            self.assertEqual(fmt["fill_ref"][1], "1000G_EUR")
            self.assertAlmostEqual(float(fmt["fill_ref"][0]), 1 - 0.87449801, places=6)
            self.assertEqual(eur.qc["reference_fills"]["EAF"], 4)
            self.assertEqual(eur.sample["AFFilled"], "4")

        with self.subTest("a filled AF is written with 6 significant digits"):
            self.assertEqual(fmt["fill_ref"][0], "0.125502")  # float32 0.87449801 would otherwise carry noise digits

        with self.subTest("the study's own frequencies are kept unchanged, with AFS study"):
            self.assertEqual(fmt["study"], ("0.1234567", "study"))

        # --- frequency check -------------------------------------------------------------------------------------
        rng = np.random.default_rng(821)
        anti = np.round(rng.uniform(0.05, 0.95, 1000), 4)
        palindromes = np.round(rng.uniform(0.05, 0.95, 200), 4)
        filled = np.round(rng.uniform(0.05, 0.95, 200), 4)
        check_ref = write_grch38_references(
            tmp / "check_reference",
            [("1", 10000 + 10 * i, 0, "G", "A", f) for i, f in enumerate(anti)]
            + [("1", 30000 + 10 * j, 0, "A", "T", f) for j, f in enumerate(palindromes)]
            + [("1", 50000 + 10 * k, 0, "G", "A", f) for k, f in enumerate(filled)])
        # the effect allele is labelled G, but the frequency column describes A: the swapped-columns error
        mislabelled = [[1, 10000 + 10 * i, f"anti{i}", "G", "A", f, 0.1, 0.02, 0.5, 5000] for i, f in enumerate(anti)]
        pal_rows = [[1, 30000 + 10 * j, f"pal{j}", "T", "A", f, 0.1, 0.02, 0.5, 5000] for j, f in enumerate(palindromes)]
        fill_rows = [[1, 50000 + 10 * k, f"fill{k}", "A", "G", None, 0.1, 0.02, 0.5, 5000] for k in range(len(filled))]

        with self.subTest("effect-allele frequencies anti-correlated with 1000 Genomes over at least 1,000 SNPs stop the row"):
            raw = write_table(tmp / "anti.tsv", HEADER, mislabelled)
            failed = run_standardise(tmp, raw, check_ref, input_build="hg38", ancestry="EUR", name="anti")
            self.assertEqual(failed.returncode, 65)
            self.assertIn("anti-correlate", failed.stderr)

        with self.subTest("the frequency check excludes filled and palindromic rows, and a weak correlation is a warning"):
            # 999 mislabelled SNPs: 1 short of the threshold, unless the palindromic or filled rows were counted
            raw = write_table(tmp / "excluded.tsv", HEADER, mislabelled[:999] + pal_rows + fill_rows)
            warned = run_standardise(tmp, raw, check_ref, input_build="hg38", ancestry="EUR", name="excluded")
            self.assertEqual(warned.returncode, 0, warned.stderr)
            self.assertEqual(warned.qc["af_concordance"]["n"], 999)
            self.assertLess(warned.qc["af_concordance"]["pearson_r"], -0.99)
            self.assertEqual(warned.qc["reference_fills"]["EAF"], 200)
            self.assertTrue(any("weak agreement" in w for w in warned.qc["warnings"]), warned.qc["warnings"])


if __name__ == "__main__":
    unittest.main()
