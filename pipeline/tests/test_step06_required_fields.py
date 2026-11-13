"""README step 6, required fields: rows without a GRCh38 position, 2 different alleles made of A, C, G and T, an effect,
a positive SE or a P value are removed, each counted against the first cause that applies; nothing else is filtered.

The table is declared GRCh38, so no liftover is involved, and the reference has no FASTA, so alignment leaves the alleles
as reported.
"""

import tempfile
import unittest
from pathlib import Path

from testlib import EXAMPLES, example_reference, make_reference, run_standardise, write_table

HEADER = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "INFO", "N"]


def row(chrom, pos, snp, ea="A", oa="G", eaf=0.3, beta=0.1, se=0.02, p=1e-6, info=0.95):
    return [chrom, pos, snp, ea, oa, eaf, beta, se, p, info, 5000]


class Step06RequiredFieldsTest(unittest.TestCase):
    def test_step06_required_fields(self):
        tmp = Path(tempfile.mkdtemp())
        reference = make_reference(tmp / "reference", lookup_rows=[("1", "rs10", 900, "A", "G", 0.3)])  # GRCh38 1000
        kept = [
            row(1, 1000, "plain"),
            row(1, 1000, "second_alt", ea="C"),  # a second ALT allele at the same position
            row(1, 1100, "rare", eaf=0.004),
            row(1, 1200, "no_eaf", eaf=None),
            row(1, 1300, "indel", ea="AT", oa="A"),
            row(1, 1400, "palindrome", ea="A", oa="T"),
            row(1, 1500, "low_info", info=0.1),
            row(1, 1600, "lower_case", ea="a", oa="g"),
            row(1, 1700, "p_missing", p=None),  # P is recovered from BETA/SE (step 3), so the row is complete
            row("X", 1000, "chrX"), row("Y", 1000, "chrY"), row("MT", 100, "chrMT"),
        ]
        removed = [
            row("", "", "rs999999"),  # rsid_not_in_coordinate_index
            row("", "", "marker_without_position"),  # no_position_or_rsid
            row(1, 2000, "missing_ea", ea=""),  # missing_allele
            row(1, 2100, "id_codes", ea="D", oa="I"),  # invalid_allele
            row(1, 2200, "dash", ea="-", oa="A"),  # invalid_allele
            row(1, 2300, "same_alleles", ea="A", oa="A"),  # invalid_allele
            row(1, 2400, "n_allele", ea="A", oa="N"),  # invalid_allele
            row("6_ssto_hap7", 1000, "alt_contig"),  # chromosome_not_primary
            row(1, 0, "position_zero"),  # invalid_position
            row(1, 2500, "no_beta", beta=None),  # missing_beta
            row(1, 2600, "no_beta_zero_se", beta=None, se=0),  # missing_beta only: the first cause that applies
            row(1, 2700, "zero_se", se=0),  # missing_or_invalid_se
            row(1, 2800, "negative_se", se=-0.02),  # missing_or_invalid_se
            row(1, 2900, "p_above_one", p=1.5),  # missing_or_invalid_p
            row(1, 3000, "p_negative", p=-0.1),  # missing_or_invalid_p
            row(1, 3100, "id_codes_no_beta", ea="D", oa="I", beta=None),  # invalid_allele, not missing_beta
        ]
        raw = write_table(tmp / "gwas.tsv", HEADER, kept + removed)
        result = run_standardise(tmp, raw, reference, input_build="hg38", ancestry="EUR")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = result.by_snpid()
        causes = {name: detail["rows"] for name, detail in result.qc["incomplete_by_cause"].items()}

        with self.subTest("rows without a GRCh38 position are removed"):
            self.assertEqual((causes["rsid_not_in_coordinate_index"], causes["no_position_or_rsid"], causes["invalid_position"]),
                             (1, 1, 1))

        with self.subTest("rows without 2 different alleles made of A, C, G and T are removed"):
            self.assertEqual((causes["missing_allele"], causes["invalid_allele"]), (1, 5))

        with self.subTest("rows without an effect are removed"):
            self.assertEqual(causes["missing_beta"], 2)

        with self.subTest("rows without a positive SE are removed"):
            self.assertEqual(causes["missing_or_invalid_se"], 2)

        with self.subTest("rows without a valid P value are removed"):
            self.assertEqual(causes["missing_or_invalid_p"], 2)

        with self.subTest("non-primary contigs are dropped"):
            self.assertEqual(causes["chromosome_not_primary"], 1)
            self.assertNotIn("alt_contig", out)

        with self.subTest("each removed row is counted once, against the first cause that applies"):
            self.assertEqual(sum(causes.values()), len(removed))
            self.assertEqual(result.qc["rows_removed_missing_required"], len(removed))
            self.assertFalse({r[2] for r in removed} & set(out))
            steps = {s["title"]: s for s in result.qc["steps"]}
            self.assertEqual(steps["Alleles and GRCh38 coordinates"]["removed"] + steps["Statistics"]["removed"], len(removed))
            for step in result.qc["steps"]:
                self.assertEqual(step["rows_in"] - step["removed"], step["rows_out"])
                self.assertEqual(sum(step["removed_by"].values()), step["removed"])

        with self.subTest("nothing else is filtered: rare variants, indels, palindromes, low imputation quality, "
                          "missing frequency, sex and mitochondrial chromosomes and a second ALT allele are kept"):
            self.assertEqual(sorted(out), sorted(r[2] for r in kept))
            self.assertEqual(len(result.rows), len(kept))
            self.assertEqual((out["indel"]["REF"], out["indel"]["ALT"]), ("A", "AT"))
            self.assertEqual((out["lower_case"]["REF"], out["lower_case"]["ALT"]), ("G", "A"))
            self.assertEqual(out["rare"]["FORMAT"]["AF"], "0.004")

        with self.subTest("low retention is reported as a warning"):
            self.assertTrue(any(f"retained {len(kept)} of {len(kept) + len(removed)}" in w for w in result.qc["warnings"]),
                            result.qc["warnings"])

        with self.subTest("no usable rows is a data error"):
            # "Estimate" is not a recognised effect-size column, so no row has an effect
            raw = write_table(tmp / "no_beta.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "Estimate", "SE", "P", "N"],
                              [[1, 9000, "rs90", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
            empty = run_standardise(tmp, raw, reference, input_build="hg38", ancestry="EUR", name="empty")
            self.assertEqual(empty.returncode, 65)
            self.assertIn("No variants remain", empty.stderr)
            self.assertIn("BETA", empty.stderr)

        with self.subTest("missing ancestry fails loudly"):
            missing = run_standardise(tmp, EXAMPLES / "tables/generic.tsv", example_reference(), ancestry="", name="ancestry")
            self.assertEqual(missing.returncode, 65)
            self.assertIn("ancestry", missing.stderr.lower())


if __name__ == "__main__":
    unittest.main()
