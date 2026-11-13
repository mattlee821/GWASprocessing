"""README step 9, Reference alignment: REF must match the GRCh38 genome (keep, swap, complement or remove), the effect
size and frequency follow a swap, indels are left-aligned and trimmed, and INFO/HM records what was done relative to the
alleles the study reported, also when liftover changed the strand. GWAS-VCF inputs are aligned the same way.

The fixture genome is seeded (821); only the bases named in each fixture matter.
"""

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import (bcftools_liftover_available, chain_block, make_reference, run_standardise, write_chain,  # noqa: I001
                     write_fasta, write_table)

import gwas_alleles as alleles
from standardise_gwas import harmonisation_action

COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}
HEADER = ["CHR", "POS38", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
GRCH38_CONTIG = "##contig=<ID=1,length=248956422,assembly=GRCh38>"


def sequence(length: int, fixed: dict[int, str], seed: int = 821) -> str:
    """A random sequence with the given bases written at 1-based positions."""
    rng = np.random.default_rng(seed)
    bases = list("".join(rng.choice(list("ACGT"), size=length)))
    for pos, text in fixed.items():
        bases[pos - 1:pos - 1 + len(text)] = list(text)
    return "".join(bases)


def write_vcf(path: Path, records: list[tuple], fmt: str = "ES:SE:LP:AF:EZ:SS") -> Path:
    """A GRCh38 GWAS-VCF with records (chrom, pos, id, REF, ALT, sample values)."""
    formats = [f'##FORMAT=<ID={key},Number=A,Type=Float,Description="{key}">' for key in fmt.split(":")]
    lines = ["##fileformat=VCFv4.2", GRCH38_CONTIG, *formats, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBMI"]
    lines += [f"{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t.\t{fmt}\t{v}" for c, p, i, r, a, v in records]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def variants(result) -> dict[tuple, dict]:
    return {(r["CHROM"], r["POS"]): r for r in result.rows}


class ReferenceAlignmentStepTest(unittest.TestCase):
    def test_reference_alignment(self):
        tmp = Path(tempfile.mkdtemp())

        # Chromosome 1 of the GRCh38 FASTA (3,000 bases); chromosome 2 is absent from it.
        genome = sequence(3000, {100: "G", 200: "G", 300: "G", 400: "G", 500: "G", 1001: "GCACACATT", 2001: "GCACACATT",
                                 2501: "TAG", 2700: "CCC", 600: "N", 700: "GNN"})
        reference = make_reference(tmp / "reference", lookup_rows=[], fasta={"1": genome})
        rows = [
            [1, 100, "keep", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],     # OA is the GRCh38 base
            [1, 200, "swap", "G", "A", 0.3, 0.2, 0.02, 1e-6, 5000],     # EA is the GRCh38 base
            [1, 300, "comp", "T", "C", 0.3, 0.3, 0.02, 1e-6, 5000],     # reverse strand: complement of OA is the base
            [1, 400, "compswap", "C", "T", 0.3, 0.4, 0.02, 1e-6, 5000],  # reverse strand and EA is the base
            [1, 500, "neither", "A", "T", 0.3, 0.5, 0.02, 1e-6, 5000],  # matches neither strand of G
            [1, 1005, "del", "A", "ACA", 0.3, 0.6, 0.02, 1e-6, 5000],   # CA deletion written at the right of the repeat
            [1, 2003, "ins", "ACAC", "AC", 0.3, 0.7, 0.02, 1e-6, 5000],  # over-long CA insertion
            [1, 2502, "indelswap", "A", "AT", 0.3, 0.8, 0.02, 1e-6, 5000],  # EA is the GRCh38 sequence
            [1, 2700, "indelbad", "A", "GTT", 0.3, 0.9, 0.02, 1e-6, 5000],  # neither indel allele matches the sequence
            [1, 5000, "beyond", "A", "G", 0.3, 1.0, 0.02, 1e-6, 5000],  # past the end of the FASTA
            [1, 600, "nbase", "A", "G", 0.3, 1.2, 0.02, 1e-6, 5000],    # the GRCh38 base is N
            [1, 700, "nindel", "A", "GTT", 0.3, 1.3, 0.02, 1e-6, 5000],  # the GRCh38 sequence holds an N
            [2, 100, "nocontig", "A", "G", 0.3, 1.1, 0.02, 1e-6, 5000],  # chromosome absent from the FASTA
        ]
        raw = write_table(tmp / "gwas.tsv", HEADER, rows)
        result = run_standardise(tmp, raw, reference, name="aligned")
        self.assertEqual(result.returncode, 0, result.stderr)
        out = variants(result)
        written = lambda key: (out[key]["REF"], out[key]["ALT"], out[key]["FORMAT"]["ES"], out[key]["FORMAT"]["AF"],  # noqa: E731
                               out[key]["INFO"]["HM"])
        summary = result.qc["reference_alignment"]

        with self.subTest("REF is the GRCh38 base and ALT the other allele: kept as reported (HM forward)"):
            self.assertEqual(written(("1", 100)), ("G", "A", "0.1", "0.3", "forward"))
            self.assertTrue(summary["performed"])
            self.assertEqual(summary["kept"], 3)  # the SNV and the 2 indels whose other allele matched

        with self.subTest("effect allele is the reference base: alleles swapped, ES negated, AF to 1 - AF (HM swapped)"):
            self.assertEqual(written(("1", 200)), ("G", "A", "-0.2", "0.7", "swapped"))
            self.assertEqual(summary["swapped"], 3)
            self.assertEqual(result.sample["SwitchedAlleles"], "3")

        with self.subTest("non-palindromic SNP reported on the reverse strand is complemented (HM complemented)"):
            self.assertEqual(written(("1", 300)), ("G", "A", "0.3", "0.3", "complemented"))

        with self.subTest("reverse strand and effect allele the reference base: complemented and swapped"):
            self.assertEqual(written(("1", 400)), ("G", "A", "-0.4", "0.7", "complemented_swapped"))
            self.assertEqual(summary["strand_flipped"], 2)

        with self.subTest("a SNV matching neither strand is removed"):
            self.assertNotIn(("1", 500), out)
            self.assertEqual(summary["dropped_mismatch"], 1)
            self.assertNotIn("neither", [r["INFO"]["SNPID"] for r in result.rows])

        with self.subTest("indels are checked against the reference, left-aligned and trimmed"):
            self.assertEqual(written(("1", 1001))[:2], ("GCA", "G"))
            self.assertEqual((out[("1", 1001)]["INFO"]["SNPID"], out[("1", 1001)]["INFO"]["HM"]), ("del", "forward"))
            self.assertEqual(written(("1", 2001))[:2], ("G", "GCA"))
            self.assertEqual(out[("1", 2001)]["INFO"]["SNPID"], "ins")
            self.assertEqual(result.qc["indels_left_aligned"], 2)
            self.assertEqual(written(("1", 2502)), ("A", "AT", "-0.8", "0.7", "swapped"))
            self.assertNotIn(("1", 2700), out)
            self.assertEqual(summary["indels_unmatched_dropped"], 1)

        with self.subTest("no reference sequence at the position: kept as reported (HM unverified)"):
            self.assertEqual(written(("1", 5000)), ("G", "A", "1", "0.3", "unverified"))
            self.assertEqual(written(("2", 100)), ("G", "A", "1.1", "0.3", "unverified"))
            self.assertEqual(summary["unalignable_kept"], 4)

        with self.subTest("an N in the GRCh38 sequence counts as no reference sequence (HM unverified)"):
            self.assertEqual(written(("1", 600)), ("G", "A", "1.2", "0.3", "unverified"))
            self.assertEqual(written(("1", 700)), ("GTT", "A", "1.3", "0.3", "unverified"))

        with self.subTest("the step log counts the rows alignment removed"):
            step = next(s for s in result.qc["steps"] if s["title"] == "GRCh38 reference alignment")
            self.assertEqual(step["removed"], 2)
            self.assertEqual(sum(step["removed_by"].values()), 2)
            self.assertEqual(len(result.rows), len(rows) - 2)

        with self.subTest("alignment is reported as skipped without a FASTA, and every variant is unverified"):
            bare = run_standardise(tmp, raw, make_reference(tmp / "no_fasta", lookup_rows=[]), name="no_fasta")
            self.assertEqual(bare.returncode, 0, bare.stderr)
            self.assertFalse(bare.qc["reference_alignment"]["performed"])
            self.assertTrue(any("no FASTA" in w for w in bare.qc["warnings"]))
            self.assertNotIn("GRCh38 reference alignment", [s["title"] for s in bare.qc["steps"]])
            self.assertEqual({r["INFO"]["HM"] for r in bare.rows if r["INFO"]["SNPID"] != "neither"}, {"unverified"})
            as_reported = variants(bare)[("1", 200)]
            self.assertEqual((as_reported["REF"], as_reported["ALT"], as_reported["FORMAT"]["ES"]), ("A", "G", "0.2"))
            self.assertIn(("1", 500), variants(bare))  # nothing to check it against, so nothing is removed

        with self.subTest("GWAS-VCF input is aligned the same way, ES, EZ and AF follow a swap, counts add up over pieces"):
            vcf_genome = sequence(3500, {1000: "A", 2000: "A", 3000: "A"})
            vcf_reference = make_reference(tmp / "vcf_reference", lookup_rows=[], fasta={"1": vcf_genome})
            vcf = write_vcf(tmp / "align.gwas.vcf", [("1", 1000, "rs1", "G", "A", "0.1:0.02:6:0.3:5:5000"),
                                                     ("1", 2000, "rs2", "G", "A", "0.2:0.02:6:0.4:10:5000"),
                                                     ("1", 3000, "rs3", "T", "C", "0.3:0.02:6:0.3:15:5000")])
            aligned = run_standardise(tmp, vcf, vcf_reference, name="vcf", env={**os.environ, "GWAS_TABLE_CHUNK": "1"})
            self.assertEqual(aligned.returncode, 0, aligned.stderr)
            self.assertEqual(aligned.qc["table_pieces"]["pieces"], 3)
            got = {r["POS"]: (r["REF"], r["ALT"], r["FORMAT"]["ES"], r["FORMAT"]["AF"], r["FORMAT"]["EZ"], r["FORMAT"]["SE"],
                              r["INFO"]["HM"]) for r in aligned.rows}
            self.assertEqual(got, {1000: ("A", "G", "-0.1", "0.7", "-5", "0.02", "swapped"),
                                   2000: ("A", "G", "-0.2", "0.6", "-10", "0.02", "swapped"),
                                   3000: ("A", "G", "0.3", "0.3", "15", "0.02", "complemented")})
            self.assertEqual((aligned.qc["reference_alignment"]["swapped"], aligned.qc["reference_alignment"]["strand_flipped"]), (2, 1))

        with self.subTest("HM is relative to the alleles the study reported, when combined with a liftover strand change"):
            # every combination of palindromic decision, alignment action and liftover complement
            pal = np.array(["", "", "", "", "", "", "", "forward", "forward", "complemented", "ambiguous"], dtype=object)
            action = np.array(["keep", "complement", "keep", "swap", "complement", "drop", "indel_unmatched", "keep", "swap",
                               "keep", "keep"], dtype=object)
            changed = np.array([False, False, True, True, True, False, False, False, True, True, True])
            self.assertEqual(harmonisation_action(pal, action, changed).tolist(), [
                "forward", "complemented", "complemented", "complemented_swapped", "forward", "dropped", "unverified",
                "palindrome_forward", "palindrome_complemented_swapped", "palindrome_forward",
                "palindrome_ambiguous_complemented"])
            self.assertEqual(harmonisation_action(np.array(["ambiguous"], dtype=object), np.array(["swap"], dtype=object),
                                                  np.array([False])).tolist(), ["palindrome_ambiguous_swapped"])
            self.assertEqual(harmonisation_action(np.array(["", ""], dtype=object),
                                                  np.array(["complement_swap", "swap"], dtype=object)).tolist(),
                             ["complemented_swapped", "swapped"])

        with self.subTest("a minus-strand liftover block complements the alleles and HM says so; 2 complements cancel"):
            if not bcftools_liftover_available():
                self.skipTest("bcftools with the +liftover plugin is required")
            self.check_liftover_strand(tmp)

        with self.subTest("the GRCh38 FASTA is found under the reference root and read by position"):
            root = Path(tempfile.mkdtemp())
            self.assertIsNone(alleles.Fasta.find(root))
            (root / "genome").mkdir()
            (root / "genome" / "._hg38.fa").write_bytes(b"\x00\x05\x16\x07")  # macOS AppleDouble file
            self.assertIsNone(alleles.Fasta.find(root))
            write_fasta(root / "genome", {"chr1": "ACGTACGTAC" * 13, "chrM": "GATTACA"}, width=7)
            fasta = alleles.Fasta.find(root)
            self.assertEqual(fasta.path.name, "GRCh38.fa")
            self.assertTrue((root / "genome" / "GRCh38.fa.fai").exists())
            self.assertEqual(fasta.bases("1", np.array([1, 2, 8, 9, 130])).tolist(), ["A", "C", "T", "A", "C"])
            self.assertEqual(fasta.bases("1", np.array([3]), length=4).tolist(), ["GTAC"])
            self.assertEqual(fasta.bases("MT", np.array([4])).tolist(), ["T"])
            self.assertEqual(fasta.bases("1", np.array([131, 0])).tolist(), ["", ""])
            self.assertEqual(fasta.bases("2", np.array([1])).tolist(), [""])

        with self.subTest("alignment actions for SNVs, palindromes and indels"):
            write_fasta(tmp / "unit", {"1": "GGGGGGGGAAGTCC"})
            fasta = alleles.Fasta(tmp / "unit" / "GRCh38.fa")
            cases = [(1, "A", "G", "keep"), (2, "G", "A", "swap"), (3, "T", "C", "complement"), (4, "C", "T", "complement_swap"),
                     (5, "A", "T", "drop"), (9, "T", "A", "keep"), (10, "A", "T", "swap"), (11, "C", "G", "keep"),
                     (12, "AT", "T", "keep"), (13, "C", "CA", "swap"), (6, "GTT", "A", "indel_unmatched")]
            result_actions = alleles.reference_alignment(pd.Series(["1"] * len(cases)), pd.Series([c[0] for c in cases]),
                                                         pd.Series([c[1] for c in cases]), pd.Series([c[2] for c in cases]), fasta)
            self.assertEqual(list(result_actions["action"]), [c[3] for c in cases])

    def check_liftover_strand(self, tmp: Path) -> None:
        """GRCh37 chr1 [1000, 2000) maps to the minus strand of GRCh38 chr1 (GRCh37 1501 is GRCh38 2000); chr2 is
        GRCh38 = GRCh37 + 100 on the plus strand."""
        rng = np.random.default_rng(821)
        genome37 = list("".join(rng.choice(list("ACGT"), size=3000)))
        genome38 = list("".join(rng.choice(list("ACGT"), size=3000)))
        for o in range(1000):
            genome38[3000 - 1 - (500 + o)] = COMPLEMENT[genome37[1000 + o]]
        chr2_37 = "".join(rng.choice(list("ACGT"), size=3000))
        chr2_38 = "".join(rng.choice(list("ACGT"), size=100)) + chr2_37
        reference = make_reference(tmp / "minus_reference", lookup_rows=[], offset_chains=False)
        write_chain(reference / "genome" / "hg19ToHg38.over.chain.gz",
                    [chain_block("chr1", "chr1", 3000, 3000, 1000, 1000, 500, 1, dst_strand="-"),
                     chain_block("chr2", "chr2", 3000, 3100, 0, 2900, 100, 2)])
        write_fasta(reference / "genome", {"1": "".join(genome37), "2": chr2_37}, name="GRCh37.fa")
        write_fasta(reference / "genome", {"1": "".join(genome38), "2": chr2_38}, name="GRCh38.fa")

        def alt_of(ref: str) -> str:
            return next(b for b in "ACGT" if b not in (ref, COMPLEMENT[ref]))

        r1, r3, r4 = genome37[1501 - 1], genome37[1511 - 1], genome37[1521 - 1]
        r2 = chr2_37[1000 - 1]
        rows = [[1, 1501, "minus", alt_of(r1), r1, 0.3, 0.1, 0.02, 1e-6],  # as reported, lifted across a minus block
                [2, 1000, "plus_reverse", COMPLEMENT[alt_of(r2)], COMPLEMENT[r2], 0.3, 0.1, 0.02, 1e-6],
                [1, 1511, "minus_reverse", COMPLEMENT[alt_of(r3)], COMPLEMENT[r3], 0.3, 0.1, 0.02, 1e-6],
                [1, 1521, "minus_swap", r4, alt_of(r4), 0.3, 0.1, 0.02, 1e-6]]  # effect allele is the reference base
        raw = write_table(tmp / "minus.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P"], rows)
        result = run_standardise(tmp, raw, reference, input_build="hg19", ancestry="EUR", name="minus")
        self.assertEqual(result.returncode, 0, result.stderr)
        got = {r["INFO"]["SNPID"]: (r["POS"], r["REF"], r["ALT"], r["FORMAT"]["ES"], r["INFO"]["HM"]) for r in result.rows}
        self.assertEqual(got["minus"], (2000, COMPLEMENT[r1], COMPLEMENT[alt_of(r1)], "0.1", "complemented"))
        self.assertEqual(got["plus_reverse"], (1100, r2, alt_of(r2), "0.1", "complemented"))
        self.assertEqual(got["minus_reverse"], (1990, COMPLEMENT[r3], COMPLEMENT[alt_of(r3)], "0.1", "forward"))
        self.assertEqual(got["minus_swap"], (1980, COMPLEMENT[r4], COMPLEMENT[alt_of(r4)], "-0.1", "complemented_swapped"))


if __name__ == "__main__":
    unittest.main()
