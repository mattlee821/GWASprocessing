"""README step 11, rsIDs: the ID column (and FORMAT/ID) is the dbSNP rsID matched on GRCh38 position and both alleles, or
`.` when there is none; variants without an rsID are kept, and INFO/SNPID keeps the study's own identifier.

The fixture genome is seeded (821); only the bases named below matter.
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import make_reference, run_standardise, write_grch38_references, write_table  # noqa: I001

import gwas_reference

HEADER = ["CHR", "POS38", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]


def sequence(length: int, fixed: dict[int, str], seed: int = 821) -> str:
    rng = np.random.default_rng(seed)
    bases = list("".join(rng.choice(list("ACGT"), size=length)))
    for pos, text in fixed.items():
        bases[pos - 1:pos - 1 + len(text)] = list(text)
    return "".join(bases)


class RsidStepTest(unittest.TestCase):
    def test_rsids(self):
        tmp = Path(tempfile.mkdtemp())
        genome = sequence(3000, {1000: "G", 1100: "G", 1200: "G", 1300: "G", 1400: "G", 1500: "G", 1700: "GC", 1800: "G"})
        reference = make_reference(tmp / "reference", lookup_rows=[], fasta={"1": genome})
        # dbSNP on GRCh38: (chrom, pos38, rs number, REF, ALT, ALT frequency)
        write_grch38_references(reference, [("1", 1000, 101, "G", "A", 0.3), ("1", 1100, 102, "G", "A", 0.3),
                                            ("1", 1200, 103, "G", "A", 0.3), ("1", 1300, 104, "G", "A", 0.3),
                                            ("1", 1400, 106, "G", "A", 0.3), ("1", 1600, 107, "G", "A", 0.3),
                                            ("1", 1700, 108, "GC", "G", 0.3), ("1", 1800, 109, "G", "A", 0.3),
                                            ("1", 1800, 110, "G", "T", 0.3)])
        rows = [
            [1, 1000, "1:1000:A:G", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],  # as dbSNP has it
            [1, 1100, "1:1100", "G", "A", 0.7, 0.1, 0.02, 1e-6, 5000],      # effect allele is the reference base
            [1, 1200, "1:1200", "T", "C", 0.3, 0.1, 0.02, 1e-6, 5000],      # reported on the reverse strand
            [1, 1300, "rs104", "T", "G", 0.3, 0.1, 0.02, 1e-6, 5000],       # dbSNP has another allele at the position
            [1, 1400, "rs555", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],       # the study's rsID differs from dbSNP's
            [1, 1500, "1:1500", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],      # no dbSNP record at the position
            [1, 1700, "1:1700:GC:G", "G", "GC", 0.3, 0.1, 0.02, 1e-6, 5000],  # deletion
            [1, 1800, "1:1800:T:G", "T", "G", 0.3, 0.1, 0.02, 1e-6, 5000],  # second ALT of a multi-allelic dbSNP site
        ]
        result = run_standardise(tmp, write_table(tmp / "gwas.tsv", HEADER, rows), reference, ancestry="EUR")
        self.assertEqual(result.returncode, 0, result.stderr)
        by_pos = {(r["POS"], r["ALT"]): r for r in result.rows}
        ids = {key: r["ID"] for key, r in by_pos.items()}

        with self.subTest("the ID is the dbSNP rsID matched on GRCh38 position and both alleles"):
            self.assertEqual(ids[(1000, "A")], "rs101")
            self.assertEqual(ids[(1100, "A")], "rs102")  # matched after alignment, so a swapped study still matches
            self.assertEqual(ids[(1200, "A")], "rs103")  # and a complemented one
            self.assertEqual(ids[(1700, "G")], "rs108")
            self.assertEqual(ids[(1800, "T")], "rs110")

        with self.subTest("FORMAT/ID is the rsID, as in the ID column"):
            self.assertEqual({r["ID"] for r in result.rows}, {r["FORMAT"]["ID"] for r in result.rows})
            self.assertTrue(all(r["ID"] == r["FORMAT"]["ID"] for r in result.rows))

        with self.subTest("no match on position and both alleles gives '.', and the variant is kept"):
            self.assertEqual(ids[(1300, "T")], ".")
            self.assertEqual(ids[(1500, "A")], ".")
            self.assertEqual(len(result.rows), len(rows))
            self.assertEqual((result.qc["rsid_dbsnp"]["matched"], result.qc["rsid_dbsnp"]["unmatched"]), (6, 2))

        with self.subTest("the study's identifier is kept in INFO/SNPID, never copied to the ID"):
            self.assertEqual((ids[(1400, "A")], by_pos[(1400, "A")]["INFO"]["SNPID"]), ("rs106", "rs555"))
            self.assertEqual((ids[(1300, "T")], by_pos[(1300, "T")]["INFO"]["SNPID"]), (".", "rs104"))
            self.assertEqual(by_pos[(1000, "A")]["INFO"]["SNPID"], "1:1000:A:G")

        with self.subTest("a GWAS-VCF record's ID is replaced by the dbSNP rsID and kept in INFO/SNPID"):
            vcf = tmp / "ids.gwas.vcf"
            vcf.write_text("\n".join([
                "##fileformat=VCFv4.2", "##contig=<ID=1,length=248956422,assembly=GRCh38>",
                *[f'##FORMAT=<ID={k},Number=A,Type=Float,Description="{k}">' for k in ("ES", "SE", "LP", "AF", "SS")],
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBMI",
                "1\t1000\trs999\tG\tA\t.\tPASS\t.\tES:SE:LP:AF:SS\t0.1:0.02:6:0.3:5000",
                "1\t1500\trs998\tG\tA\t.\tPASS\t.\tES:SE:LP:AF:SS\t0.1:0.02:6:0.3:5000"]) + "\n", encoding="utf-8")
            from_vcf = run_standardise(tmp, vcf, reference, ancestry="EUR", name="vcf")
            self.assertEqual(from_vcf.returncode, 0, from_vcf.stderr)
            self.assertEqual([(r["POS"], r["ID"], r["FORMAT"]["ID"], r["INFO"]["SNPID"]) for r in from_vcf.rows],
                             [(1000, "rs101", "rs101", "rs999"), (1500, ".", ".", "rs998")])

        with self.subTest("the allele key in the dbSNP index ignores case and keeps allele order"):
            # the key is written into the dbSNP index; a change would silently break every rsID match
            self.assertEqual(gwas_reference.allele_hash(["A", "AT"], ["G", "A"]).tolist(),
                             gwas_reference.allele_hash(["a", "AT"], ["g", "A"]).tolist())
            self.assertNotEqual(gwas_reference.allele_hash(["A"], ["G"])[0], gwas_reference.allele_hash(["G"], ["A"])[0])


if __name__ == "__main__":
    unittest.main()
