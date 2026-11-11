"""End-to-end tests with the GRCh38 references (dbSNP index, 1000 Genomes GRCh38 frequencies, paired build sample).

The fixture genome is seeded (821). testlib's offset chains place GRCh38 at GRCh37 + 100 on chromosomes 1 and 2.
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import make_reference, run_standardise, write_grch38_references, write_table  # noqa: I001

import gwas_reference

OFFSET = 100
HEADER = ["SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]


def fixture(root: Path, n: int = 30):
    rng = np.random.default_rng(821)
    genome = "".join(rng.choice(list("ACGT"), size=7000))
    variants = []
    for i in range(n):
        pos38 = 1000 + 33 * i  # never 100 apart, so no GRCh37 position equals another GRCh38 one
        ref = genome[pos38 - 1]
        alt = next(b for b in "ACGT" if b not in (ref, {"A": "T", "C": "G", "G": "C", "T": "A"}[ref]))
        variants.append(("1", pos38, 5000 + i, ref, alt, 0.3))
    reference = make_reference(root / "reference", lookup_rows=[], fasta={"1": genome, "2": genome})
    write_grch38_references(reference, variants, paired=[(v[0], v[1] - OFFSET, v[1], v[3], v[4]) for v in variants])
    return reference, variants


class Grch38ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.reference, self.variants = fixture(self.tmp)

    def test_rsid_only_table_takes_grch38_positions_from_dbsnp_without_liftover(self):
        rows = [[f"rs{v[2]}", v[4], v[3], 0.3 if i else "", 0.1, 0.02, 1e-6, 5000] for i, v in enumerate(self.variants)]
        rows.append(["rs999999", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000])  # not in dbSNP
        result = run_standardise(self.tmp, write_table(self.tmp / "rs.tsv", HEADER, rows), self.reference, ancestry="EUR")
        self.assertEqual(result.returncode, 0, result.stderr)
        by_id = result.by_id()
        self.assertEqual(len(by_id), 30)
        first = by_id["rs5000"]
        self.assertEqual((first["CHROM"], first["POS"], first["REF"], first["ALT"]), ("1", 1000, self.variants[0][3], self.variants[0][4]))
        self.assertEqual((first["FORMAT"]["AF"], first["FORMAT"]["AFS"], first["FORMAT"]["ID"]), ("0.30000001", "1000G_EUR", "rs5000"))
        self.assertEqual(by_id["rs5001"]["FORMAT"]["AFS"], "study")
        self.assertEqual(result.qc["coordinate_source"], "dbsnp")
        self.assertEqual(result.qc["reference_fills"]["POS38"], 30)
        self.assertNotIn("hg19_to_hg38", result.qc.get("liftover", {}))
        self.assertEqual(result.qc["incomplete_by_cause"]["rsid_not_in_coordinate_index"]["rows"], 1)
        self.assertEqual(result.sample["AFFilled"], "1")

    def test_ids_come_from_dbsnp_matched_on_position_and_alleles(self):
        header = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P"]
        rows = [[1, v[1] - OFFSET, f"1:{v[1] - OFFSET}", v[4], v[3], 0.3, 0.1, 0.02, 1e-6] for v in self.variants]
        wrong = self.variants[3]
        other = next(b for b in "ACGT" if b not in (wrong[3], wrong[4], {"A": "T", "C": "G", "G": "C", "T": "A"}[wrong[3]]))
        rows[3] = [1, wrong[1] - OFFSET, "rs5003", other, wrong[3], 0.3, 0.1, 0.02, 1e-6]  # alleles differ from dbSNP's
        result = run_standardise(self.tmp, write_table(self.tmp / "pos.tsv", header, rows), self.reference, ancestry="EUR")
        self.assertEqual(result.returncode, 0, result.stderr)
        by_pos = {r["POS"]: r for r in result.rows}
        self.assertEqual(by_pos[1000]["ID"], "rs5000")
        self.assertEqual((by_pos[wrong[1]]["ID"], by_pos[wrong[1]]["INFO"]["SNPID"]), (".", "rs5003"))
        self.assertEqual(result.qc["rsid_dbsnp"], {"matched": 29, "unmatched": 1})

    def test_build_is_detected_from_the_paired_sample(self):
        header = ["CHR", "POS", "EA", "OA", "EAF", "BETA", "SE", "P"]
        grch37 = [[1, v[1] - OFFSET, v[4], v[3], 0.3, 0.1, 0.02, 1e-6] for v in self.variants]
        result = run_standardise(self.tmp, write_table(self.tmp / "b37.tsv", header, grch37), self.reference, ancestry="EUR",
                                 name="b37")
        self.assertEqual(result.returncode, 0, result.stderr)
        detection = result.qc["build_detection"]
        self.assertEqual((detection["build"], detection["method"], detection["hits_hg19"], detection["hits_hg38"]), ("hg19", "detected", 30, 0))
        self.assertEqual(sorted(r["POS"] for r in result.rows)[:2], [1000, 1033])
        grch38 = [[1, v[1], v[4], v[3], 0.3, 0.1, 0.02, 1e-6] for v in self.variants]
        result = run_standardise(self.tmp, write_table(self.tmp / "b38.tsv", header, grch38), self.reference, ancestry="EUR",
                                 name="b38")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((result.qc["build_detection"]["build"], result.qc["build_detection"]["hits_hg38"]), ("hg38", 30))

    def test_allele_hash_is_stable(self):
        # The hash is written into the dbSNP index; a change would silently break every rsID match.
        self.assertEqual(gwas_reference.allele_hash(["A", "AT"], ["G", "A"]).tolist(), gwas_reference.allele_hash(["a", "AT"], ["g", "A"]).tolist())
        self.assertNotEqual(gwas_reference.allele_hash(["A"], ["G"])[0], gwas_reference.allele_hash(["G"], ["A"])[0])


if __name__ == "__main__":
    unittest.main()
