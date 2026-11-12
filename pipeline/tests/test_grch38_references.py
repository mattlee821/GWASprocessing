"""End-to-end tests with the GRCh38 references (dbSNP index, 1000 Genomes GRCh38 frequencies, paired build sample).

The fixture genome is seeded (821). testlib's offset chains place GRCh38 at GRCh37 + 100 on chromosomes 1 and 2.
"""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import chain_block, make_reference, needs_liftover, run_standardise, write_chain, write_fasta, write_grch38_references, write_table  # noqa: I001

import gwas_reference
from standardise_gwas import harmonisation_action

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
    write_fasta(reference / "genome", {"1": genome[OFFSET:], "2": genome[OFFSET:]}, name="GRCh37.fa")  # GRCh38 = GRCh37 + 100
    write_grch38_references(reference, variants, paired=[(v[0], v[1] - OFFSET, v[1], v[3], v[4]) for v in variants])
    return reference, variants


@needs_liftover
class BcftoolsLiftoverTests(unittest.TestCase):
    def test_grch37_inputs_are_lifted_with_the_effect_allele_kept_through_a_reference_change(self):
        tmp = Path(tempfile.mkdtemp())
        reference, variants = fixture(tmp)
        genome38 = (reference / "genome" / "GRCh38.fa").read_text().split("\n", 1)[1].replace("\n", "")
        swap = variants[5]  # at this site GRCh37 carries the ALT base: the reference allele changed between builds
        genome37 = list(genome38[OFFSET:])
        genome37[swap[1] - OFFSET - 1] = swap[4]
        write_fasta(reference / "genome", {"1": "".join(genome37), "2": genome38[OFFSET:]}, name="GRCh37.fa")
        header = ["CHR", "POS", "EA", "OA", "EAF", "BETA", "SE", "P"]
        rows = [[1, v[1] - OFFSET, v[4], v[3], 0.3, 0.1, 0.02, 1e-6] for v in variants]
        rows[5] = [1, swap[1] - OFFSET, swap[3], swap[4], 0.3, 0.1, 0.02, 1e-6]  # effect = the GRCh38 base, other = the GRCh37 base
        result = run_standardise(tmp, write_table(tmp / "b37.tsv", header, rows), reference, ancestry="EUR")
        self.assertEqual(result.returncode, 0, result.stderr)
        lift = result.qc["liftover"]["hg19_to_hg38"]
        self.assertEqual((lift["tool"], lift["ok"], lift["swapped_between_builds"]), ("bcftools +liftover", 30, 1))
        moved = next(r for r in result.rows if r["POS"] == swap[1])
        # the effect allele is the GRCh38 reference base, so it is written as REF and ES is negated relative to ALT
        self.assertEqual((moved["REF"], moved["ALT"], moved["FORMAT"]["ES"]), (swap[3], swap[4], "-0.1"))
        self.assertTrue(moved["INFO"]["HM"].startswith("swapped"))


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
        self.assertEqual((first["FORMAT"]["AF"], first["FORMAT"]["AFS"], first["FORMAT"]["ID"]), ("0.3", "1000G_EUR", "rs5000"))
        self.assertEqual(by_id["rs5001"]["FORMAT"]["AFS"], "study")
        self.assertEqual(result.qc["coordinate_source"], "dbsnp")
        self.assertEqual(result.qc["reference_fills"]["POS38"], 30)
        self.assertNotIn("hg19_to_hg38", result.qc.get("liftover", {}))
        self.assertEqual(result.qc["incomplete_by_cause"]["rsid_not_in_coordinate_index"]["rows"], 1)
        self.assertEqual(result.sample["AFFilled"], "1")

    def test_rsids_merged_in_dbsnp_take_the_position_of_the_current_rsid(self):
        write_grch38_references(self.reference, self.variants, merges={9000: 5002, 9001: 5003, 9002: 999998, 9003: 5004})
        rows = [[f"rs{v[2]}", v[4], v[3], 0.3, 0.1, 0.02, 1e-6, 5000] for v in self.variants]
        rows[2][0], rows[3][0] = "rs9000", "rs9001"  # the study used retired rsIDs for 2 variants
        rows.append(["rs9002", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000])  # merged into an rsID that is not in the index
        v4 = self.variants[4]
        rows.insert(0, ["rs9003", v4[4], v4[3], 0.3, 0.5, 0.02, 1e-6, 4000])  # also reported under its current rsID, later
        result = run_standardise(self.tmp, write_table(self.tmp / "merged.tsv", HEADER, rows), self.reference, ancestry="EUR")
        self.assertEqual(result.returncode, 0, result.stderr)
        by_id = result.by_id()
        self.assertEqual(len(by_id), 30)
        moved = by_id["rs5002"]
        self.assertEqual((moved["CHROM"], moved["POS"], moved["INFO"]["SNPID"], moved["FORMAT"]["ID"]),
                         ("1", self.variants[2][1], "rs9000", "rs5002"))
        self.assertEqual(by_id["rs5003"]["INFO"]["SNPID"], "rs9001")
        # the row under the current rsID is kept over the earlier one under the retired rsID
        self.assertEqual((by_id["rs5004"]["INFO"]["SNPID"], abs(float(by_id["rs5004"]["FORMAT"]["ES"]))), ("rs5004", 0.1))
        self.assertEqual(result.qc["duplicates_removed"], 1)
        self.assertEqual(result.qc["reference_fills"]["POS38"], 31)
        self.assertEqual(result.qc["reference_fills"]["POS38_merged_rsid"], 3)
        self.assertEqual(result.qc["incomplete_by_cause"]["rsid_not_in_coordinate_index"]["rows"], 1)

    @needs_liftover
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
        self.assertEqual((result.qc["rsid_dbsnp"]["matched"], result.qc["rsid_dbsnp"]["unmatched"]), (29, 1))

    @needs_liftover
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


COMPLEMENT = {"A": "T", "C": "G", "G": "C", "T": "A"}


class HarmonisationActionTests(unittest.TestCase):
    def test_hm_is_relative_to_the_alleles_the_study_reported(self):
        pal = np.array(["", "", "", "", "", "", "", "forward", "forward", "complemented", "ambiguous"], dtype=object)
        action = np.array(["keep", "complement", "keep", "swap", "complement", "drop", "indel_unmatched", "keep", "swap", "keep",
                           "keep"], dtype=object)
        # strand_changed: liftover already complemented the alleles (orientation to GRCh37 or a minus-strand chain block)
        changed = np.array([False, False, True, True, True, False, False, False, True, True, True])
        self.assertEqual(harmonisation_action(pal, action, changed).tolist(), [
            "forward", "complemented", "complemented", "complemented_swapped", "forward", "dropped", "unverified",
            "palindrome_forward", "palindrome_complemented_swapped", "palindrome_forward", "palindrome_ambiguous_complemented"])
        # an ambiguous palindrome that liftover did not complement stays ambiguous
        self.assertEqual(harmonisation_action(np.array(["ambiguous"], dtype=object), np.array(["swap"], dtype=object),
                                              np.array([False])).tolist(), ["palindrome_ambiguous_swapped"])

    def test_without_liftover_hm_follows_the_alignment_alone(self):
        self.assertEqual(harmonisation_action(np.array(["", ""], dtype=object), np.array(["complement_swap", "swap"], dtype=object)).tolist(),
                         ["complemented_swapped", "swapped"])


@needs_liftover
class LiftoverStrandTests(unittest.TestCase):
    """INFO/HM records a strand change made during liftover: GRCh37 chr1 [1000, 2000) maps to the minus strand of GRCh38
    chr1, so GRCh37 position 1501 is GRCh38 position 2000 and its alleles are complemented."""

    SITE37, SITE38 = 1501, 2000

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        rng = np.random.default_rng(821)
        self.genome37 = list("".join(rng.choice(list("ACGT"), size=3000)))
        genome38 = list("".join(rng.choice(list("ACGT"), size=3000)))
        for o in range(1000):
            genome38[3000 - 1 - (500 + o)] = COMPLEMENT[self.genome37[1000 + o]]
        self.genome38 = genome38
        self.reference = make_reference(self.tmp / "reference", lookup_rows=[], offset_chains=False)
        write_chain(self.reference / "genome" / "hg19ToHg38.over.chain.gz",
                    [chain_block("chr1", "chr1", 3000, 3000, 1000, 1000, 500, 1, dst_strand="-"),
                     chain_block("chr2", "chr2", 3000, 3100, 0, 2900, 100, 2)])
        self.genome37_2 = "".join(rng.choice(list("ACGT"), size=3000))
        self.genome38_2 = "".join(rng.choice(list("ACGT"), size=100)) + self.genome37_2

    def write_genomes(self):
        write_fasta(self.reference / "genome", {"1": "".join(self.genome37), "2": self.genome37_2}, name="GRCh37.fa")
        write_fasta(self.reference / "genome", {"1": "".join(self.genome38), "2": self.genome38_2}, name="GRCh38.fa")

    def run_rows(self, rows, name):
        header = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P"]
        result = run_standardise(self.tmp, write_table(self.tmp / f"{name}.tsv", header, rows), self.reference, input_build="hg19",
                                 ancestry="EUR", name=name)
        self.assertEqual(result.returncode, 0, result.stderr)
        return result.by_snpid()

    def test_a_minus_strand_block_complements_and_hm_says_so(self):
        self.write_genomes()
        ref = self.genome37[self.SITE37 - 1]
        alt = next(b for b in "ACGT" if b not in (ref, COMPLEMENT[ref]))
        row = self.run_rows([[1, self.SITE37, "rs1", alt, ref, 0.3, 0.1, 0.02, 1e-6]], "minus")["rs1"]
        self.assertEqual((row["POS"], row["REF"], row["ALT"], row["FORMAT"]["ES"]), (self.SITE38, COMPLEMENT[ref], COMPLEMENT[alt], "0.1"))
        self.assertEqual(row["INFO"]["HM"], "complemented")

    def test_a_reverse_strand_report_lifted_on_a_plus_block_reads_complemented(self):
        self.write_genomes()
        ref = self.genome37_2[999]
        alt = next(b for b in "ACGT" if b not in (ref, COMPLEMENT[ref]))
        row = self.run_rows([[2, 1000, "rs2", COMPLEMENT[alt], COMPLEMENT[ref], 0.3, 0.1, 0.02, 1e-6]], "plus")["rs2"]
        self.assertEqual((row["POS"], row["REF"], row["ALT"], row["INFO"]["HM"]), (1100, ref, alt, "complemented"))

    def test_a_reverse_strand_report_on_a_minus_block_cancels_to_forward(self):
        self.write_genomes()
        ref = self.genome37[self.SITE37 - 1]
        alt = next(b for b in "ACGT" if b not in (ref, COMPLEMENT[ref]))
        row = self.run_rows([[1, self.SITE37, "rs3", COMPLEMENT[alt], COMPLEMENT[ref], 0.3, 0.1, 0.02, 1e-6]], "both")["rs3"]
        self.assertEqual((row["REF"], row["ALT"], row["INFO"]["HM"]), (COMPLEMENT[ref], COMPLEMENT[alt], "forward"))

    def test_a_palindrome_resolved_after_a_minus_strand_liftover_reads_complemented(self):
        self.genome37[self.SITE37 - 1] = "A"
        self.genome38[self.SITE38 - 1] = "T"
        self.write_genomes()
        # 1000 Genomes on GRCh38: REF T, ALT A at frequency 0.8, so T (the complement of the study's A) has 0.2
        write_grch38_references(self.reference, [("1", self.SITE38, 0, "T", "A", 0.8)])
        row = self.run_rows([[1, self.SITE37, "rs4", "A", "T", 0.2, 0.1, 0.02, 1e-6]], "palindrome")["rs4"]
        self.assertEqual((row["POS"], row["REF"], row["ALT"]), (self.SITE38, "T", "A"))
        self.assertEqual(row["INFO"]["HM"], "palindrome_complemented_swapped")
        self.assertEqual(row["FORMAT"]["ES"], "-0.1")  # the study's A is GRCh38 T, the reference base


if __name__ == "__main__":
    unittest.main()
