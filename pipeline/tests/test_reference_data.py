"""Tests for src/reference_data: the reference validator, the dbSNP merge map builder and the build-detection
sample builder."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

from testlib import (SRC, bcftools_liftover_available, write_chain, write_genomes, write_grch38_references,
                     write_offset_chains)

import gwas_alleles as alleles

POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")
VALIDATOR = SRC / "reference_data" / "validate_references.py"
DETECTION_BUILDER = SRC / "reference_data" / "build_build_detection_sample.py"


def run(script, *args):
    return subprocess.run([sys.executable, str(script), *map(str, args)], capture_output=True, text=True, timeout=300)


class ValidatorTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())

    def complete(self):
        write_grch38_references(self.root, [("1", 1100, 1, "G", "A", 0.2)], populations=POPULATIONS,
                                paired=[("1", 1000, 1100, "G", "A")])
        write_genomes(self.root / "genome", [("1", 1000, "G")])
        for name in ("GRCh38.fa", "GRCh37.fa"):
            alleles.Fasta(self.root / "genome" / name)  # writes the .fai
        write_chain(self.root / "genome" / "hg19ToHg38.over.chain.gz", [])

    def test_empty_reference_root_fails_listing_what_is_missing(self):
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 1)
        for item in ("GRCh38.fa", "GRCh37.fa", "hg19ToHg38", "dbSNP GRCh38 index", "1000 Genomes GRCh38 EUR", "build-detection"):
            self.assertIn(item, proc.stderr)

    def test_a_missing_population_fails(self):
        self.complete()
        (self.root / "1000genomes" / "GRCh38" / "AFR.af.f4").unlink()
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("AFR frequencies (af.f4)", proc.stderr)

    def test_complete_references_pass_and_optional_items_are_reported(self):
        self.complete()
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertIn("merge map missing", proc.stdout)
        self.assertIn("hg38ToHg19", proc.stdout)
        write_grch38_references(self.root, [("1", 1100, 1, "G", "A", 0.2)], populations=POPULATIONS,
                                paired=[("1", 1000, 1100, "G", "A")], merges={7: 1})
        (self.root / "dbsnp" / "GRCh38" / "merged_meta.json").write_text(json.dumps({"retired_rsids": 1}), encoding="utf-8")
        write_chain(self.root / "genome" / "hg38ToHg19.over.chain.gz", [])
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 0)
        self.assertNotIn("WARNING", proc.stdout)
        self.assertIn("1 retired rsIDs", proc.stdout)


@unittest.skipUnless(bcftools_liftover_available(), "bcftools with the +liftover plugin is required (GWAS_BCFTOOLS, BCFTOOLS_PLUGINS)")
class BuildDetectionSampleTests(unittest.TestCase):
    def test_common_grch38_snps_are_placed_in_grch37_by_liftover(self):
        root = Path(tempfile.mkdtemp())
        rng = np.random.default_rng(821)
        genome37 = "".join(rng.choice(list("ACGT"), size=6000))
        sites = [(1000 + 50 * i, genome37[999 + 50 * i]) for i in range(40)]  # GRCh37 position, REF
        alt = {"A": "G", "C": "T", "G": "A", "T": "C"}
        variants = [("1", p + 100, 0, ref, alt[ref], 0.02 if i < 5 else 0.3) for i, (p, ref) in enumerate(sites)]  # 5 rare
        write_grch38_references(root, variants, populations=("EUR",))
        write_offset_chains(root / "genome")
        write_genomes(root / "genome", [("1", p, ref) for p, ref in sites])
        proc = run(DETECTION_BUILDER, "--reference-root", root, "--size", "30", "--force")  # the fixture wrote an empty sample
        self.assertEqual(proc.returncode, 0, proc.stderr)
        meta = json.loads((root / "build_detection" / "meta.json").read_text(encoding="utf-8"))
        self.assertEqual((meta["snps"], meta["common_candidates"]), (30, 35))
        key37 = np.fromfile(root / "build_detection" / "grch37.key.i64", dtype="<i8")
        key38 = np.fromfile(root / "build_detection" / "grch38.key.i64", dtype="<i8")
        np.testing.assert_array_equal(key38 - key37, 100)
        np.testing.assert_array_equal(np.fromfile(root / "build_detection" / "grch37.ref.u1", dtype="u1"),
                                      np.fromfile(root / "build_detection" / "grch38.ref.u1", dtype="u1"))


class DbsnpMergeBuilderTests(unittest.TestCase):
    def test_retired_rsids_map_to_the_end_of_their_merge_chain(self):
        import bz2

        root = Path(tempfile.mkdtemp())
        records = [
            {"refsnp_id": "10", "dbsnp1_merges": [{"merged_rsid": "5"}], "merged_snapshot_data": {"merged_into": ["20"]}},
            {"refsnp_id": "20", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["30"]}},  # 10 -> 20 -> 30
            {"refsnp_id": "40", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["41", "42"]}},  # ambiguous
            {"refsnp_id": "50", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["50"]}},  # self-merge
            # as dbSNP has rs33969899: its own record says 60 -> 61; 61's record lists 60 as merged earlier and 61 -> 62
            {"refsnp_id": "61", "dbsnp1_merges": [{"merged_rsid": "60"}], "merged_snapshot_data": {"merged_into": ["62"]}},
            {"refsnp_id": "60", "dbsnp1_merges": [], "merged_snapshot_data": {"merged_into": ["61"]}},
            # 70 appears only as an earlier merge of 2 records whose chains end differently: dropped
            {"refsnp_id": "71", "dbsnp1_merges": [{"merged_rsid": "70"}], "merged_snapshot_data": {"merged_into": ["80"]}},
            {"refsnp_id": "72", "dbsnp1_merges": [{"merged_rsid": "70"}], "merged_snapshot_data": {"merged_into": ["81"]}},
        ]
        source = root / "refsnp-merged.json.bz2"
        source.write_bytes(bz2.compress("".join(json.dumps(r, separators=(",", ":")) + "\n" for r in records).encode()))
        out = root / "dbsnp" / "GRCh38"
        result = run(SRC / "reference_data" / "build_dbsnp_merges.py", "--merged", source, "--out", out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(np.fromfile(out / "merged_from.u32", dtype="<u4").tolist(), [5, 10, 20, 60, 61, 71, 72])
        self.assertEqual(np.fromfile(out / "merged_into.u32", dtype="<u4").tolist(), [30, 30, 30, 62, 62, 80, 81])
        meta = json.loads((out / "merged_meta.json").read_text())
        self.assertEqual((meta["retired_rsids"], meta["retired_with_more_than_one_target"]), (7, 1))
        self.assertEqual((meta["retired_with_conflicting_targets"], meta["earlier_merge_disagreeing_with_own_record"]), (1, 0))


if __name__ == "__main__":
    unittest.main()
