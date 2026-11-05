"""Tests for src/reference_data: lookup-table builder, rsID coordinate index builder and validator."""

import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from testlib import SRC, write_chain, write_fasta, write_lookup

import gwas_coordinates as coords

POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]
BUILDER = SRC / "reference_data" / "make_1000genomes_lookup_tables.py"
VALIDATOR = SRC / "reference_data" / "validate_references.py"


def run(script, *args):
    return subprocess.run([sys.executable, str(script), *map(str, args)], capture_output=True, text=True, timeout=300)


class LookupBuilderTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        self.processed = self.root / "1000genomes" / "phase3" / "processed"
        for pop in POPULATIONS:
            folder = self.processed / pop
            folder.mkdir(parents=True)
            (folder / f"{pop}.bim").write_text("1\trs1\t0\t1000\tA\tG\nX\trs2\t0\t500\tT\tC\n", encoding="utf-8")
            self.write_afreq(pop, 0.8)

    def write_afreq(self, pop, af):
        (self.processed / pop / "stats.afreq").write_text(
            f"#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n1\trs1\tG\tA\t{af}\t1000\nX\trs2\tC\tT\t0.1\t900\n", encoding="utf-8")

    def lookup(self, pop="EUR"):
        return pd.read_csv(self.root / "1000genomes" / "phase3" / "lookup" / f"{pop}.tsv.gz", sep="\t", dtype={"CHR": str})

    def test_lookup_records_alleles_and_the_frequency_of_a1(self):
        proc = run(BUILDER, "--reference-root", self.root, "--skip-coordinates")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        table = self.lookup()
        self.assertEqual(list(table.columns), ["CHR", "SNP", "POS19", "A1", "A2", "AF1"])
        self.assertEqual(table.iloc[0].tolist(), ["1", "rs1", 1000, "A", "G", 0.8])
        self.assertEqual(table.iloc[1].tolist(), ["X", "rs2", 500, "T", "C", 0.1])

    def test_existing_tables_are_kept_unless_forced(self):
        run(BUILDER, "--reference-root", self.root, "--skip-coordinates")
        self.write_afreq("EUR", 0.6)
        run(BUILDER, "--reference-root", self.root, "--skip-coordinates")
        self.assertEqual(self.lookup()["AF1"].iloc[0], 0.8)
        run(BUILDER, "--reference-root", self.root, "--skip-coordinates", "--force")
        self.assertEqual(self.lookup()["AF1"].iloc[0], 0.6)

    def pvar_text(self):
        return ("##fileformat=PVARv1.0\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\n"
                "1\t1000\trs1\tA\tG\nX\t500\trs99\tC\tT\n2\t300\t.\tA\tC\n3\t10\trs7\tA\tG\n4\t20\trs7\tA\tG\n")

    def check_index(self):
        index = coords.CoordinateIndex.load(self.root / "1000genomes" / "phase3" / "coordinates")
        found = index.lookup(pd.Series(["rs1", "rs99", "rs7"]))
        self.assertEqual(found["CHR"].tolist()[:2], ["1", "X"])
        self.assertEqual(found["POS19"].tolist()[:2], [1000, 500])
        self.assertTrue(pd.isna(found["CHR"].iloc[2]), "an rsID at two positions must not be used")

    def test_coordinate_index_is_built_from_the_decompressed_variant_list(self):
        (self.processed / "all_phase3.pvar").write_text(self.pvar_text(), encoding="utf-8")
        proc = run(BUILDER, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.check_index()

    @unittest.skipUnless(shutil.which("zstd"), "zstd is required")
    def test_coordinate_index_is_built_from_the_compressed_raw_variant_list(self):
        raw = self.root / "1000genomes" / "phase3" / "raw"
        raw.mkdir(parents=True)
        (raw / "all_phase3.pvar").write_text(self.pvar_text(), encoding="utf-8")
        subprocess.run(["zstd", "-q", "--rm", str(raw / "all_phase3.pvar")], check=True)
        proc = run(BUILDER, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.check_index()


class ValidatorTests(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())

    def complete_lookups(self, columns=None, rows=None):
        for pop in POPULATIONS:
            write_lookup(self.root, rows or [("1", "rs1", 1000, "A", "G", 0.2)], population=pop, columns=columns)

    def test_empty_reference_root_fails_listing_what_is_missing(self):
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("lookup", proc.stdout + proc.stderr)
        self.assertIn("hg19ToHg38", proc.stdout + proc.stderr)

    def test_lookups_without_directional_frequencies_fail(self):
        self.complete_lookups(columns=["CHR", "SNP", "SNPID", "POS19", "POS38", "EAF_ALL"], rows=[("1", "rs1", "1:1000", 1000, 1000, 0.2)])
        write_chain(self.root / "genome" / "hg19ToHg38.over.chain.gz", [])
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 1)
        self.assertIn("A1, A2, AF1", proc.stdout + proc.stderr)

    def test_complete_references_pass_and_optional_items_are_reported(self):
        self.complete_lookups()
        write_chain(self.root / "genome" / "hg19ToHg38.over.chain.gz", [])
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
        self.assertIn("coordinate index", proc.stdout)
        self.assertIn("FASTA", proc.stdout)
        write_fasta(self.root / "genome", {"1": "ACGT"})
        coords.write_coordinate_index(self.root / "1000genomes" / "phase3" / "coordinates", rs=[1], chrom=["1"], pos=[1000])
        write_chain(self.root / "genome" / "hg38ToHg19.over.chain.gz", [])
        proc = run(VALIDATOR, "--reference-root", self.root)
        self.assertEqual(proc.returncode, 0)
        self.assertNotIn("WARNING", proc.stdout)


if __name__ == "__main__":
    unittest.main()
