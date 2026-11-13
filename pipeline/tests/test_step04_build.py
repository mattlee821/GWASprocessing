"""Step 4, Genome build (README "Pipeline in detail", 4. Genome build).

The build is detected from the data in every case: a fixed-seed sample of up to 20,000 rows from the first 200,000 is
matched on chromosome, position and alleles (either order, either strand) against common SNPs placed in GRCh37 and
GRCh38. A build is called when it has at least 20 matches and at least 10 times the other build's matches. Only when
detection fails is a declared build used: the manifest or YAML input_build, then the source metadata, then the VCF
header. A declared build that disagrees with the detected one is reported and not used; with neither, the row stops
(exit 65: that row fails, the other rows run). Files with rsIDs and no positions need no build.

The fixture places 30 common SNPs at GRCh37 10,000 + 1,000 i and GRCh38 = GRCh37 + 100 (testlib's offset chains).
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import bcftools_liftover_available, make_reference, run_standardise, write_table

import gwas_coordinates as coords
import gwas_reference
import standardise_gwas
from standardise_gwas import declared_build, detect_build_paired

LOOKUP = [("1", f"rs{100 + i}", 10_000 + 1000 * i, "A", "G", 0.3) for i in range(30)]  # CHR, SNP, POS19, ALT, REF, AF
GRCH37 = [row[2] for row in LOOKUP]
GRCH38 = [pos + 100 for pos in GRCH37]
NO_MATCH = [500_001 + 7 * i for i in range(30)]  # in neither build of the sample
HEADER = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]


def write_vcf(path: Path, records: list[tuple], contig: str) -> Path:
    lines = ["##fileformat=VCFv4.2", *([contig] if contig else []),
             *[f'##FORMAT=<ID={k},Number=A,Type=Float,Description="{k}">' for k in ("ES", "SE", "LP", "AF", "SS")],
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGWAS"]
    lines += [f"1\t{pos}\t{rsid}\tG\tA\t.\tPASS\t.\tES:SE:LP:AF:SS\t0.1:0.02:6:0.3:5000" for pos, rsid in records]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class Step04BuildTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="step04_"))
        self.reference = make_reference(self.tmp / "reference", lookup_rows=LOOKUP, genomes=True)
        self.paired = gwas_reference.paired_build_sample(self.reference)

    def table(self, name: str, positions: list[int], ea: str = "A", oa: str = "G") -> Path:
        return write_table(self.tmp / f"{name}.tsv", HEADER, [[1, pos, f"1:{pos}", ea, oa, 0.3, 0.1, 0.02, 1e-6, 5000]
                                                             for pos in positions])

    def run_table(self, name: str, positions: list[int], **kwargs):
        return run_standardise(self.tmp, self.table(name, positions), self.reference, name=name, **kwargs)

    def detect(self, positions: list[int], ea: str = "A", oa: str = "G") -> dict:
        n = len(positions)
        return detect_build_paired(pd.DataFrame({"CHR": ["1"] * n, "POS": positions, "EA": [ea] * n, "OA": [oa] * n}), self.paired)

    def test_genome_build(self):
        with self.subTest("GRCh37 is detected from the data"):
            result = self.detect(GRCH37)
            self.assertEqual((result["build"], result["hits_hg19"], result["hits_hg38"]), ("hg19", 30, 0))
            if not bcftools_liftover_available():
                self.skipTest("bcftools with the +liftover plugin is needed to standardise a GRCh37 file end to end")
            lifted = self.run_table("grch37", GRCH37)
            self.assertEqual(lifted.returncode, 0, lifted.stderr)
            detection = lifted.qc["build_detection"]
            self.assertEqual((detection["build"], detection["method"], detection["hits_hg19"], detection["hits_hg38"]),
                             ("hg19", "detected", 30, 0))
            self.assertEqual(sorted(r["POS"] for r in lifted.rows), GRCH38)

        with self.subTest("GRCh38 is detected from the data"):
            result = self.run_table("grch38", GRCH38)
            self.assertEqual(result.returncode, 0, result.stderr)
            detection = result.qc["build_detection"]
            self.assertEqual((detection["build"], detection["method"], detection["hits_hg19"], detection["hits_hg38"]),
                             ("hg38", "detected", 0, 30))
            self.assertEqual(sorted(r["POS"] for r in result.rows), GRCH38)

        with self.subTest("Alleles match in either order and on either strand; other alleles do not count"):
            self.assertEqual(self.detect(GRCH37, ea="G", oa="A")["hits_hg19"], 30)
            self.assertEqual(self.detect(GRCH37, ea="T", oa="C")["hits_hg19"], 30)
            self.assertEqual((self.detect(GRCH37, ea="A", oa="C")["build"], self.detect(GRCH37, ea="A", oa="C")["hits_hg19"]), ("", 0))

        with self.subTest("A build is called with at least 20 matches"):
            enough = self.run_table("hits_20", GRCH38[:20] + NO_MATCH[:5])
            self.assertEqual(enough.returncode, 0, enough.stderr)
            self.assertEqual((enough.qc["build_detection"]["build"], enough.qc["build_detection"]["method"]), ("hg38", "detected"))
            too_few = self.run_table("hits_19", GRCH38[:19] + NO_MATCH[:5])
            self.assertEqual(too_few.returncode, 65)
            self.assertIn("input_build", too_few.stderr)

        with self.subTest("A build is called with at least 10 times the other build's matches"):
            self.assertEqual(self.detect(GRCH37[:20] + GRCH38[25:27])["build"], "hg19")  # 20 against 2
            self.assertEqual(self.detect(GRCH37[:20] + GRCH38[25:28])["build"], "")  # 20 against 3

        with self.subTest("A fixed-seed sample of up to 20,000 rows from the first 200,000 is matched"):
            self.assertEqual((coords.DETECTION_SAMPLE, standardise_gwas.BUILD_HEAD_ROWS), (20_000, 200_000))
            positions = GRCH38 + [600_001 + i for i in range(25_000)]
            first, second = self.detect(positions), self.detect(positions)
            self.assertEqual(first["sampled"], 20_000)
            self.assertEqual((first["hits_hg38"], first["hits_hg19"]), (second["hits_hg38"], second["hits_hg19"]))

        with self.subTest("A declared build is used only when detection fails"):
            result = self.run_table("declared", NO_MATCH, input_build="hg38")
            self.assertEqual(result.returncode, 0, result.stderr)
            detection = result.qc["build_detection"]
            self.assertEqual((detection["build"], detection["method"], detection["declared_source"]), ("hg38", "declared", "manifest_or_yaml"))
            self.assertEqual(sorted(r["POS"] for r in result.rows), NO_MATCH)

        with self.subTest("Declared builds are taken from the manifest or YAML, then the source metadata, then the VCF header"):
            metadata = {"metadata": {"gwasinfo": [{"build": "HG19/GRCh37"}]}}
            self.assertEqual(declared_build({"input_build": "hg38", **metadata}, {"_vcf_header_build": "hg19"}), ("hg38", "manifest_or_yaml"))
            self.assertEqual(declared_build({"input_build": "", **metadata}, {"input_build": "hg38"}), ("hg38", "manifest_or_yaml"))
            self.assertEqual(declared_build({"input_build": "", **metadata}, {"_vcf_header_build": "hg38"}), ("hg19", "source_metadata"))
            self.assertEqual(declared_build({"input_build": "", "metadata": {}}, {"_vcf_header_build": "hg38"}), ("hg38", "vcf_header"))
            from_metadata = self.run_table("source_metadata", NO_MATCH, metadata={"genome_assembly": "GRCh38"})
            self.assertEqual(from_metadata.returncode, 0, from_metadata.stderr)
            self.assertEqual((from_metadata.qc["build_detection"]["method"], from_metadata.qc["build_detection"]["declared_source"]),
                             ("declared", "source_metadata"))
            vcf = write_vcf(self.tmp / "header.gwas.vcf", [(123_457, "rsX")], contig="##contig=<ID=1,length=248956422,assembly=GRCh38>")
            from_header = run_standardise(self.tmp, vcf, self.reference, name="vcf_header")
            self.assertEqual(from_header.returncode, 0, from_header.stderr)
            self.assertEqual((from_header.qc["build_detection"]["build"], from_header.qc["build_detection"]["method"],
                              from_header.qc["build_detection"]["declared_source"]), ("hg38", "declared", "vcf_header"))

        with self.subTest("A declared build that disagrees with the detected one is reported and not used"):
            result = self.run_table("contradiction", GRCH38, input_build="hg19")
            self.assertEqual(result.returncode, 0, result.stderr)
            detection = result.qc["build_detection"]
            self.assertEqual((detection["build"], detection["method"], detection["declared"]), ("hg38", "detected", "hg19"))
            self.assertTrue(any("disagrees" in w for w in result.qc["warnings"]), result.qc["warnings"])
            self.assertEqual(sorted(r["POS"] for r in result.rows), GRCH38)  # not lifted as GRCh37

        with self.subTest("With neither a detected nor a declared build the row fails on its own (exit 65)"):
            table = self.run_table("undeterminable", NO_MATCH)
            self.assertEqual(table.returncode, 65)
            self.assertIn("input_build", table.stderr)
            vcf = write_vcf(self.tmp / "nobuild.gwas.vcf", [(123_457, "rsX")], contig="")
            gwas_vcf = run_standardise(self.tmp, vcf, self.reference, name="undeterminable_vcf")
            self.assertEqual(gwas_vcf.returncode, 65)
            self.assertIn("input_build", gwas_vcf.stderr)

        with self.subTest("Files with rsIDs and no positions need no build"):
            raw = write_table(self.tmp / "rsids.tsv", ["SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"],
                              [[snp, alt, ref, 0.3, 0.1, 0.02, 1e-6, 5000] for _, snp, _, alt, ref, _ in LOOKUP[:5]])
            result = run_standardise(self.tmp, raw, self.reference, name="rsids")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("build_detection", result.qc)
            self.assertEqual(result.qc["coordinate_source"], "dbsnp")
            self.assertEqual(sorted(r["POS"] for r in result.rows), GRCH38[:5])


if __name__ == "__main__":
    unittest.main()
