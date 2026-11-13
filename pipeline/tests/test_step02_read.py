"""Step 2, Read (README "Pipeline in detail", 2. Read).

The delimiter comes from the manifest or YAML, or is detected from the header (tab, comma or whitespace); `.` reads as
missing and leading `##` lines are skipped. Columns are matched against built-in aliases, and a YAML `columns:` block
names any the aliases miss. Allele-dependent columns come from 1 family only (original or GWAS Catalog `hm_*`).
Compound identifiers keep their text in INFO/SNPID and supply positions when a file has none. GWAS-VCF records are read
into the same columns, a multi-allelic record becoming 1 row per ALT allele.

Positions are given in GRCh38 with input_build hg38 and an empty reference, so no liftover is involved.
"""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import pandas as pd

from testlib import EXAMPLES, example_reference, make_reference, run_standardise, write_table

import gwas_pipeline_common as common
from standardise_gwas import select_columns

HEADER = ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"]
GRCH38_CONTIG = "##contig=<ID=1,length=248956422,assembly=GRCh38>"
VCF_FORMATS = [f'##FORMAT=<ID={key},Number=A,Type=Float,Description="{key}">' for key in ("ES", "SE", "LP", "AF", "SS")]


def write_vcf(path: Path, records: list[tuple], fmt: str = "ES:SE:LP:AF:SS") -> Path:
    lines = ["##fileformat=VCFv4.2", GRCH38_CONTIG, *VCF_FORMATS, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGWAS"]
    lines += [f"{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t.\t{fmt}\t{v}" for c, p, i, r, a, v in records]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


class Step02ReadTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="step02_"))
        self.reference = make_reference(self.tmp / "reference", lookup_rows=[])

    def run_hg38(self, raw: Path, name: str, **kwargs):
        kwargs.setdefault("input_build", "hg38")
        return run_standardise(self.tmp, raw, kwargs.pop("reference", self.reference), name=name, **kwargs)

    def write_text(self, name: str, text: str) -> Path:
        path = self.tmp / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_read(self):
        # --- Delimiter ------------------------------------------------------------------------------------------
        with self.subTest("The delimiter is detected from the header: tab, comma or whitespace"):
            self.assertEqual(common.infer_delimiter(self.write_text("t.txt", "CHR\tPOS\tSNP\n1\t100\trs1\n")), "\t")
            self.assertEqual(common.infer_delimiter(self.write_text("crlf.txt", "CHR\tPOS\r\n1\t100\r\n")), "\t")
            self.assertEqual(common.infer_delimiter(self.write_text("c.txt", "CHR,POS,SNP\n1,100,rs1\n")), ",")
            self.assertEqual(common.infer_delimiter(self.write_text("s.txt", "SNP  A1 A2   b\nrs1  A  G    0.1\n")), r"\s+")

        tab = self.write_text("tab.tsv", "## produced by tool x, version 2\n## another note\n"
                              "CHR\tPOS\tSNP\tEA\tOA\tEAF\tBETA\tSE\tP\tN\n"
                              "1\t1100\trs1\tA\tG\t.\t0.10\t0.02\t1e-6\t5000\n"
                              "1\t1200\trs2\tC\tT\t\t0.05\t0.01\t0.2\t5000\n")
        tab_result = self.run_hg38(tab, "tab")
        with self.subTest("A tab-delimited table is read with the detected delimiter"):
            self.assertEqual(tab_result.returncode, 0, tab_result.stderr)
            self.assertEqual(tab_result.qc["input_delimiter"], "\t")
            self.assertEqual(sorted(tab_result.by_snpid()), ["rs1", "rs2"])

        with self.subTest("Leading ## lines are skipped"):
            table = common.read_table(tab)
            self.assertEqual(list(table.columns), HEADER)
            self.assertEqual(table.loc[0, "P"], 1e-6)

        with self.subTest("`.` reads as missing, and an empty cell does not shift later columns"):
            rows = tab_result.by_snpid()
            self.assertEqual(rows["rs1"]["FORMAT"], {"ES": "0.1", "SE": "0.02", "LP": "6", "AF": ".", "AFS": ".", "SS": "5000", "ID": "."})
            self.assertEqual((rows["rs2"]["FORMAT"]["ES"], rows["rs2"]["FORMAT"]["AF"], rows["rs2"]["FORMAT"]["SS"]), ("0.05", ".", "5000"))

        comma = write_table(self.tmp / "comma.csv", HEADER, [[1, 1100, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]], sep=",")
        with self.subTest("A comma-delimited table is read with the detected delimiter"):
            result = self.run_hg38(comma, "comma")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.qc["input_delimiter"], ",")
            self.assertEqual((result.rows[0]["FORMAT"]["ES"], result.rows[0]["FORMAT"]["AF"]), ("0.1", "0.3"))

        with self.subTest("A whitespace-aligned table is read with the detected delimiter"):
            result = run_standardise(self.tmp, EXAMPLES / "tables/locke_direct.tsv", example_reference(), name="whitespace")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.qc["input_delimiter"], r"\s+")
            self.assertEqual({k: result.qc["source_columns"][k] for k in ("EA", "OA", "EAF", "BETA", "SE", "P")},
                             {"EA": "A1", "OA": "A2", "EAF": "Freq1.Hapmap", "BETA": "b", "SE": "se", "P": "p"})
            self.assertEqual(len(result.rows), 2)

        with self.subTest("The delimiter comes from the manifest"):
            piped = write_table(self.tmp / "piped.txt", HEADER, [[1, 1100, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]], sep="|")
            result = self.run_hg38(piped, "manifest_delimiter", delimiter="|")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((result.qc["input_delimiter"], result.rows[0]["FORMAT"]["ES"]), ("|", "0.1"))

        with self.subTest("The delimiter comes from the study YAML"):
            semicolon = write_table(self.tmp / "semicolon.txt", HEADER, [[1, 1100, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]], sep=";")
            yaml_path = self.write_text("semicolon.yaml", 'delimiter: ";"\n')
            result = self.run_hg38(semicolon, "yaml_delimiter", study_yaml=str(yaml_path))
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((result.qc["input_delimiter"], result.rows[0]["FORMAT"]["ES"]), (";", "0.1"))

        # --- Column aliases -------------------------------------------------------------------------------------
        with self.subTest("Columns are matched against built-in aliases (EA, effect_allele, ALLELE1, A1; NCASE; INFO, Rsq)"):
            def chosen(columns: list[str]) -> dict:
                return select_columns(pd.DataFrame({c: ["1"] for c in columns}), {})

            for alias in ("EA", "effect_allele", "ALLELE1", "A1"):
                self.assertEqual(chosen([alias, "OA", "BETA"])["EA"], alias)
            for alias in ("NCASE", "N_CASES", "Ncases"):
                self.assertEqual(chosen(["EA", "OA", "BETA", alias])["NCASE"], alias)
            for alias in ("INFO", "Rsq"):
                self.assertEqual(chosen(["EA", "OA", "BETA", alias])["IMPINFO"], alias)
            self.assertEqual(chosen(["ea", "oa", "beta"])["EA"], "ea")  # case does not matter
            raw = write_table(self.tmp / "aliases.tsv", ["CHR", "POS", "SNP", "ALLELE1", "ALLELE0", "Freq_Tested_Allele_in_HRS",
                                                         "BETA", "SE", "P", "NCASE", "Rsq", "N"],
                              [[1, 1100, "rs1", "A", "G", 0.31, 0.1, 0.02, 1e-6, 500, 0.9, 1000]])
            result = self.run_hg38(raw, "aliases")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual({k: result.qc["source_columns"][k] for k in ("EA", "OA", "EAF", "NCASE", "IMPINFO")},
                             {"EA": "ALLELE1", "OA": "ALLELE0", "EAF": "Freq_Tested_Allele_in_HRS", "NCASE": "NCASE", "IMPINFO": "Rsq"})
            fmt = result.rows[0]["FORMAT"]
            self.assertEqual((fmt["AF"], fmt["AFS"], fmt["NC"], fmt["SI"]), ("0.31", "study", "500", "0.9"))

        with self.subTest("A YAML columns: block names columns the aliases miss"):
            raw = write_table(self.tmp / "estimate.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "Estimate", "SE", "P", "N"],
                              [[1, 1100, "rs1", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
            unmapped = self.run_hg38(raw, "estimate_unmapped")
            self.assertEqual(unmapped.returncode, 65)  # "Estimate" is no alias, so no row has an effect
            yaml_path = self.write_text("estimate.yaml", "columns:\n  BETA: Estimate\n")
            mapped = self.run_hg38(raw, "estimate_mapped", study_yaml=str(yaml_path))
            self.assertEqual(mapped.returncode, 0, mapped.stderr)
            self.assertEqual((mapped.qc["source_columns"]["BETA"], mapped.rows[0]["FORMAT"]["ES"]), ("Estimate", "0.1"))

        # --- Allele families ------------------------------------------------------------------------------------
        family_header = ["chromosome", "base_pair_location", "variant_id", "effect_allele", "other_allele", "beta",
                         "hm_effect_allele", "hm_other_allele", "hm_beta", "standard_error", "p_value"]
        with self.subTest("Allele-dependent columns come from 1 family: the hm_* family when the original is incomplete"):
            raw = write_table(self.tmp / "hm.tsv", family_header, [[1, 1100, "rs1", "A", "G", None, "G", "A", -0.1, 0.02, 1e-6]])
            result = self.run_hg38(raw, "family_hm")
            self.assertEqual(result.returncode, 0, result.stderr)
            row = result.rows[0]
            self.assertEqual((row["REF"], row["ALT"], row["FORMAT"]["ES"]), ("A", "G", "-0.1"))
            self.assertEqual({k: result.qc["source_columns"][k] for k in ("EA", "OA", "BETA")},
                             {"EA": "hm_effect_allele", "OA": "hm_other_allele", "BETA": "hm_beta"})

        with self.subTest("Allele-dependent columns come from 1 family: the original family when it is complete"):
            raw = write_table(self.tmp / "orig.tsv", family_header, [[1, 1100, "rs1", "A", "G", 0.2, "G", "A", -0.1, 0.02, 1e-6]])
            result = self.run_hg38(raw, "family_original")
            self.assertEqual(result.returncode, 0, result.stderr)
            row = result.rows[0]
            self.assertEqual((row["REF"], row["ALT"], row["FORMAT"]["ES"]), ("G", "A", "0.2"))
            self.assertEqual({k: result.qc["source_columns"][k] for k in ("EA", "OA", "BETA")},
                             {"EA": "effect_allele", "OA": "other_allele", "BETA": "beta"})

        # --- Compound identifiers -------------------------------------------------------------------------------
        with self.subTest("Compound identifiers (rs123:A:G, 1:12345:A:G) keep their text in INFO/SNPID"):
            raw = write_table(self.tmp / "compound.tsv", HEADER, [[1, 1100, "rs123:A:G", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000],
                                                                  [1, 12345, "1:12345:A:G", "A", "G", 0.3, 0.1, 0.02, 1e-6, 5000]])
            result = self.run_hg38(raw, "compound")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(sorted(result.by_snpid()), ["1:12345:A:G", "rs123:A:G"])

        with self.subTest("Compound identifiers supply positions when a file has none"):
            raw = write_table(self.tmp / "ids_only.tsv", ["MarkerName", "EA", "OA", "BETA", "SE", "P", "N"],
                              [["1:1100:A:G", "A", "G", 0.1, 0.02, 1e-6, 5000], ["chr1_1200_C_T", "C", "T", 0.1, 0.02, 1e-6, 5000]])
            result = self.run_hg38(raw, "ids_only")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual({r["INFO"]["SNPID"]: (r["CHROM"], r["POS"]) for r in result.rows},
                             {"1:1100:A:G": ("1", 1100), "chr1_1200_C_T": ("1", 1200)})
            self.assertEqual(result.qc["positions_from_variant_ids"], 2)

        # --- GWAS-VCF -------------------------------------------------------------------------------------------
        with self.subTest("GWAS-VCF records are read into the same columns, the format detected from the file"):
            raw = write_vcf(self.tmp / "download.txt", [("1", 1100, "rs1", "G", "A", "0.1:0.02:6:0.3:5000"),
                                                        ("1", 1200, "1:1200:C:T", "C", "T", "-0.2:0.04:2:0.4:4000")])
            result = run_standardise(self.tmp, raw, self.reference, name="vcf")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.qc["input_format"], "gwas_vcf")
            rows = result.by_snpid()
            self.assertEqual((rows["rs1"]["REF"], rows["rs1"]["ALT"]), ("G", "A"))
            self.assertEqual({k: rows["rs1"]["FORMAT"][k] for k in ("ES", "SE", "LP", "AF", "SS")},
                             {"ES": "0.1", "SE": "0.02", "LP": "6", "AF": "0.3", "SS": "5000"})
            self.assertEqual({k: rows["1:1200:C:T"]["FORMAT"][k] for k in ("ES", "SE", "LP", "AF", "SS")},
                             {"ES": "-0.2", "SE": "0.04", "LP": "2", "AF": "0.4", "SS": "4000"})

        with self.subTest("A multi-allelic GWAS-VCF record becomes 1 row per ALT allele"):
            raw = write_vcf(self.tmp / "multi.gwas.vcf", [("1", 1100, "rs1", "G", "A,T", "0.1,0.2:0.02,0.03:6,4:0.3,0.1:5000"),
                                                          ("1", 1200, "rs2", "C", "T", "0.1:0.02:6:0.3:5000")])
            result = run_standardise(self.tmp, raw, self.reference, name="multi")
            self.assertEqual(result.returncode, 0, result.stderr)
            split = sorted((r["POS"], r["REF"], r["ALT"], r["FORMAT"]["ES"], r["FORMAT"]["SS"]) for r in result.rows if r["POS"] == 1100)
            # SS has 1 value for the record, so it applies to both ALT alleles
            self.assertEqual(split, [(1100, "G", "A", "0.1", "5000"), (1100, "G", "T", "0.2", "5000")])

        # --- PLINK ----------------------------------------------------------------------------------------------
        with self.subTest("A PLINK table with REF, ALT and A1 derives the other allele from REF and ALT"):
            result = run_standardise(self.tmp, EXAMPLES / "tables/zhang_plink.tsv", example_reference(),
                                     config_root=EXAMPLES / "study_yaml", phenotype="protein", ancestry="EUR",
                                     author="Zhang", year="2022", PMID="35501419", study_id="zhang_2022_35501419", name="plink")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(result.qc["derived_oa_from_ref_alt"])
            rows = result.by_snpid()
            self.assertEqual((rows["rs5"]["POS"], rows["rs5"]["REF"], rows["rs5"]["ALT"]), (5100, "G", "A"))
            self.assertEqual((rows["rs6"]["REF"], rows["rs6"]["ALT"]), ("C", "T"))
            self.assertEqual(rows["rs5"]["INFO"]["ERRCODE"], ".")


if __name__ == "__main__":
    unittest.main()
