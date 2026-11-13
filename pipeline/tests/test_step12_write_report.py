"""README step 12, Write and report: the GWAS-VCF header and body, bgzip and tabix, the step log in standardised.qc.json
(also printed to the run log) whose counts reconcile from input rows to written rows, the OpenGWAS QC flags, the
self-contained HTML report with QQ plots, and the audit scripts that read an output in Python and in R (gwasvcf).

The fixture genome is seeded (821); only the bases named below matter.
"""

import contextlib
import csv
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from testlib import (PIPELINE, bcftools_liftover_available, make_reference, parse_vcf_rows, read_vcf_text,  # noqa: I001
                     run_script, run_standardise, write_table)

import gwas_vcf
import report
from standardise_gwas import FLAG_RULES, qc_flags, study_totals

PY_AUDIT = PIPELINE / "tests" / "gwasvcf_guide_audit.py"
R_AUDIT = PIPELINE / "tests" / "gwasvcf_guide_audit.R"
GRCH38_ORDER = [*map(str, range(1, 23)), "X", "Y", "MT"]
SAMPLE_COUNTS = ("TotalVariants", "VariantsNotRead", "VariantsNotHarmonised", "HarmonisedVariants", "SwitchedAlleles")


def sequence(length: int, fixed: dict[int, str], seed: int = 821) -> str:
    rng = np.random.default_rng(seed)
    bases = list("".join(rng.choice(list("ACGT"), size=length)))
    for pos, text in fixed.items():
        bases[pos - 1:pos - 1 + len(text)] = list(text)
    return "".join(bases)


def r_has_gwasvcf() -> bool:
    if not shutil.which("Rscript"):
        return False
    probe = subprocess.run(["Rscript", "-e", "quit(status = !requireNamespace('gwasvcf', quietly = TRUE))"], capture_output=True)
    return probe.returncode == 0


def read_flat(path: Path) -> list[dict[str, str]]:
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def without_htslib() -> dict[str, str]:
    path = os.pathsep.join(p for p in os.environ["PATH"].split(os.pathsep) if not Path(p, "bgzip").exists())
    return {**os.environ, "PATH": path}


class WriteAndReportStepTest(unittest.TestCase):
    def assert_reconciles(self, qc: dict) -> None:
        steps = qc["steps"]
        self.assertEqual(steps[0]["rows_in"], qc["input_rows"])
        self.assertEqual(steps[-1]["rows_out"], qc["output_rows"])
        for before, after in zip(steps, steps[1:]):
            self.assertEqual(after["rows_in"], before["rows_out"])
        for step in steps:
            self.assertEqual(step["rows_in"] - step["removed"], step["rows_out"])
            self.assertEqual(sum(step["removed_by"].values()), step["removed"])
        self.assertNotIn("does not reconcile", " ".join(qc["warnings"]))

    def assert_sample_counts_add_up(self, result) -> None:
        sample = {k: int(result.sample[k]) for k in SAMPLE_COUNTS}
        self.assertEqual(sample["TotalVariants"],
                         sample["VariantsNotRead"] + sample["VariantsNotHarmonised"] + sample["HarmonisedVariants"])
        self.assertEqual((sample["TotalVariants"], sample["HarmonisedVariants"]), (result.qc["input_rows"], result.qc["output_rows"]))
        alignment = result.qc["reference_alignment"]
        self.assertEqual(sample["VariantsNotHarmonised"], result.qc["duplicates_removed"] + alignment.get("dropped_mismatch", 0)
                         + alignment.get("indels_unmatched_dropped", 0))
        self.assertEqual(sample["SwitchedAlleles"], sum(r["INFO"]["HM"].endswith("swapped") for r in result.rows))

    def test_write_and_report(self):
        tmp = Path(tempfile.mkdtemp())

        # A table with every kind of removed row, rows out of genome order, an imputation column (FORMAT/SI) and an extra
        # column named in the study YAML (INFO/note).
        genome = sequence(3000, {100: "G", 200: "G", 300: "G", 400: "G", 500: "G"})
        reference = make_reference(tmp / "reference", lookup_rows=[], fasta={"1": genome})
        config = tmp / "manifests"
        config.mkdir()
        (config / "test_2026_1.yaml").write_text("extras:\n  - note\nsample_size: 12345\n", encoding="utf-8")
        header = ["CHR", "POS38", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N", "INFO", "note"]
        rows = [
            ["chr10", 500, "v1", "A", "G", 0.3, 0.01, 0.02, 0.3, 5000, 0.9, "a b"],
            [2, 900, "v2", "A", "G", 0.3, 0.02, 0.02, 0.01, 5000, 0.9, "x"],
            ["X", 50, "v3", "C", "T", 0.3, 0.03, 0.02, 0.5, 5000, 0.9, "x"],
            [1, 300, "v4", "A", "G", 0.3, 0.04, 0.02, 0.2, 5000, 0.9, "x"],
            ["MT", 70, "v5", "A", "G", 0.3, 0.05, 0.02, 0.7, 5000, 0.9, "x"],
            [2, 100, "v6", "A", "G", 0.3, 0.06, 0.02, 0.04, 5000, 0.9, "x"],
            [1, 200, "v7", "G", "A", 0.3, 0.07, 0.02, 0.6, 5000, 0.9, "x"],  # swapped against the FASTA
            [1, 100, "v8", "A", "G", 0.3, 0.08, None, None, 5000, 0.9, "x"],  # no standard error, and no P to derive one
            [1, 400, "v9", "D", "I", 0.3, 0.09, 0.02, 0.6, 5000, 0.9, "x"],  # alleles not made of A, C, G and T
            [1, 500, "v10", "A", "T", 0.3, 0.10, 0.02, 0.6, 5000, 0.9, "x"],  # matches neither strand of the FASTA
            ["", "", "rs99999999", "A", "G", 0.3, 0.11, 0.02, 0.6, 5000, 0.9, "x"],  # rsID absent from dbSNP
            [2, 900, "v2_again", "A", "G", 0.3, 0.12, 0.02, 0.6, 5000, 0.9, "x"],  # duplicate
        ]
        raw = write_table(tmp / "gwas.tsv", header, rows)
        full = run_standardise(tmp, raw, reference, config_root=config, name="full", sex="female",
                               output_id="GWAScatalog_GCST000001", phenotype="Body mass index")
        self.assertEqual(full.returncode, 0, full.stderr)
        vcf = Path(full.record["standard_vcf"])
        keys = ["ES", "SE", "LP", "AF", "AFS", "SS", "SI", "ID"]

        with self.subTest("the header declares GRCh38 contigs, the reference genome and the settings fingerprint"):
            contigs = [h for h in full.header if h.startswith("##contig=")]
            self.assertEqual([c.split("ID=")[1].split(",")[0] for c in contigs], GRCH38_ORDER)
            self.assertTrue(all(c.endswith(",assembly=GRCh38>") for c in contigs))
            self.assertIn("##contig=<ID=1,length=248956422,assembly=GRCh38>", contigs)
            self.assertIn("##reference=GRCh38.fa", full.header)
            self.assertIn(f"##settings_hash={full.qc['settings_hash']}", full.header)
            self.assertEqual(full.header[0], "##fileformat=VCFv4.2")

        with self.subTest("the sample is the output ID, with a ##SAMPLE line and a ##META definition for each key"):
            self.assertTrue(full.header[-1].endswith("\tFORMAT\tGWAScatalog_GCST000001"))
            self.assertEqual(sum(h.startswith("##SAMPLE=") for h in full.header), 1)
            sample = full.sample
            self.assertEqual((sample["ID"], sample["Sex"], sample["Phenotype"], sample["StudyID"], sample["Ancestry"], sample["PMID"]),
                             ("GWAScatalog_GCST000001", "female", "Body mass index", "test_2026_1", "ALL", "1"))
            self.assertEqual((sample["StudyType"], sample["TotalControls"], sample["AFFilled"]), ("Continuous", "12345", "0"))
            meta = {h.split("ID=")[1].split(",")[0] for h in full.header if h.startswith("##META=")}
            self.assertEqual(meta, set(sample) - {"ID"})
            # trait text is put on 1 line and quoted; an empty sex is written NA
            unit = gwas_vcf.table_header({"study_id": "s", "output_id": "local_test 1", "phenotype": "Body mass index\t(BMI)\nx",
                                          "ancestry": "EUR", "sex": "", "PMID": "1"}, [], gwas_vcf.format_keys([]),
                                         {"TotalVariants": 1, "Phenotype": "Body mass index\t(BMI)\nx", "Sex": ""})
            line = next(h for h in unit if h.startswith("##SAMPLE="))
            self.assertEqual(line, '##SAMPLE=<ID=local_test_1,TotalVariants=1,Phenotype="Body mass index (BMI) x",Sex=NA>')
            self.assertEqual(study_totals({"metadata": {}}, {"ncase": 100, "ncontrol": 200}, False), ("CaseControl", 100.0, 200.0))

        with self.subTest("every FORMAT and INFO field is declared; optional FORMAT keys follow the columns present"):
            formats = [h.split("ID=")[1].split(",")[0] for h in full.header if h.startswith("##FORMAT=")]
            self.assertEqual(formats, keys)
            self.assertTrue(all(list(r["FORMAT"]) == keys for r in full.rows))
            infos = [h.split("ID=")[1].split(",")[0] for h in full.header if h.startswith("##INFO=")]
            self.assertEqual(infos, ["SNPID", "HM", "note"])
            self.assertFalse(any("POS19" in h or "EAF_SOURCE" in h for h in full.header))
            self.assertEqual(gwas_vcf.format_keys([]), ["ES", "SE", "LP", "AF", "AFS", "ID"])
            self.assertEqual(gwas_vcf.format_keys(["N", "NCASE", "EZ", "IMPINFO"]),
                             ["ES", "SE", "LP", "AF", "AFS", "SS", "NC", "EZ", "SI", "ID"])

        with self.subTest("extra columns are written as INFO strings"):
            self.assertIn('##INFO=<ID=note,Number=1,Type=String,Description="Source column note">', full.header)
            self.assertEqual(full.by_snpid()["v1"]["INFO"]["note"], "a_b")
            extras = gwas_vcf.table_header({"output_id": "x"}, ["9score", "HM"], gwas_vcf.format_keys([]), {})
            self.assertIn('##INFO=<ID=X_9score,Number=1,Type=String,Description="Source column 9score">', extras)
            self.assertIn('##INFO=<ID=SRC_HM,Number=1,Type=String,Description="Source column HM">', extras)

        with self.subTest("record values: AF from 1000 Genomes to 6 significant digits, the study's as given, '.' when missing"):
            # 1000 Genomes frequencies arrive as float32, so 0.874498 reads back as 0.87449801
            body = pd.DataFrame({"CHR": ["1"] * 4, "POS38": [10, 20, 30, 40], "SNP": ["rs1", pd.NA, "1:30", "rs4"],
                                 "SNPID": ["1:10;x", pd.NA, "", "rs4"], "EA": ["A"] * 4, "OA": ["G"] * 4,
                                 "EAF": [float(np.float32(0.874498)), float(np.float32(0.0592369)), np.nan, 0.123456789],
                                 "EAF_SOURCE": ["reference", "reference", "", "input"], "BETA": [0.1] * 4, "SE": [0.02] * 4,
                                 "LP": [544.665305866, 6.0, np.nan, 2.0], "N": [5000.4, np.nan, 10, 10], "HM": ["forward"] * 4})
            lines = parse_vcf_rows("\n".join(gwas_vcf.table_body_lines(body, [], gwas_vcf.format_keys(["N"]), "1000G_EUR")))
            self.assertEqual([(r["FORMAT"]["AF"], r["FORMAT"]["AFS"]) for r in lines],
                             [("0.874498", "1000G_EUR"), ("0.0592369", "1000G_EUR"), (".", "."), ("0.12345679", "study")])
            self.assertEqual([r["ID"] for r in lines], ["rs1", ".", ".", "rs4"])
            self.assertEqual([r["INFO"]["SNPID"] for r in lines], ["1:10_x", ".", ".", "rs4"])
            self.assertEqual([r["FORMAT"]["LP"] for r in lines], ["544.66531", "6", ".", "2"])
            self.assertEqual([r["FORMAT"]["SS"] for r in lines], ["5000.4", ".", "10", "10"])
            self.assertEqual((lines[0]["REF"], lines[0]["ALT"]), ("G", "A"))

        with self.subTest("records are sorted in numeric chromosome order and by position, without a chr prefix"):
            self.assertEqual([(r["CHROM"], r["POS"]) for r in full.rows],
                             [("1", 200), ("1", 300), ("2", 100), ("2", 900), ("10", 500), ("X", 50), ("MT", 70)])
            self.assertTrue(all(r["line"].split("\t")[5:7] == [".", "PASS"] for r in full.rows))

        with self.subTest("the output is compressed with bgzip and indexed with tabix"):
            index = Path(str(vcf) + ".tbi")
            self.assertTrue(index.exists())
            self.assertEqual(subprocess.run(["bgzip", "-t", str(vcf)], capture_output=True).returncode, 0)
            listed = subprocess.run(["tabix", "-l", str(vcf)], capture_output=True, text=True, check=True).stdout.split()
            self.assertEqual(listed, ["1", "2", "10", "X", "MT"])
            region = subprocess.run(["tabix", str(vcf), "2:1-1000"], capture_output=True, text=True, check=True).stdout
            self.assertEqual([r["POS"] for r in parse_vcf_rows(region)], [100, 900])
            if shutil.which("bcftools"):
                check = subprocess.run(["bcftools", "view", "-Ou", "-o", os.devnull, str(vcf)], capture_output=True, text=True)
                self.assertEqual((check.returncode, check.stderr.strip()), (0, ""))

        with self.subTest("missing htslib is a configuration error"):
            path = os.environ["PATH"]
            try:
                os.environ["PATH"] = without_htslib()["PATH"]
                with self.assertRaises(SystemExit) as caught, contextlib.redirect_stderr(io.StringIO()):
                    gwas_vcf.require_htslib()
                self.assertEqual(caught.exception.code, 78)
            finally:
                os.environ["PATH"] = path
            missing = run_standardise(tmp, raw, reference, config_root=config, name="no_htslib", env=without_htslib())
            self.assertEqual(missing.returncode, 78)
            self.assertIn("bgzip", missing.stderr)

        with self.subTest("the step log reconciles for a table: every removed row counted once with its cause"):
            self.assert_reconciles(full.qc)
            self.assertEqual((full.qc["input_rows"], full.qc["output_rows"]), (12, 7))
            causes = full.qc["incomplete_by_cause"]
            self.assertEqual({k: v["rows"] for k, v in causes.items()},
                             {"rsid_not_in_coordinate_index": 1, "invalid_allele": 1, "missing_or_invalid_se": 1})
            self.assertEqual(sum(c["rows"] for c in causes.values()), full.qc["rows_removed_missing_required"])
            removed = {s["title"]: s["removed"] for s in full.qc["steps"] if s["removed"]}
            self.assertEqual(removed, {"Alleles and GRCh38 coordinates": 2, "Statistics": 1, "Duplicates": 1,
                                       "GRCh38 reference alignment": 1})

        with self.subTest("the step log is printed to the run log"):
            self.assertIn("Alleles and GRCh38 coordinates", full.stderr)
            self.assertIn("Written: 7 of 12 input rows", full.stderr)
            self.assertNotIn("UNRECONCILED", full.stderr)

        with self.subTest("the ##SAMPLE header counts add up and match the step log"):
            self.assert_sample_counts_add_up(full)
            self.assertEqual(full.sample["SwitchedAlleles"], "1")

        with self.subTest("the step log reconciles for a lifted GWAS-VCF"):
            if not bcftools_liftover_available():
                self.skipTest("bcftools with the +liftover plugin is required")
            lifted_reference = make_reference(tmp / "vcf_reference", genomes=True,
                                              lookup_rows=[("1", "rs1", 1000, "G", "A", 0.3), ("1", "rs2", 2000, "C", "T", 0.004)])
            records = [("1", 1000, "rs1", "A", "G"), ("1", 2000, "rs2", "T", "C"), ("1", 1000, "rs1", "A", "G"),
                       ("1", 2_000_000, "rs3", "A", "G"), ("GL000192.1", 10, "rs4", "A", "G"), ("2", 2_000_000, "rs5", "A", "G")]
            text = ["##fileformat=VCFv4.2", "##contig=<ID=1,length=249250621,assembly=HG19/GRCh37>",
                    *[f'##FORMAT=<ID={k},Number=A,Type=Float,Description="{k}">' for k in ("ES", "SE", "LP", "AF", "SS")],
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tBMI"]
            text += [f"{c}\t{p}\t{i}\t{r}\t{a}\t.\tPASS\t.\tES:SE:LP:AF:SS\t0.01:0.01:1:0.3:1000" for c, p, i, r, a in records]
            lifted_raw = tmp / "grch37.gwas.vcf"
            lifted_raw.write_text("\n".join(text) + "\n", encoding="utf-8")
            lifted = run_standardise(tmp, lifted_raw, lifted_reference, name="lifted", input_build="hg19")
            self.assertEqual(lifted.returncode, 0, lifted.stderr)
            self.assert_reconciles(lifted.qc)
            self.assertEqual((lifted.qc["input_rows"], lifted.qc["input_format"]), (6, "gwas_vcf"))
            steps = {s["title"]: s["removed"] for s in lifted.qc["steps"]}
            self.assertEqual((steps["Alleles and GRCh38 coordinates"], steps["Duplicates"]), (3, 1))
            self.assert_sample_counts_add_up(lifted)

        with self.subTest("the Python audit flattens an output"):
            flat = tmp / "flat.python.tsv"
            proc = subprocess.run([sys.executable, str(PY_AUDIT), str(vcf), str(flat)], capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            flat_rows = read_flat(flat)
            self.assertEqual([(r["CHR"], int(r["POS"]), r["SNP"], r["REF"], r["ALT"], r["BETA"], r["N"]) for r in flat_rows],
                             [(r["CHROM"], r["POS"], r["ID"], r["REF"], r["ALT"], r["FORMAT"]["ES"], r["FORMAT"]["SS"])
                              for r in full.rows])
            self.assertAlmostEqual(float(flat_rows[2]["P"]), 0.04, places=7)  # chromosome 2 position 100

        with self.subTest("the Python audit rejects records without the required FORMAT keys"):
            bad = tmp / "noss.vcf"
            bad.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
                           "1\t100\trs1\tA\tG\t.\tPASS\t.\tES:SE:LP\t0.1:0.01:2\n", encoding="utf-8")
            proc = subprocess.run([sys.executable, str(PY_AUDIT), str(bad), str(tmp / "x.tsv")], capture_output=True, text=True)
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("SS", proc.stderr)

        with self.subTest("R gwasvcf reads the output identically"):
            if not r_has_gwasvcf():
                self.skipTest("R with the gwasvcf package is required")
            r_flat = tmp / "flat.r.tsv"
            proc = subprocess.run(["Rscript", str(R_AUDIT), str(vcf), str(r_flat)], capture_output=True, text=True, timeout=600)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            r_rows, py_rows = read_flat(r_flat), read_flat(tmp / "flat.python.tsv")
            key = lambda r: (r["CHR"], r["POS"], r["REF"], r["ALT"])  # noqa: E731
            self.assertEqual([key(r) for r in r_rows], [key(r) for r in py_rows])
            for column in ("BETA", "SE", "EAF"):
                self.assertEqual([float(r[column]) for r in r_rows], [float(r[column]) for r in py_rows], column)
            np.testing.assert_allclose([float(r["P"]) for r in r_rows], [float(r["P"]) for r in py_rows], rtol=1e-6)

        # A small table without a FASTA: 3 P values give lambda GC 8.44394650653 (R: qchisq(1 - 0.05, 1, lower.tail = TRUE)
        # / qchisq(0.5, 1)), and N of 1,000 and a mean effect of 0.6 raise QC flags.
        scratch = tmp / "scratch"
        scratch.mkdir()
        small_raw = write_table(tmp / "small.tsv", ["CHR", "POS", "SNP", "EA", "OA", "EAF", "BETA", "SE", "P", "N"],
                                [[1, 1000 + i, f"rs{i}", "A", "G", 0.3, 0.6, 0.02, p, 1000] for i, p in enumerate([0.01, 0.05, 0.2])])
        small = run_standardise(tmp, small_raw, make_reference(tmp / "small_reference", lookup_rows=[]), name="small",
                                input_build="hg38", env={**os.environ, "GWAS_TMPDIR": str(scratch)})
        self.assertEqual(small.returncode, 0, small.stderr)

        with self.subTest("QC flags follow the OpenGWAS thresholds and never stop a row"):
            self.assertEqual(set(small.qc["qc_flags"]), {"inflation_factor", "n", "mean_EFFECT_05", "mean_EFFECT_01"})
            self.assertEqual(len(small.rows), 3)
            clean = {"af_concordance": {"pearson_r": 0.95}, "lambda_gc": 1.05,
                     "qc_metrics": {"max_sample_size": 50000, "mean_effect": 0.01, "n_p_below_5e-8": 1000}}
            self.assertEqual(qc_flags(clean, {}), {})
            raised = {"af_concordance": {"pearson_r": 0.69}, "lambda_gc": 1.21,
                      "qc_metrics": {"max_sample_size": 9999, "mean_effect": -0.51, "n_p_below_5e-8": 1001}}
            self.assertEqual(qc_flags(raised, {}), FLAG_RULES)
            between = {"qc_metrics": {"mean_effect": 0.2}, "lambda_gc": 1.2}
            self.assertEqual(set(qc_flags(between, {"TotalControls": 20000})), {"mean_EFFECT_01"})
            self.assertIn("n", qc_flags({"qc_metrics": {}}, {"TotalControls": 4000, "TotalCases": 1000}))

        with self.subTest("genomic inflation and provenance are recorded"):
            self.assertAlmostEqual(small.qc["lambda_gc"], 8.44394650653, places=5)
            self.assertEqual(small.qc["provenance"]["schema_version"], small.metadata["pipeline_version"])
            self.assertEqual(len(small.qc["settings_hash"]), 16)
            self.assertEqual(small.metadata["settings_hash"], small.qc["settings_hash"])

        with self.subTest("no intermediate files are left in the output or temporary folders"):
            outputs = sorted(p.name for p in vcf.parent.iterdir())
            self.assertEqual(outputs, ["standardised.gwas.vcf.gz", "standardised.gwas.vcf.gz.tbi", "standardised.metadata.json",
                                       "standardised.qc.json"])
            self.assertEqual(sorted(p.name for p in Path(small.record["standard_vcf"]).parent.iterdir()), outputs)
            self.assertEqual(list(scratch.iterdir()), [])

        record = tmp / "small" / "record.json"

        def plot(*extra):
            out = tmp / "plot.json"
            run_script("report.py", "--standardise-record", record, "--output", out, *extra)
            return json.loads(out.read_text(encoding="utf-8"))

        with self.subTest("the HTML report holds the steps, flags, header counts, QQ plots, warnings and timings in 1 file"):
            made = plot()
            self.assertTrue(made["report_created"])
            text = Path(made["report_html"]).read_text(encoding="utf-8")
            for section in ("<h2>Steps</h2>", "<h2>QC flags</h2>", "<h2>##SAMPLE header</h2>", "<h2>QQ plots</h2>",
                            "<h2>Warnings</h2>", "<h2>Time (seconds)</h2>"):
                self.assertIn(section, text)
            self.assertIn("Palindromic (A/T, C/G) SNPs", text)
            self.assertIn("TotalVariants", text)
            self.assertIn("lambda GC above 1.2", text)
            self.assertIn("data:image/png;base64,", text)
            self.assertNotIn("http://", text.replace("http://www.w3.org", ""))  # no external links or files
            self.assertNotIn("<script src", text)
            self.assertAlmostEqual(made["lambda_gc_raw"], 8.44394650653, places=5)  # the raw panel reads the raw table
            self.assertAlmostEqual(made["lambda_gc_standardised"], 8.44394650653, places=5)

        with self.subTest("QQ plots keep every point above LP 2 and thin the bulk"):
            rng = np.random.default_rng(821)
            lp = -np.log10(rng.uniform(size=100_000))
            expected, observed = report.qq_points(lp)
            self.assertEqual(int((observed >= 2).sum()), int((lp >= 2).sum()))
            self.assertLess(len(observed), 5_000)
            self.assertAlmostEqual(float(observed.max()), float(lp.max()))
            self.assertAlmostEqual(float(expected.max()), -np.log10(0.5 / 100_000))

        with self.subTest("an existing report is kept only while it is newer than the VCF"):
            self.assertFalse(plot()["report_created"])
            later = time.time() + 5
            os.utime(small.record["standard_vcf"], (later, later))
            self.assertTrue(plot()["report_created"])

        with self.subTest("a report can be made after processing from the output folder"):
            folder = Path(small.record["standard_vcf"]).parent
            (folder / "standardised.report.html").unlink()
            run_script("report.py", "--dir", folder)
            text = (folder / "standardised.report.html").read_text(encoding="utf-8")
            self.assertIn("<h2>Steps</h2>", text)
            self.assertIn("data:image/png;base64,", text)


if __name__ == "__main__":
    unittest.main()
