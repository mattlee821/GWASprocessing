"""The manifest: reading and validating it (prepare_manifest.py, gwas_pipeline_common), deriving IDs and output folders,
choosing which rows run, and the state file of completed rows (update_state.py).

1 test; each README rule is a subTest. Synthetic data only, no network.
"""

from __future__ import annotations

import contextlib
import csv
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from testlib import EXAMPLES, make_reference, run_script

import gwas_pipeline_common as common

REQUIRED = ["GWAS_location", "phenotype", "ancestry", "author", "year", "PMID"]
OPTIONAL = ["sex", "source_type", "input_build", "population", "delimiter", "format", "study_yaml"]
HEADER = "\t".join(REQUIRED) + "\n"
FULL_HEADER = "\t".join(REQUIRED[:3] + ["sex"] + REQUIRED[3:] + OPTIONAL[1:]) + "\n"
OUTPUT_FILES = {"standard_vcf": "standardised.gwas.vcf.gz", "standard_index": "standardised.gwas.vcf.gz.tbi",
                "metadata_json": "standardised.metadata.json", "qc_json": "standardised.qc.json",
                "report_html": "standardised.report.html"}
RAW = EXAMPLES / "tables" / "generic.tsv"


class Case:
    """1 isolated run of prepare_manifest.py: a manifests folder, a work root, a state file and a rows folder."""

    def __init__(self, base: Path, name: str, reference_root: Path):
        self.dir = base / name
        self.config_root = self.dir / "manifests"
        self.config_root.mkdir(parents=True)
        self.reference_root = reference_root
        self.state_file = self.dir / "work" / ".standardise_state.tsv"
        self.path = self.config_root / "manifest.tsv"

    def manifest(self, body: str, header: str = HEADER) -> Path:
        self.path.write_text(header + body, encoding="utf-8")
        return self.path

    def prepare(self, *extra, check: bool = True):
        args = ["--manifest", self.path, "--work-root", self.dir / "work", "--state-file", self.state_file,
                "--output", self.dir / "run.tsv", "--summary", self.dir / "summary.json", "--rows-dir", self.dir / "rows",
                "--config-root", self.config_root, "--reference-root", self.reference_root, *extra]
        return run_script("prepare_manifest.py", *args, check=check)

    def rows(self) -> list[dict]:
        with open(self.dir / "run.tsv", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def row_jsons(self) -> list[dict]:
        return [json.loads(p.read_text(encoding="utf-8")) for p in sorted((self.dir / "rows").glob("*.json"))]

    def summary(self) -> dict:
        return json.loads((self.dir / "summary.json").read_text(encoding="utf-8"))

    def complete(self, row: dict, report: bool = True) -> dict[str, Path]:
        """Write a row's outputs and record it with update_state.py, as the end of a successful run does."""
        prefix = Path(row["output_prefix"])
        prefix.parent.mkdir(parents=True, exist_ok=True)
        files = {key: prefix.parent / name for key, name in OUTPUT_FILES.items()}
        for key, path in files.items():
            if key != "report_html" or report:
                path.write_text("ok", encoding="utf-8")
        record = {"row_hash": row["row_hash"], "study_id": row["study_id"], "output_id": row["output_id"],
                  "GWAS_location": row["GWAS_location"], **{k: str(v) for k, v in files.items()},
                  "qc": {"settings_hash": row["settings_hash"]}}
        record_path = self.dir / f"{row['output_id']}.row_record.json"
        record_path.write_text(json.dumps(record), encoding="utf-8")
        run_script("update_state.py", "--state-file", self.state_file, "--report", str(report).lower(),
                   "--output", self.dir / "state_update_summary.json", "--records", record_path)
        return files


class ManifestTest(unittest.TestCase):
    def assert_rejected(self, c: Case, proc, run: int) -> None:
        """The manifest was prepared, `run` rows will run, and the rows failing a check were set aside as failed."""
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(len(c.row_jsons()), run)
        self.assertGreater(c.summary()["rows_invalid"], 0)
        self.assertIn("not run (the other rows run)", proc.stderr)

    def test_manifest(self):
        base = Path(tempfile.mkdtemp(prefix="gwas_test_manifest_"))
        reference = make_reference(base / "references")
        counter = iter(range(1000))

        def case(reference_root: Path = reference) -> Case:
            return Case(base, f"case{next(counter)}", reference_root)

        # ---- Columns ------------------------------------------------------------------------------------------------
        with self.subTest("a manifest without a required column is rejected, naming the column"):
            path = base / "missing_pmid.tsv"
            path.write_text("GWAS_location\tphenotype\tancestry\tauthor\tyear\nx\tBMI\tEUR\tA\t2020\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "PMID"):
                common.read_manifest(path)

        with self.subTest("a row with an empty GWAS_location is rejected"):
            path = base / "no_location.tsv"
            path.write_text(HEADER + "\tBMI\tEUR\tLocke\t2015\t25673413\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "empty GWAS_location"):
                common.read_manifest(path)

        with self.subTest("blank rows are skipped, values are trimmed and absent optional columns read as empty"):
            path = base / "blank.tsv"
            path.write_text(HEADER + " ieu-a-2 \tBMI\tALL\tLocke\t2015\t25673413\n\t\t\t\t\t\n", encoding="utf-8")
            rows = common.read_manifest(path)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["GWAS_location"], "ieu-a-2")
            self.assertEqual({col: rows[0][col] for col in OPTIONAL}, {col: "" for col in OPTIONAL})

        with self.subTest("sex must be combined, male or female (any case), or empty"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tFemale\tLocke\t2015\t25673413\n"
                       f"{RAW}\tWHR\tEUR\twomen\tPulit\t2018\t30239722\n",
                       "GWAS_location\tphenotype\tancestry\tsex\tauthor\tyear\tPMID\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=1)
            self.assertIn("row 3: sex 'women' must be combined, male or female", proc.stderr)
            self.assertNotIn("row 2", proc.stderr)

        with self.subTest("format must be table or vcf, or empty to detect it"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\tvcf\n{RAW}\tWHR\tEUR\tPulit\t2018\t30239722\tparquet\n",
                       HEADER.rstrip("\n") + "\tformat\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=1)
            self.assertIn("row 3: format 'parquet' must be table or vcf", proc.stderr)
            self.assertNotIn("row 2", proc.stderr)

        with self.subTest("input_build must be hg19 or hg38 (GRCh37 and GRCh38 labels accepted), or empty"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\tGRCh37\n{RAW}\tWHR\tEUR\tPulit\t2018\t30239722\thg18\n",
                       HEADER.rstrip("\n") + "\tinput_build\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=1)
            self.assertIn("row 3: input_build 'hg18' must be hg19 (GRCh37) or hg38 (GRCh38)", proc.stderr)
            self.assertNotIn("row 2", proc.stderr)

        with self.subTest("population must be a 1000 Genomes super-population code, or empty"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\teur\n{RAW}\tWHR\tEUR\tPulit\t2018\t30239722\tCEU\n",
                       HEADER.rstrip("\n") + "\tpopulation\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=1)
            self.assertIn("row 3: population 'CEU' must be ALL, AFR, AMR, EAS, EUR or SAS", proc.stderr)
            self.assertNotIn("row 2", proc.stderr)

        with self.subTest("source_type must be gwas_catalog, opengwas, direct_url or local"):
            for given, expected in {"gwas_catalog": "gwas_catalog", "GWAS Catalog": "gwas_catalog", "opengwas": "opengwas",
                                    "direct_url": "direct_url", "url": "direct_url", "local": "local"}.items():
                self.assertEqual(common.normalize_source_type(given), expected, given)
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\tdropbox\n", HEADER.rstrip("\n") + "\tsource_type\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=0)
            self.assertIn("unsupported source_type 'dropbox'", proc.stderr)

        with self.subTest("a GWAS Catalog or OpenGWAS row must carry its accession or dataset ID"):
            c = case()
            c.manifest("https://example.org/no-accession\tBMI\tEUR\tLocke\t2015\t25673413\tgwas_catalog\n"
                       "https://example.org/no-dataset\tWHR\tEUR\tPulit\t2018\t30239722\topengwas\n",
                       HEADER.rstrip("\n") + "\tsource_type\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=0)
            self.assertIn("row 2: no GWAS Catalog accession", proc.stderr)
            self.assertIn("row 3: no OpenGWAS dataset ID", proc.stderr)

        with self.subTest("every optional column is preserved in the run manifest and the row JSON"):
            c = case()
            (c.config_root / "custom.yaml").write_text("description: custom\n", encoding="utf-8")
            c.manifest(f"{RAW}\tBMI\tall\tfemale\tLocke\t2015\t25673413\tlocal\thg19\tEUR\tspace\ttable\tcustom.yaml\n",
                       FULL_HEADER)
            c.prepare()
            expected = {"sex": "female", "source_type": "local", "input_build": "hg19", "population": "EUR",
                        "delimiter": "space", "format": "table", "study_yaml": "custom.yaml"}
            self.assertEqual({col: c.rows()[0][col] for col in OPTIONAL}, expected)
            self.assertEqual({col: c.row_jsons()[0][col] for col in OPTIONAL}, expected)

        with self.subTest("every optional column is part of the row hash and the settings hash"):
            plain = {"GWAS_location": str(RAW), "phenotype": "BMI", "ancestry": "EUR", "author": "Locke", "year": "2015",
                     "PMID": "25673413", "study_id": "locke_2015_25673413"}
            for col, value in {"sex": "male", "source_type": "local", "input_build": "hg38", "population": "AFR",
                               "delimiter": "comma", "format": "vcf", "study_yaml": "x.yaml"}.items():
                changed = {**plain, col: value}
                self.assertNotEqual(common.row_hash(plain), common.row_hash(changed), col)
                self.assertNotEqual(common.settings_fingerprint(plain, {}, reference),
                                    common.settings_fingerprint(changed, {}, reference), col)
            self.assertEqual(common.row_hash(plain), common.row_hash({**plain, "sex": ""}),
                             "an empty optional column must hash as an absent one")

        with self.subTest("manifest input_build, population, delimiter and format override the study YAML"):
            c = case()
            (c.config_root / "locke_2015_25673413.yaml").write_text(
                "input_build: hg38\nreference_population: AFR\ndelimiter: ','\nformat: vcf\npalindromic_tolerance: 0.1\n",
                encoding="utf-8")
            record = {"study_id": "locke_2015_25673413", "input_build": "hg19", "population": "eur", "delimiter": "tab",
                      "format": "table", "study_yaml": ""}
            config = common.resolve_row_config(c.config_root, record)
            self.assertEqual((config["input_build"], config["reference_population"], config["delimiter"], config["format"],
                              config["palindromic_tolerance"]), ("hg19", "EUR", "\t", "table", 0.1))
            delimiters = {"tab": "\t", "comma": ",", "space": r"\s+", ";": ";", "|": "|", "": None}
            for given, expected in delimiters.items():
                self.assertEqual(common.normalize_delimiter_value(given), expected, given)

        # ---- Ancestry ----------------------------------------------------------------------------------------------
        with self.subTest("ancestry codes and unambiguous labels resolve to a 1000 Genomes super-population"):
            cases = {"AFR": "AFR", "amr": "AMR", "EAS": "EAS", "EUR": "EUR", "SAS": "SAS", "ALL": "ALL", "all": "ALL",
                     "European": "EUR", "East Asian": "EAS", "south asian": "SAS", "African American": "AFR",
                     "Hispanic": "AMR", "Mixed": "ALL", "multi-ancestry": "ALL"}
            for raw, expected in cases.items():
                self.assertEqual(common.resolve_ancestry(raw), expected, raw)

        with self.subTest("blank, unknown and ambiguous ancestry labels do not resolve"):
            for raw in ["", "unknown", "asian", "american", "Middle East", "european ancestry"]:
                self.assertIsNone(common.resolve_ancestry(raw), raw)

        with self.subTest("a row that fails a check is listed as failed with its reason, and the other rows run"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\n{RAW}\tWHR\tunknown\tPulit\t2018\t30239722\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=1)
            self.assertIn("row 3", proc.stderr)
            self.assertIn("ancestry", proc.stderr)
            self.assertEqual([r["study_id"] for r in c.row_jsons()], ["locke_2015_25673413"])
            invalid = c.summary()["invalid"]
            self.assertEqual([r["study_id"] for r in invalid], ["pulit_2018_30239722"])
            self.assertIn("ancestry 'unknown'", invalid[0]["errors"][0])
            self.assertEqual(sorted(r["study_id"] for r in c.rows()), ["locke_2015_25673413", "pulit_2018_30239722"],
                             "rejected rows are in the run manifest, so the run lists them as failed")

        with self.subTest("a problem with the manifest as a whole (2 rows writing the same outputs) stops the run before any row starts"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\n{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\n")
            proc = c.prepare(check=False)
            self.assertEqual(proc.returncode, 78)
            self.assertIn("write to the same output directory", proc.stderr)
            self.assertFalse((c.dir / "run.tsv").exists())

        # ---- Sources and IDs ---------------------------------------------------------------------------------------
        with self.subTest("source_type is worked out from GWAS_location"):
            cases = {"https://www.ebi.ac.uk/gwas/studies/GCST002783": "gwas_catalog", "GCST002783": "gwas_catalog",
                     "gcst002783": "gwas_catalog", "https://opengwas.io/datasets/ebi-a-GCST002783": "opengwas",
                     "https://gwas.mrcieu.ac.uk/datasets/ieu-a-2/": "opengwas", "ebi-a-GCST002783": "opengwas",
                     "https://example.org/sumstats/gwas.tsv.gz": "direct_url", "ftp://example.org/gwas.txt.gz": "direct_url",
                     str(RAW): "local", "data/gwas.tsv": "local"}
            for location, expected in cases.items():
                self.assertEqual(common.infer_source_type(location), expected, location)

        with self.subTest("OpenGWAS IDs with text or suffixes are recognised and keep their case"):
            for dataset in ["ieu-b-5129", "ukb-d-30780_irnt", "ebi-a-GCST90000123", "finn-b-I9_HEARTFAIL", "met-d-LDL_C"]:
                self.assertEqual(common.infer_source_type(dataset), "opengwas", dataset)
                self.assertEqual(common.opengwas_dataset_id(f"https://opengwas.io/datasets/{dataset}"), dataset)
                self.assertEqual(common.output_id_from_row({"GWAS_location": dataset}, "opengwas"), f"openGWAS_{dataset}")

        with self.subTest("an existing local file named like an OpenGWAS ID stays local"):
            here = base / "named_like_opengwas"
            here.mkdir()
            (here / "ieu-a-2").write_text(RAW.read_text(encoding="utf-8"), encoding="utf-8")
            cwd = os.getcwd()
            os.chdir(here)
            try:
                self.assertEqual(common.infer_source_type("ieu-a-2"), "local")
            finally:
                os.chdir(cwd)
            c = case()  # beside the manifest, not in the working directory: source_type local says so
            (c.config_root / "ieu-a-2").write_text(RAW.read_text(encoding="utf-8"), encoding="utf-8")
            c.manifest("ieu-a-2\tBMI\tEUR\tLocke\t2015\t25673413\tlocal\n", HEADER.rstrip("\n") + "\tsource_type\n")
            c.prepare()
            (row,) = c.row_jsons()
            self.assertEqual((row["source_type"], row["GWAS_location"]), ("local", str((c.config_root / "ieu-a-2").resolve())))
            self.assertTrue(row["output_id"].startswith("local_ieu-a-2_"))

        with self.subTest("output IDs name the source: GWAScatalog_, openGWAS_, directURL_ or local_"):
            c = case()
            c.manifest("https://www.ebi.ac.uk/gwas/studies/gcst002783\tBMI\tALL\tLocke\t2015\t25673413\n"
                       "https://opengwas.io/datasets/ebi-a-GCST002783\tBMI\tALL\tLocke\t2015\t25673413\n"
                       "https://example.org/sumstats/bmi.tsv.gz\tBMI\tALL\tLocke\t2015\t25673413\n"
                       f"{RAW}\tBMI\tALL\tLocke\t2015\t25673413\n")
            c.prepare()
            rows = c.rows()
            self.assertEqual([r["source_type"] for r in rows], ["gwas_catalog", "opengwas", "direct_url", "local"])
            self.assertEqual(rows[0]["output_id"], "GWAScatalog_GCST002783")
            self.assertEqual(rows[1]["output_id"], "openGWAS_ebi-a-GCST002783")
            self.assertRegex(rows[2]["output_id"], r"^directURL_bmi_[0-9a-f]{8}$")
            self.assertRegex(rows[3]["output_id"], r"^local_generic_[0-9a-f]{8}$")

        with self.subTest("author, year and PMID name the study folder, which holds standardise/<output ID>/"):
            self.assertEqual(common.study_id_from_row({"author": "Locke ", "year": "2015", "PMID": "PMID:25673413"}),
                             "locke_2015_25673413")
            self.assertEqual(common.study_id_from_row({"author": "O'Mara", "year": "2018", "PMID": "30093612"}),
                             "o_mara_2018_30093612")
            c = case()
            c.manifest("GCST002783\tBMI\tALL\tLocke\t2015\t25673413\n")
            c.prepare()
            (row,) = c.rows()
            self.assertEqual(row["study_id"], "locke_2015_25673413")
            self.assertEqual(Path(row["standardise_dir"]), c.dir / "work" / "locke_2015_25673413" / "standardise" / "GWAScatalog_GCST002783")
            self.assertEqual(row["output_prefix"], str(Path(row["standardise_dir"]) / "standardised"))

        with self.subTest("relative local paths resolve against the manifest directory"):
            c = case()
            (c.config_root / "data").mkdir()
            local = c.config_root / "data" / "gwas.tsv"
            local.write_text(RAW.read_text(encoding="utf-8"), encoding="utf-8")
            c.manifest("data/gwas.tsv\tBMI\tEUR\tLocke\t2015\t25673413\n")
            c.prepare()
            self.assertEqual(c.row_jsons()[0]["GWAS_location"], str(local.resolve()))

        with self.subTest("row JSON files carry every field, manifest text unchanged"):
            c = case()
            c.manifest(f'{RAW}\tBMI "adjusted" $HOME `x`\tEUR\tLocke\t2015\t25673413\n')
            c.prepare()
            (row,) = c.row_jsons()
            self.assertEqual(row["phenotype"], 'BMI "adjusted" $HOME `x`')
            self.assertEqual(row["GWAS_location"], str(RAW))
            self.assertEqual(len(row["settings_hash"]), 16)
            for key in REQUIRED + OPTIONAL + ["study_id", "output_id", "row_hash", "raw_dir", "standardise_dir", "output_prefix"]:
                self.assertIn(key, row)

        # ---- Study YAML --------------------------------------------------------------------------------------------
        with self.subTest("a study YAML is found in the manifests folder or the repository root"):
            repo = base / "repo"
            (repo / "manifests").mkdir(parents=True)
            (repo / "manifests" / "x.yaml").write_text("study_type: casecontrol\n", encoding="utf-8")
            (repo / "manifests" / "example_2020_1.yaml").write_text("input_build: hg19\n", encoding="utf-8")
            manifests = repo / "manifests"
            self.assertEqual(common.resolve_study_config(manifests, "s", None, "x.yaml")["study_type"], "casecontrol")
            self.assertEqual(common.resolve_study_config(manifests, "s", None, "manifests/x.yaml")["study_type"], "casecontrol")
            self.assertEqual(common.resolve_study_config(manifests, "example_2020_1")["input_build"], "hg19",
                             "<study ID>.yaml in the manifests folder is used without being named")
            self.assertEqual(common.resolve_study_config(manifests, "locke_2015_25673413"), {})

        with self.subTest("a named study YAML must exist"):
            with self.assertRaises(SystemExit) as caught, contextlib.redirect_stderr(io.StringIO()):
                common.resolve_study_config(base / "repo" / "manifests", "s", None, "manifests/missing.yaml")
            self.assertEqual(caught.exception.code, 65)
            c = case()
            c.manifest(f"{RAW}\tBMI\tEUR\tLocke\t2015\t25673413\tmissing.yaml\n", HEADER.rstrip("\n") + "\tstudy_yaml\n")
            proc = c.prepare(check=False)
            self.assert_rejected(c, proc, run=0)
            self.assertIn("missing.yaml", proc.stderr)

        # ---- Selecting rows ----------------------------------------------------------------------------------------
        two_rows = f"{RAW}\tBMI\tall\tLocke\t2015\t25673413\n{RAW}\tWHR\tall\tPulit\t2018\t30239722\n"

        with self.subTest("--row-limit runs the first N eligible rows and defers the rest; 0 runs every row"):
            c = case()
            c.manifest(two_rows)
            c.prepare("--row-limit", "1")
            summary = c.summary()
            self.assertEqual((len(c.rows()), len(c.row_jsons())), (1, 1))
            self.assertEqual((summary["rows_eligible"], summary["rows_to_run"], summary["rows_deferred_due_to_limit"]), (2, 1, 1))
            c.prepare("--row-limit", "0")
            self.assertEqual((len(c.rows()), len(c.row_jsons()), c.summary()["rows_deferred_due_to_limit"]), (2, 2, 0))

        with self.subTest("--only-study runs 1 study"):
            c = case()
            c.manifest(two_rows)
            c.prepare("--only-study", "pulit_2018_30239722")
            self.assertEqual([r["study_id"] for r in c.rows()], ["pulit_2018_30239722"])
            self.assertEqual(c.summary()["rows_filtered"], 1)

        with self.subTest("a completed row is skipped while its outputs exist and its settings are unchanged"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tall\tLocke\t2015\t25673413\n")
            c.prepare()
            c.complete(c.row_jsons()[0])
            c.prepare()
            self.assertEqual((c.rows(), c.row_jsons()), ([], []))
            self.assertEqual(c.summary()["rows_skipped_complete"], 1)

        with self.subTest("--force reruns completed rows"):
            c.prepare("--force", "true")
            self.assertEqual(len(c.rows()), 1)

        with self.subTest("a completed row whose output was deleted or emptied runs again"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tall\tLocke\t2015\t25673413\n")
            c.prepare()
            files = c.complete(c.row_jsons()[0])
            files["standard_index"].write_text("", encoding="utf-8")
            c.prepare()
            self.assertEqual(len(c.rows()), 1)
            files["standard_index"].unlink()
            c.prepare()
            self.assertEqual(len(c.rows()), 1)

        with self.subTest("a change to the manifest row, study YAML, reference files or output version reruns a completed row"):
            c = case(make_reference(base / "own_references"))
            c.manifest(f"{RAW}\tBMI\tall\tLocke\t2015\t25673413\n")
            c.prepare()
            first = c.row_jsons()[0]
            c.complete(first)
            c.prepare()
            self.assertEqual(len(c.rows()), 0, "unchanged settings must not rerun a complete row")
            (c.config_root / "locke_2015_25673413.yaml").write_text("palindromic_tolerance: 0.1\n", encoding="utf-8")
            c.prepare()
            self.assertEqual(len(c.rows()), 1, "a new study YAML must rerun the row")
            c.complete(c.row_jsons()[0])
            c.prepare()
            self.assertEqual(len(c.rows()), 0)
            with open(c.reference_root / "dbsnp" / "GRCh38" / "meta.json", "a", encoding="utf-8") as handle:
                handle.write("\n")
            c.prepare()
            self.assertEqual(len(c.rows()), 1, "a changed reference file must rerun the row")
            c.complete(c.row_jsons()[0])
            c.manifest(f"{RAW}\tBMI\tall\tLocke\t2015\t25673413\thg19\n", HEADER.rstrip("\n") + "\tinput_build\n")
            c.prepare()
            self.assertEqual(len(c.rows()), 1, "a changed manifest row must rerun")
            record = {**first, "study_id": first["study_id"]}
            fingerprint = common.settings_fingerprint(record, {}, c.reference_root)
            with mock.patch.object(common, "OUTPUT_SCHEMA_VERSION", "999.0.0"):
                self.assertNotEqual(common.settings_fingerprint(record, {}, c.reference_root), fingerprint,
                                    "a new pipeline output version must change the settings hash")

        with self.subTest("with --report false a completed row needs no HTML report"):
            c = case()
            c.manifest(f"{RAW}\tBMI\tall\tLocke\t2015\t25673413\n")
            c.prepare("--report", "false")
            c.complete(c.row_jsons()[0], report=False)
            c.prepare("--report", "false")
            self.assertEqual(len(c.rows()), 0)
            c.prepare("--report", "true")
            self.assertEqual(len(c.rows()), 1, "a report run must redo a row that has no report")

        # ---- State -------------------------------------------------------------------------------------------------
        with self.subTest("complete means every output file exists and is non-empty"):
            folder = base / "state_complete"
            folder.mkdir()
            paths = {key: folder / name for key, name in OUTPUT_FILES.items()}
            for path in paths.values():
                path.write_text("x", encoding="utf-8")
            row = {key: str(path) for key, path in paths.items()}
            self.assertTrue(common.state_row_complete(row, report=True))
            paths["report_html"].write_text("", encoding="utf-8")
            self.assertFalse(common.state_row_complete(row, report=True))
            self.assertTrue(common.state_row_complete(row, report=False))
            paths["standard_vcf"].unlink()
            self.assertFalse(common.state_row_complete(row, report=False))

        with self.subTest("a failed JSON write leaves the previous file intact"):
            folder = base / "json_write"
            folder.mkdir()
            path = folder / "qc.json"
            common.write_json(path, {"ok": 1})
            with self.assertRaises(TypeError):
                common.write_json(path, {"bad": object()})
            self.assertEqual(json.loads(path.read_text(encoding="utf-8")), {"ok": 1})
            self.assertEqual(sorted(p.name for p in folder.iterdir()), ["qc.json"])

        with self.subTest("a failed state write leaves the previous state intact"):
            folder = base / "state_write"
            folder.mkdir()
            path = folder / ".standardise_state.tsv"
            common.write_state(path, [{"row_hash": "h1", "study_id": "s"}])
            before = path.read_text(encoding="utf-8")

            def rows():
                yield {"row_hash": "h2"}
                raise RuntimeError("interrupted")

            with self.assertRaises(RuntimeError):
                common.write_state(path, rows())
            self.assertEqual(path.read_text(encoding="utf-8"), before)
            self.assertEqual(list(common.read_state(path)), ["h1"])
            self.assertEqual(sorted(p.name for p in folder.iterdir()), [".standardise_state.tsv"])

        with self.subTest("update_state records complete rows with their settings hash and keeps earlier rows"):
            c = case()
            c.manifest(two_rows)
            c.prepare("--row-limit", "0")
            first, second = c.row_jsons()
            c.complete(first)
            summary = json.loads((c.dir / "state_update_summary.json").read_text(encoding="utf-8"))
            self.assertEqual((summary["updated"], summary["updated_row_hashes"]), (1, [first["row_hash"]]))
            c.complete(second)
            state = common.read_state(c.state_file)
            self.assertEqual(sorted(state), sorted([first["row_hash"], second["row_hash"]]))
            self.assertEqual(state[first["row_hash"]]["settings_hash"], first["settings_hash"])

        with self.subTest("update_state ignores incomplete and missing records, giving reasons"):
            folder = base / "update_state"
            folder.mkdir()
            out = folder / "noreport"
            out.mkdir()
            paths = {key: out / name for key, name in OUTPUT_FILES.items()}
            for key, path in paths.items():
                if key != "report_html":
                    path.write_text("x", encoding="utf-8")
            record = folder / "noreport.row_record.json"
            record.write_text(json.dumps({"row_hash": "h", **{k: str(v) for k, v in paths.items()}}), encoding="utf-8")
            absent = folder / "absent.row_record.json"
            state, summary_path = folder / ".standardise_state.tsv", folder / "summary.json"
            run_script("update_state.py", "--state-file", state, "--report", "true", "--output", summary_path,
                       "--records", record, absent)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual((summary["updated"], summary["ignored"]), (0, 2))
            reasons = {item["record"]: item["reason"] for item in summary["ignored_records"]}
            self.assertEqual((reasons[str(absent)], reasons[str(record)]), ("record file not found", "outputs incomplete"))
            self.assertEqual(common.read_state(state), {})

        with self.subTest("update_state with no records still writes a summary"):
            summary_path = base / "update_state" / "empty_summary.json"
            run_script("update_state.py", "--state-file", base / "update_state" / "empty.tsv", "--report", "true",
                       "--output", summary_path)
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertEqual((summary["updated"], summary["ignored"]), (0, 0))


if __name__ == "__main__":
    unittest.main()
