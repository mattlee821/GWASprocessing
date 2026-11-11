"""Tests for gwas_pipeline_common: parsing, identifiers, normalisation and persistent state."""

import json
import math
import os
import tempfile
import unittest
from pathlib import Path

import testlib  # noqa: F401  (puts pipeline/bin on sys.path)

import gwas_pipeline_common as common


class DelimiterAndTableReadingTests(unittest.TestCase):
    def write(self, text, name="table.txt"):
        tmp = Path(tempfile.mkdtemp())
        path = tmp / name
        path.write_text(text, encoding="utf-8")
        return path

    def test_tab_file_is_detected_as_tab_even_though_line_ends_in_newline(self):
        self.assertEqual(common.infer_delimiter(self.write("CHR\tPOS\tSNP\n1\t100\trs1\n")), "\t")

    def test_tab_file_with_windows_line_endings_is_tab(self):
        self.assertEqual(common.infer_delimiter(self.write("CHR\tPOS\r\n1\t100\r\n")), "\t")

    def test_comma_file_is_detected_as_comma(self):
        self.assertEqual(common.infer_delimiter(self.write("CHR,POS,SNP\n1,100,rs1\n")), ",")

    def test_space_aligned_file_is_detected_as_whitespace(self):
        self.assertEqual(common.infer_delimiter(self.write("SNP  A1 A2   b\nrs1  A  G    0.1\n")), r"\s+")

    def test_leading_double_hash_metadata_lines_are_ignored(self):
        path = self.write("## produced by tool x, version 2\n## another note\nCHR\tPOS\tP\n1\t100\t0.5\n")
        self.assertEqual(common.infer_delimiter(path), "\t")
        table = common.read_table(path)
        self.assertEqual(list(table.columns), ["CHR", "POS", "P"])
        self.assertEqual(table.loc[0, "P"], 0.5)

    def test_empty_tab_cell_does_not_shift_later_columns(self):
        path = self.write("SNP\tEAF\tBETA\tSE\tP\tN\nrs1\t\t0.10\t0.02\t1e-6\t5000\n")
        row = common.read_table(path).iloc[0]
        self.assertTrue(math.isnan(row["EAF"]))
        self.assertEqual((row["BETA"], row["SE"], row["P"], row["N"]), (0.10, 0.02, 1e-6, 5000))

    def test_dot_is_read_as_missing(self):
        path = self.write("SNP\tEAF\nrs1\t.\n")
        self.assertTrue(math.isnan(common.read_table(path).loc[0, "EAF"]))


class ChromosomeAndSourceTests(unittest.TestCase):
    def test_normalize_chromosome_aliases(self):
        self.assertEqual(common.normalize_chromosome("chr23"), "X")
        self.assertEqual(common.normalize_chromosome("24"), "Y")
        self.assertEqual(common.normalize_chromosome("chrM"), "MT")

    def test_plink_codes_25_and_26_are_pseudo_autosomal_x_and_mitochondrial(self):
        self.assertEqual(common.normalize_chromosome("25"), "X")
        self.assertEqual(common.normalize_chromosome("26"), "MT")
        self.assertEqual(common.normalize_chromosome("XY"), "X")

    def test_synapse_requires_a_syn_accession(self):
        self.assertEqual(common.infer_source_type("syn26123456"), "synapse")
        self.assertEqual(common.infer_source_type("synovitis_gwas.tsv.gz"), "local")

    def test_opengwas_ids_with_text_or_suffixes_are_recognised(self):
        for dataset in ["ieu-b-5129", "ukb-d-30780_irnt", "ebi-a-GCST90000123", "finn-b-I9_HEARTFAIL", "met-d-LDL_C"]:
            with self.subTest(dataset=dataset):
                self.assertEqual(common.infer_source_type(dataset), "opengwas")
                self.assertEqual(common.opengwas_dataset_id(f"https://opengwas.io/datasets/{dataset}"), dataset)

    def test_existing_local_file_named_like_an_opengwas_id_stays_local(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "ieu-a-2").write_text("SNP\tP\n", encoding="utf-8")
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            self.assertEqual(common.infer_source_type("ieu-a-2"), "local")
        finally:
            os.chdir(cwd)

    def test_opengwas_output_ids_keep_case_and_suffix(self):
        self.assertEqual(common.output_id_from_row({"GWAS_location": "ebi-a-GCST90000123"}, "opengwas"), "openGWAS_ebi-a-GCST90000123")
        self.assertEqual(
            common.output_id_from_row({"GWAS_location": "https://opengwas.io/datasets/ukb-d-30780_irnt"}, "opengwas"),
            "openGWAS_ukb-d-30780_irnt",
        )

    def test_study_id_strips_pmid_prefix_and_punctuation(self):
        self.assertEqual(common.study_id_from_row({"author": "Locke ", "year": "2015", "PMID": "PMID:25673413"}), "locke_2015_25673413")


class AncestryAndBuildTests(unittest.TestCase):
    def test_codes_and_unambiguous_synonyms_resolve(self):
        cases = {"EUR": "EUR", "european": "EUR", "East Asian": "EAS", "south asian": "SAS",
                 "African American": "AFR", "all": "ALL", "Mixed": "ALL", "Hispanic": "AMR"}
        for raw, expected in cases.items():
            with self.subTest(raw=raw):
                self.assertEqual(common.resolve_ancestry(raw), expected)

    def test_blank_unknown_and_ambiguous_labels_do_not_resolve(self):
        for raw in ["", "unknown", "asian", "american", "Middle East", "european ancestry"]:
            with self.subTest(raw=raw):
                self.assertIsNone(common.resolve_ancestry(raw))

    def test_build_labels_from_sources_are_normalised(self):
        cases = {"HG19/GRCh37": "hg19", "GRCh38.p13": "hg38", "b37": "hg19", "hg38": "hg38", "37": "hg19",
                 "38": "hg38", "NCBI36": "", "": "", "hg19 or hg38": ""}
        for raw, expected in cases.items():
            with self.subTest(raw=raw):
                self.assertEqual(common.normalize_build_value(raw), expected)


class ConfigAndStateTests(unittest.TestCase):
    def test_config_resolution_from_manifests_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "example_2020_1.yaml").write_text("input_build: hg19\n", encoding="utf-8")
            (tmp / "explicit.yaml").write_text("input_build: hg38\n", encoding="utf-8")
            self.assertEqual(common.resolve_study_config(tmp, "example_2020_1")["input_build"], "hg19")
            self.assertEqual(common.resolve_study_config(tmp, "locke_2015_25673413"), {})
            self.assertEqual(common.resolve_study_config(tmp, "example_2020_1", str(tmp / "explicit.yaml"))["input_build"], "hg38")

    def test_merge_dicts_merges_nested_mappings_without_mutating_inputs(self):
        base = {"analysis_qc": {"maf_threshold": 0.01, "deduplicate": True}, "delimiter": "\t"}
        merged = common.merge_dicts(base, {"analysis_qc": {"maf_threshold": 0.05}})
        self.assertEqual(merged["analysis_qc"], {"maf_threshold": 0.05, "deduplicate": True})
        self.assertEqual(base["analysis_qc"]["maf_threshold"], 0.01)

    def test_failed_json_write_leaves_previous_file_intact(self):
        tmp = Path(tempfile.mkdtemp())
        path = tmp / "qc.json"
        common.write_json(path, {"ok": 1})
        with self.assertRaises(TypeError):
            common.write_json(path, {"bad": object()})
        self.assertEqual(json.loads(path.read_text(encoding="utf-8")), {"ok": 1})
        self.assertEqual(sorted(p.name for p in tmp.iterdir()), ["qc.json"])

    def test_failed_state_write_leaves_previous_state_intact(self):
        tmp = Path(tempfile.mkdtemp())
        path = tmp / ".standardise_state.tsv"
        common.write_state(path, [{"row_hash": "h1", "study_id": "s"}])
        before = path.read_text(encoding="utf-8")

        def rows():
            yield {"row_hash": "h2"}
            raise RuntimeError("interrupted")

        with self.assertRaises(RuntimeError):
            common.write_state(path, rows())
        self.assertEqual(path.read_text(encoding="utf-8"), before)
        self.assertEqual(list(common.read_state(path)), ["h1"])

    def test_state_row_complete_requires_every_file_to_exist_and_be_non_empty(self):
        tmp = Path(tempfile.mkdtemp())
        paths = {key: tmp / f"{key}.out" for key in ["standard_vcf", "standard_index", "metadata_json", "qc_json", "report_html"]}
        for path in paths.values():
            path.write_text("x", encoding="utf-8")
        row = {key: str(path) for key, path in paths.items()}
        self.assertTrue(common.state_row_complete(row, report=True))
        paths["report_html"].write_text("", encoding="utf-8")
        self.assertFalse(common.state_row_complete(row, report=True))
        self.assertTrue(common.state_row_complete(row, report=False))


class ReferenceIdentityTests(unittest.TestCase):
    def test_appledouble_files_do_not_change_the_identity(self):
        root = Path(tempfile.mkdtemp())
        (root / "genome").mkdir()
        (root / "genome" / "hg19ToHg38.over.chain.gz").write_bytes(b"chain")
        before = common.reference_identity(root)
        (root / "genome" / "._hg19ToHg38.over.chain.gz").write_bytes(b"\x00\x05\x16\x07")
        self.assertEqual(common.reference_identity(root), before)


class ManifestReadingTests(unittest.TestCase):
    def test_missing_required_column_is_rejected(self):
        tmp = Path(tempfile.mkdtemp())
        path = tmp / "m.tsv"
        path.write_text("GWAS_location\tphenotype\tancestry\tauthor\tyear\nx\tBMI\tEUR\tA\t2020\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "PMID"):
            common.read_manifest(path)

    def test_blank_rows_are_skipped_and_values_trimmed(self):
        tmp = Path(tempfile.mkdtemp())
        path = tmp / "m.tsv"
        path.write_text(
            "GWAS_location\tphenotype\tancestry\tauthor\tyear\tPMID\n"
            " ieu-a-2 \tBMI\tALL\tLocke\t2015\t25673413\n"
            "\t\t\t\t\t\n",
            encoding="utf-8",
        )
        rows = common.read_manifest(path)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["GWAS_location"], "ieu-a-2")
        self.assertEqual(rows[0]["sex"], "")


if __name__ == "__main__":
    unittest.main()


class StudyYamlResolutionTests(unittest.TestCase):
    def test_named_study_yaml_resolves_from_the_repository_root_and_must_exist(self):
        import tempfile
        repo = Path(tempfile.mkdtemp())
        (repo / "manifests").mkdir()
        (repo / "manifests" / "x.yaml").write_text("study_type: casecontrol\n", encoding="utf-8")
        config = common.resolve_study_config(repo / "manifests", "s", None, "manifests/x.yaml")
        self.assertEqual(config["study_type"], "casecontrol")
        with self.assertRaises(SystemExit) as caught:
            common.resolve_study_config(repo / "manifests", "s", None, "manifests/missing.yaml")
        self.assertEqual(caught.exception.code, 78)
