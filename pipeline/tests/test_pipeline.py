import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PIPELINE = ROOT / "pipeline"
BIN = PIPELINE / "bin"
FIXTURES = PIPELINE / "tests" / "fixtures"


class PipelineTests(unittest.TestCase):
    def run_py(self, script, *args):
        cmd = [sys.executable, str(BIN / script), *map(str, args)]
        return subprocess.run(cmd, cwd=ROOT, text=True, capture_output=True, check=True)

    def test_prepare_manifest_derives_ids_and_skip_state(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            manifest = tmp / "manifest.tsv"
            raw = FIXTURES / "tables/generic.tsv"
            manifest.write_text(
                "GWAS_location\tphenotype\tancestry\tauthor\tyear\tPMID\n"
                f"{raw}\tBMI\tall\tLocke\t2015\t25673413\n",
                encoding="utf-8",
            )
            run_manifest = tmp / "run.tsv"
            summary = tmp / "summary.json"
            state = tmp / ".standardise_state.tsv"
            self.run_py(
                "prepare_manifest.py",
                "--manifest",
                manifest,
                "--work-root",
                tmp / "work",
                "--state-file",
                state,
                "--output",
                run_manifest,
                "--summary",
                summary,
            )
            with open(run_manifest, encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["study_id"], "locke_2015_25673413")
            self.assertEqual(rows[0]["source_type"], "local")
            self.assertEqual(rows[0]["sex"], "")
            self.assertEqual(rows[0]["input_build"], "")
            self.assertEqual(rows[0]["study_yaml"], "")
            self.assertTrue(rows[0]["output_id"].startswith("local_generic_"))
            self.assertTrue(rows[0]["raw_dir"].endswith(f"raw/{rows[0]['output_id']}"))
            self.assertTrue(rows[0]["standardise_dir"].endswith(f"standardise/{rows[0]['output_id']}"))
            self.assertTrue(rows[0]["output_prefix"].endswith(f"standardise/{rows[0]['output_id']}/standardised"))

            standardise_dir = Path(rows[0]["standardise_dir"])
            standardise_dir.mkdir(parents=True)
            files = {
                "standard_vcf": standardise_dir / "standardised.gwas.vcf.gz",
                "standard_index": standardise_dir / "standardised.gwas.vcf.gz.tbi",
                "metadata_json": standardise_dir / "standardised.metadata.json",
                "qc_json": standardise_dir / "standardised.qc.json",
                "plot_tiff": standardise_dir / "standardised.tiff",
            }
            for path in files.values():
                path.write_text("ok", encoding="utf-8")
            state.write_text(
                "row_hash\tstudy_id\toutput_id\tgwas_location\tstandard_vcf\tstandard_index\tmetadata_json\tqc_json\tplot_tiff\tupdated_at\n"
                f"{rows[0]['row_hash']}\t{rows[0]['study_id']}\t{rows[0]['output_id']}\t{raw}\t"
                f"{files['standard_vcf']}\t{files['standard_index']}\t{files['metadata_json']}\t{files['qc_json']}\t{files['plot_tiff']}\tnow\n",
                encoding="utf-8",
            )
            self.run_py(
                "prepare_manifest.py",
                "--manifest",
                manifest,
                "--work-root",
                tmp / "work",
                "--state-file",
                state,
                "--output",
                run_manifest,
                "--summary",
                summary,
            )
            with open(run_manifest, encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows, [])

    def test_manifest_optional_fields_are_preserved_and_hashed(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = FIXTURES / "tables/generic.tsv"
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "GWAS_location\tphenotype\tancestry\tsex\tauthor\tyear\tPMID\tsource_type\tinput_build\tpopulation\tdelimiter\tformat\tstudy_yaml\n"
                f"{raw}\tBMI\tall\tfemale\tLocke\t2015\t25673413\tlocal\thg19\tALL\t\\s+\ttable\tcustom.yaml\n",
                encoding="utf-8",
            )
            run_manifest = tmp / "run.tsv"
            summary = tmp / "summary.json"
            state = tmp / ".standardise_state.tsv"
            self.run_py(
                "prepare_manifest.py",
                "--manifest",
                manifest,
                "--work-root",
                tmp / "work",
                "--state-file",
                state,
                "--output",
                run_manifest,
                "--summary",
                summary,
            )
            with open(run_manifest, encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["sex"], "female")
            self.assertEqual(rows[0]["input_build"], "hg19")
            self.assertEqual(rows[0]["population"], "ALL")
            self.assertEqual(rows[0]["delimiter"], "\\s+")
            self.assertEqual(rows[0]["format"], "table")
            self.assertEqual(rows[0]["study_yaml"], "custom.yaml")

            manifest.write_text(
                "GWAS_location\tphenotype\tancestry\tauthor\tyear\tPMID\tsource_type\tinput_build\n"
                f"{raw}\tBMI\tall\tLocke\t2015\t25673413\tlocal\thg38\n",
                encoding="utf-8",
            )
            self.run_py(
                "prepare_manifest.py",
                "--manifest",
                manifest,
                "--work-root",
                tmp / "work",
                "--state-file",
                state,
                "--output",
                run_manifest,
                "--summary",
                summary,
            )
            with open(run_manifest, encoding="utf-8") as handle:
                changed_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertNotEqual(rows[0]["row_hash"], changed_rows[0]["row_hash"])

    def test_prepare_manifest_row_limit_defers_remaining_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = FIXTURES / "tables/generic.tsv"
            manifest = tmp / "manifest.tsv"
            manifest.write_text(
                "GWAS_location\tphenotype\tancestry\tauthor\tyear\tPMID\n"
                f"{raw}\tBMI\tall\tLocke\t2015\t25673413\n"
                f"{raw}\tWHR\tall\tPulit\t2018\t30239722\n",
                encoding="utf-8",
            )
            run_manifest = tmp / "run.tsv"
            summary = tmp / "summary.json"
            self.run_py(
                "prepare_manifest.py",
                "--manifest",
                manifest,
                "--work-root",
                tmp / "work",
                "--state-file",
                tmp / ".standardise_state.tsv",
                "--output",
                run_manifest,
                "--summary",
                summary,
                "--row-limit",
                "1",
            )
            with open(run_manifest, encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            data = json.loads(summary.read_text(encoding="utf-8"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(data["rows_eligible"], 2)
            self.assertEqual(data["rows_to_run"], 1)
            self.assertEqual(data["rows_deferred_due_to_limit"], 1)
            self.assertEqual(data["row_limit"], 1)

    def test_config_resolution_from_manifests_root(self):
        sys.path.insert(0, str(BIN))
        from gwas_pipeline_common import resolve_study_config

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "example_2020_1.yaml").write_text("input_build: hg19\n", encoding="utf-8")
            (tmp / "explicit.yaml").write_text("input_build: hg38\n", encoding="utf-8")
            self.assertEqual(resolve_study_config(tmp, "example_2020_1")["input_build"], "hg19")
            self.assertEqual(resolve_study_config(tmp, "locke_2015_25673413"), {})
            self.assertEqual(resolve_study_config(tmp, "example_2020_1", str(tmp / "explicit.yaml"))["input_build"], "hg38")

    def test_gwas_catalog_selects_legacy_tbl_gz(self):
        sys.path.insert(0, str(BIN))
        import download_gwas_catalog

        original = download_gwas_catalog.list_absolute

        def fake_list_absolute(url):
            if url.endswith("/GCST008996/"):
                return [
                    url.rsplit("/", 2)[0] + "/",
                    url + "fat-distn.giant.ukbb.meta-analysis.whr.combined.tbl.gz",
                    url + "md5sum.txt",
                    url + "readme.txt",
                ]
            return []

        try:
            download_gwas_catalog.list_absolute = fake_list_absolute
            selected = download_gwas_catalog.choose_summary_file("https://example.org/GCST008996/", "GCST008996", "30239722")
            self.assertTrue(selected.endswith(".tbl.gz"))
        finally:
            download_gwas_catalog.list_absolute = original

    def test_gwas_catalog_md5_uses_exact_basename(self):
        sys.path.insert(0, str(BIN))
        from download_gwas_catalog import md5, verify_md5

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            downloaded = tmp / "file.tbl.gz"
            downloaded.write_text("actual", encoding="utf-8")
            sidecar = tmp / "md5sum.txt"
            sidecar.write_text(
                "00000000000000000000000000000000 ._file.tbl.gz\n"
                f"{md5(downloaded)} file.tbl.gz\n",
                encoding="utf-8",
            )
            result = verify_md5(downloaded, [str(sidecar)])
            self.assertTrue(result["ok"])

    def test_standardise_maps_snpid_and_reference_eaf(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = FIXTURES / "tables/snpid_missing_eaf.tsv"
            stage_record = tmp / "stage.json"
            output_prefix = tmp / "standardise/test"
            record = {
                "GWAS_location": str(raw),
                "phenotype": "BMI",
                "ancestry": "all",
                "sex": "female",
                "author": "Locke",
                "year": "2015",
                "PMID": "25673413",
                "study_id": "locke_2015_25673413",
                "output_id": "test",
                "row_hash": "abc",
                "source_type": "local",
                "raw_path": str(raw),
                "output_prefix": str(output_prefix),
                "standardise_dir": str(output_prefix.parent),
                "metadata": {},
            }
            stage_record.write_text(json.dumps(record), encoding="utf-8")
            out_record = tmp / "record.json"
            self.run_py(
                "standardise_gwas.py",
                "--stage-record",
                stage_record,
                "--config-root",
                PIPELINE / "config",
                "--reference-root",
                FIXTURES / "reference",
                "--output",
                out_record,
            )
            result = json.loads(out_record.read_text(encoding="utf-8"))
            vcf = Path(result["standard_vcf"])
            self.assertTrue(vcf.exists())
            text = subprocess.check_output(["gzip", "-cd", str(vcf)], text=True)
            self.assertIn("rs3", text)
            self.assertIn(":0.33:", text)
            qc = json.loads(Path(result["qc_json"]).read_text(encoding="utf-8"))
            self.assertGreaterEqual(qc["reference_fills"].get("rsid_from_snpid", 0), 1)
            self.assertGreaterEqual(qc["reference_fills"].get("EAF", 0), 1)

    def test_locke_like_direct_url_standardises_without_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = FIXTURES / "tables/locke_direct.tsv"
            stage_record = tmp / "stage.json"
            output_prefix = tmp / "standardise/test"
            record = {
                "GWAS_location": str(raw),
                "phenotype": "BMI",
                "ancestry": "all",
                "sex": "female",
                "author": "Locke",
                "year": "2015",
                "PMID": "25673413",
                "study_id": "locke_2015_25673413",
                "output_id": "test",
                "row_hash": "abc",
                "source_type": "local",
                "raw_path": str(raw),
                "output_prefix": str(output_prefix),
                "standardise_dir": str(output_prefix.parent),
                "metadata": {},
            }
            stage_record.write_text(json.dumps(record), encoding="utf-8")
            out_record = tmp / "record.json"
            self.run_py(
                "standardise_gwas.py",
                "--stage-record",
                stage_record,
                "--config-root",
                tmp / "empty_manifests",
                "--reference-root",
                FIXTURES / "reference",
                "--output",
                out_record,
            )
            result = json.loads(out_record.read_text(encoding="utf-8"))
            qc = json.loads(Path(result["qc_json"]).read_text(encoding="utf-8"))
            self.assertEqual(qc["config"], "")
            self.assertEqual(qc["source_columns"]["SNP"], "SNP")
            self.assertEqual(qc["source_columns"]["EA"], "A1")
            self.assertEqual(qc["source_columns"]["OA"], "A2")
            self.assertEqual(qc["source_columns"]["EAF"], "Freq1.Hapmap")
            self.assertGreaterEqual(qc["reference_fills"].get("CHR", 0), 2)
            self.assertGreaterEqual(qc["reference_fills"].get("POS38", 0), 2)
            text = subprocess.check_output(["gzip", "-cd", result["standard_vcf"]], text=True)
            metadata = json.loads(Path(result["metadata_json"]).read_text(encoding="utf-8"))
            self.assertEqual(metadata["sex"], "female")
            self.assertIn("##sex=female", text)
            self.assertIn("rs1", text)
            self.assertIn("\t1000\trs1\tG\tA\t", text)

    def test_existing_gwas_vcf_gets_manifest_sex_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = FIXTURES / "tables/existing_gwas.vcf"
            stage_record = tmp / "stage.json"
            output_prefix = tmp / "standardise/test"
            record = {
                "GWAS_location": str(raw),
                "phenotype": "BMI",
                "ancestry": "EUR",
                "sex": "male",
                "author": "Example",
                "year": "2026",
                "PMID": "1",
                "study_id": "example_2026_1",
                "output_id": "test",
                "row_hash": "abc",
                "source_type": "local",
                "raw_path": str(raw),
                "output_prefix": str(output_prefix),
                "standardise_dir": str(output_prefix.parent),
                "metadata": {},
            }
            stage_record.write_text(json.dumps(record), encoding="utf-8")
            out_record = tmp / "record.json"
            self.run_py(
                "standardise_gwas.py",
                "--stage-record",
                stage_record,
                "--config-root",
                tmp / "empty_manifests",
                "--reference-root",
                FIXTURES / "reference",
                "--output",
                out_record,
            )
            result = json.loads(out_record.read_text(encoding="utf-8"))
            text = subprocess.check_output(["gzip", "-cd", result["standard_vcf"]], text=True)
            metadata = json.loads(Path(result["metadata_json"]).read_text(encoding="utf-8"))
            self.assertEqual(metadata["sex"], "male")
            self.assertIn("##sex=male", text)
            self.assertIn("#CHROM", text)

    def test_qc_plot_generates_tiff(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = FIXTURES / "tables/generic.tsv"
            stage_record = tmp / "stage.json"
            output_prefix = tmp / "standardise/test"
            record = {
                "GWAS_location": str(raw),
                "phenotype": "BMI",
                "ancestry": "all",
                "author": "Locke",
                "year": "2015",
                "PMID": "25673413",
                "study_id": "locke_2015_25673413",
                "output_id": "test",
                "row_hash": "abc",
                "source_type": "local",
                "raw_path": str(raw),
                "output_prefix": str(output_prefix),
                "standardise_dir": str(output_prefix.parent),
                "metadata": {},
            }
            stage_record.write_text(json.dumps(record), encoding="utf-8")
            out_record = tmp / "record.json"
            self.run_py(
                "standardise_gwas.py",
                "--stage-record",
                stage_record,
                "--config-root",
                PIPELINE / "config",
                "--reference-root",
                FIXTURES / "reference",
                "--output",
                out_record,
            )
            plot_record = tmp / "plot.json"
            self.run_py("qc_plot.py", "--standardise-record", out_record, "--output", plot_record)
            result = json.loads(plot_record.read_text(encoding="utf-8"))
            tiff = Path(result["plot_tiff"])
            self.assertTrue(tiff.exists())
            self.assertGreater(tiff.stat().st_size, 100)

    def test_qc_plot_fills_raw_positions_from_reference(self):
        sys.path.insert(0, str(BIN))
        from qc_plot import clean_manhattan_df, clean_qq_df, extract_plot_df
        from gwas_pipeline_common import normalize_chromosome

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = tmp / "rsid_only.tsv"
            raw.write_text("SNP\tp\nrs1:A:G\t0.001\nrs2\t0.5\n", encoding="utf-8")
            record = {"reference_root": str(FIXTURES / "reference"), "ancestry": "all"}
            df = extract_plot_df(raw, record)
            self.assertEqual(len(clean_qq_df(df)), 2)
            self.assertEqual(len(clean_manhattan_df(df)), 2)
            self.assertTrue(df["POS"].notna().all())
            self.assertEqual(df.loc[0, "SNP"], "rs1")
            self.assertEqual(df.loc[0, "SNPID"], "rs1:A:G")
            self.assertEqual(normalize_chromosome("chr23"), "X")
            self.assertEqual(normalize_chromosome("24"), "Y")
            self.assertEqual(normalize_chromosome("chrM"), "MT")


if __name__ == "__main__":
    unittest.main()
