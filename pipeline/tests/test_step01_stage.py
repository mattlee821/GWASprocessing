"""Step 1, Stage (README "Pipeline in detail", 1. Stage).

GWAS Catalog sources: the raw author-submitted file first, the harmonised file URL recorded, the harmonised file as
fallback when the raw file cannot be standardised, and `gwas_catalog_file: harmonised` choosing it from the start.
OpenGWAS sources through the authenticated API. URLs downloaded with retries; local files copied and copied again when
they change; downloads checked by checksum or gzip integrity.

Every remote service is mimicked by a server bound to 127.0.0.1; testlib sends any other request to a dead proxy.
"""

from __future__ import annotations

import gzip
import hashlib
import http.server
import json
import os
import shutil
import tempfile
import threading
import unittest
from pathlib import Path
from unittest import mock

from testlib import EXAMPLES, LocalHttpServer, make_reference, run_script, run_standardise, write_table

import download_gwas_catalog
import gwas_pipeline_common as common

ACCESSION, PMID = "GCST008996", "30239722"
RAW_NAME = f"{PMID}-{ACCESSION}-EFO_0004343.tsv.gz"
HARMONISED_NAME = f"{PMID}-{ACCESSION}-EFO_0004343.h.tsv.gz"
TEST_TOKEN = "synthetic-test-token"  # not a real OpenGWAS token


class ScriptedHandler(http.server.BaseHTTPRequestHandler):
    """Answers from a script {path: [(status, body, headers), ...]}, consumed in order (the last entry repeats)."""

    script: dict = {}
    requests: list = []

    def log_message(self, *args):
        pass

    def respond(self):
        body_in = self.rfile.read(int(self.headers.get("Content-Length", 0) or 0)).decode()
        self.requests.append((self.command, self.path, self.headers.get("Authorization", ""), body_in))
        queue = self.script.get(self.path.split("?")[0], [(404, b"missing", {})])
        status, body, headers = queue.pop(0) if len(queue) > 1 else queue[0]
        self.send_response(status)
        for key, value in headers.items():
            self.send_header(key, value)
        if "Content-Length" not in headers:
            self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    do_GET = do_POST = respond


class ScriptedServer:
    def __init__(self, script: dict):
        self.handler = type("Handler", (ScriptedHandler,), {"script": script, "requests": []})
        self.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), self.handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)

    @property
    def url(self) -> str:
        return f"http://127.0.0.1:{self.server.server_address[1]}"

    def __enter__(self) -> "ScriptedServer":
        self.thread.start()
        return self

    def __exit__(self, *exc):
        self.server.shutdown()
        self.server.server_close()


def gwas_catalog_mirror(root: Path, raw: bytes | None = b"SNP\tA1\tA2\nrs1\tA\tG\n", md5_text: str | None = None,
                        harmonised: bytes = b"hm_rsid\thm_pos\nrs1\t100\n") -> Path:
    """A GWAS Catalog FTP study folder: a raw file, harmonised/ with its file, -meta.yaml and md5sum.txt."""
    study = root / download_gwas_catalog.accession_range(ACCESSION) / ACCESSION
    (study / "harmonised").mkdir(parents=True)
    payload = gzip.compress(harmonised)
    (study / "harmonised" / HARMONISED_NAME).write_bytes(payload)
    (study / "harmonised" / f"{HARMONISED_NAME}-meta.yaml").write_text("genome_assembly: GRCh38\n", encoding="utf-8")
    (study / "harmonised" / "md5sum.txt").write_text(md5_text or f"{hashlib.md5(payload).hexdigest()}  {HARMONISED_NAME}\n",
                                                     encoding="utf-8")
    (study / "readme.txt").write_text("readme\n", encoding="utf-8")
    if raw is not None:
        (study / RAW_NAME).write_bytes(raw)
    return study


class Step01StageTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="step01_"))

    def row(self, location, source_type, name="row", **extra) -> Path:
        case = self.tmp / name
        row = {"GWAS_location": str(location), "phenotype": 'BMI "adj"', "ancestry": "EUR", "author": "Locke", "year": "2015",
               "PMID": PMID, "sex": "", "source_type": source_type, "input_build": "", "population": "", "delimiter": "",
               "format": "", "study_yaml": "", "study_id": "locke_2015_25673413", "output_id": "out1", "row_hash": "h1",
               "raw_dir": str(case / "raw"), "standardise_dir": str(case / "std"),
               "output_prefix": str(case / "std" / "standardised"), "settings_hash": "s" * 16}
        row.update(extra)
        case.mkdir(parents=True, exist_ok=True)
        path = case / "row.json"
        path.write_text(json.dumps(row), encoding="utf-8")
        return path

    def stage(self, row_json: Path, *extra, env=None):
        out = row_json.parent / "stage_record.json"
        out.unlink(missing_ok=True)
        proc = run_script("stage_source.py", "--row-json", row_json, "--config-root", self.tmp, "--reference-root",
                          self.tmp / "references", "--output", out, *extra, env=env, check=False)
        return proc, (json.loads(out.read_text(encoding="utf-8")) if out.exists() else None)

    def test_stage(self):
        # --- GWAS Catalog ---------------------------------------------------------------------------------------
        mirror = self.tmp / "ftp"
        raw_payload = gzip.compress(b"SNP\tA1\tA2\nrs1\tA\tG\n")
        gwas_catalog_mirror(mirror, raw=raw_payload)
        yaml_path = self.tmp / "harmonised.yaml"
        yaml_path.write_text("gwas_catalog_file: harmonised\n", encoding="utf-8")
        with LocalHttpServer(mirror) as server:
            env = {**os.environ, "GWAS_CATALOG_FTP_ROOT": server.url}
            raw_proc, raw_record = self.stage(self.row(ACCESSION, "gwas_catalog", name="catalog_raw"), env=env)
            hm_proc, hm_record = self.stage(self.row(f"https://www.ebi.ac.uk/gwas/studies/{ACCESSION}", "gwas_catalog",
                                                     name="catalog_hm", study_yaml=str(yaml_path)), env=env)

        with self.subTest("GWAS Catalog sources are downloaded as the raw author-submitted file"):
            self.assertEqual(raw_proc.returncode, 0, raw_proc.stderr)
            meta = raw_record["metadata"]
            self.assertEqual((meta["source_type"], meta["file_kind"]), ("gwas_catalog", "raw"))
            self.assertEqual(Path(raw_record["raw_path"]).read_bytes(), raw_payload)
            self.assertTrue(meta["summary_url"].endswith(RAW_NAME))

        with self.subTest("The URL of the Catalog's harmonised file is recorded"):
            self.assertTrue(raw_record["metadata"]["harmonised_url"].endswith(f"/harmonised/{HARMONISED_NAME}"))

        with self.subTest("gwas_catalog_file: harmonised in a study YAML chooses the harmonised file from the start"):
            self.assertEqual(hm_proc.returncode, 0, hm_proc.stderr)
            meta = hm_record["metadata"]
            self.assertEqual(meta["file_kind"], "harmonised")
            self.assertTrue(hm_record["raw_path"].endswith(HARMONISED_NAME))
            self.assertTrue(meta["raw_url"].endswith(RAW_NAME))
            self.assertEqual(meta["genome_assembly"], "GRCh38")

        with self.subTest("Downloads are checked by checksum: a file matching its md5sum sidecar passes"):
            self.assertEqual((hm_record["metadata"]["md5"]["checked"], hm_record["metadata"]["md5"]["ok"]), (True, True))
            # a sidecar line for an AppleDouble file of the same name is not taken for the file's own line
            downloaded = self.tmp / "file.tbl.gz"
            downloaded.write_text("actual", encoding="utf-8")
            sidecar = self.tmp / "md5sum.txt"
            sidecar.write_text("0" * 32 + " ._file.tbl.gz\n" + f"{download_gwas_catalog.md5(downloaded)} file.tbl.gz\n",
                               encoding="utf-8")
            self.assertTrue(download_gwas_catalog.verify_md5(downloaded, [str(sidecar)])["ok"])

        with self.subTest("Downloads are checked by gzip integrity when there is no checksum"):
            self.assertEqual(raw_record["metadata"]["md5"]["checked"], False)
            self.assertEqual((raw_record["metadata"]["gzip"]["checked"], raw_record["metadata"]["gzip"]["ok"]), (True, True))

        with self.subTest("A download that does not match its checksum is rejected"):
            bad = self.tmp / "ftp_bad_md5"
            gwas_catalog_mirror(bad, raw=None, md5_text="0" * 32 + f"  {HARMONISED_NAME}\n")
            with LocalHttpServer(bad) as server, mock.patch.dict(os.environ, {"GWAS_CATALOG_FTP_ROOT": server.url}):
                with self.assertRaisesRegex(ValueError, "MD5"):
                    download_gwas_catalog.download_gwas_catalog(ACCESSION, self.tmp / "bad_md5", PMID)

        with self.subTest("A download that fails the gzip integrity check is rejected"):
            corrupt = self.tmp / "ftp_corrupt"
            gwas_catalog_mirror(corrupt, raw=gzip.compress(b"SNP\tP\n" * 1000)[:-20])  # truncated gzip stream
            with LocalHttpServer(corrupt) as server, mock.patch.dict(os.environ, {"GWAS_CATALOG_FTP_ROOT": server.url}):
                with self.assertRaisesRegex(ValueError, "Gzip integrity"):
                    download_gwas_catalog.download_gwas_catalog(ACCESSION, self.tmp / "corrupt", PMID)

        with self.subTest("GWAS Catalog study folders are found in the FTP accession ranges, legacy file names included"):
            self.assertEqual(download_gwas_catalog.accession_range("GCST002783"), "GCST002001-GCST003000")
            self.assertEqual(download_gwas_catalog.accession_range("GCST90000123"), "GCST90000001-GCST90001000")
            legacy = self.tmp / "ftp_legacy"
            study = legacy / "GCST90000001-GCST90001000" / "GCST90000123"
            study.mkdir(parents=True)
            payload = gzip.compress(b"SNP\tP\nrs1\t0.1\n")
            (study / "fat-distn.giant.ukbb.meta-analysis.whr.combined.tbl.gz").write_bytes(payload)
            (study / "readme.txt").write_text("readme\n", encoding="utf-8")
            with LocalHttpServer(legacy) as server, mock.patch.dict(os.environ, {"GWAS_CATALOG_FTP_ROOT": server.url}):
                meta = download_gwas_catalog.download_gwas_catalog("GCST90000123", self.tmp / "legacy", "")
            self.assertEqual((meta["file_kind"], Path(meta["raw_path"]).read_bytes()), ("raw", payload))

        # Fallback to the harmonised file: the row is standardised from the raw file and, when that fails, from the
        # harmonised file named in the stage metadata. 25 common SNPs place the harmonised file in GRCh38.
        lookup = [("1", f"rs{100 + i}", 10_000 + 1000 * i, "A", "G", 0.3) for i in range(25)]
        reference = make_reference(self.tmp / "fallback_reference", lookup_rows=lookup)
        served = self.tmp / "served_harmonised"
        served.mkdir()
        write_table(served / "x.h.tsv", ["chromosome", "base_pair_location", "rsid", "effect_allele", "other_allele",
                                         "effect_allele_frequency", "beta", "standard_error", "p_value"],
                    [[1, pos19 + 100, rs, alt, ref, 0.3, 0.1, 0.02, 1e-6] for _, rs, pos19, alt, ref, _ in lookup])
        unmappable = write_table(self.tmp / "raw_unmappable.tsv", ["SNP", "A1", "A2", "P"], [["rs100", "A", "G", 0.1]])
        no_build = write_table(self.tmp / "raw_no_build.tsv", ["CHR", "POS", "SNP", "EA", "OA", "BETA", "SE", "P"],
                               [[1, 500_001 + 7 * i, f"rs{900 + i}", "A", "G", 0.1, 0.02, 1e-6] for i in range(25)])
        with LocalHttpServer(served) as server:
            metadata = {"source_type": "gwas_catalog", "file_kind": "raw", "harmonised_url": f"{server.url}/x.h.tsv"}
            columns = run_standardise(self.tmp, unmappable, reference, metadata=metadata, name="fallback_columns")
            build = run_standardise(self.tmp, no_build, reference, metadata=metadata, name="fallback_build")

        with self.subTest("Raw file with columns that cannot be mapped: the harmonised file is used and the fallback recorded"):
            self.assertEqual(columns.returncode, 0, columns.stderr)
            self.assertEqual(columns.qc["gwas_catalog_fallback"]["harmonised_file"], "x.h.tsv")
            self.assertEqual(columns.qc["gwas_catalog_fallback"]["raw_file"], "raw_unmappable.tsv")
            self.assertEqual(len(columns.rows), 25)
            self.assertTrue(any("harmonised file x.h.tsv" in w for w in columns.qc["warnings"]))

        with self.subTest("Raw file whose build cannot be detected: the harmonised file is used and the fallback recorded"):
            self.assertEqual(build.returncode, 0, build.stderr)
            self.assertEqual((build.qc["gwas_catalog_fallback"]["harmonised_file"], build.qc["gwas_catalog_fallback"]["exit_code"]),
                             ("x.h.tsv", 78))
            self.assertEqual(build.qc["build_detection"]["build"], "hg38")
            self.assertEqual(len(build.rows), 25)

        # --- OpenGWAS -------------------------------------------------------------------------------------------
        dataset = "ebi-a-GCST90000123"
        vcf = gzip.compress(b"##fileformat=VCFv4.2\n")
        script: dict = {}
        with ScriptedServer(script) as server:
            files = {dataset: [f"{server.url}/files/{dataset}.vcf.gz", f"{server.url}/files/{dataset}.vcf.gz.tbi"]}
            script.update({
                "/api/gwasinfo": [(200, json.dumps([{"id": dataset, "build": "HG19/GRCh37", "sample_size": 1000}]).encode(), {})],
                "/api/gwasinfo/files": [(200, json.dumps(files).encode(), {})],
                f"/files/{dataset}.vcf.gz": [(200, vcf, {})],
                f"/files/{dataset}.vcf.gz.tbi": [(200, b"tbi", {})],
            })
            env = {**os.environ, "OPENGWAS_API_ROOT": f"{server.url}/api", "OPENGWAS_JWT": TEST_TOKEN}
            og_proc, og_record = self.stage(self.row(f"https://opengwas.io/datasets/{dataset}", "opengwas", name="opengwas"), env=env)
            requests = list(server.handler.requests)

        with self.subTest("OpenGWAS sources use the authenticated API"):
            self.assertEqual(og_proc.returncode, 0, og_proc.stderr)
            self.assertEqual(Path(og_record["raw_path"]).read_bytes(), vcf)
            self.assertTrue(og_record["metadata"]["index_path"].endswith(".vcf.gz.tbi"))
            self.assertEqual(og_record["metadata"]["dataset_id"], dataset)  # ID case preserved
            info_call = next(r for r in requests if r[1] == "/api/gwasinfo")
            self.assertEqual((info_call[0], info_call[2], json.loads(info_call[3])), ("POST", f"Bearer {TEST_TOKEN}", {"id": [dataset]}))
            self.assertTrue(all(r[2] == f"Bearer {TEST_TOKEN}" for r in requests if r[1].startswith("/api/")))

        with self.subTest("A missing OpenGWAS token is a configuration error (exit 78)"):
            env = {**os.environ, "OPENGWAS_API_ROOT": "http://127.0.0.1:9/api", "OPENGWAS_JWT": ""}
            proc, _ = self.stage(self.row(dataset, "opengwas", name="no_token"), env=env)
            self.assertEqual(proc.returncode, 78, proc.stderr)
            self.assertIn("OPENGWAS_JWT", proc.stderr)

        for status in (401, 403):
            with self.subTest(f"A rejected OpenGWAS token (HTTP {status}) is a configuration error (exit 78)"):
                with ScriptedServer({"/api/gwasinfo": [(status, b'{"message": "expired"}', {})]}) as server:
                    env = {**os.environ, "OPENGWAS_API_ROOT": f"{server.url}/api", "OPENGWAS_JWT": TEST_TOKEN}
                    proc, _ = self.stage(self.row(dataset, "opengwas", name=f"rejected_{status}"), env=env)
                self.assertEqual(proc.returncode, 78, proc.stderr)
                self.assertIn("expired", proc.stderr)
                self.assertNotIn(TEST_TOKEN, proc.stderr + proc.stdout)  # the token is never echoed

        # --- URLs ---------------------------------------------------------------------------------------------
        with self.subTest("URLs are downloaded"):
            served = self.tmp / "served_url"
            served.mkdir()
            shutil.copy(EXAMPLES / "tables/generic.tsv", served / "generic.tsv")
            with LocalHttpServer(served) as server:
                proc, record = self.stage(self.row(f"{server.url}/generic.tsv", "direct_url", name="url"))
            self.assertEqual(proc.returncode, 0, proc.stderr)
            self.assertEqual(Path(record["raw_path"]).read_bytes(), (EXAMPLES / "tables/generic.tsv").read_bytes())

        with self.subTest("URLs are downloaded with retries"):
            with ScriptedServer({"/f.tsv": [(500, b"busy", {}), (200, b"SNP\tP\n", {})]}) as server:
                path = common.download_url(f"{server.url}/f.tsv", self.tmp / "retry", retries=2, backoff=0)
                self.assertEqual(path.read_bytes(), b"SNP\tP\n")
                self.assertEqual(len(server.handler.requests), 2)

        with self.subTest("A truncated transfer is rejected and leaves no file"):
            target = self.tmp / "truncated"
            with ScriptedServer({"/f.tsv": [(200, b"short", {"Content-Length": "50"})]}) as server:
                with self.assertRaisesRegex(OSError, "Incomplete download"):
                    common.download_url(f"{server.url}/f.tsv", target, retries=1, backoff=0, timeout=5)
            self.assertEqual(list(target.iterdir()), [])

        with self.subTest("A complete earlier download is not fetched again"):
            target = self.tmp / "cached"
            target.mkdir()
            (target / "f.tsv").write_text("cached", encoding="utf-8")
            with ScriptedServer({"/f.tsv": [(200, b"new", {})]}) as server:
                path = common.download_url(f"{server.url}/f.tsv", target)
                self.assertEqual((path.read_text(encoding="utf-8"), len(server.handler.requests)), ("cached", 0))

        # --- Local files --------------------------------------------------------------------------------------
        with self.subTest("Local files are copied, and copied again when they change"):
            source = self.tmp / "gwas.tsv"
            source.write_text("SNP\tP\nrs1\t0.1\n", encoding="utf-8")
            row = self.row(source, "local", name="local")
            _, record = self.stage(row)
            self.assertNotEqual(Path(record["raw_path"]), source)
            self.assertEqual(Path(record["raw_path"]).read_text(encoding="utf-8"), "SNP\tP\nrs1\t0.1\n")
            source.write_text("SNP\tP\nrs1\t0.1\nrs2\t0.2\n", encoding="utf-8")  # a different size
            _, record = self.stage(row)
            self.assertEqual(Path(record["raw_path"]).read_text(encoding="utf-8"), "SNP\tP\nrs1\t0.1\nrs2\t0.2\n")
            source.write_text("SNP\tP\nrs1\t0.1\nrs2\t0.3\n", encoding="utf-8")  # the same size, a later time
            later = Path(record["raw_path"]).stat().st_mtime + 10
            os.utime(source, (later, later))
            _, record = self.stage(row)
            self.assertEqual(Path(record["raw_path"]).read_text(encoding="utf-8"), "SNP\tP\nrs1\t0.1\nrs2\t0.3\n")

        with self.subTest("A missing local file fails the row"):
            proc, _ = self.stage(self.row(self.tmp / "absent.tsv", "local", name="absent"))
            self.assertNotEqual(proc.returncode, 0)
            self.assertIn("not found", proc.stderr)

        with self.subTest("An unsupported source type is a configuration error (exit 78)"):
            proc, _ = self.stage(self.row(self.tmp / "gwas.tsv", "unsupported", name="unsupported"))
            self.assertEqual(proc.returncode, 78, proc.stderr)

        # --- Stage record -------------------------------------------------------------------------------------
        with self.subTest("The stage record carries the manifest fields and the settings hash"):
            _, record = self.stage(self.row(source, "local", name="fields"))
            self.assertEqual((record["phenotype"], record["settings_hash"], record["study_id"], record["source_type"]),
                             ('BMI "adj"', "s" * 16, "locke_2015_25673413", "local"))

        with self.subTest("A dry run records the staged path without staging"):
            _, record = self.stage(self.row(self.tmp / "absent.tsv", "local", name="dry"), "--dry-run", "true")
            self.assertEqual(Path(record["raw_path"]).parent, self.tmp / "dry" / "raw")
            self.assertFalse(Path(record["raw_path"]).exists())


if __name__ == "__main__":
    unittest.main()
