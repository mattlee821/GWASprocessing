"""Tests for the download layer: download_url, the GWAS Catalog FTP client and the OpenGWAS API client.

Real HTTP is used throughout, against servers bound to 127.0.0.1 that mimic the remote services.
"""

import gzip
import hashlib
import http.server
import json
import os
import tempfile
import threading
import unittest
from pathlib import Path

from testlib import LocalHttpServer, run_script

import download_gwas_catalog
import gwas_pipeline_common as common


class ScriptedHandler(http.server.BaseHTTPRequestHandler):
    """Serves responses from a class-level script: {path: [(status, body bytes, headers), ...]} consumed in order."""

    script: dict = {}
    requests: list = []

    def log_message(self, *args):
        pass

    def respond(self):
        self.requests.append((self.command, self.path, self.headers.get("Authorization", ""),
                              self.rfile.read(int(self.headers.get("Content-Length", 0) or 0)).decode()))
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
    def __init__(self, script):
        handler = type("Handler", (ScriptedHandler,), {"script": script, "requests": []})
        self.handler = handler
        self.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)

    @property
    def url(self):
        return f"http://127.0.0.1:{self.server.server_address[1]}"

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *exc):
        self.server.shutdown()
        self.server.server_close()


class DownloadUrlTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_transient_failure_is_retried(self):
        with ScriptedServer({"/f.tsv": [(500, b"busy", {}), (200, b"SNP\tP\n", {})]}) as server:
            path = common.download_url(f"{server.url}/f.tsv", self.tmp, retries=2, backoff=0)
            self.assertEqual(path.read_bytes(), b"SNP\tP\n")
            self.assertEqual(len(server.handler.requests), 2)

    def test_truncated_transfer_is_rejected_and_leaves_no_file(self):
        with ScriptedServer({"/f.tsv": [(200, b"short", {"Content-Length": "50"})]}) as server:
            with self.assertRaisesRegex(OSError, "Incomplete download"):
                common.download_url(f"{server.url}/f.tsv", self.tmp, retries=1, backoff=0, timeout=5)
        self.assertEqual(sorted(p.name for p in self.tmp.iterdir()), [])

    def test_existing_complete_file_is_not_downloaded_again(self):
        (self.tmp / "f.tsv").write_text("cached", encoding="utf-8")
        with ScriptedServer({"/f.tsv": [(200, b"new", {})]}) as server:
            path = common.download_url(f"{server.url}/f.tsv", self.tmp)
            self.assertEqual((path.read_text(encoding="utf-8"), len(server.handler.requests)), ("cached", 0))


class GwasCatalogTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def test_gwas_catalog_selects_legacy_tbl_gz(self):
        original = download_gwas_catalog.list_absolute

        def fake_list_absolute(url):
            if url.endswith("/GCST008996/"):
                return [url.rsplit("/", 2)[0] + "/", url + "fat-distn.giant.ukbb.meta-analysis.whr.combined.tbl.gz",
                        url + "md5sum.txt", url + "readme.txt"]
            return []

        try:
            download_gwas_catalog.list_absolute = fake_list_absolute
            selected = download_gwas_catalog.choose_summary_file("https://example.org/GCST008996/", "GCST008996", "30239722")
            self.assertTrue(selected.endswith(".tbl.gz"))
        finally:
            download_gwas_catalog.list_absolute = original

    def test_gwas_catalog_md5_uses_exact_basename(self):
        downloaded = self.tmp / "file.tbl.gz"
        downloaded.write_text("actual", encoding="utf-8")
        sidecar = self.tmp / "md5sum.txt"
        sidecar.write_text("00000000000000000000000000000000 ._file.tbl.gz\n"
                           f"{download_gwas_catalog.md5(downloaded)} file.tbl.gz\n", encoding="utf-8")
        self.assertTrue(download_gwas_catalog.verify_md5(downloaded, [str(sidecar)])["ok"])

    def test_accession_ranges_follow_the_ftp_layout(self):
        self.assertEqual(download_gwas_catalog.accession_range("GCST002783"), "GCST002001-GCST003000")
        self.assertEqual(download_gwas_catalog.accession_range("GCST90000123"), "GCST90000001-GCST90001000")

    def serve_catalog(self, md5_text=None):
        root = self.tmp / "ftp"
        study = root / "GCST008001-GCST009000" / "GCST008996"
        (study / "harmonised").mkdir(parents=True)
        payload = gzip.compress(b"hm_rsid\thm_pos\nrs1\t100\n")
        name = "30239722-GCST008996-EFO_0004343.h.tsv.gz"
        (study / "harmonised" / name).write_bytes(payload)
        (study / "harmonised" / f"{name}-meta.yaml").write_text("genome_assembly: GRCh38\n", encoding="utf-8")
        (study / "harmonised" / "md5sum.txt").write_text(md5_text or f"{hashlib.md5(payload).hexdigest()}  {name}\n", encoding="utf-8")
        (study / "readme.txt").write_text("readme\n", encoding="utf-8")
        return root, payload

    def test_harmonised_file_is_downloaded_verified_and_its_assembly_recorded(self):
        root, payload = self.serve_catalog()
        with LocalHttpServer(root) as server:
            os.environ["GWAS_CATALOG_FTP_ROOT"] = server.url
            try:
                meta = download_gwas_catalog.download_gwas_catalog("https://www.ebi.ac.uk/gwas/studies/GCST008996", self.tmp / "raw", "30239722")
            finally:
                del os.environ["GWAS_CATALOG_FTP_ROOT"]
        self.assertTrue(meta["summary_url"].endswith(".h.tsv.gz"))
        self.assertEqual(Path(meta["raw_path"]).read_bytes(), payload)
        self.assertEqual((meta["md5"]["ok"], meta["genome_assembly"]), (True, "GRCh38"))

    def test_checksum_mismatch_is_an_error(self):
        root, _ = self.serve_catalog(md5_text="0" * 32 + "  30239722-GCST008996-EFO_0004343.h.tsv.gz\n")
        with LocalHttpServer(root) as server:
            os.environ["GWAS_CATALOG_FTP_ROOT"] = server.url
            try:
                with self.assertRaisesRegex(ValueError, "MD5"):
                    download_gwas_catalog.download_gwas_catalog("GCST008996", self.tmp / "raw", "30239722")
            finally:
                del os.environ["GWAS_CATALOG_FTP_ROOT"]


class OpenGwasTests(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def run_download(self, server_url, token="valid-token", dataset="ebi-a-GCST90000123"):
        env = {**os.environ, "OPENGWAS_API_ROOT": f"{server_url}/api", "OPENGWAS_JWT": token}
        return run_script("download_opengwas.py", "--gwas-location", f"https://opengwas.io/datasets/{dataset}",
                          "--out-dir", self.tmp / "raw", "--metadata", self.tmp / "meta.json", env=env, check=False)

    def test_dataset_is_downloaded_with_the_bearer_token_and_id_case_preserved(self):
        vcf = gzip.compress(b"##fileformat=VCFv4.2\n")
        script = {}
        with ScriptedServer(script) as server:
            files = {"ebi-a-GCST90000123": [f"{server.url}/files/ebi-a-GCST90000123.vcf.gz", f"{server.url}/files/ebi-a-GCST90000123.vcf.gz.tbi"]}
            script.update({
                "/api/gwasinfo": [(200, json.dumps([{"id": "ebi-a-GCST90000123", "build": "HG19/GRCh37", "sample_size": 1000}]).encode(), {})],
                "/api/gwasinfo/files": [(200, json.dumps(files).encode(), {})],
                "/files/ebi-a-GCST90000123.vcf.gz": [(200, vcf, {})],
                "/files/ebi-a-GCST90000123.vcf.gz.tbi": [(200, b"tbi", {})],
            })
            proc = self.run_download(server.url)
            requests = list(server.handler.requests)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        meta = json.loads((self.tmp / "meta.json").read_text(encoding="utf-8"))
        self.assertEqual(Path(meta["raw_path"]).read_bytes(), vcf)
        info_call = next(r for r in requests if r[1] == "/api/gwasinfo")
        self.assertEqual((info_call[0], info_call[2], json.loads(info_call[3])), ("POST", "Bearer valid-token", {"id": ["ebi-a-GCST90000123"]}))

    def test_missing_token_is_a_configuration_error(self):
        proc = self.run_download("http://127.0.0.1:9", token="")
        self.assertEqual(proc.returncode, 78)
        self.assertIn("OPENGWAS_JWT", proc.stderr)

    def test_rejected_token_is_a_configuration_error(self):
        with ScriptedServer({"/api/gwasinfo": [(401, b'{"message": "expired"}', {})]}) as server:
            proc = self.run_download(server.url)
        self.assertEqual(proc.returncode, 78)
        self.assertIn("expired", proc.stderr)


if __name__ == "__main__":
    unittest.main()
