#!/usr/bin/env python3
"""Download GWAS Catalog summary statistics from a study URL or GCST accession."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import urllib.parse
from pathlib import Path

from gwas_pipeline_common import download_url, ensure_dir, url_links, write_json


FTP_ROOT = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"


def accession_from_location(location: str) -> str:
    match = re.search(r"GCST\d+", location, flags=re.I)
    if not match:
        raise ValueError(f"Could not find GWAS Catalog accession in: {location}")
    return match.group(0).upper()


def accession_range(accession: str) -> str:
    number = int(re.search(r"\d+", accession).group(0))
    width = max(6, len(str(number)))
    start = ((number - 1) // 1000) * 1000 + 1
    end = start + 999
    return f"GCST{start:0{width}d}-GCST{end:0{width}d}"


def list_absolute(url: str) -> list[str]:
    links = []
    for link in url_links(url.rstrip("/") + "/"):
        if link.startswith("?") or link.startswith("../"):
            continue
        links.append(urllib.parse.urljoin(url.rstrip("/") + "/", link))
    return links


def choose_summary_file(study_url: str, accession: str, pmid: str = "") -> str:
    top_links = list_absolute(study_url)
    harmonised = [u for u in top_links if u.rstrip("/").endswith("/harmonised")]
    candidates: list[str] = []
    search_dirs = harmonised or [study_url]
    for directory in search_dirs:
        for link in list_absolute(directory):
            name = Path(urllib.parse.urlparse(link).path).name
            if not name.endswith((".h.tsv.gz", ".tsv.gz", ".txt.gz", ".tbl.gz")):
                continue
            candidates.append(link)
    if not candidates:
        raise FileNotFoundError(f"No summary statistics file found in {study_url}")
    if pmid:
        pmid_matches = [u for u in candidates if Path(urllib.parse.urlparse(u).path).name.startswith(f"{pmid}-{accession}")]
        if pmid_matches:
            pmid_harmonised = [u for u in pmid_matches if Path(urllib.parse.urlparse(u).path).name.endswith(".h.tsv.gz")]
            return sorted(pmid_harmonised or pmid_matches)[0]
    harmonised_matches = [u for u in candidates if ".h.tsv.gz" in Path(urllib.parse.urlparse(u).path).name]
    return sorted(harmonised_matches or candidates)[0]


def download_sidecars(study_url: str, out_dir: Path) -> list[str]:
    sidecars = []
    directories = [study_url]
    directories.extend([u for u in list_absolute(study_url) if u.rstrip("/").endswith("/harmonised")])
    for directory in directories:
        for link in list_absolute(directory):
            name = Path(urllib.parse.urlparse(link).path).name
            lower = name.lower()
            if lower in {"readme.txt", "md5sum", "md5sums", "md5sum.txt"} or lower.endswith((".md5", ".md5sum")):
                try:
                    sidecars.append(str(download_url(link, out_dir, name)))
                except Exception:
                    pass
    return sidecars


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_md5(downloaded: Path, sidecars: list[str]) -> dict[str, str | bool]:
    result: dict[str, str | bool] = {"checked": False, "ok": False}
    for sidecar in sidecars:
        text = Path(sidecar).read_text(errors="replace")
        expected = None
        for line in text.splitlines():
            parts = line.split()
            if len(parts) < 2:
                continue
            if Path(parts[-1]).name != downloaded.name:
                continue
            match = re.search(r"\b[a-fA-F0-9]{32}\b", line)
            if match:
                expected = match.group(0).lower()
                break
        if expected:
            observed = md5(downloaded)
            result.update({"checked": True, "expected": expected, "observed": observed, "ok": expected == observed})
            return result
    return result


def verify_gzip(path: Path) -> dict[str, str | bool]:
    if not str(path).endswith(".gz"):
        return {"checked": False, "ok": False}
    try:
        with gzip.open(path, "rb") as handle:
            for _ in iter(lambda: handle.read(1024 * 1024), b""):
                pass
    except Exception as exc:
        return {"checked": True, "ok": False, "error": str(exc)}
    return {"checked": True, "ok": True}


def download_gwas_catalog(location: str, out_dir: str | Path, pmid: str = "") -> dict:
    out_dir = ensure_dir(out_dir)
    accession = accession_from_location(location)
    study_url = f"{FTP_ROOT}/{accession_range(accession)}/{accession}/"
    summary_url = choose_summary_file(study_url, accession, pmid)
    summary = download_url(summary_url, out_dir)
    sidecars = download_sidecars(study_url, out_dir)
    md5_result = verify_md5(summary, sidecars)
    if md5_result.get("checked") and not md5_result.get("ok"):
        summary.unlink(missing_ok=True)
        summary = download_url(summary_url, out_dir)
        md5_result = verify_md5(summary, sidecars)
        if md5_result.get("checked") and not md5_result.get("ok"):
            raise ValueError(f"MD5 verification failed for {summary}: {md5_result}")
    gzip_result = {"checked": False, "ok": False}
    if not md5_result.get("checked"):
        gzip_result = verify_gzip(summary)
        if gzip_result.get("checked") and not gzip_result.get("ok"):
            summary.unlink(missing_ok=True)
            summary = download_url(summary_url, out_dir)
            gzip_result = verify_gzip(summary)
            if gzip_result.get("checked") and not gzip_result.get("ok"):
                raise ValueError(f"Gzip integrity check failed for {summary}: {gzip_result}")
    metadata = {
        "source_type": "gwas_catalog",
        "accession": accession,
        "study_url": study_url,
        "summary_url": summary_url,
        "raw_path": str(summary),
        "sidecars": sidecars,
        "md5": md5_result,
        "gzip": gzip_result,
    }
    write_json(out_dir / f"{accession}.gwas_catalog_metadata.json", metadata)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas-location", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--pmid", default="")
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()
    write_json(args.metadata, download_gwas_catalog(args.gwas_location, args.out_dir, args.pmid))


if __name__ == "__main__":
    main()
