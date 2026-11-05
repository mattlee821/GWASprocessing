#!/usr/bin/env python3
"""Download OpenGWAS GWAS-VCF files via the authenticated OpenGWAS API.

A valid JWT in the OPENGWAS_JWT environment variable is required. If the token
is missing or rejected (expired/invalid) the download fails with an explicit
error telling the user to provide a valid token; there is no unauthenticated
fallback.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import requests

from gwas_pipeline_common import download_url, ensure_dir, opengwas_dataset_id, write_json

DEFAULT_API_ROOT = "https://api.opengwas.io/api"
# sysexits.h EX_CONFIG: a configuration/authentication error the user must fix.
# The Nextflow errorStrategy terminates the whole run on this code (rather than
# ignoring the row) so the prompt to supply a valid token is not missed.
EXIT_CONFIG = 78
TOKEN_HELP = (
    "Provide a valid token in server.env as OPENGWAS_JWT=<token>. "
    "Generate or renew one at https://api.opengwas.io/profile/ ."
)


def config_error(message: str) -> "SystemExit":
    print(message, file=sys.stderr)
    return SystemExit(EXIT_CONFIG)


def api_root() -> str:
    """OpenGWAS API root; OPENGWAS_API_ROOT points it at a mirror (or a test server)."""
    return os.environ.get("OPENGWAS_API_ROOT", DEFAULT_API_ROOT).rstrip("/")


def opengwas_id(location: str) -> str:
    """Dataset ID such as ieu-b-5129, ukb-d-30780_irnt or ebi-a-GCST90000123, with its case preserved."""
    dataset = opengwas_dataset_id(location)
    if not dataset:
        raise config_error(f"Could not find an OpenGWAS dataset ID in: {location}")
    return dataset


def collect_urls(payload) -> list[str]:
    urls: list[str] = []
    stack = [payload]
    while stack:
        item = stack.pop()
        if isinstance(item, str) and item.startswith(("http://", "https://")):
            urls.append(item)
        elif isinstance(item, dict):
            stack.extend(item.values())
        elif isinstance(item, list):
            stack.extend(item)
    return sorted(set(urls))


def opengwas_jwt() -> str:
    token = os.environ.get("OPENGWAS_JWT", "").strip()
    if not token:
        raise config_error(
            "OpenGWAS authentication is required but OPENGWAS_JWT is not set. " + TOKEN_HELP
        )
    return token


def opengwas_api(method: str, path: str, token: str, **kwargs) -> dict | list:
    headers = {"Authorization": f"Bearer {token}", "Accept": "application/json"}
    try:
        response = requests.request(method, f"{api_root()}{path}", headers=headers, timeout=120, **kwargs)
    except requests.RequestException as exc:
        raise SystemExit(f"OpenGWAS API request failed ({method} {path}): {exc}") from exc
    if response.status_code in (401, 403):
        try:
            detail = str(response.json().get("message", response.text))
        except ValueError:
            detail = response.text
        raise config_error(
            f"OpenGWAS rejected the JWT in OPENGWAS_JWT (HTTP {response.status_code}): "
            f"{detail[:200]} The token is missing, expired, or invalid. " + TOKEN_HELP
        )
    response.raise_for_status()
    return response.json()


def download_opengwas(location: str, out_dir: str | Path) -> dict:
    out_dir = ensure_dir(out_dir)
    dataset_id = opengwas_id(location)
    token = opengwas_jwt()

    info = opengwas_api("POST", "/gwasinfo", token, json={"id": [dataset_id]})
    files = opengwas_api("POST", "/gwasinfo/files", token, json={"id": dataset_id, "type": "vcf"})

    downloaded: list[str] = []
    for url in collect_urls(files):
        name = Path(url.split("?", 1)[0]).name
        if any(token_str in name for token_str in [".vcf.gz", ".tbi"]):
            downloaded.append(str(download_url(url, out_dir, name)))

    vcf = next((p for p in downloaded if p.endswith(".vcf.gz")), "")
    index = next((p for p in downloaded if p.endswith(".tbi")), "")
    if not vcf:
        raise SystemExit(
            f"OpenGWAS returned no downloadable VCF for {dataset_id}: {json.dumps(files)[:300]}"
        )
    metadata = {
        "source_type": "opengwas",
        "dataset_id": dataset_id,
        "gwasinfo": info,
        "gwasinfo_files": files,
        "raw_path": vcf,
        "index_path": index,
        "downloaded": downloaded,
    }
    write_json(out_dir / f"{dataset_id}.opengwas_metadata.json", metadata)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas-location", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()
    write_json(args.metadata, download_opengwas(args.gwas_location, args.out_dir))


if __name__ == "__main__":
    main()
