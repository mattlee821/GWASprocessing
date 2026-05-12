#!/usr/bin/env python3
"""Download OpenGWAS GWAS-VCF files using ieugwaspy gwasinfo_files()."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import requests

from gwas_pipeline_common import download_url, ensure_dir, write_json

MAIN_FIELDS_ORDER = [
    "trait",
    "build",
    "category",
    "subcategory",
    "population",
    "sex",
    "author",
    "year",
    "ontology",
    "unit",
    "sample_size",
    "consortium",
    "mr",
    "priority",
]


def opengwas_id(location: str) -> str:
    match = re.search(r"([a-z]+-[a-z]-\d+)", location, flags=re.I)
    if not match:
        raise ValueError(f"Could not find OpenGWAS dataset ID in: {location}")
    return match.group(1).lower()


def first_url(value):
    if isinstance(value, str) and value.startswith(("http://", "https://")):
        return value
    if isinstance(value, dict):
        for item in value.values():
            result = first_url(item)
            if result:
                return result
    if isinstance(value, list):
        for item in value:
            result = first_url(item)
            if result:
                return result
    return None


def api_error(payload) -> str:
    if isinstance(payload, dict):
        message = str(payload.get("message", ""))
        if message.upper().startswith("ERROR"):
            return message
    return ""


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


def trial_files(dataset_id: str, file_type: str):
    response = requests.post(
        "https://api.opengwas.io/api/gwasinfo/files/trial",
        json={"id": dataset_id, "type": file_type},
        timeout=60,
    )
    response.raise_for_status()
    return response.json()


def public_gwasinfo(dataset_id: str) -> dict:
    response = requests.get("https://opengwas.io/data/gwasinfo/gwasinfo.json", timeout=120)
    response.raise_for_status()
    payload = response.json()
    compressed = payload.get("datasets_compressed", {}).get(dataset_id)
    if not compressed:
        return {}
    main_values, additional = compressed
    coding = payload.get("majority_fields_and_coding", {})
    info = {}
    for index, field in enumerate(MAIN_FIELDS_ORDER):
        value = main_values[index] if index < len(main_values) else None
        values = coding.get(field)
        if isinstance(values, list) and isinstance(value, int) and 0 <= value < len(values):
            value = values[value]
        info[field] = value
    if isinstance(additional, dict):
        info.update(additional)
    info["id"] = dataset_id
    return info


def download_opengwas(location: str, out_dir: str | Path) -> dict:
    try:
        from ieugwaspy import query
    except ImportError as exc:
        raise SystemExit(
            "ieugwaspy is required for OpenGWAS downloads. Install it or provide a local/direct GWAS_location."
        ) from exc

    out_dir = ensure_dir(out_dir)
    dataset_id = opengwas_id(location)
    try:
        info = query.gwasinfo([dataset_id])
    except TypeError:
        info = query.gwasinfo(dataset_id)
    if api_error(info):
        info = public_gwasinfo(dataset_id)

    try:
        files = query.gwasinfo_files([dataset_id])
    except TypeError:
        files = query.gwasinfo_files(dataset_id)
    files_error = api_error(files)
    if files_error:
        files = {
            "source": "opengwas_trial_endpoint",
            "api_error": files_error,
            "vcf": trial_files(dataset_id, "vcf"),
            "report": trial_files(dataset_id, "report"),
        }

    downloaded: list[str] = []
    for payload in [files]:
        urls = collect_urls(payload)
        for url in sorted(set(urls)):
            name = Path(url.split("?", 1)[0]).name
            if any(token in name for token in [".vcf.gz", ".tbi", "report"]):
                downloaded.append(str(download_url(url, out_dir, name)))

    vcf = next((p for p in downloaded if p.endswith(".vcf.gz")), "")
    index = next((p for p in downloaded if p.endswith(".tbi")), "")
    if not vcf:
        message = api_error(files) or json.dumps(files)[:500]
        raise SystemExit(
            f"OpenGWAS did not provide a downloadable VCF for {dataset_id}. "
            f"{message} Set OPENGWAS_JWT or provide a direct/local GWAS_location."
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
