#!/usr/bin/env python3
"""Download summary statistics from Synapse via the authenticated Synapse REST API.

A personal access token in SYNAPSE_TOKEN (or SYNAPSE_AUTH_TOKEN) is required, with
View and Download permissions and any dataset Terms of Use accepted on synapse.org.
Missing or rejected tokens exit 78 (EX_CONFIG) so Nextflow stops the run immediately;
a missing entity or non-file entity exits 65 (EX_DATAERR) so only that row fails.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from pathlib import Path
from typing import Any

import requests

from download_gwas_catalog import md5, verify_gzip
from gwas_pipeline_common import (
    atomic_write_text,
    config_error,
    data_error,
    download_url,
    ensure_dir,
    safe_path_component,
    sanitize_filename,
    synapse_entity_id,
    write_json,
)

DEFAULT_REPO_ROOT = "https://repo-prod.prod.sagebase.org/repo/v1"
DEFAULT_FILE_ROOT = "https://file-prod.prod.sagebase.org/file/v1"
TOKEN_HELP = (
    "Provide a personal access token in .env as SYNAPSE_TOKEN=<token> "
    "(generate one at https://www.synapse.org/ under Account Settings -> Personal Access Tokens "
    "with View and Download permissions, and ensure dataset Terms of Use are accepted)."
)


def repo_root() -> str:
    """Synapse repository service root; SYNAPSE_REPO_ROOT or SYNAPSE_API_ROOT overrides it."""
    if os.environ.get("SYNAPSE_REPO_ROOT"):
        return os.environ["SYNAPSE_REPO_ROOT"].rstrip("/")
    if os.environ.get("SYNAPSE_API_ROOT"):
        return os.environ["SYNAPSE_API_ROOT"].rstrip("/") + "/repo/v1"
    return DEFAULT_REPO_ROOT


def file_root() -> str:
    """Synapse file-handle service root; SYNAPSE_FILE_ROOT or SYNAPSE_API_ROOT overrides it."""
    if os.environ.get("SYNAPSE_FILE_ROOT"):
        return os.environ["SYNAPSE_FILE_ROOT"].rstrip("/")
    if os.environ.get("SYNAPSE_API_ROOT"):
        return os.environ["SYNAPSE_API_ROOT"].rstrip("/") + "/file/v1"
    return DEFAULT_FILE_ROOT


def synapse_id(location: str) -> str:
    """Synapse entity ID such as syn51365303 or syn51365303.1, in lowercase."""
    entity = synapse_entity_id(location)
    if not entity:
        raise data_error(f"Could not find a Synapse ID (syn...) in: {location}")
    return entity


def synapse_token() -> str:
    token = (os.environ.get("SYNAPSE_TOKEN") or os.environ.get("SYNAPSE_AUTH_TOKEN") or "").strip()
    if not token:
        raise config_error("Synapse authentication is required but SYNAPSE_TOKEN is not set. " + TOKEN_HELP)
    return token


def _error_detail(response: requests.Response) -> str:
    try:
        payload = response.json()
        if isinstance(payload, dict):
            return str(payload.get("reason") or payload.get("message") or response.text)[:240]
    except ValueError:
        pass
    return response.text[:240]


def synapse_api(method: str, url: str, token: str, retries: int = 3, **kwargs: Any) -> dict[str, Any]:
    headers = {"Authorization": f"Bearer {token}", "Accept": "application/json"}
    last_exc: Exception | None = None
    for attempt in range(1, max(1, retries) + 1):
        try:
            response = requests.request(method, url, headers=headers, timeout=120, **kwargs)
        except requests.RequestException as exc:
            last_exc = exc
            if attempt < retries:
                time.sleep(min(10.0, 1.5 ** attempt))
                continue
            raise SystemExit(f"Synapse API request failed ({method} {url}): {exc}") from exc

        if response.status_code in (401, 403):
            detail = _error_detail(response)
            raise config_error(
                f"Synapse rejected the request with HTTP {response.status_code}: {detail} "
                "The token in SYNAPSE_TOKEN is expired, invalid, lacks Download permission, "
                "or dataset Terms of Use have not been accepted. " + TOKEN_HELP
            )
        if response.status_code == 404:
            raise data_error(f"Synapse entity not found (HTTP 404 at {url}): {_error_detail(response)}")
        if response.status_code in (429, 502, 503, 504) and attempt < retries:
            wait = float(response.headers.get("Retry-After") or min(15.0, 2.0 ** attempt))
            time.sleep(min(30.0, max(0.1, wait)))
            continue
        response.raise_for_status()
        data = response.json()
        return data if isinstance(data, dict) else {"data": data}
    raise SystemExit(f"Synapse API request failed ({method} {url}): {last_exc}")


def request_file_handle(base_id: str, file_handle_id: str, token: str) -> dict[str, Any]:
    payload = {
        "includeFileHandles": True,
        "includePreSignedURLs": True,
        "requestedFiles": [
            {
                "fileHandleId": file_handle_id,
                "associateObjectId": base_id,
                "associateObjectType": "FileEntity",
            }
        ],
    }
    batch = synapse_api("POST", f"{file_root()}/fileHandle/batch", token, json=payload)
    items = batch.get("requestedFiles") or []
    if not items:
        raise SystemExit(f"Synapse returned no file handle entry for {base_id} ({file_handle_id})")
    item = items[0]
    failure = item.get("failureCode")
    if failure in ("UNAUTHORIZED", "FORBIDDEN"):
        raise config_error(
            f"Synapse denied download access to {base_id} ({failure}). Ensure SYNAPSE_TOKEN has "
            "Download permission and that your Synapse account has accepted any Terms of Use / "
            "Data Access conditions for this dataset. " + TOKEN_HELP
        )
    if failure == "NOT_FOUND":
        raise data_error(f"Synapse file handle {file_handle_id} for {base_id} was not found")
    if failure:
        raise SystemExit(f"Synapse fileHandle/batch failed for {base_id}: {failure}")
    return item


def download_synapse(location: str, out_dir: str | Path) -> dict[str, Any]:
    out_dir = ensure_dir(out_dir)
    syn_id = synapse_id(location)
    token = synapse_token()
    base_id, _, version = syn_id.partition(".")

    entity_url = f"{repo_root()}/entity/{base_id}/version/{version}" if version else f"{repo_root()}/entity/{base_id}"
    entity = synapse_api("GET", entity_url, token)
    file_handle_id = str(entity.get("dataFileHandleId") or "").strip()
    if not file_handle_id:
        raise data_error(
            f"Synapse entity {syn_id} ({entity.get('concreteType', 'unknown')}) is not a downloadable file "
            "(missing dataFileHandleId)"
        )

    item = request_file_handle(base_id, file_handle_id, token)
    fh = item.get("fileHandle") or {}
    presigned_url = str(item.get("preSignedURL") or "").strip()
    if not presigned_url:
        raise SystemExit(f"Synapse returned no preSignedURL for {syn_id} (concreteType={fh.get('concreteType')})")

    filename = sanitize_filename(fh.get("fileName") or entity.get("name") or f"{base_id}.dat")
    expected_md5 = str(fh.get("contentMd5") or "").strip().lower()
    summary = download_url(presigned_url, out_dir, filename)

    md5_result: dict[str, Any] = {"checked": False, "ok": False}
    gzip_result: dict[str, Any] = {"checked": False, "ok": False}
    if expected_md5:
        observed = md5(summary)
        md5_result = {"checked": True, "expected": expected_md5, "observed": observed, "ok": observed == expected_md5}
        if not md5_result["ok"]:
            summary.unlink(missing_ok=True)
            fresh_url = str(request_file_handle(base_id, file_handle_id, token).get("preSignedURL") or presigned_url)
            summary = download_url(fresh_url, out_dir, filename)
            observed = md5(summary)
            md5_result = {"checked": True, "expected": expected_md5, "observed": observed, "ok": observed == expected_md5}
            if not md5_result["ok"]:
                summary.unlink(missing_ok=True)
                raise ValueError(f"MD5 verification failed for {summary}: {md5_result}")
    else:
        gzip_result = verify_gzip(summary)
        if gzip_result.get("checked") and not gzip_result.get("ok"):
            summary.unlink(missing_ok=True)
            fresh_url = str(request_file_handle(base_id, file_handle_id, token).get("preSignedURL") or presigned_url)
            summary = download_url(fresh_url, out_dir, filename)
            gzip_result = verify_gzip(summary)
            if gzip_result.get("checked") and not gzip_result.get("ok"):
                summary.unlink(missing_ok=True)
                raise ValueError(f"Gzip integrity check failed for {summary}: {gzip_result}")

    metadata = {
        "source_type": "synapse",
        "synapse_id": syn_id,
        "entity_id": base_id,
        "version_number": entity.get("versionNumber"),
        "entity_name": entity.get("name", ""),
        "file_handle_id": file_handle_id,
        "file_name": filename,
        "content_size": fh.get("contentSize"),
        "raw_path": str(summary),
        "md5": md5_result,
        "gzip": gzip_result,
    }
    write_json(out_dir / f"{safe_path_component(syn_id)}.synapse_metadata.json", metadata)
    return metadata


def list_synapse_children(parent_location: str) -> list[dict[str, Any]]:
    """List all file entities directly inside a Synapse folder or project (paginated)."""
    parent_id = synapse_id(parent_location).partition(".")[0]
    token = synapse_token()
    children: list[dict[str, Any]] = []
    next_token: str | None = None
    while True:
        body: dict[str, Any] = {
            "parentId": parent_id,
            "includeTypes": ["file"],
            "sortBy": "NAME",
            "sortDirection": "ASC",
        }
        if next_token:
            body["nextPageToken"] = next_token
        resp = synapse_api("POST", f"{repo_root()}/entity/children", token, json=body)
        for entry in resp.get("page") or []:
            syn = str(entry.get("id") or "").strip().lower()
            name = str(entry.get("name") or "").strip()
            if syn and name:
                children.append({
                    "synapse_id": syn,
                    "file": name,
                    "version_number": entry.get("versionNumber", ""),
                })
        next_token = resp.get("nextPageToken")
        if not next_token:
            break
    return children


def write_children_tsv(children: list[dict[str, Any]], out_path: str = "") -> None:
    lines = ["synapse_id\tfile\tversion_number\n"]
    for row in children:
        lines.append(f"{row['synapse_id']}\t{row['file']}\t{row.get('version_number', '')}\n")
    text = "".join(lines)
    if out_path:
        atomic_write_text(out_path, text)
    else:
        sys.stdout.write(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas-location", default="", help="Synapse file ID or URL to download")
    parser.add_argument("--out-dir", default="", help="Output directory for downloaded file")
    parser.add_argument("--metadata", default="", help="Output JSON metadata path")
    parser.add_argument("--list-children", default="", help="List all file entities inside a Synapse folder ID")
    parser.add_argument("--out", default="", help="Output TSV path for --list-children (default: stdout)")
    args = parser.parse_args()

    if args.list_children:
        write_children_tsv(list_synapse_children(args.list_children), args.out)
        return
    if not args.gwas_location or not args.out_dir or not args.metadata:
        parser.error("--gwas-location, --out-dir and --metadata are required (or use --list-children)")
    write_json(args.metadata, download_synapse(args.gwas_location, args.out_dir))


if __name__ == "__main__":
    main()
