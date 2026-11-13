#!/usr/bin/env python3
"""Stage the raw GWAS file for one manifest row and write its stage record.

Usage:
    stage_source.py --row-json rows/0001_x.json --config-root manifests --reference-root references \
        [--study-yaml extra.yaml] [--dry-run false] --output stage_record.json
The row JSON comes from prepare_manifest.py, so manifest text never passes through a shell.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from download_gwas_catalog import download_gwas_catalog, verify_gzip
from download_opengwas import download_opengwas
from gwas_pipeline_common import (
    MANIFEST_COLUMNS,
    MANIFEST_OPTIONAL_COLUMNS,
    data_error,
    download_url,
    ensure_dir,
    guess_filename_from_url,
    parse_bool,
    read_json,
    resolve_row_config,
    sanitize_filename,
    write_json,
)


def copy_local(location: str, raw_dir: Path) -> Path:
    """Copy a local file into raw_dir, refreshing the copy when the source size or mtime changed."""
    src = Path(location.replace("file://", ""))
    if not src.exists():
        raise FileNotFoundError(f"Local GWAS file not found: {src}")
    dest = raw_dir / sanitize_filename(src.name)
    stat = src.stat()
    if not dest.exists() or dest.stat().st_size != stat.st_size or int(dest.stat().st_mtime) != int(stat.st_mtime):
        tmp = dest.with_name(dest.name + ".part")
        shutil.copy2(src, tmp)
        tmp.replace(dest)
    return dest


def checked_download(fetch, metadata: dict) -> Path:
    """Download with fetch(), then check a .gz file's gzip integrity: a file that fails (including one left by an
    earlier interrupted run) is deleted and downloaded once more, and a second failure stops the row."""
    path = fetch()
    result = verify_gzip(path)
    if result.get("checked") and not result.get("ok"):
        path.unlink(missing_ok=True)
        path = fetch()
        result = verify_gzip(path)
        if result.get("checked") and not result.get("ok"):
            path.unlink(missing_ok=True)
            raise ValueError(f"gzip integrity check failed for {path}: {result.get('error', '')}")
    metadata["gzip"] = result
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--row-json", required=True)
    parser.add_argument("--config-root", required=True)
    parser.add_argument("--study-yaml", default="")
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    row = read_json(args.row_json)
    raw_dir = ensure_dir(row["raw_dir"])
    ensure_dir(row["standardise_dir"])
    dry_run = parse_bool(args.dry_run)
    source_type = row["source_type"]
    location = row["GWAS_location"]
    config = resolve_row_config(args.config_root, row, args.study_yaml)
    metadata: dict = {"source_type": source_type, "dry_run": dry_run, "config": config.get("_config_path", "")}

    if dry_run:
        raw_path = raw_dir / sanitize_filename(guess_filename_from_url(location, row["output_id"]))
    elif source_type == "gwas_catalog":
        prefer = str(config.get("gwas_catalog_file", "raw") or "raw").strip().lower()
        metadata.update(download_gwas_catalog(location, raw_dir, row.get("PMID", ""), "harmonised" if prefer == "harmonised" else "raw"))
        raw_path = Path(metadata["raw_path"])
    elif source_type == "opengwas":
        info: dict = {}

        def fetch_opengwas() -> Path:
            info.update(download_opengwas(location, raw_dir))
            return Path(info["raw_path"])
        raw_path = checked_download(fetch_opengwas, metadata)
        metadata.update(info)
    elif source_type == "direct_url":
        raw_path = checked_download(lambda: download_url(location, raw_dir), metadata)
    elif source_type == "local":
        raw_path = copy_local(location, raw_dir)
    else:
        raise data_error(f"Unsupported source_type '{source_type}' for {location}")

    record = {col: row.get(col, "") for col in MANIFEST_COLUMNS + MANIFEST_OPTIONAL_COLUMNS}
    record.update({
        "study_id": row["study_id"], "output_id": row["output_id"], "row_hash": row["row_hash"], "source_type": source_type,
        "raw_dir": str(raw_dir), "standardise_dir": row["standardise_dir"], "output_prefix": row["output_prefix"],
        "settings_hash": row.get("settings_hash", ""), "raw_path": str(raw_path), "reference_root": args.reference_root,
        "metadata": metadata,
    })
    write_json(args.output, record)


if __name__ == "__main__":
    main()
