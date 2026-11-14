#!/usr/bin/env python3
"""Stage the raw GWAS file for one manifest row and write its stage record.

Usage:
    stage_source.py --row-json rows/0001_x.json --config-root manifests --reference-root references \
        [--study-yaml extra.yaml] [--dry-run false] --output stage_record.json
The row JSON comes from prepare_manifest.py, so manifest text never passes through a shell.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import shutil
from pathlib import Path

from download_gwas_catalog import download_gwas_catalog, md5, verify_gzip
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


def md5_sidecar(path: Path) -> Path | None:
    """The checksum file beside a local file: <name>.md5sum or <name>.md5, or the same with .gz taken off the name
    (deCODE ships X.txt.md5sum for X.txt.gz)."""
    names = [path.name] + ([path.name[:-3]] if path.name.endswith(".gz") else [])
    for name in names:
        for suffix in (".md5sum", ".md5"):
            candidate = path.with_name(name + suffix)
            if candidate.is_file():
                return candidate
    return None


def verify_local_md5(path: Path) -> dict:
    """Check a local file against its checksum file, when it has one. The checksum may be of the file as it is or, for
    a .gz, of the text inside it (deCODE's md5sum is of the uncompressed file). A file without a checksum file is
    not checked here (the gzip check covers it)."""
    sidecar = md5_sidecar(path)
    if sidecar is None:
        return {"checked": False, "ok": False}
    stems = {path.name, path.name[:-3] if path.name.endswith(".gz") else path.name}
    lines = [line.split() for line in sidecar.read_text(errors="replace").splitlines() if line.strip()]
    expected = None
    for parts in lines:
        found = re.fullmatch(r"[a-fA-F0-9]{32}", parts[0])
        if found and (len(lines) == 1 or len(parts) == 1 or Path(parts[-1].lstrip("*")).name in stems):
            expected = parts[0].lower()
            break
    result: dict = {"checked": True, "ok": False, "sidecar": sidecar.name, "expected": expected}
    if expected is None:
        result["error"] = f"no md5 for {path.name} found in {sidecar.name}"
        return result
    if md5(path) == expected:
        return {**result, "ok": True, "matched": "file"}
    if path.name.endswith(".gz"):
        digest = hashlib.md5()
        try:
            with gzip.open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
        except Exception as exc:
            return {**result, "error": f"gzip error: {exc}"}
        if digest.hexdigest() == expected:
            return {**result, "ok": True, "matched": "uncompressed"}
    return result


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
        source = Path(location.replace("file://", ""))
        if not source.exists():
            raise FileNotFoundError(f"Local GWAS file not found: {source}")
        metadata["md5"] = verify_local_md5(source)
        if metadata["md5"]["checked"] and not metadata["md5"]["ok"]:
            raise data_error(f"MD5 check failed for {source}: {metadata['md5'].get('error') or 'the file does not match ' + metadata['md5']['sidecar']}")
        study_raw = raw_dir.parent.resolve()  # <work root>/<study>/raw: a file already here is the study's raw file
        if source.resolve().is_relative_to(study_raw):
            raw_path = source.resolve()
            metadata["in_place"] = True
            try:
                raw_dir.rmdir()  # nothing is staged, so the row's own raw folder stays empty
            except OSError:
                pass
        else:
            raw_path = copy_local(location, raw_dir)
        if not metadata["md5"]["checked"]:
            metadata["gzip"] = verify_gzip(raw_path)
            if metadata["gzip"].get("checked") and not metadata["gzip"].get("ok"):
                raise data_error(f"gzip integrity check failed for {source}: {metadata['gzip'].get('error', '')}")
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
