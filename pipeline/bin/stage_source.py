#!/usr/bin/env python3
"""Stage raw GWAS data from source-specific or generic locations."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from download_gwas_catalog import download_gwas_catalog
from download_opengwas import download_opengwas
from gwas_pipeline_common import (
    download_url,
    ensure_dir,
    guess_filename_from_url,
    parse_bool,
    read_json,
    resolve_study_config,
    sanitize_filename,
    shell_copy_or_download,
    write_json,
)


def copy_local(location: str, raw_dir: Path) -> Path:
    src = Path(location.replace("file://", ""))
    if not src.exists():
        raise FileNotFoundError(f"Local GWAS file not found: {src}")
    dest = raw_dir / sanitize_filename(src.name)
    if not dest.exists():
        shutil.copy2(src, dest)
    return dest


def stage_s3_or_gs(location: str, raw_dir: Path, source_type: str) -> Path:
    raw_dir.mkdir(parents=True, exist_ok=True)
    command = ["aws", "s3", "cp", location, str(raw_dir)] if source_type == "s3" else ["gsutil", "cp", location, str(raw_dir)]
    subprocess.run(command, check=True)
    return raw_dir / Path(location).name


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gwas-location", required=True)
    parser.add_argument("--phenotype", required=True)
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--author", required=True)
    parser.add_argument("--year", required=True)
    parser.add_argument("--pmid", required=True)
    parser.add_argument("--sex", default="")
    parser.add_argument("--input-build", default="")
    parser.add_argument("--population", default="")
    parser.add_argument("--delimiter", default="")
    parser.add_argument("--format", default="")
    parser.add_argument("--manifest-study-yaml", default="")
    parser.add_argument("--study-id", required=True)
    parser.add_argument("--output-id", required=True)
    parser.add_argument("--row-hash", required=True)
    parser.add_argument("--source-type", required=True)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--standardise-dir", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--config-root", required=True)
    parser.add_argument("--study-yaml", default="")
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    raw_dir = ensure_dir(args.raw_dir)
    ensure_dir(args.standardise_dir)
    dry_run = parse_bool(args.dry_run)
    config = resolve_study_config(args.config_root, args.study_id, args.study_yaml or None, args.manifest_study_yaml or None)

    metadata = {
        "source_type": args.source_type,
        "dry_run": dry_run,
        "config": config.get("_config_path", ""),
    }

    if dry_run:
        raw_path = raw_dir / sanitize_filename(guess_filename_from_url(args.gwas_location, args.output_id))
    elif args.source_type == "gwas_catalog":
        metadata.update(download_gwas_catalog(args.gwas_location, raw_dir, args.pmid))
        raw_path = Path(metadata["raw_path"])
    elif args.source_type == "opengwas":
        metadata.update(download_opengwas(args.gwas_location, raw_dir))
        raw_path = Path(metadata["raw_path"])
    elif args.source_type == "direct_url":
        raw_path = download_url(args.gwas_location, raw_dir)
    elif args.source_type == "local":
        raw_path = copy_local(args.gwas_location, raw_dir)
    elif args.source_type in {"s3", "gs"}:
        raw_path = stage_s3_or_gs(args.gwas_location, raw_dir, args.source_type)
    elif args.source_type == "synapse":
        token = os.environ.get("SYNAPSE_AUTH_TOKEN", "")
        if not token:
            raise SystemExit("SYNAPSE_AUTH_TOKEN is required for Synapse staging")
        raise NotImplementedError("Synapse staging is reserved for the authenticated environment")
    else:
        raw_path = shell_copy_or_download(args.gwas_location, raw_dir)

    record = {
        "GWAS_location": args.gwas_location,
        "phenotype": args.phenotype,
        "ancestry": args.ancestry,
        "author": args.author,
        "year": args.year,
        "PMID": args.pmid,
        "sex": args.sex,
        "input_build": args.input_build,
        "population": args.population,
        "delimiter": args.delimiter,
        "format": args.format,
        "study_yaml": args.manifest_study_yaml,
        "study_id": args.study_id,
        "output_id": args.output_id,
        "row_hash": args.row_hash,
        "source_type": args.source_type,
        "raw_dir": str(raw_dir),
        "standardise_dir": args.standardise_dir,
        "output_prefix": args.output_prefix,
        "raw_path": str(raw_path),
        "reference_root": args.reference_root,
        "metadata": metadata,
    }
    write_json(args.output, record)


if __name__ == "__main__":
    main()
