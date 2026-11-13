#!/usr/bin/env python3
"""Prepare a run: validate manifest rows, derive IDs and output paths, skip rows whose outputs are complete
for the current settings, and write one JSON file per row to run (no manifest text reaches a shell).

Usage:
    prepare_manifest.py --manifest manifests/test.tsv --work-root GWAS --state-file GWAS/.standardise_state.tsv \
        --output run_manifest.tsv --summary run_manifest_summary.json --rows-dir rows \
        [--config-root manifests] [--reference-root references] [--study-yaml extra.yaml] \
        [--only-study STUDYID] [--force false] [--report true] [--row-limit 0] [--repo-root .]

Exits 78 (configuration error) when any row is invalid, before anything is downloaded.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import sys
from pathlib import Path
from typing import Any

from gwas_pipeline_common import (
    EXIT_CONFIG,
    MANIFEST_COLUMNS,
    MANIFEST_OPTIONAL_COLUMNS,
    SUPER_POPULATIONS,
    atomic_write_text,
    infer_source_type,
    normalize_build_value,
    normalize_source_type,
    opengwas_dataset_id,
    output_id_from_row,
    parse_bool,
    read_manifest,
    read_state,
    resolve_ancestry,
    resolve_row_config,
    row_hash,
    safe_path_component,
    sanitize,
    settings_fingerprint,
    state_row_complete,
    study_id_from_row,
    write_json,
)

DERIVED_COLUMNS = ["study_id", "output_id", "row_hash", "raw_dir", "standardise_dir", "output_prefix", "settings_hash"]
RUN_COLUMNS = MANIFEST_COLUMNS + MANIFEST_OPTIONAL_COLUMNS + DERIVED_COLUMNS
SUPPORTED_SOURCES = {"gwas_catalog", "opengwas", "direct_url", "local"}
SEX_VALUES = {"combined", "male", "female"}
FORMAT_VALUES = {"table", "vcf"}


def resolve_local_location(location: str, repo_root: Path, manifest_dir: Path) -> str:
    candidate = Path(location.replace("file://", ""))
    if candidate.is_absolute():
        return str(candidate)
    for base in (repo_root, manifest_dir):
        if (base / candidate).exists():
            return str((base / candidate).resolve())
    return location


def validate_row(row: dict[str, str], line: int, source_type: str) -> list[str]:
    errors = []
    if resolve_ancestry(row.get("ancestry")) is None:
        errors.append(f"row {line}: ancestry '{row.get('ancestry', '')}' is not a 1000 Genomes super-population "
                      "(ALL, AFR, AMR, EAS, EUR, SAS) or an unambiguous synonym")
    if source_type not in SUPPORTED_SOURCES:
        errors.append(f"row {line}: unsupported source_type '{source_type}'")
    if source_type == "opengwas" and not opengwas_dataset_id(row["GWAS_location"]):
        errors.append(f"row {line}: no OpenGWAS dataset ID found in {row['GWAS_location']}")
    if source_type == "gwas_catalog" and not re.search(r"GCST\d+", row["GWAS_location"], flags=re.I):
        errors.append(f"row {line}: no GWAS Catalog accession (GCST...) found in {row['GWAS_location']}")
    sex = str(row.get("sex", "") or "").strip().lower()
    if sex and sex not in SEX_VALUES:
        errors.append(f"row {line}: sex '{row.get('sex')}' must be combined, male or female (or empty)")
    build = str(row.get("input_build", "") or "").strip()
    if build and not normalize_build_value(build):
        errors.append(f"row {line}: input_build '{build}' must be hg19 (GRCh37) or hg38 (GRCh38) (or empty)")
    population = str(row.get("population", "") or "").strip()
    if population and population.upper() not in SUPER_POPULATIONS:
        errors.append(f"row {line}: population '{population}' must be ALL, AFR, AMR, EAS, EUR or SAS (or empty)")
    data_format = str(row.get("format", "") or "").strip().lower()
    if data_format and data_format not in FORMAT_VALUES:
        errors.append(f"row {line}: format '{row.get('format')}' must be table or vcf (or empty to detect it)")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--state-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--rows-dir", default="rows")
    parser.add_argument("--config-root", default="")
    parser.add_argument("--reference-root", default="")
    parser.add_argument("--study-yaml", default="")
    parser.add_argument("--only-study", default="")
    parser.add_argument("--force", default="false")
    parser.add_argument("--report", default="true")
    parser.add_argument("--row-limit", type=int, default=0)
    args = parser.parse_args()

    manifest = Path(args.manifest)
    force, report = parse_bool(args.force), parse_bool(args.report)
    work_root, repo_root = Path(args.work_root), Path(args.repo_root)
    config_root = Path(args.config_root) if args.config_root else manifest.parent
    reference_root = Path(args.reference_root) if args.reference_root else repo_root / "references"
    state = read_state(args.state_file)

    errors: list[str] = []  # problems with the manifest as a whole: they stop the run before any row starts
    invalid: list[dict[str, Any]] = []  # rows that fail a check: listed as failed while the other rows run
    seen_outputs: dict[str, int] = {}
    eligible: list[dict[str, Any]] = []
    skipped = filtered = 0
    for row in read_manifest(manifest):
        line = int(row.pop("_line", 0))
        study_id = study_id_from_row(row)
        if args.only_study and sanitize(args.only_study) != sanitize(study_id):
            filtered += 1
            continue
        source_type = normalize_source_type(row["source_type"]) if row.get("source_type") else infer_source_type(row["GWAS_location"])
        row_errors = validate_row(row, line, source_type)
        original_location = row["GWAS_location"]
        if source_type == "local":
            row["GWAS_location"] = resolve_local_location(original_location, repo_root, manifest.parent)
        output_id = output_id_from_row({**row, "GWAS_location": original_location}, source_type)
        standardise_dir = work_root / study_id / "standardise" / output_id
        output_prefix = standardise_dir / "standardised"
        if str(output_prefix) in seen_outputs:
            errors.append(f"rows {seen_outputs[str(output_prefix)]} and {line} write to the same output directory "
                          f"{standardise_dir}; give them distinct GWAS_location values or split the study")
            continue
        seen_outputs[str(output_prefix)] = line
        out_row: dict[str, Any] = {col: row.get(col, "") for col in MANIFEST_COLUMNS + MANIFEST_OPTIONAL_COLUMNS}
        out_row.update({
            "source_type": source_type, "study_id": study_id, "output_id": output_id,
            "row_hash": row_hash({**row, "GWAS_location": original_location}),
            "raw_dir": str(work_root / study_id / "raw" / output_id), "standardise_dir": str(standardise_dir),
            "output_prefix": str(output_prefix),
        })
        config = None
        if not row_errors:
            try:
                config = resolve_row_config(config_root, out_row, args.study_yaml)
            except SystemExit:
                row_errors.append(f"row {line}: study YAML '{out_row.get('study_yaml') or args.study_yaml}' not found")
        if row_errors:
            invalid.append({**out_row, "errors": row_errors})
            continue
        out_row["settings_hash"] = settings_fingerprint(out_row, config, reference_root)
        expected = {
            "standard_vcf": str(output_prefix.with_suffix(".gwas.vcf.gz")),
            "standard_index": str(output_prefix.with_suffix(".gwas.vcf.gz")) + ".tbi",
            "metadata_json": str(output_prefix.with_suffix(".metadata.json")),
            "qc_json": str(output_prefix.with_suffix(".qc.json")),
            "report_html": str(output_prefix.with_suffix(".report.html")),
            "settings_hash": out_row["settings_hash"],
        }
        previous = state.get(out_row["row_hash"])
        if not force and previous and all(previous.get(k, "") == v for k, v in expected.items()) and state_row_complete(previous, report):
            skipped += 1
            continue
        eligible.append(out_row)

    if errors:
        print("Manifest validation failed:\n  " + "\n  ".join(errors + [e for r in invalid for e in r["errors"]]), file=sys.stderr)
        raise SystemExit(EXIT_CONFIG)
    if invalid:
        print("These rows fail the manifest checks and are not run (the other rows run):\n  "
              + "\n  ".join(e for r in invalid for e in r["errors"]), file=sys.stderr)

    row_limit = max(args.row_limit, 0)
    to_run = eligible[:row_limit] if row_limit else eligible
    rows_dir = Path(args.rows_dir)
    rows_dir.mkdir(parents=True, exist_ok=True)
    for stale in rows_dir.glob("*.json"):
        stale.unlink()
    for index, out_row in enumerate(to_run, start=1):
        name = f"{index:04d}_{safe_path_component(out_row['study_id'])}__{safe_path_component(out_row['output_id'])}.json"
        write_json(rows_dir / name, out_row)

    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=RUN_COLUMNS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(to_run)
    writer.writerows(invalid)  # listed so the run reports them as failed and manifest.html keeps them
    atomic_write_text(args.output, buffer.getvalue())
    write_json(args.summary, {
        "manifest": str(manifest.resolve()), "state_file": str(Path(args.state_file).resolve()),
        "rows_to_run": len(to_run), "rows_eligible": len(eligible), "rows_deferred_due_to_limit": len(eligible) - len(to_run),
        "rows_skipped_complete": skipped, "rows_filtered": filtered, "force": force, "report": report,
        "row_limit": row_limit, "rows_invalid": len(invalid),
        "invalid": [{"study_id": r["study_id"], "output_id": r["output_id"], "errors": r["errors"]} for r in invalid],
    })


if __name__ == "__main__":
    main()
