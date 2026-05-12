#!/usr/bin/env python3
"""Prepare a persistent run manifest and skip already completed rows."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from gwas_pipeline_common import (
    MANIFEST_COLUMNS,
    MANIFEST_OPTIONAL_COLUMNS,
    infer_source_type,
    normalize_source_type,
    output_id_from_row,
    read_manifest,
    read_state,
    row_hash,
    sanitize,
    state_row_complete,
    study_id_from_row,
    write_json,
)


RUN_COLUMNS = MANIFEST_COLUMNS + [
    "sex",
    "source_type",
    "input_build",
    "population",
    "delimiter",
    "format",
    "study_yaml",
    "study_id",
    "output_id",
    "row_hash",
    "raw_dir",
    "standardise_dir",
    "output_prefix",
]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--repo-root", default=".")
    parser.add_argument("--state-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--only-study", default="")
    parser.add_argument("--force", default="false")
    parser.add_argument("--qcplot", default="true")
    parser.add_argument("--row-limit", type=int, default=0)
    args = parser.parse_args()

    force = str(args.force).lower() == "true"
    qcplot = str(args.qcplot).lower() == "true"
    work_root = Path(args.work_root)
    repo_root = Path(args.repo_root)
    state = read_state(args.state_file)

    eligible_rows: list[dict[str, str]] = []
    skipped = 0
    filtered = 0

    for row in read_manifest(args.manifest):
        study_id = study_id_from_row(row)
        if args.only_study and sanitize(args.only_study) != sanitize(study_id):
            filtered += 1
            continue

        source_type = normalize_source_type(row.get("source_type", "")) if row.get("source_type") else infer_source_type(row["GWAS_location"])
        original_location = row["GWAS_location"]
        if source_type == "local":
            candidate = Path(original_location.replace("file://", ""))
            if not candidate.is_absolute():
                repo_candidate = repo_root / candidate
                manifest_candidate = Path(args.manifest).parent / candidate
                if repo_candidate.exists():
                    row["GWAS_location"] = str(repo_candidate.resolve())
                elif manifest_candidate.exists():
                    row["GWAS_location"] = str(manifest_candidate.resolve())
        output_row = dict(row)
        output_row["GWAS_location"] = original_location
        output_id = output_id_from_row(output_row, source_type)
        hash_row = dict(row)
        hash_row["GWAS_location"] = original_location
        rhash = row_hash(hash_row)
        raw_dir = work_root / study_id / "raw" / output_id
        standardise_dir = work_root / study_id / "standardise" / output_id
        output_prefix = standardise_dir / "standardised"

        expected_state_paths = {
            "standard_vcf": str(output_prefix.with_suffix(".gwas.vcf.gz")),
            "standard_index": str(output_prefix.with_suffix(".gwas.vcf.gz")) + ".tbi",
            "metadata_json": str(output_prefix.with_suffix(".metadata.json")),
            "qc_json": str(output_prefix.with_suffix(".qc.json")),
            "plot_tiff": str(output_prefix.with_suffix(".tiff")),
        }
        state_matches_expected = rhash in state and all(state[rhash].get(key, "") == value for key, value in expected_state_paths.items())

        if not force and state_matches_expected and state_row_complete(state[rhash], qcplot):
            skipped += 1
            continue

        out_row = dict(row)
        for col in MANIFEST_OPTIONAL_COLUMNS:
            out_row.setdefault(col, "")
        out_row.update(
            {
                "study_id": study_id,
                "output_id": output_id,
                "row_hash": rhash,
                "source_type": source_type,
                "raw_dir": str(raw_dir),
                "standardise_dir": str(standardise_dir),
                "output_prefix": str(output_prefix),
            }
        )
        eligible_rows.append(out_row)

    row_limit = max(args.row_limit, 0)
    rows_out = eligible_rows[:row_limit] if row_limit else eligible_rows
    deferred = max(len(eligible_rows) - len(rows_out), 0)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RUN_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    write_json(
        args.summary,
        {
            "manifest": str(Path(args.manifest).resolve()),
            "state_file": str(Path(args.state_file).resolve()),
            "rows_to_run": len(rows_out),
            "rows_eligible": len(eligible_rows),
            "rows_deferred_due_to_limit": deferred,
            "rows_skipped_complete": skipped,
            "rows_filtered": filtered,
            "force": force,
            "qcplot": qcplot,
            "row_limit": row_limit,
        },
    )


if __name__ == "__main__":
    main()
