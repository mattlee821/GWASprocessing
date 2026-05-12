#!/usr/bin/env python3
"""Persist successful standardisation records."""

from __future__ import annotations

import argparse
from pathlib import Path

from gwas_pipeline_common import STATE_COLUMNS, read_json, read_state, state_row_complete, utc_now, write_json, write_state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-file", required=True)
    parser.add_argument("--records", nargs="*", default=[])
    parser.add_argument("--qcplot", default="true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    qcplot = str(args.qcplot).lower() == "true"
    state = read_state(args.state_file)
    updated = 0
    ignored = 0

    for record_path in args.records:
        if not Path(record_path).exists():
            ignored += 1
            continue
        record = read_json(record_path)
        row = {
            "row_hash": record.get("row_hash", ""),
            "study_id": record.get("study_id", ""),
            "output_id": record.get("output_id", ""),
            "gwas_location": record.get("GWAS_location", ""),
            "standard_vcf": record.get("standard_vcf", ""),
            "standard_index": record.get("standard_index", ""),
            "metadata_json": record.get("metadata_json", ""),
            "qc_json": record.get("qc_json", ""),
            "plot_tiff": record.get("plot_tiff", ""),
            "updated_at": utc_now(),
        }
        if row["row_hash"] and state_row_complete(row, qcplot):
            state[row["row_hash"]] = row
            updated += 1
        else:
            ignored += 1

    write_state(args.state_file, state.values())
    write_json(args.output, {"state_file": str(Path(args.state_file).resolve()), "updated": updated, "ignored": ignored})


if __name__ == "__main__":
    main()
