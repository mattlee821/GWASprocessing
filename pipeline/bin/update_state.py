#!/usr/bin/env python3
"""Record completed rows in the persistent state file.

Usage:
    update_state.py --state-file GWAS/.standardise_state.tsv --qcplot true --output state_update_summary.json \
        [--records a.row_record.json b.row_record.json ...]
A row is recorded only when all of its outputs exist and are non-empty; the settings hash stored with it
lets prepare_manifest.py regenerate the row when settings change.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from gwas_pipeline_common import parse_bool, read_json, read_state, state_row_complete, utc_now, write_json, write_state


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-file", required=True)
    parser.add_argument("--records", nargs="*", default=[])
    parser.add_argument("--qcplot", default="true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    qcplot = parse_bool(args.qcplot)
    state = read_state(args.state_file)
    updated: list[str] = []
    ignored: list[dict[str, str]] = []
    for record_path in args.records:
        if not Path(record_path).exists():
            ignored.append({"record": record_path, "reason": "record file not found"})
            continue
        record = read_json(record_path)
        row = {
            "row_hash": record.get("row_hash", ""), "study_id": record.get("study_id", ""), "output_id": record.get("output_id", ""),
            "gwas_location": record.get("GWAS_location", ""), "standard_vcf": record.get("standard_vcf", ""),
            "standard_index": record.get("standard_index", ""), "metadata_json": record.get("metadata_json", ""),
            "qc_json": record.get("qc_json", ""), "plot_tiff": record.get("plot_tiff", ""),
            "settings_hash": (record.get("qc") or {}).get("settings_hash", "") or record.get("settings_hash", ""),
            "updated_at": utc_now(),
        }
        if not row["row_hash"]:
            ignored.append({"record": record_path, "reason": "no row_hash"})
        elif not state_row_complete(row, qcplot):
            ignored.append({"record": record_path, "reason": "outputs incomplete"})
        else:
            state[row["row_hash"]] = row
            updated.append(row["row_hash"])
    write_state(args.state_file, state.values())
    write_json(args.output, {"state_file": str(Path(args.state_file).resolve()), "updated": len(updated),
                             "updated_row_hashes": updated, "ignored": len(ignored), "ignored_records": ignored})


if __name__ == "__main__":
    main()
