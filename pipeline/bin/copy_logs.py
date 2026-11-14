#!/usr/bin/env python3
"""Copy logs into <work_root>/<study_id>/logs/.

  --row-json ROW.json   (run by Nextflow after every row, also when it failed) the row's task output as
                        logs/<output_id>.<time>.log and .err, taken from .command.log and .command.err of the task folder
  --run-manifest TSV    (run by src/standardise.sh when the run ends) the Nextflow log and the SLURM output of the run,
                        into the logs/ of every study in the run, as run.<time>.nextflow.log and run.<time>.slurm.out
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import time
from pathlib import Path


def stamp() -> str:
    return time.strftime("%Y%m%d-%H%M%S")


def copy(source: Path, target: Path) -> None:
    if source.is_file():
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--row-json")
    parser.add_argument("--run-manifest")
    parser.add_argument("--nextflow-log")
    parser.add_argument("--slurm-out")
    args = parser.parse_args()
    work_root = Path(args.work_root)
    now = stamp()
    if args.row_json:
        row = json.loads(Path(args.row_json).read_text(encoding="utf-8"))
        logs = work_root / row["study_id"] / "logs"
        copy(Path(".command.log"), logs / f"{row['output_id']}.{now}.log")
        copy(Path(".command.err"), logs / f"{row['output_id']}.{now}.err")
    if args.run_manifest and Path(args.run_manifest).is_file():
        with open(args.run_manifest, encoding="utf-8") as handle:
            studies = sorted({r["study_id"] for r in csv.DictReader(handle, delimiter="\t")})
        for study in studies:
            logs = work_root / study / "logs"
            if args.nextflow_log:
                copy(Path(args.nextflow_log), logs / f"run.{now}.nextflow.log")
            if args.slurm_out:
                copy(Path(args.slurm_out), logs / f"run.{now}.slurm.out")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
