#!/usr/bin/env python3
"""Validate reference lookup tables required by the pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", default="references")
    args = parser.parse_args()
    lookup = Path(args.reference_root) / "1000genomes" / "phase3" / "lookup"
    missing = []
    for pop in ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]:
        if not (lookup / f"{pop}.tsv.gz").exists() and not (lookup / f"{pop}.tsv").exists():
            missing.append(pop)
    if missing:
        raise SystemExit(f"Missing lookup tables for: {', '.join(missing)}")
    print(f"Reference lookup tables found in {lookup}")


if __name__ == "__main__":
    main()
