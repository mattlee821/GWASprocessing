#!/usr/bin/env python3
"""Check that the reference data GWASprocessing needs is present and in the current format.

Usage: python3 src/reference_data/validate_references.py [--reference-root references]

Required (exit 1 when missing): lookup tables for AFR, AMR, EAS, EUR, SAS and ALL with columns
CHR, SNP, POS19, A1, A2, AF1, and the hg19ToHg38 chain. Recommended (reported as WARNING): the
hg38ToHg19 chain, the rsID coordinate index and the GRCh38 FASTA.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]
REQUIRED_COLUMNS = ["CHR", "SNP", "POS19", "A1", "A2", "AF1"]


def check(root: Path) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    lookup = root / "1000genomes" / "phase3" / "lookup"
    for pop in POPULATIONS:
        path = next((p for p in (lookup / f"{pop}.tsv.gz", lookup / f"{pop}.tsv") if p.exists()), None)
        if path is None:
            errors.append(f"missing lookup table {lookup / (pop + '.tsv.gz')}")
            continue
        columns = list(pd.read_csv(path, sep="\t", nrows=0).columns)
        missing = [c for c in REQUIRED_COLUMNS if c not in columns]
        if missing:
            errors.append(f"{path} lacks {', '.join(missing)}; rebuild it with make_1000genomes_lookup_tables.py --force")
    genome = root / "genome"
    if not (genome / "hg19ToHg38.over.chain.gz").exists():
        errors.append(f"missing {genome / 'hg19ToHg38.over.chain.gz'} (run download_fasta_and_chains.sh)")
    if not (genome / "hg38ToHg19.over.chain.gz").exists():
        warnings.append("hg38ToHg19 chain missing: GRCh38 inputs will have no POS19")
    if not (root / "1000genomes" / "phase3" / "coordinates" / "rs.npy").exists():
        warnings.append("rsID coordinate index missing: files without positions are limited to common SNPs "
                        "(run make_1000genomes_lookup_tables.py)")
    if not any((genome / name).exists() for name in ("GRCh38.fa", "GRCh38.fasta")):
        warnings.append("GRCh38 FASTA missing: REF/ALT will not be aligned to the reference genome (run download_fasta_and_chains.sh)")
    return errors, warnings


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", default="references")
    args = parser.parse_args()
    errors, warnings = check(Path(args.reference_root))
    for warning in warnings:
        print(f"WARNING: {warning}")
    if errors:
        print("Reference data incomplete:\n  " + "\n  ".join(errors), file=sys.stderr)
        raise SystemExit(1)
    print(f"Reference data OK under {args.reference_root}")


if __name__ == "__main__":
    main()
