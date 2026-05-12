#!/usr/bin/env python3
"""Create lightweight rsID/SNPID/EAF lookup tables from 1000 Genomes outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def make_lookup(processed_root: Path, out_root: Path, population: str) -> None:
    bim = processed_root / population / f"{population}.bim"
    stats = processed_root / population / "stats.frq"
    if not bim.exists():
        raise FileNotFoundError(bim)
    df = pd.read_csv(
        bim,
        sep=r"\s+",
        names=["CHR", "SNP", "CM", "POS19", "A1", "A2"],
        dtype={"CHR": str, "SNP": str},
    )
    df["SNPID"] = df["CHR"].astype(str) + ":" + df["POS19"].astype(str)
    df["POS38"] = df["POS19"]
    if stats.exists():
        frq = pd.read_csv(stats, sep=r"\s+")
        if "SNP" in frq.columns and "MAF" in frq.columns:
            df = df.merge(frq[["SNP", "MAF"]], on="SNP", how="left")
            df[f"EAF_{population}"] = df["MAF"]
            df = df.drop(columns=["MAF"])
    out_root.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_root / f"{population}.tsv.gz", sep="\t", index=False, compression="gzip")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", default="references")
    args = parser.parse_args()
    root = Path(args.reference_root) / "1000genomes" / "phase3"
    processed = root / "processed"
    lookup = root / "lookup"
    for pop in ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]:
        make_lookup(processed, lookup, pop)


if __name__ == "__main__":
    main()
