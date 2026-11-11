#!/usr/bin/env python3
"""Check a GWASprocessing install against stored benchmark results.

Usage:
    python3 src/benchmark/check_install.py --expected src/benchmark/expected/locke_2015.tsv --output standardised.gwas.vcf.gz

Run the benchmark GWAS through the pipeline on the new machine (see src/benchmark/README.md), then point this at its
output. Every expected variant must be present with the same GRCh38 position, alleles and values (ES, SE and LP to 6
significant digits, AF and AFS exactly). Exits 1 and lists the first differences if anything differs.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

from compare import read_vcf


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--expected", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    expected = pd.read_csv(args.expected, sep="\t", dtype={"chrom": str, "afs": str}).fillna({"afs": ""})
    observed = read_vcf(args.output)
    merged = expected.merge(observed[observed["rsid"].isin(expected["rsid"])].drop_duplicates("rsid"), on="rsid", how="left",
                            suffixes=("", "_new"))
    problems = []
    for _, row in merged.iterrows():
        if pd.isna(row["pos_new"]):
            problems.append(f"{row.rsid}: missing")
            continue
        for col in ("chrom", "ea", "oa", "afs"):
            if str(row[col]) != str(row[f"{col}_new"]):
                problems.append(f"{row.rsid}: {col} {row[col]} != {row[col + '_new']}")
        if int(row["pos"]) != int(row["pos_new"]):
            problems.append(f"{row.rsid}: pos {row['pos']} != {row['pos_new']}")
        for col in ("effect", "se", "lp", "af"):
            a, b = row[col], row[f"{col}_new"]
            if not (pd.isna(a) and pd.isna(b)) and not np.isclose(a, b, rtol=1e-6, atol=1e-12):
                problems.append(f"{row.rsid}: {col} {a} != {b}")
    print(f"{len(expected)} expected variants, {len(problems)} differences")
    for problem in problems[:20]:
        print("  " + problem)
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
