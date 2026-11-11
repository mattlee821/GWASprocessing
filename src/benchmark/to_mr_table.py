#!/usr/bin/env python3
"""Convert 1 pipeline's output to the common table the MR script reads (src/benchmark/mr.R).

Usage: python3 src/benchmark/to_mr_table.py --format vcf|catalog|mungesumstats --input FILE --output TABLE.tsv.gz

Columns: SNP (rsID), effect_allele, other_allele, eaf, beta, se, pval (10^-LP for GWAS-VCF; 0 where LP is beyond double
precision), chr, pos. Rows without an rsID are left out, because TwoSampleMR joins exposure and outcome by rsID.
"""

from __future__ import annotations

import argparse

import numpy as np

from compare import read_catalog, read_mungesumstats, read_vcf

READERS = {"vcf": read_vcf, "catalog": read_catalog, "mungesumstats": read_mungesumstats}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--format", required=True, choices=sorted(READERS))
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    df = READERS[args.format](args.input)
    df = df[df["rsid"].astype(str).str.startswith("rs")]
    out = df.rename(columns={"rsid": "SNP", "ea": "effect_allele", "oa": "other_allele", "af": "eaf", "effect": "beta",
                             "chrom": "chr"})
    out["pval"] = np.power(10.0, -out["lp"].astype(float))
    out[["SNP", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "chr", "pos"]].to_csv(
        args.output, sep="\t", index=False, float_format="%.6g")
    print(f"{args.output}: {len(out)} rows")


if __name__ == "__main__":
    main()
