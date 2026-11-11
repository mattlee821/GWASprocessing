#!/usr/bin/env python3
"""Cross-pipeline benchmark: compare standardised outputs of the same GWAS with each other and with the raw study file.

Usage:
    python3 src/benchmark/compare.py --raw RAW.tsv.gz --raw-columns SNP,A1,A2,Freq1.Hapmap,b,se,p \
        --ours ours.gwas.vcf.gz [--opengwas X.vcf.gz] [--gwas2vcf X.vcf.gz] [--catalog X.h.tsv.gz] [--mungesumstats X.tsv.gz] \
        --leads rs1558902,rs13021737 --out-dir DIR [--expected src/benchmark/expected/NAME.tsv]

--raw-columns names the raw file's rsID, effect allele, other allele, effect-allele frequency, effect, SE and P columns.
Variants are matched by rsID (the outputs may be in different builds). For known lead SNPs plus 10,000 rsIDs drawn
from our output (seed 821), each pipeline is compared with the raw study values after aligning to the same effect
allele (either allele order, either strand): present, allele pair, effect sign, |effect|, SE, -log10 P and allele
frequency. Writes summary.tsv (1 row per pipeline), variants.tsv (every sampled value) and, with --expected, our
sampled values as the expected results a new install must reproduce (see check_install.py).
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
TOL = {"effect": 1e-4, "se": 1e-4, "lp": 1e-2, "af": 0.02}


def comp(allele: str) -> str:
    return "".join(COMP.get(b, b) for b in str(allele).upper()[::-1]) if len(str(allele)) > 1 else COMP.get(str(allele).upper(), "")


def read_vcf(path: str) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            rsid = f[2] if f[2].startswith("rs") else fmt.get("ID", ".")
            rows.append((rsid, f[0], int(f[1]), f[4], f[3], fmt.get("ES"), fmt.get("SE"), fmt.get("LP"), fmt.get("AF"), fmt.get("AFS", "")))
    out = pd.DataFrame(rows, columns=["rsid", "chrom", "pos", "ea", "oa", "effect", "se", "lp", "af", "afs"])
    for c in ("effect", "se", "lp", "af"):
        out[c] = pd.to_numeric(out[c].replace(".", np.nan), errors="coerce")
    return out


def lp_from_p(p: pd.Series) -> pd.Series:
    p = pd.to_numeric(p, errors="coerce")
    return -np.log10(p.where(p > 0))


def read_catalog(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, na_values=["NA", ""])
    if "hm_rsid" in df.columns:  # pre-2023 layout with hm_* columns
        out = pd.DataFrame({"rsid": df["hm_rsid"], "chrom": df["hm_chrom"], "pos": pd.to_numeric(df["hm_pos"], errors="coerce"),
                            "ea": df["hm_effect_allele"], "oa": df["hm_other_allele"], "effect": df["hm_beta"], "se": df["standard_error"],
                            "lp": lp_from_p(df["p_value"]), "af": df["hm_effect_allele_frequency"]})
    else:
        out = pd.DataFrame({"rsid": df["rsid"], "chrom": df["chromosome"], "pos": pd.to_numeric(df["base_pair_location"], errors="coerce"),
                            "ea": df["effect_allele"], "oa": df["other_allele"], "effect": df.get("beta"), "se": df["standard_error"],
                            "lp": lp_from_p(df["p_value"]), "af": df.get("effect_allele_frequency")})
    for c in ("effect", "se", "af"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out["afs"] = ""
    return out


def read_mungesumstats(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    out = pd.DataFrame({"rsid": df["SNP"], "chrom": df["CHR"], "pos": pd.to_numeric(df["BP"], errors="coerce"),
                        "ea": df["A2"], "oa": df["A1"], "effect": pd.to_numeric(df["BETA"], errors="coerce"),
                        "se": pd.to_numeric(df["SE"], errors="coerce"), "lp": lp_from_p(df["P"]),
                        "af": pd.to_numeric(df["FRQ"], errors="coerce") if "FRQ" in df.columns else np.nan})
    out["afs"] = ""
    return out


def read_raw(path: str, columns: list[str]) -> pd.DataFrame:
    snp, ea, oa, eaf, beta, se, p = columns
    df = pd.read_csv(path, sep=None, engine="python", dtype=str, usecols=[c for c in columns if c])
    return pd.DataFrame({"rsid": df[snp].str.split(":").str[0], "ea": df[ea].str.upper(), "oa": df[oa].str.upper(),
                         "af": pd.to_numeric(df[eaf], errors="coerce") if eaf else np.nan,
                         "effect": pd.to_numeric(df[beta], errors="coerce"), "se": pd.to_numeric(df[se], errors="coerce"),
                         "lp": lp_from_p(df[p])})


def align(raw: pd.Series, other: pd.Series) -> dict:
    """Compare 1 pipeline's record with the raw record for the same rsID, after matching the effect allele."""
    if other is None:
        return {"present": False}
    r_ea, r_oa, o_ea, o_oa = str(raw.ea), str(raw.oa), str(other.ea).upper(), str(other.oa).upper()
    if (o_ea, o_oa) == (r_ea, r_oa) or (o_ea, o_oa) == (comp(r_ea), comp(r_oa)):
        sign, af = 1.0, other.af
    elif (o_ea, o_oa) == (r_oa, r_ea) or (o_ea, o_oa) == (comp(r_oa), comp(r_ea)):
        sign, af = -1.0, 1 - other.af if pd.notna(other.af) else np.nan
    else:
        return {"present": True, "alleles_match": False}
    effect = sign * other.effect
    return {"present": True, "alleles_match": True,
            "sign_agrees": bool(np.sign(effect) == np.sign(raw.effect)) if raw.effect != 0 else True,
            "effect_equal": bool(abs(effect - raw.effect) <= TOL["effect"]), "se_equal": bool(abs(other.se - raw.se) <= TOL["se"]),
            "lp_equal": bool(abs(other.lp - raw.lp) <= TOL["lp"]) if pd.notna(raw.lp) and pd.notna(other.lp) else None,
            "af_close": bool(abs(af - raw.af) <= TOL["af"]) if pd.notna(af) and pd.notna(raw.af) else None,
            "afs": getattr(other, "afs", "")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw", required=True)
    parser.add_argument("--raw-columns", required=True)
    parser.add_argument("--ours", required=True)
    parser.add_argument("--opengwas", default="")
    parser.add_argument("--gwas2vcf", default="")
    parser.add_argument("--catalog", default="")
    parser.add_argument("--mungesumstats", default="")
    parser.add_argument("--leads", default="")
    parser.add_argument("--sample", type=int, default=10_000)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--expected", default="")
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw = read_raw(args.raw, args.raw_columns.split(",")).drop_duplicates("rsid").set_index("rsid")
    pipelines = {"GWASprocessing": read_vcf(args.ours)}
    for name, path, reader in (("OpenGWAS", args.opengwas, read_vcf), ("gwas2vcf_GRCh38", args.gwas2vcf, read_vcf),
                               ("GWAS_Catalog_harmoniser", args.catalog, read_catalog), ("MungeSumstats", args.mungesumstats, read_mungesumstats)):
        if path:
            pipelines[name] = reader(path)
    ours = pipelines["GWASprocessing"]
    pool = ours.loc[ours["rsid"].str.startswith("rs") & ours["rsid"].isin(raw.index), "rsid"].drop_duplicates()
    sample = pool.sample(n=min(args.sample, len(pool)), random_state=821).tolist()
    leads = [r for r in args.leads.split(",") if r]
    chosen = list(dict.fromkeys(leads + sample))

    records, summary = [], []
    for name, frame in pipelines.items():
        indexed = frame[frame["rsid"].isin(chosen)].drop_duplicates("rsid").set_index("rsid")
        results = []
        for rsid in chosen:
            if rsid not in raw.index:
                continue
            other = indexed.loc[rsid] if rsid in indexed.index else None
            result = {"pipeline": name, "rsid": rsid, "lead": rsid in leads, **align(raw.loc[rsid], other)}
            if other is not None:
                result.update({"chrom": other.chrom, "pos": other.pos, "ea": other.ea, "oa": other.oa, "effect": other.effect,
                               "se": other.se, "lp": other.lp, "af": other.af})
            results.append(result)
        df = pd.DataFrame(results)
        present = df[df["present"]]
        aligned = present[present.get("alleles_match", pd.Series(dtype=bool)).fillna(False).astype(bool)] if len(present) else present
        rate = lambda col: float(aligned[col].dropna().astype(bool).mean()) if len(aligned) and col in aligned and aligned[col].notna().any() else math.nan  # noqa: E731
        summary.append({"pipeline": name, "sampled": len(df), "present": int(df["present"].sum()),
                        "alleles_match": int(len(aligned)), "sign_agrees": rate("sign_agrees"), "effect_equal": rate("effect_equal"),
                        "se_equal": rate("se_equal"), "lp_equal": rate("lp_equal"), "af_close": rate("af_close"),
                        "leads_present": int(df.loc[df["lead"], "present"].sum()), "leads": len(leads)})
        records.append(df)
    pd.DataFrame(summary).to_csv(out_dir / "summary.tsv", sep="\t", index=False, float_format="%.4f")
    pd.concat(records).to_csv(out_dir / "variants.tsv", sep="\t", index=False)
    print(pd.DataFrame(summary).to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    if args.expected:
        expected = ours[ours["rsid"].isin(chosen)].drop_duplicates("rsid").sort_values("rsid")
        Path(args.expected).parent.mkdir(parents=True, exist_ok=True)
        expected.to_csv(args.expected, sep="\t", index=False, float_format="%.8g")
        print(f"expected results: {args.expected} ({len(expected)} variants)")


if __name__ == "__main__":
    main()
