#!/usr/bin/env python3
"""Compare the same protein across studies, to catch a wrong effect allele or a flipped sign in the standardised output.

For each protein (output name) found in 2 or more studies, the variants that are suggestive (-log10 P >= 6) in any
study are looked up in the others and the effects compared after both are expressed for the same allele (the one that
sorts last). Reported per pair of studies: the correlation of the z scores, the share of variants genome-wide significant
in both that have the same sign, whether each study's lead variant has the same sign in the other, and the correlation of
the allele frequencies at variants that are not palindromic. A flipped sign or a swapped allele shows as a negative
correlation or sign agreement near 0; a correct pair of studies gives a positive correlation and agreement near 1.

Usage:
  python3 pipeline/bin/concordance.py --study ferkingstad=GWAS/ferkingstad_2021_34857953 --study pietzner=GWAS/pietzner_2021_34648354 \\
      [--study zhang=...] [--proteins P03951_2190-55,...] [--jobs 4] [--out concordance.tsv]
Each study folder holds standardise/<output name>/<output name>.vcf.gz. Needs bcftools on the PATH. Reads only the files; nothing is written but --out.
"""

from __future__ import annotations

import argparse
import io
import itertools
import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

LP_PASS = 6.0  # suggestive: variants this strong in any study are compared
LP_SIG = -np.log10(5e-8)  # genome-wide significant
MIN_MATCHED = 10  # matched variants needed for a verdict
PALINDROMIC = {frozenset("AT"), frozenset("CG")}


def read_vcf(path: Path, keep: set[str] | None = None, lp_min: float | None = None) -> pd.DataFrame:
    """Rows of a GWAS-VCF with the effect for the allele that sorts last: rows with -log10 P >= lp_min, or whose position
    is in `keep` (CHROM:POS:LO:HI keys). Read with bcftools, which scans a 25-million-variant file in about a minute.
    Columns: es, se, lp, af, palindromic (indexed by key)."""
    command = ["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t[%ES\t%SE\t%LP\t%AF]\n"]
    targets = None
    if keep is None:
        command += ["-i", f"FORMAT/LP>={lp_min}"]
    else:
        targets = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
        targets.write("".join(f"{k.split(':')[0]}\t{k.split(':')[1]}\n" for k in sorted(keep)))
        targets.close()
        command += ["-T", targets.name]
    try:
        text = subprocess.run([*command, str(path)], check=True, capture_output=True, text=True).stdout
    finally:
        if targets:
            os.unlink(targets.name)
    columns = ["chrom", "pos", "ref", "alt", "es", "se", "lp", "af"]
    raw = pd.read_csv(io.StringIO(text), sep="\t", header=None, names=columns, dtype=str, na_values=["."]) if text else pd.DataFrame(columns=columns)
    last = raw["alt"] >= raw["ref"]  # the effect allele (ALT) is the one that sorts last
    lo, hi = raw["ref"].where(last, raw["alt"]), raw["alt"].where(last, raw["ref"])
    sign = np.where(last, 1.0, -1.0)
    af = pd.to_numeric(raw["af"], errors="coerce")
    frame = pd.DataFrame({
        "key": raw["chrom"] + ":" + raw["pos"] + ":" + lo + ":" + hi,
        "es": pd.to_numeric(raw["es"], errors="coerce") * sign,
        "se": pd.to_numeric(raw["se"], errors="coerce"),
        "lp": pd.to_numeric(raw["lp"], errors="coerce"),
        "af": af.where(last, 1 - af),
        "palindromic": [frozenset((a, b)) in PALINDROMIC for a, b in zip(lo, hi)],
    })
    if keep is not None:
        frame = frame[frame["key"].isin(keep)]
    return frame.drop_duplicates("key").set_index("key")


def compare(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    both = a.join(b, how="inner", lsuffix="_a", rsuffix="_b")
    either = both[(both["lp_a"] >= LP_PASS) | (both["lp_b"] >= LP_PASS)]
    z_a, z_b = either["es_a"] / either["se_a"], either["es_b"] / either["se_b"]
    r_z = float(np.corrcoef(z_a, z_b)[0, 1]) if len(either) >= 3 and z_a.std() > 0 and z_b.std() > 0 else np.nan
    strong = either[(either["lp_a"] >= LP_SIG) & (either["lp_b"] >= LP_SIG)]
    agree = float((np.sign(strong["es_a"]) == np.sign(strong["es_b"])).mean()) if len(strong) else np.nan
    lead = {}
    for tag in ("a", "b"):
        if len(either):
            row = either.loc[either[f"lp_{tag}"].idxmax()]
            lead[tag] = bool(np.sign(row["es_a"]) == np.sign(row["es_b"]))
        else:
            lead[tag] = None
    plain = either[~either["palindromic_a"] & either["af_a"].notna() & either["af_b"].notna()]
    r_af = float(np.corrcoef(plain["af_a"], plain["af_b"])[0, 1]) if len(plain) >= 3 and plain["af_a"].std() > 0 and plain["af_b"].std() > 0 else np.nan
    n = len(either)
    if n < MIN_MATCHED or np.isnan(r_z):
        verdict = "too few shared variants"
    elif r_z <= -0.3 or (not np.isnan(agree) and agree <= 0.1):
        verdict = "SIGN OR ALLELE FLIPPED"
    elif r_z >= 0.5 and (np.isnan(agree) or agree >= 0.9) and lead["a"] is not False and lead["b"] is not False:
        verdict = "ok"
    else:
        verdict = "check"
    return {"n_shared_suggestive": n, "n_both_significant": len(strong), "r_z": r_z, "sign_agree_both_significant": agree,
            "lead_a_same_sign": lead["a"], "lead_b_same_sign": lead["b"], "r_af": r_af, "verdict": verdict}


def protein_pairs(job: tuple[str, dict[str, Path]]) -> list[dict]:
    protein, paths = job
    suggestive = {name: read_vcf(path, lp_min=LP_PASS) for name, path in paths.items()}
    keep = set().union(*(set(frame.index) for frame in suggestive.values()))
    full = {name: read_vcf(path, keep=keep) for name, path in paths.items()}
    return [{"protein": protein, "study_a": x, "study_b": y, **compare(full[x], full[y])} for x, y in itertools.combinations(sorted(full), 2)]


def find_files(studies: dict[str, Path]) -> dict[str, dict[str, Path]]:
    found: dict[str, dict[str, Path]] = {}
    for name, folder in studies.items():
        for path in sorted(Path(folder).glob("standardise/*/*.vcf.gz")):
            found.setdefault(path.parent.name, {})[name] = path
    return {protein: paths for protein, paths in found.items() if len(paths) >= 2}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--study", action="append", required=True, help="NAME=FOLDER, at least 2")
    parser.add_argument("--proteins", help="comma-separated output names (default: every protein in 2 or more studies)")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--out", default="concordance.tsv")
    args = parser.parse_args()
    studies = dict(item.split("=", 1) for item in args.study)
    if len(studies) < 2:
        parser.error("give at least 2 --study NAME=FOLDER")
    found = find_files({name: Path(folder) for name, folder in studies.items()})
    for name in studies:
        count = sum(name in paths for paths in found.values())
        print(f"{name}: {count} proteins shared with at least 1 other study")
        if not count:
            parser.error(f"no standardised files of {name} match another study (looked for {studies[name]}/standardise/<name>/*.vcf.gz)")
    if args.proteins:
        found = {p: found[p] for p in args.proteins.split(",") if p in found}
    if not found:
        print("No protein is in 2 or more of the studies")
        return 1
    with ProcessPoolExecutor(args.jobs) as pool:
        rows = [row for result in pool.map(protein_pairs, sorted(found.items())) for row in result]
    table = pd.DataFrame(rows)
    table.to_csv(args.out, sep="\t", index=False, float_format="%.3f")
    with pd.option_context("display.width", 250, "display.max_columns", 20):
        print(table.round(3).to_string(index=False))
    print("\nverdicts:", table["verdict"].value_counts().to_dict(), "| written", args.out)
    return 0 if not (table["verdict"] == "SIGN OR ALLELE FLIPPED").any() else 2


if __name__ == "__main__":
    raise SystemExit(main())
