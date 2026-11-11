#!/usr/bin/env python3
"""Build GRCh38 1000 Genomes allele-frequency tables for the 2,504 unrelated phase 3 samples.

Usage:
    python3 src/reference_data/build_1000genomes_grch38.py --raw <dir> --out <reference_root>/1000genomes/GRCh38

<dir> holds the plink2 resource files for the 2022 high-coverage GRCh38 callset (Byrska-Bishop et al.), split by
chromosome (the whole-genome .pvar needs more memory than an 8 GB machine has): by_chromosome/chr<C>_hg38.pgen.zst
and chr<C>_hg38.pvar.zst for C in 1-22 and X, all_hg38.psam (hg38_corrected.psam), and phase3_corrected.psam (the
2,504 phase 3 samples; only these are used, so frequencies are comparable with phase 3). Requires plink2 on PATH
(PLINK2 overrides). Pseudoautosomal regions (PAR1, PAR2) are written as X.

For each super-population (AFR, AMR, EAS, EUR, SAS) and ALL, writes <POP>.key.i64, <POP>.ref.u1, <POP>.alt.u1 and
<POP>.af.f4: SNVs (multi-allelic sites split into 1 row per ALT allele) with 0 < ALT frequency < 1 in that
population, sorted by key = chromosome code * 1e9 + GRCh38 position (codes 1-22, 23 = X, 24 = Y, 25 = MT).
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipeline" / "bin"))
from gwas_reference import BASES, CHROM_CODE, KEY_SCALE  # noqa: E402

POPULATIONS = ("AFR", "AMR", "EAS", "EUR", "SAS", "ALL")
CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X"]  # the callset's chrY/chrM have no rsIDs or are haploid-only


def run(command: list[str]) -> None:
    print(" ".join(command), file=sys.stderr, flush=True)
    subprocess.run(command, check=True)


def convert(afreqs: list[Path], prefix: Path) -> int:
    """plink2 .afreq (cols=chrom,pos,ref,alt,altfreq) to sorted key/ref/alt/af arrays; returns rows kept."""
    keys, refs, alts, afs = [], [], [], []
    for afreq in afreqs:
        _read_one(afreq, keys, refs, alts, afs)
    key, ref, alt, af = (np.concatenate(x) if x else np.zeros(0) for x in (keys, refs, alts, afs))
    order = np.argsort(key, kind="stable")
    key[order].astype("<i8").tofile(f"{prefix}.key.i64")
    ref[order].astype("u1").tofile(f"{prefix}.ref.u1")
    alt[order].astype("u1").tofile(f"{prefix}.alt.u1")
    af[order].astype("<f4").tofile(f"{prefix}.af.f4")
    return int(len(key))


def _read_one(afreq: Path, keys: list, refs: list, alts: list, afs: list) -> None:
    reader = pd.read_csv(afreq, sep="\t", dtype={"#CHROM": str, "POS": "int64", "REF": str, "ALT": str, "ALT_FREQS": "float64"},
                         usecols=["#CHROM", "POS", "REF", "ALT", "ALT_FREQS"], chunksize=5_000_000, na_filter=False)
    for chunk in reader:
        code = chunk["#CHROM"].map(lambda c: CHROM_CODE.get({"23": "X", "24": "Y", "25": "MT", "M": "MT", "PAR1": "X", "PAR2": "X",
                                                             "XY": "X"}.get(c, c), 0))
        ref, alt = chunk["REF"].map(BASES), chunk["ALT"].map(BASES)
        af = chunk["ALT_FREQS"].to_numpy()
        keep = (code > 0).to_numpy() & ref.notna().to_numpy() & alt.notna().to_numpy() & (af > 0) & (af < 1)
        keys.append(code.to_numpy()[keep].astype("int64") * KEY_SCALE + chunk["POS"].to_numpy()[keep])
        refs.append(ref.to_numpy()[keep].astype("uint8"))
        alts.append(alt.to_numpy()[keep].astype("uint8"))
        afs.append(af[keep].astype("float32"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--work", default="", help="scratch folder for the decompressed .pgen (default: <out>/build.tmp)")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    raw, out = Path(args.raw), Path(args.out)
    if (out / "meta.json").exists() and not args.force:
        print(f"Skip: {out} already built")
        return
    out.mkdir(parents=True, exist_ok=True)
    work = Path(args.work) if args.work else out / "build.tmp"
    work.mkdir(parents=True, exist_ok=True)
    plink2 = os.environ.get("PLINK2", "plink2")
    phase3 = set(pd.read_csv(raw / "phase3_corrected.psam", sep="\t", dtype=str)["#IID"])
    samples = pd.read_csv(raw / "all_hg38.psam", sep="\t", dtype=str)
    samples = samples[samples["#IID"].isin(phase3)]
    for pop in POPULATIONS:
        keep = samples if pop == "ALL" else samples[samples["SuperPop"] == pop]
        keep["#IID"].to_csv(work / f"{pop}.keep", index=False, header=False)
    for chrom in CHROMOSOMES:
        if all(Path(work / f"{pop}.chr{chrom}.afreq").exists() for pop in POPULATIONS):
            continue
        prefix = work / f"chr{chrom}_hg38"
        run([plink2, "--zst-decompress", str(raw / "by_chromosome" / f"chr{chrom}_hg38.pgen.zst"), f"{prefix}.pgen"])
        shutil.copy(raw / "by_chromosome" / f"chr{chrom}_hg38.pvar.zst", f"{prefix}.pvar.zst")
        shutil.copy(raw / "all_hg38.psam", f"{prefix}.psam")
        for pop in POPULATIONS:
            run([plink2, "--pfile", str(prefix), "vzs", "--keep", str(work / f"{pop}.keep"), "--snps-only", "just-acgt",
                 "--freq", "cols=chrom,pos,ref,alt,altfreq", "--memory", "2500", "--threads", "3",
                 "--out", str(work / f"{pop}.chr{chrom}")])
        for suffix in (".pgen", ".pvar.zst", ".psam"):
            Path(f"{prefix}{suffix}").unlink(missing_ok=True)
    counts = {}
    for pop in POPULATIONS:
        parts = [work / f"{pop}.chr{chrom}.afreq" for chrom in CHROMOSOMES]
        rows = convert(parts, out / pop)
        counts[pop] = {"samples": int((samples["SuperPop"] == pop).sum() if pop != "ALL" else len(samples)), "variants": rows}
        print(f"{pop}: {counts[pop]}", file=sys.stderr, flush=True)
    shutil.rmtree(work, ignore_errors=True)
    meta = {"callset": "1000 Genomes 30x GRCh38 (Byrska-Bishop et al. 2022), plink2 resource all_hg38",
            "samples": "2,504 unrelated phase 3 samples", "populations": counts, "version": 1}
    (out / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
