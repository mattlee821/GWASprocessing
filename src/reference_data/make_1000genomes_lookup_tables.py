#!/usr/bin/env python3
"""Build the 1000 Genomes lookup tables and the rsID coordinate index used by standardisation.

Usage:
    python3 src/reference_data/make_1000genomes_lookup_tables.py [--reference-root references] [--force] [--skip-coordinates]

Lookup tables, <root>/1000genomes/phase3/lookup/<POP>.tsv.gz, one per super-population plus ALL:
    CHR, SNP, POS19 (GRCh37), A1 (ALT), A2 (REF), AF1 (frequency of A1 in that population)
for the common biallelic SNPs written by build_1000genomes_phase3.sh. AF1 is directional, so the pipeline
can derive the frequency of whichever allele a GWAS reports as its effect allele.

Coordinate index, <root>/1000genomes/phase3/coordinates/: every rsID in the full 1000 Genomes variant list
(all frequencies, indels and chromosome X included) with its chromosome and GRCh37 position, for GWAS
files that report rsIDs without positions. rsIDs listed at more than one position are left out.
Existing outputs are kept unless --force is given.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipeline" / "bin"))
from gwas_coordinates import write_coordinate_index  # noqa: E402

POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]
LOOKUP_COLUMNS = ["CHR", "SNP", "POS19", "A1", "A2", "AF1"]
CHUNK = 2_000_000


def make_lookup(processed_root: Path, out_root: Path, population: str) -> Path:
    bim = processed_root / population / f"{population}.bim"
    stats = processed_root / population / "stats.afreq"
    for path in (bim, stats):
        if not path.exists():
            raise FileNotFoundError(f"{path} is missing; run src/reference_data/build_1000genomes_phase3.sh first")
    variants = pd.read_csv(bim, sep=r"\s+", header=None, names=["CHR", "SNP", "CM", "POS19", "BIM_A1", "BIM_A2"],
                           dtype={"CHR": str, "SNP": str})
    # plink2 --freq: #CHROM ID REF ALT ALT_FREQS OBS_CT. ALT_FREQS is the frequency of ALT, so A1=ALT, AF1=ALT_FREQS.
    freq = pd.read_csv(stats, sep=r"\s+", dtype={"ID": str, "REF": str, "ALT": str})
    freq = freq.rename(columns={"ID": "SNP", "ALT": "A1", "REF": "A2", "ALT_FREQS": "AF1"})[["SNP", "A1", "A2", "AF1"]]
    table = variants.merge(freq, on="SNP", how="inner")[LOOKUP_COLUMNS]
    out_root.mkdir(parents=True, exist_ok=True)
    out = out_root / f"{population}.tsv.gz"
    tmp = out.with_name(f".{out.name}.tmp")
    table.to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, out)
    return out


def open_pvar(path: Path):
    """Text stream of a .pvar, decompressing .zst with zstd (or plink2) on the fly."""
    if path.suffix != ".zst":
        return open(path, "r", encoding="utf-8"), None
    tool = shutil.which("zstd")
    command = [tool, "-dc", str(path)] if tool else [os.environ.get("PLINK2", "plink2"), "--zst-decompress", str(path)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    return process.stdout, process


def make_coordinate_index(pvar_path: Path, out_dir: Path) -> dict:
    handle, process = open_pvar(pvar_path)
    try:
        header = None
        for line in handle:
            if line.startswith("##"):
                continue
            header = line.lstrip("#").rstrip("\n").split("\t")
            break
        if header is None or not {"CHROM", "POS", "ID"}.issubset(header):
            raise ValueError(f"{pvar_path} has no #CHROM POS ID header")
        rs_parts, chrom_parts, pos_parts = [], [], []
        reader = pd.read_csv(handle, sep="\t", header=None, names=header, usecols=["CHROM", "POS", "ID"],
                             dtype={"CHROM": str, "ID": str}, chunksize=CHUNK)
        for chunk in reader:
            numbers = pd.to_numeric(chunk["ID"].str.extract(r"^rs(\d+)(?:$|[:;_])", flags=re.I)[0], errors="coerce")
            keep = numbers.notna().to_numpy()
            rs_parts.append(numbers[keep].to_numpy(dtype="int64"))
            chrom_parts.append(chunk.loc[keep, "CHROM"].to_numpy())
            pos_parts.append(pd.to_numeric(chunk.loc[keep, "POS"], errors="coerce").fillna(0).to_numpy(dtype="int64"))
    finally:
        handle.close()
        if process is not None and process.wait() != 0:
            raise RuntimeError(f"decompressing {pvar_path} failed")
    stats = write_coordinate_index(out_dir, np.concatenate(rs_parts) if rs_parts else np.array([], dtype="int64"),
                                   np.concatenate(chrom_parts) if chrom_parts else np.array([], dtype=object),
                                   np.concatenate(pos_parts) if pos_parts else np.array([], dtype="int64"), source=str(pvar_path))
    return stats


def find_pvar(phase3: Path) -> Path | None:
    for candidate in (phase3 / "processed" / "all_phase3.pvar", phase3 / "raw" / "all_phase3.pvar.zst"):
        if candidate.exists():
            return candidate
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", default="references")
    parser.add_argument("--force", action="store_true", help="Rebuild outputs that already exist.")
    parser.add_argument("--skip-coordinates", action="store_true", help="Do not build the rsID coordinate index.")
    args = parser.parse_args()
    phase3 = Path(args.reference_root) / "1000genomes" / "phase3"
    lookup = phase3 / "lookup"
    for pop in POPULATIONS:
        out = lookup / f"{pop}.tsv.gz"
        if out.exists() and not args.force:
            print(f"Skip: {out} already exists")
            continue
        print(f"Wrote {make_lookup(phase3 / 'processed', lookup, pop)}")
    if args.skip_coordinates:
        return
    coordinates = phase3 / "coordinates"
    if (coordinates / "rs.npy").exists() and not args.force:
        print(f"Skip: {coordinates} already exists")
        return
    pvar = find_pvar(phase3)
    if pvar is None:
        print("WARNING: no 1000 Genomes .pvar found; rsID coordinate index not built "
              "(run src/reference_data/download_1000genomes_phase3.sh)", file=sys.stderr)
        return
    stats = make_coordinate_index(pvar, coordinates)
    print(f"Wrote coordinate index for {stats['variants']} rsIDs to {coordinates}")


if __name__ == "__main__":
    main()
