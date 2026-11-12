#!/usr/bin/env python3
"""Build the paired GRCh37/GRCh38 SNP sample used to detect an input's genome build.

Usage:
    python3 src/reference_data/build_build_detection_sample.py --reference-root <root> [--size 2000000]

Takes common biallelic SNPs (ALT frequency 0.05-0.95, the only record at their position) from the 1000 Genomes
GRCh38 EUR frequencies (<root>/1000genomes/GRCh38), draws a fixed-seed (821) sample, and lifts it to GRCh37 with
BCFtools/liftover (GRCh38 and GRCh37 FASTAs in <root>/genome, UCSC hg38ToHg19 chain). SNPs that lift to the same
chromosome as 1 biallelic record are kept, so the same SNPs are placed in each build and neither build is favoured.
Up to --size SNPs are written to <root>/build_detection/ as grch37.{key.i64,ref.u1,alt.u1} and
grch38.{key.i64,ref.u1,alt.u1}, each sorted by key (chromosome code * 1e9 + position). Needs bcftools with the
+liftover plugin (GWAS_BCFTOOLS, BCFTOOLS_PLUGINS).
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "pipeline" / "bin"))
import gwas_alleles as alleles  # noqa: E402
import gwas_reference as ref  # noqa: E402
from standardise_gwas import bcftools_liftover  # noqa: E402

COMMON = (0.05, 0.95)
DRAW_MARGIN = 1.1  # draw 10% more than --size: a few SNPs do not lift cleanly
LETTERS = np.array(["", "A", "C", "G", "T"], dtype=object)


def common_snps(freqs: ref.Frequencies) -> np.ndarray:
    """Indices of common SNPs that are the only 1000 Genomes record at their position."""
    key = np.asarray(freqs.key)
    af = np.asarray(freqs.af, dtype="float64")
    single = np.r_[True, key[1:] != key[:-1]] & np.r_[key[1:] != key[:-1], True]
    return np.where(single & (af >= COMMON[0]) & (af <= COMMON[1]))[0]


def lift_to_grch37(chrom: np.ndarray, pos: np.ndarray, ref38: np.ndarray, alt38: np.ndarray, root: Path, tool: str
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """(row index, GRCh37 chromosome, GRCh37 position, REF, ALT) for each SNP lifted to 1 biallelic GRCh37 record."""
    fasta38, fasta37 = alleles.Fasta.find(root), alleles.Fasta.find(root, "GRCh37")
    chain = root / "genome" / "hg38ToHg19.over.chain.gz"
    if fasta38 is None or fasta37 is None or not chain.exists():
        raise SystemExit(f"need the GRCh38 and GRCh37 FASTAs and {chain} (run download_fasta_and_chains.sh)")
    work = Path(tempfile.mkdtemp(prefix="build_detection_"))
    try:
        vcf = work / "grch38.vcf"
        contigs = [f"##contig=<ID={name},length={length}>" for name, (length, *_rest) in fasta38.index.items()]
        with open(vcf, "w", encoding="utf-8") as handle:
            handle.write("\n".join(["##fileformat=VCFv4.2", *contigs, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]) + "\n")
            handle.writelines(f"{c}\t{p}\t{i}\t{r}\t{a}\t.\t.\t.\n" for i, (c, p, r, a) in enumerate(zip(chrom, pos, ref38, alt38)))
        command = [tool, "+liftover", "--no-version", str(vcf), "--", "-s", str(fasta38.path), "-f", str(fasta37.path),
                   "-c", str(chain)]
        lifted = subprocess.run(command, capture_output=True, text=True, check=True).stdout
    finally:
        shutil.rmtree(work, ignore_errors=True)
    rows, chroms, positions, refs, alts = [], [], [], [], []
    for line in lifted.splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t", 5)
        i = int(f[2])
        if f[0] != chrom[i] or len(f[3]) != 1 or len(f[4]) != 1:  # another chromosome, or a new reference allele added
            continue
        rows.append(i)
        chroms.append(f[0])
        positions.append(int(f[1]))
        refs.append(f[3])
        alts.append(f[4])
    return (np.array(rows, dtype="int64"), np.array(chroms, dtype=object), np.array(positions, dtype="int64"),
            np.array(refs, dtype=object), np.array(alts, dtype=object))


def write_build(out: Path, name: str, keys: np.ndarray, ref_codes: np.ndarray, alt_codes: np.ndarray) -> None:
    order = np.argsort(keys, kind="stable")
    keys[order].astype("<i8").tofile(out / f"{name}.key.i64")
    ref_codes[order].astype("u1").tofile(out / f"{name}.ref.u1")
    alt_codes[order].astype("u1").tofile(out / f"{name}.alt.u1")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--size", type=int, default=2_000_000)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    root = Path(args.reference_root)
    out = root / "build_detection"
    if (out / "meta.json").exists() and not args.force:
        print(f"Skip: {out} already built")
        return
    freqs = ref.Frequencies.find(root, "EUR")
    if freqs is None:
        raise SystemExit(f"no 1000 Genomes GRCh38 frequencies under {root}/1000genomes/GRCh38; run build_1000genomes_grch38.py first")
    tool = bcftools_liftover()
    if not tool:
        raise SystemExit("bcftools with the +liftover plugin not found (set GWAS_BCFTOOLS and BCFTOOLS_PLUGINS)")
    candidates = common_snps(freqs)
    draw = min(len(candidates), int(args.size * DRAW_MARGIN))
    drawn = np.sort(np.random.default_rng(821).choice(candidates, size=draw, replace=False))
    key38 = np.asarray(freqs.key)[drawn]
    code = key38 // ref.KEY_SCALE
    chrom = np.array([ref.CHROMS[c - 1] for c in code], dtype=object)
    pos38 = key38 % ref.KEY_SCALE
    ref38 = LETTERS[np.asarray(freqs.ref)[drawn]]
    alt38 = LETTERS[np.asarray(freqs.alt)[drawn]]
    rows, chrom37, pos37, ref37, alt37 = lift_to_grch37(chrom, pos38, ref38, alt38, root, tool)
    if len(rows) > args.size:
        keep = np.sort(np.random.default_rng(821).choice(len(rows), size=args.size, replace=False))
        rows, chrom37, pos37, ref37, alt37 = rows[keep], chrom37[keep], pos37[keep], ref37[keep], alt37[keep]
    out.mkdir(parents=True, exist_ok=True)
    write_build(out, "grch38", key38[rows].copy(), ref.base_codes(ref38[rows]), ref.base_codes(alt38[rows]))
    write_build(out, "grch37", ref.chrom_codes(chrom37) * ref.KEY_SCALE + pos37, ref.base_codes(ref37), ref.base_codes(alt37))
    meta = {"snps": int(len(rows)), "drawn": int(draw), "common_candidates": int(len(candidates)),
            "source": f"1000 Genomes GRCh38 {freqs.population} common biallelic SNPs, lifted to GRCh37 with BCFtools/liftover",
            "seed": 821, "version": 2}
    (out / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta))


if __name__ == "__main__":
    main()
