#!/usr/bin/env python3
"""Build the dbSNP indel allele table (GRCh38) that gives I/D-coded indels their sequences.

Usage:
    python3 src/reference_data/build_dbsnp_indels.py --dbsnp GCF_000001405.40.gz \
        --assembly-report GCF_000001405.40_assembly_report.txt --out <reference_root>/dbsnp/GRCh38/indels

Requires bcftools on PATH. Reads the NCBI dbSNP VCF once, start to finish (it can stay on a slow drive), keeping every
record with an rsID whose REF and ALT differ in length (a multi-ALT record gives 1 row per ALT; alleles over 50 bases
are skipped, as in the main index). Writes to --out, which should be on local disk:
    rs.u32, chrom.u8, pos.u32   1 row per distinct (rsID, chromosome, position, REF, ALT), sorted by rsID
    off.u64, alleles.bin        n + 1 byte offsets, and the text "REF>ALT" of each row joined with no separator
    meta.json                   counts, dbSNP build, source file
Chromosome codes as in build_dbsnp_index.py (1-22, 23 = X, 24 = Y, 25 = MT).
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
from gwas_reference import CHROM_CODE  # noqa: E402

MAX_ALLELE = 50
RS_BUCKETS = 64
CHUNK = 2_000_000
VERSION = 1


def refseq_to_chrom(report: Path) -> dict[str, str]:
    mapping = {}
    for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) > 6 and fields[1] == "assembled-molecule":
            chrom = {"MT": "MT"}.get(fields[0], fields[0])
            if chrom in CHROM_CODE:
                mapping[fields[6]] = chrom
    return mapping


def indel_rows(chunk: pd.DataFrame, contigs: dict[str, str]) -> pd.DataFrame:
    """The (rs, chrom, pos, ref, alt) rows of a chunk of `bcftools query` output that are indels with an rsID."""
    chunk = chunk.assign(CHR=chunk["CHROM"].map(contigs))
    chunk = chunk[chunk["CHR"].notna() & chunk["ID"].str.startswith("rs")]
    if chunk["ALT"].str.contains(",", regex=False).any():
        chunk = chunk.assign(ALT=chunk["ALT"].str.split(",")).explode("ALT")
    ref_len, alt_len = chunk["REF"].str.len(), chunk["ALT"].str.len()
    chunk = chunk[(ref_len != alt_len) & (ref_len <= MAX_ALLELE) & (alt_len <= MAX_ALLELE) & (chunk["ALT"] != ".") & ~chunk["ALT"].str.startswith("<")]
    rs = pd.to_numeric(chunk["ID"].str.slice(2), errors="coerce")
    keep = rs.notna() & (rs < 2 ** 32)
    return pd.DataFrame({"rs": rs[keep].astype("uint32"), "chrom": chunk["CHR"][keep].map(CHROM_CODE).astype("uint8"),
                         "pos": chunk["POS"][keep].astype("uint32"), "ref": chunk["REF"][keep], "alt": chunk["ALT"][keep]})


def write_index(buckets: list[Path], out: Path) -> dict:
    """Sort each bucket of TSV rows (rs, chrom, pos, ref, alt) by rsID, drop repeats and append to the index files."""
    rows = 0
    offset = 0
    with open(out / "rs.u32", "wb") as rs_out, open(out / "chrom.u8", "wb") as chrom_out, open(out / "pos.u32", "wb") as pos_out, \
            open(out / "off.u64", "wb") as off_out, open(out / "alleles.bin", "wb") as blob_out:
        for bucket in buckets:
            if not bucket.exists() or bucket.stat().st_size == 0:
                continue
            frame = pd.read_csv(bucket, sep="\t", header=None, names=["rs", "chrom", "pos", "ref", "alt"],
                                dtype={"rs": "uint32", "chrom": "uint8", "pos": "uint32", "ref": str, "alt": str}, na_filter=False)
            frame = frame.drop_duplicates().sort_values(["rs", "chrom", "pos", "ref", "alt"], kind="stable")
            text = (frame["ref"] + ">" + frame["alt"]).to_numpy(dtype=object)
            lengths = np.fromiter((len(t) for t in text), dtype="int64", count=len(text))
            off_out.write((offset + np.cumsum(lengths) - lengths).astype("<u8").tobytes())
            blob_out.write("".join(text).encode("ascii"))
            offset += int(lengths.sum())
            rs_out.write(frame["rs"].to_numpy().astype("<u4").tobytes())
            chrom_out.write(frame["chrom"].to_numpy().astype("u1").tobytes())
            pos_out.write(frame["pos"].to_numpy().astype("<u4").tobytes())
            rows += len(frame)
            os.remove(bucket)
        off_out.write(np.array([offset], dtype="<u8").tobytes())
    return {"rows": rows}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dbsnp", required=True)
    parser.add_argument("--assembly-report", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    out = Path(args.out)
    if (out / "meta.json").exists() and not args.force:
        print(f"Skip: {out} already built")
        return
    shutil.rmtree(out, ignore_errors=True)
    work = out / "build.tmp"
    work.mkdir(parents=True)
    contigs = refseq_to_chrom(Path(args.assembly_report))
    buckets = [work / f"bucket_{i:03d}.tsv" for i in range(RS_BUCKETS)]
    handles = [open(b, "w", encoding="ascii") for b in buckets]
    span = 2 ** 32 // RS_BUCKETS
    command = ["bcftools", "query", "-i", 'TYPE~"indel"', "-f", "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n", args.dbsnp]
    seen = 0
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as proc:
        reader = pd.read_csv(proc.stdout, sep="\t", header=None, names=["CHROM", "POS", "ID", "REF", "ALT"],
                             dtype={"CHROM": str, "POS": "uint32", "ID": str, "REF": str, "ALT": str}, chunksize=CHUNK,
                             na_filter=False, engine="c")
        for chunk in reader:
            rows = indel_rows(chunk, contigs)
            bucket = (rows["rs"].to_numpy(dtype="int64") // span)
            for b in np.unique(bucket):
                rows[bucket == b].to_csv(handles[b], sep="\t", header=False, index=False)
            seen += len(rows)
            print(f"{seen:,} indel rows", file=sys.stderr, flush=True)
    for handle in handles:
        handle.close()
    if proc.returncode != 0:
        raise SystemExit(f"bcftools query exited {proc.returncode}")
    counts = write_index(buckets, out)
    shutil.rmtree(work, ignore_errors=True)
    header = subprocess.run(["bcftools", "view", "-h", args.dbsnp], capture_output=True, text=True, check=True).stdout
    build = next((l.split("=", 1)[1] for l in header.splitlines() if l.startswith("##dbSNP_BUILD_ID=")), "")
    meta = {"source": Path(args.dbsnp).name, "dbsnp_build": build, "assembly": "GRCh38", "rows": counts["rows"],
            "rows_read": seen, "max_allele": MAX_ALLELE, "version": VERSION}
    (out / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta), file=sys.stderr)


if __name__ == "__main__":
    main()
