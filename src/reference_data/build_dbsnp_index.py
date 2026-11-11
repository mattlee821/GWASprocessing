#!/usr/bin/env python3
"""Build the GRCh38 dbSNP index used for rsID-to-position lookups and rsID annotation.

Usage:
    python3 src/reference_data/build_dbsnp_index.py --dbsnp GCF_000001405.40.gz \
        --assembly-report GCF_000001405.40_assembly_report.txt --out <reference_root>/dbsnp/GRCh38

Requires bcftools on PATH. Reads the NCBI dbSNP VCF once, start to finish (it can stay on a slow drive), and
writes to --out, which should be on local disk because the lookups read it at random positions:
    chr<C>.pos.u32, chr<C>.hash.u32, chr<C>.rs.u32   1 row per (variant, ALT allele), sorted by position, for
                                                       C in 1-22, X, Y, MT; hash = allele_hash(REF, ALT)
    rs.u32, rs_chrom.u8, rs_pos.u32                    rsIDs sorted, with chromosome code and position; rsIDs at more
                                                       than 1 position are left out
    meta.json                                          counts, dbSNP build, source file
Chromosome codes: 1-22, 23 = X, 24 = Y, 25 = MT. Alleles longer than 50 bases are skipped.
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
from gwas_reference import CHROM_CODE, allele_hash  # noqa: E402

MAX_ALLELE = 50
RS_BUCKETS = 64
CHUNK = 5_000_000


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


def dbsnp_build(path: Path) -> str:
    header = subprocess.run(["bcftools", "view", "-h", str(path)], capture_output=True, text=True, check=True).stdout
    for line in header.splitlines():
        if line.startswith("##dbSNP_BUILD_ID="):
            return line.split("=", 1)[1]
    return ""


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
    work = out / "build.tmp"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    contigs = refseq_to_chrom(Path(args.assembly_report))
    chrom_files = {c: {kind: open(work / f"chr{c}.{kind}", "wb") for kind in ("pos", "hash", "rs")} for c in CHROM_CODE}
    bucket_files = [open(work / f"rs_bucket_{i:03d}.bin", "wb") for i in range(RS_BUCKETS)]
    bucket_dtype = np.dtype([("rs", "<u4"), ("chrom", "u1"), ("pos", "<u4")])
    rs_span = 2 ** 32 // RS_BUCKETS
    command = ["bcftools", "query", "-f", "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n", args.dbsnp]
    rows = skipped = 0
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as proc:
        reader = pd.read_csv(proc.stdout, sep="\t", header=None, names=["CHROM", "POS", "ID", "REF", "ALT"],
                             dtype={"CHROM": str, "POS": "uint32", "ID": str, "REF": str, "ALT": str}, chunksize=CHUNK,
                             na_filter=False, engine="c")
        for chunk in reader:
            chunk["CHR"] = chunk["CHROM"].map(contigs)
            chunk = chunk[chunk["CHR"].notna() & chunk["ID"].str.startswith("rs")]
            if chunk["ALT"].str.contains(",", regex=False).any():
                chunk = chunk.assign(ALT=chunk["ALT"].str.split(",")).explode("ALT")
            keep = (chunk["REF"].str.len() <= MAX_ALLELE) & (chunk["ALT"].str.len() <= MAX_ALLELE) & (chunk["ALT"] != ".")
            skipped += int((~keep).sum())
            chunk = chunk[keep]
            rs = pd.to_numeric(chunk["ID"].str.slice(2), errors="coerce")
            chunk = chunk[rs.notna() & (rs < 2 ** 32)]
            rs = rs[rs.notna() & (rs < 2 ** 32)].to_numpy(dtype="uint32")
            hashes = allele_hash(chunk["REF"].to_numpy(dtype=object), chunk["ALT"].to_numpy(dtype=object))
            pos = chunk["POS"].to_numpy(dtype="uint32")
            chrom = chunk["CHR"].to_numpy(dtype=object)
            for c in pd.unique(chrom):
                sel = chrom == c
                files = chrom_files[c]
                files["pos"].write(pos[sel].astype("<u4").tobytes())
                files["hash"].write(hashes[sel].astype("<u4").tobytes())
                files["rs"].write(rs[sel].astype("<u4").tobytes())
            codes = np.array([CHROM_CODE[c] for c in chrom], dtype="u1")
            records = np.empty(len(rs), dtype=bucket_dtype)
            records["rs"], records["chrom"], records["pos"] = rs, codes, pos
            bucket = (rs // rs_span).astype(int)
            for b in np.unique(bucket):
                bucket_files[b].write(records[bucket == b].tobytes())
            rows += len(rs)
            print(f"{rows:,} rows", file=sys.stderr, flush=True)
    if proc.returncode != 0:
        raise SystemExit(f"bcftools query exited {proc.returncode}")
    for files in chrom_files.values():
        for handle in files.values():
            handle.close()
    for handle in bucket_files:
        handle.close()
    counts = {}
    for c in CHROM_CODE:
        pos = np.fromfile(work / f"chr{c}.pos", dtype="<u4")
        order = None if np.all(pos[1:] >= pos[:-1]) else np.argsort(pos, kind="stable")
        for kind in ("pos", "hash", "rs"):
            data = np.fromfile(work / f"chr{c}.{kind}", dtype="<u4")
            (data if order is None else data[order]).tofile(out / f"chr{c}.{kind}.u32")
            os.remove(work / f"chr{c}.{kind}")
        counts[c] = int(len(pos))
    total_rs = ambiguous = 0
    with open(out / "rs.u32", "wb") as rs_out, open(out / "rs_chrom.u8", "wb") as chrom_out, open(out / "rs_pos.u32", "wb") as pos_out:
        for i in range(RS_BUCKETS):
            records = np.fromfile(work / f"rs_bucket_{i:03d}.bin", dtype=bucket_dtype)
            os.remove(work / f"rs_bucket_{i:03d}.bin")
            if not len(records):
                continue
            records = np.unique(records)  # sorts by rs, chrom, pos; drops multi-allelic repeats
            first = np.r_[True, records["rs"][1:] != records["rs"][:-1]]
            last = np.r_[records["rs"][1:] != records["rs"][:-1], True]
            single = first & last
            ambiguous += int(first.sum() - single.sum())
            records = records[single]
            rs_out.write(records["rs"].astype("<u4").tobytes())
            chrom_out.write(records["chrom"].tobytes())
            pos_out.write(records["pos"].astype("<u4").tobytes())
            total_rs += len(records)
    shutil.rmtree(work, ignore_errors=True)
    meta = {"source": Path(args.dbsnp).name, "dbsnp_build": dbsnp_build(Path(args.dbsnp)), "assembly": "GRCh38",
            "rows_by_chromosome": counts, "rows": rows, "skipped_long_alleles": skipped, "rsids": total_rs,
            "rsids_at_more_than_one_position": ambiguous, "version": 1}
    (out / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta), file=sys.stderr)


if __name__ == "__main__":
    main()
