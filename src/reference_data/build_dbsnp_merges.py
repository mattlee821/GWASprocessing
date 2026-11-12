#!/usr/bin/env python3
"""Build the dbSNP merge map: retired rsIDs to the current rsID they were merged into.

Usage:
    python3 src/reference_data/build_dbsnp_merges.py --merged refsnp-merged.json.bz2 --out <reference_root>/dbsnp/GRCh38

Studies often report rsIDs that dbSNP has since merged into another (rs12016871 into rs9581854 in build 131). The
GRCh38 index (build_dbsnp_index.py) holds only current rsIDs, so without this map a file that gives rsIDs and no
positions loses those variants. Reads NCBI's refsnp-merged.json.bz2 once, start to finish, and writes to --out:
    merged_from.u32, merged_into.u32   retired rsID number (sorted) and its current rsID number
    merged_meta.json                   counts and source file
Each record gives a retired rsID (refsnp_id), what it was merged into (merged_snapshot_data.merged_into) and rsIDs
merged into it earlier (dbsnp1_merges); all of them map to the same target. Chains (A into B, B later into C) are
followed to the end, and a retired rsID's own record wins over another record's list of earlier merges. Retired
rsIDs merged into more than 1 rsID are left out.
"""

from __future__ import annotations

import argparse
from array import array
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REFSNP = re.compile(rb'"refsnp_id":"(\d+)"')
MERGED_INTO = re.compile(rb'"merged_into":\[([^\]]*)\]')
EARLIER = re.compile(rb'"merged_rsid":"(\d+)"')
CHAIN_LIMIT = 20


def parse(stream) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """(retired rsID, target, primary) per pair: primary = the record's own merge, else an rsID merged into it earlier."""
    sources = array("Q")
    targets = array("Q")  # 8 bytes per entry; Python lists of ints would need about 5 times more
    primary = array("B")
    multiple = 0
    for line in stream:
        rs = REFSNP.search(line)
        into = MERGED_INTO.search(line)
        if not rs or not into:
            continue
        ids = re.findall(rb"\d+", into.group(1))
        if len(ids) != 1:
            multiple += 1
            continue
        target = int(ids[0])
        sources.append(int(rs.group(1)))
        targets.append(target)
        primary.append(1)
        for earlier in EARLIER.findall(line):
            sources.append(int(earlier))
            targets.append(target)
            primary.append(0)
    as_np = lambda a, t: np.frombuffer(a, dtype=t).copy() if len(a) else np.zeros(0, dtype=t)  # noqa: E731
    return as_np(sources, "uint64"), as_np(targets, "uint64"), as_np(primary, "uint8").astype(bool), multiple


def follow(keys: np.ndarray, values: np.ndarray, query: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Follow query through the sorted map keys -> values until it leaves the map; (result, still in the map)."""
    out = query.copy()
    if not len(keys):
        return out, np.zeros(len(out), dtype=bool)
    for _ in range(CHAIN_LIMIT):
        idx = np.minimum(np.searchsorted(keys, out), len(keys) - 1)
        hit = keys[idx] == out
        if not hit.any():
            break
        out = np.where(hit, values[idx], out)
    idx = np.minimum(np.searchsorted(keys, out), len(keys) - 1)
    return out, keys[idx] == out


def resolve(sources: np.ndarray, targets: np.ndarray, primary: np.ndarray | None = None
            ) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """Unique retired rsIDs (sorted) with their current rsID.

    Each record's own merge (primary) is followed through the other records' own merges to the end of its chain
    (A into B, B later into C gives A into C). rsIDs listed only as merged into a record earlier take that record's
    final target. A retired rsID's own record wins over earlier-merge listings; an rsID listed only as an earlier
    merge of 2 records whose chains end differently is dropped."""
    if primary is None:
        primary = np.ones(len(sources), dtype=bool)
    keep = (sources != targets) & (sources != 0) & (targets != 0)  # self-merges carry no information
    sources, targets, primary = sources[keep], targets[keep], primary[keep]
    stats = {"retired_with_conflicting_targets": 0, "earlier_merge_disagreeing_with_own_record": 0,
             "unresolved_chains": 0, "outside_uint32": 0}
    if not len(sources):
        return sources.astype("uint32"), targets.astype("uint32"), stats
    order = np.argsort(sources[primary], kind="stable")
    own_keys, own_values = sources[primary][order], targets[primary][order]
    final, looped = follow(own_keys, own_values, targets)
    stats["unresolved_chains"] = int(looped.sum())
    keep = ~looped & (final != sources)
    sources, final, primary = sources[keep], final[keep], primary[keep]
    order = np.lexsort((final, ~primary, sources))  # per rsID: its own record first
    sources, final, primary = sources[order], final[order], primary[order]
    first = np.r_[True, sources[1:] != sources[:-1]]
    group = np.cumsum(first) - 1
    group_primary = primary[first]
    distinct = np.r_[True, (sources[1:] != sources[:-1]) | (final[1:] != final[:-1])]
    targets_per_group = np.bincount(group, weights=distinct).astype(int)
    conflict = ~group_primary & (targets_per_group > 1)
    stats["retired_with_conflicting_targets"] = int(conflict.sum())
    stats["earlier_merge_disagreeing_with_own_record"] = int((group_primary & (targets_per_group > 1)).sum())
    sources, final = sources[first][~conflict], final[first][~conflict]
    fits = (sources < 2 ** 32) & (final < 2 ** 32)
    stats["outside_uint32"] = int((~fits).sum())
    return sources[fits].astype("uint32"), final[fits].astype("uint32"), stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--merged", required=True, help="refsnp-merged.json.bz2 (or an uncompressed .json)")
    parser.add_argument("--out", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    out = Path(args.out)
    if (out / "merged_meta.json").exists() and not args.force:
        print(f"Skip: {out / 'merged_meta.json'} already built")
        return
    out.mkdir(parents=True, exist_ok=True)
    source = Path(args.merged)
    if source.suffix == ".bz2":
        proc = subprocess.Popen(["bzip2", "-dc", str(source)], stdout=subprocess.PIPE, bufsize=1 << 24)
        sources, targets, primary, multiple = parse(proc.stdout)
        if proc.wait() != 0:
            raise SystemExit(f"bzip2 exited {proc.returncode} reading {source}")
    else:
        with open(source, "rb") as handle:
            sources, targets, primary, multiple = parse(handle)
    records = len(sources)
    merged_from, merged_into, stats = resolve(sources, targets, primary)
    merged_from.astype("<u4").tofile(out / "merged_from.u32")
    merged_into.astype("<u4").tofile(out / "merged_into.u32")
    meta = {"source": source.name, "pairs_read": records, "retired_with_more_than_one_target": multiple,
            "retired_rsids": int(len(merged_from)), **stats, "version": 1}
    (out / "merged_meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(json.dumps(meta), file=sys.stderr)


if __name__ == "__main__":
    main()
