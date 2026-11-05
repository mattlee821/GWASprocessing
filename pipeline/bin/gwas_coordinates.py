#!/usr/bin/env python3
"""Coordinates: strand- and chromosome-aware liftover, genome-build detection, rsID coordinate index,
and chr:pos variant-ID parsing.

Liftover outcomes are explicit: a variant that lands on another chromosome or contig, maps to more
than one place, or does not map is reported as such and never given a fabricated position. Minus-strand
mappings are reported so callers can complement alleles.
"""

from __future__ import annotations

import json
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gwas_pipeline_common import (
    PRIMARY_CHROMS,
    alleles_match,
    normalize_chromosome,
    normalize_chromosome_series,
    numeric_series,
    upper_alleles,
    without_appledouble,
)

CHROM_CODES = {chrom: code for code, chrom in enumerate(PRIMARY_CHROMS, start=1)}
CODE_TO_CHROM = {code: chrom for chrom, code in CHROM_CODES.items()}
KEY_FACTOR = 10_000_000_000  # position keys: chromosome code * 1e10 + position (fits int64)


def position_key(chrom: Any, pos: Any) -> np.ndarray:
    """int64 key per (chromosome, position); -1 where either is missing or not a primary chromosome."""
    codes = normalize_chromosome_series(pd.Series(chrom).reset_index(drop=True)).map(CHROM_CODES).astype("float64").to_numpy()
    positions = numeric_series(pd.Series(pos).reset_index(drop=True)).to_numpy(dtype="float64")
    valid = np.isfinite(codes) & np.isfinite(positions) & (positions > 0)
    key = np.full(len(positions), -1, dtype="int64")
    key[valid] = codes[valid].astype("int64") * KEY_FACTOR + positions[valid].astype("int64")
    return key
DETECTION_SAMPLE = 20_000
DETECTION_SEED = 821
_LIFTER_CACHE: dict[str, Any] = {}


def ucsc_name(chrom: str) -> str:
    return "chrM" if chrom == "MT" else f"chr{chrom}"


def find_chain_file(reference_root: str | Path, source: str, target: str) -> Path | None:
    """Locate e.g. <reference_root>/genome/hg19ToHg38.over.chain.gz (case-insensitive name match)."""
    genome = Path(reference_root) / "genome"
    if not genome.is_dir():
        return None
    wanted = f"{source}to{target}".lower()
    for path in without_appledouble(sorted(genome.glob("*.chain*"))):
        if wanted in path.name.lower().replace("_", "").replace("-", ""):
            return path
    return None


def load_liftover(reference_root: str | Path, source: str, target: str) -> tuple[Any, str]:
    """Return (pyliftover.LiftOver or None, reason it is unavailable)."""
    chain = find_chain_file(reference_root, source, target)
    if chain is None:
        return None, f"liftover chain {source}->{target} not found under {Path(reference_root) / 'genome'}"
    key = str(chain.resolve())
    if key not in _LIFTER_CACHE:
        from pyliftover import LiftOver

        _LIFTER_CACHE[key] = LiftOver(key)
    return _LIFTER_CACHE[key], ""


@dataclass
class LiftResult:
    pos: np.ndarray      # 1-based lifted position, NaN unless status == "ok"
    strand: np.ndarray   # "+", "-", or "" when not lifted
    status: np.ndarray   # "ok", "unmapped", "multiple", "other_chrom", "missing"

    @property
    def ok(self) -> np.ndarray:
        return self.status == "ok"

    def counts(self) -> dict[str, int]:
        values, counts = np.unique(self.status, return_counts=True)
        return {str(v): int(c) for v, c in zip(values, counts)}


def lift(chrom: Any, pos: Any, lifter: Any) -> LiftResult:
    """Lift 1-based positions; only same-chromosome, single, primary-contig mappings count as ok."""
    chroms = normalize_chromosome_series(pd.Series(chrom).reset_index(drop=True))
    positions = numeric_series(pd.Series(pos).reset_index(drop=True)).to_numpy(dtype="float64")
    n = len(positions)
    out_pos = np.full(n, np.nan)
    out_strand = np.full(n, "", dtype=object)
    out_status = np.full(n, "missing", dtype=object)
    cache: dict[tuple[str, int], tuple[str, float, str]] = {}
    for i in range(n):
        c, p = chroms.iat[i], positions[i]
        if pd.isna(c) or not math.isfinite(p) or p < 1 or str(c) not in CHROM_CODES:
            continue
        key = (str(c), int(p))
        if key not in cache:
            hits = lifter.convert_coordinate(ucsc_name(key[0]), key[1] - 1)
            if not hits:
                cache[key] = ("unmapped", math.nan, "")
            elif len(hits) > 1:
                cache[key] = ("multiple", math.nan, "")
            elif normalize_chromosome(hits[0][0]) != key[0]:
                cache[key] = ("other_chrom", math.nan, "")
            else:
                cache[key] = ("ok", float(hits[0][1] + 1), hits[0][2])
        out_status[i], out_pos[i], out_strand[i] = cache[key]
    return LiftResult(pos=out_pos, strand=out_strand, status=out_status)


def lookup_keys(lookup: pd.DataFrame) -> pd.DataFrame:
    """KEY/A1/A2 view of a reference lookup (uses a precomputed KEY column when present)."""
    key = lookup["KEY"].to_numpy(dtype="int64") if "KEY" in lookup.columns else position_key(lookup["CHR"], lookup["POS19"])
    ref = pd.DataFrame({"KEY": key, "A1": lookup["A1"].to_numpy(dtype=object), "A2": lookup["A2"].to_numpy(dtype=object)})
    return ref[ref["KEY"] >= 0]


def _positional_hits(chrom: pd.Series, pos: Any, a: pd.Series, b: pd.Series, reference: pd.DataFrame) -> int:
    probe = pd.DataFrame({"_row": np.arange(len(chrom)), "KEY": position_key(chrom, pos),
                          "A": upper_alleles(pd.Series(a)).to_numpy(), "B": upper_alleles(pd.Series(b)).to_numpy()})
    probe = probe[probe["KEY"] >= 0]
    merged = probe.merge(reference[reference["KEY"].isin(probe["KEY"])], on="KEY", how="inner")
    if merged.empty:
        return 0
    matched = alleles_match(merged["A"], merged["B"], merged["A1"], merged["A2"])
    return int(merged.loc[matched, "_row"].nunique())


def detect_build(chrom: pd.Series, pos: pd.Series, a: pd.Series, b: pd.Series, lookup: pd.DataFrame,
                 lifter_38_to_19: Any, sample_size: int = DETECTION_SAMPLE, seed: int = DETECTION_SEED) -> dict[str, Any]:
    """Decide GRCh37 vs GRCh38 by matching (chromosome, position, alleles) to the GRCh37 reference.

    The GRCh38 hypothesis is tested by lifting input positions to GRCh37 first. A build is called when
    its hits reach max(min(5, n), 5% of n) and at least ten times the other build's hits.
    """
    frame = pd.DataFrame({"CHR": normalize_chromosome_series(pd.Series(chrom).reset_index(drop=True)),
                          "POS": numeric_series(pd.Series(pos).reset_index(drop=True)),
                          "A": upper_alleles(pd.Series(a).reset_index(drop=True)),
                          "B": upper_alleles(pd.Series(b).reset_index(drop=True))})
    valid = frame["CHR"].astype(str).isin(CHROM_CODES) & np.isfinite(frame["POS"]) & (frame["POS"] > 0) \
        & (frame["A"] != "") & (frame["B"] != "")
    frame = frame[valid.to_numpy()].reset_index(drop=True)
    if len(frame) > sample_size:
        frame = frame.sample(n=sample_size, random_state=seed).reset_index(drop=True)
    n = len(frame)
    result: dict[str, Any] = {"build": "", "hits_hg19": 0, "hits_hg38": 0, "sampled": n,
                              "tested_hg38": lifter_38_to_19 is not None}
    if n == 0 or lookup is None or lookup.empty or not {"CHR", "POS19", "A1", "A2"}.issubset(lookup.columns):
        result["reason"] = "no positions or no reference with alleles"
        return result
    reference = lookup_keys(lookup)
    result["hits_hg19"] = _positional_hits(frame["CHR"], frame["POS"], frame["A"], frame["B"], reference)
    if lifter_38_to_19 is not None:
        lifted = lift(frame["CHR"], frame["POS"], lifter_38_to_19)
        ok = lifted.ok
        result["hits_hg38"] = _positional_hits(frame.loc[ok, "CHR"], lifted.pos[ok], frame.loc[ok, "A"],
                                               frame.loc[ok, "B"], reference)
    min_hits = max(min(5, n), math.ceil(0.05 * n))
    h19, h38 = result["hits_hg19"], result["hits_hg38"]
    if h19 >= min_hits and h19 >= 10 * h38:
        result["build"] = "hg19"
    elif lifter_38_to_19 is not None and h38 >= min_hits and h38 >= 10 * h19:
        result["build"] = "hg38"
    result["min_hits"] = min_hits
    return result


_VARIANT_ID = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|X|Y|MT|M|XY)[:_](?P<pos>[0-9]+)(?:[:_](?P<a1>[ACGTN]+)[:_](?P<a2>[ACGTN]+))?$",
    flags=re.I,
)


def parse_variant_ids(ids: pd.Series) -> pd.DataFrame:
    """Parse chr:pos[:a1:a2] / chr_pos_a1_a2 identifiers into CHR and POS (missing where unparseable)."""
    text = pd.Series(ids).astype("string").str.strip()
    parts = text.str.extract(_VARIANT_ID)
    chrom = parts["chrom"].map(lambda v: normalize_chromosome(v) if isinstance(v, str) else pd.NA)
    return pd.DataFrame({"CHR": chrom.astype(object), "POS": numeric_series(parts["pos"]).astype("float64")}, index=text.index)


def rs_number(ids: pd.Series) -> np.ndarray:
    """Numeric part of bare rsIDs (rs123 -> 123); -1 for anything else."""
    digits = pd.Series(ids).astype("string").str.strip().str.extract(r"^[rR][sS](\d+)$")[0]
    return numeric_series(digits).fillna(-1).astype("int64").to_numpy()


def coordinate_index_dir(reference_root: str | Path) -> Path:
    return Path(reference_root) / "1000genomes" / "phase3" / "coordinates"


def _save_npy(path: Path, array: np.ndarray) -> None:
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with open(tmp, "wb") as handle:
        np.save(handle, array)
    os.replace(tmp, path)


def write_coordinate_index(out_dir: str | Path, rs: Any, chrom: Any, pos: Any, source: str = "") -> dict[str, int]:
    """Write a sorted rsID -> (chromosome, GRCh37 position) index; rsIDs with conflicting positions are dropped."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rs = np.asarray(rs, dtype="int64")
    codes = pd.Series(chrom).map(lambda c: CHROM_CODES.get(str(normalize_chromosome(c)), 0)).to_numpy(dtype="int64")
    pos = np.asarray(pos, dtype="int64")
    keep = (rs >= 0) & (codes > 0) & (pos > 0)
    rs, codes, pos = rs[keep], codes[keep], pos[keep]
    order = np.lexsort((pos, codes, rs))
    rs, codes, pos = rs[order], codes[order], pos[order]
    if len(rs):
        exact_dup = np.zeros(len(rs), dtype=bool)
        exact_dup[1:] = (rs[1:] == rs[:-1]) & (codes[1:] == codes[:-1]) & (pos[1:] == pos[:-1])
        rs, codes, pos = rs[~exact_dup], codes[~exact_dup], pos[~exact_dup]
    conflict = np.zeros(len(rs), dtype=bool)
    if len(rs) > 1:
        same = rs[1:] == rs[:-1]
        conflict[1:] |= same
        conflict[:-1] |= same
    stats = {"variants": int((~conflict).sum()), "rsids_with_conflicting_positions": int(len(np.unique(rs[conflict])))}
    _save_npy(out_dir / "rs.npy", rs[~conflict])
    _save_npy(out_dir / "chrom.npy", codes[~conflict].astype("int8"))
    _save_npy(out_dir / "pos.npy", pos[~conflict].astype("int32"))
    (out_dir / "meta.json").write_text(json.dumps({**stats, "source": source}, indent=2) + "\n", encoding="utf-8")
    return stats


class CoordinateIndex:
    """Memory-mapped rsID -> (CHR, GRCh37 POS) index built from the full 1000 Genomes variant list."""

    def __init__(self, rs: np.ndarray, chrom: np.ndarray, pos: np.ndarray, path: Path):
        self.rs, self.chrom, self.pos, self.path = rs, chrom, pos, path

    @classmethod
    def load(cls, directory: str | Path) -> "CoordinateIndex | None":
        directory = Path(directory)
        if not all((directory / name).exists() for name in ("rs.npy", "chrom.npy", "pos.npy")):
            return None
        load = lambda name: np.load(directory / name, mmap_mode="r")  # noqa: E731
        return cls(load("rs.npy"), load("chrom.npy"), load("pos.npy"), directory)

    @property
    def size(self) -> int:
        return int(len(self.rs))

    def lookup(self, rsids: pd.Series) -> pd.DataFrame:
        numbers = rs_number(rsids)
        idx = np.searchsorted(self.rs, numbers)
        idx_clipped = np.minimum(idx, max(self.size - 1, 0))
        found = (numbers >= 0) & (idx < self.size) & (np.asarray(self.rs)[idx_clipped] == numbers) if self.size else np.zeros(len(numbers), bool)
        chrom = np.full(len(numbers), None, dtype=object)
        pos = np.full(len(numbers), np.nan)
        if found.any():
            chrom[found] = [CODE_TO_CHROM[int(c)] for c in np.asarray(self.chrom)[idx_clipped[found]]]
            pos[found] = np.asarray(self.pos)[idx_clipped[found]].astype("float64")
        index = pd.Series(rsids).index
        return pd.DataFrame({"CHR": pd.Series(chrom, index=index, dtype=object), "POS19": pd.Series(pos, index=index)})
