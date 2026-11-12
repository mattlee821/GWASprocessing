#!/usr/bin/env python3
"""GRCh38 reference lookups: dbSNP (rsID <-> position) and 1000 Genomes allele frequencies.

Both are read from memory-mapped arrays built once by src/reference_data (build_dbsnp_index.py and
build_1000genomes_grch38.py), so a process holds only the pages it touches. Keep them on local disk: the lookups
read at random positions.

    <root>/dbsnp/GRCh38/          chr<C>.{pos,hash,rs}.u32 sorted by position; rs.u32 + rs_chrom.u8 + rs_pos.u32
    <root>/1000genomes/GRCh38/    <POP>.{key.i64,ref.u1,alt.u1,af.f4}: biallelic-split SNVs, sorted by key
    <root>/build_detection/       grch37.key.i64 / grch38.key.i64 (+ .ref.u1/.alt.u1): the same SNPs in both builds

Chromosome codes: 1-22, 23 = X, 24 = Y, 25 = MT; key = code * 1e9 + position.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]
CHROM_CODE = {c: i for i, c in enumerate(CHROMS, start=1)}
BASES = {"A": 1, "C": 2, "G": 3, "T": 4}
KEY_SCALE = 1_000_000_000


def allele_hash(ref: Any, alt: Any) -> np.ndarray:
    """Deterministic 32-bit hash of the allele pair 'REF>ALT' (uppercase), for matching dbSNP records."""
    text = pd.Series(ref, dtype=object).astype(str).str.upper() + ">" + pd.Series(alt, dtype=object).astype(str).str.upper()
    return (pd.util.hash_array(text.to_numpy(dtype=object)) & np.uint64(0xFFFFFFFF)).astype("uint32")


def chrom_codes(chrom: Any) -> np.ndarray:
    values = pd.Series(chrom, dtype=object).astype(str).str.upper().str.replace("^CHR", "", regex=True)
    values = values.replace({"23": "X", "24": "Y", "25": "MT", "M": "MT"})
    return values.map(CHROM_CODE).fillna(0).to_numpy(dtype="int64")


def position_keys(chrom: Any, pos: Any) -> np.ndarray:
    """int64 key per row (code * 1e9 + position); -1 where the chromosome or position is unusable."""
    codes = chrom_codes(chrom)
    positions = pd.to_numeric(pd.Series(pos), errors="coerce").to_numpy(dtype="float64")
    ok = (codes > 0) & np.isfinite(positions) & (positions > 0)
    keys = np.full(len(codes), -1, dtype="int64")
    keys[ok] = codes[ok] * KEY_SCALE + positions[ok].astype("int64")
    return keys


def base_codes(alleles: Any) -> np.ndarray:
    """1-4 for single A/C/G/T bases, 0 otherwise."""
    return pd.Series(alleles, dtype=object).astype(str).str.upper().map(BASES).fillna(0).to_numpy(dtype="uint8")


BLOCK = 4096  # coarse-index stride: a lookup reads about 4 pages of a memory-mapped array instead of about 30


def coarse_index(sorted_array: np.ndarray) -> np.ndarray:
    """Every BLOCK-th value of a sorted array, held in memory, to narrow binary searches to 1 block."""
    return np.asarray(sorted_array[::BLOCK])


def searchsorted_blocked(sorted_array: np.ndarray, coarse: np.ndarray, query: np.ndarray, side: str = "left") -> np.ndarray:
    """np.searchsorted on a large memory-mapped array, with the top levels of the search done in memory (coarse)."""
    n = len(sorted_array)
    if n == 0:
        return np.zeros(len(query), dtype="int64")
    block = np.searchsorted(coarse, query, side=side) - 1
    lo = np.clip(block, 0, None).astype("int64") * BLOCK
    hi = np.minimum(lo + BLOCK + 1, n).astype("int64")
    lo = np.where(block < 0, 0, lo)
    while True:
        active = lo < hi
        if not active.any():
            return lo
        mid = (lo + hi) // 2
        value = np.asarray(sorted_array[np.minimum(mid, n - 1)])
        go_right = (value < query) if side == "left" else (value <= query)
        go_right &= active
        lo = np.where(go_right, mid + 1, lo)
        hi = np.where(active & ~go_right, mid, hi)


def _memmap(path: Path, dtype: str) -> np.ndarray:
    return np.memmap(path, dtype=dtype, mode="r") if path.exists() and path.stat().st_size else np.zeros(0, dtype=dtype)


class Dbsnp:
    """dbSNP GRCh38 index: positions for rsIDs, and rsIDs for (position, REF, ALT)."""

    def __init__(self, root: Path):
        self.root = Path(root)
        self.meta = json.loads((self.root / "meta.json").read_text(encoding="utf-8"))
        self.rs = _memmap(self.root / "rs.u32", "<u4")
        self.rs_chrom = _memmap(self.root / "rs_chrom.u8", "u1")
        self.rs_pos = _memmap(self.root / "rs_pos.u32", "<u4")
        self.merged_from = _memmap(self.root / "merged_from.u32", "<u4")
        self.merged_into = _memmap(self.root / "merged_into.u32", "<u4")
        self._chrom: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    @property
    def has_merges(self) -> bool:
        return len(self.merged_from) > 0

    def current_rsids(self, rs_numbers: np.ndarray) -> np.ndarray:
        """The current rsID number for each retired (merged) rsID number, from build_dbsnp_merges.py; 0 when the rsID
        was never merged or the merge map is not built."""
        rs_numbers = np.asarray(rs_numbers, dtype="int64")
        out = np.zeros(len(rs_numbers), dtype="int64")
        valid = (rs_numbers > 0) & (rs_numbers < 2 ** 32)
        if not valid.any() or not self.has_merges:
            return out
        query = rs_numbers[valid].astype("uint32")
        order = np.argsort(query, kind="stable")
        if not hasattr(self, "_merged_coarse"):
            self._merged_coarse = coarse_index(self.merged_from)
        idx = np.minimum(searchsorted_blocked(self.merged_from, self._merged_coarse, query[order]), len(self.merged_from) - 1)
        found = np.asarray(self.merged_from[idx]) == query[order]
        rows = np.where(valid)[0][order]
        out[rows[found]] = np.asarray(self.merged_into[idx[found]], dtype="int64")
        return out

    @classmethod
    def find(cls, reference_root: str | Path) -> "Dbsnp | None":
        root = Path(reference_root) / "dbsnp" / "GRCh38"
        return cls(root) if (root / "meta.json").exists() else None

    @property
    def label(self) -> str:
        return f"dbSNP {self.meta.get('dbsnp_build', '')} GRCh38".strip()

    def positions(self, rs_numbers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """(chromosome code, position) per rsID number; (0, 0) when the rsID is absent or at more than 1 position."""
        rs_numbers = np.asarray(rs_numbers, dtype="int64")
        code = np.zeros(len(rs_numbers), dtype="int64")
        pos = np.zeros(len(rs_numbers), dtype="int64")
        valid = (rs_numbers > 0) & (rs_numbers < 2 ** 32)
        if not valid.any() or not len(self.rs):
            return code, pos
        query = rs_numbers[valid].astype("uint32")
        order = np.argsort(query, kind="stable")
        if not hasattr(self, "_rs_coarse"):
            self._rs_coarse = coarse_index(self.rs)
        idx = searchsorted_blocked(self.rs, self._rs_coarse, query[order])
        idx_clipped = np.minimum(idx, len(self.rs) - 1)
        found = np.asarray(self.rs[idx_clipped]) == query[order]
        rows = np.where(valid)[0][order]
        code[rows[found]] = np.asarray(self.rs_chrom[idx_clipped[found]], dtype="int64")
        pos[rows[found]] = np.asarray(self.rs_pos[idx_clipped[found]], dtype="int64")
        return code, pos

    def _chromosome(self, chrom: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        if chrom not in self._chrom:
            arrays = tuple(_memmap(self.root / f"chr{chrom}.{kind}.u32", "<u4") for kind in ("pos", "hash", "rs"))
            self._chrom[chrom] = (*arrays, coarse_index(arrays[0]))
        return self._chrom[chrom]

    def rsids(self, chrom: Any, pos: Any, ref: Any, alt: Any) -> np.ndarray:
        """dbSNP rsID number per row matching (chromosome, position, REF, ALT) exactly; 0 when none."""
        codes = chrom_codes(chrom)
        positions = pd.to_numeric(pd.Series(pos), errors="coerce").fillna(0).to_numpy(dtype="int64")
        wanted = allele_hash(ref, alt)
        out = np.zeros(len(codes), dtype="int64")
        for code in np.unique(codes[codes > 0]):
            rows = np.where(codes == code)[0]
            db_pos, db_hash, db_rs, coarse = self._chromosome(CHROMS[code - 1])
            if not len(db_pos):
                continue
            rows = rows[np.argsort(positions[rows], kind="stable")]  # sorted queries read the memory map nearly sequentially
            query = positions[rows].astype("uint32")
            lo = searchsorted_blocked(db_pos, coarse, query, "left")
            hi = searchsorted_blocked(db_pos, coarse, query, "right")
            span = hi - lo
            for k in range(int(span.max()) if len(span) else 0):
                todo = (span > k) & (out[rows] == 0)
                if not todo.any():
                    continue
                at = lo[todo] + k
                match = np.asarray(db_hash[at]) == wanted[rows[todo]]
                out[rows[todo][match]] = np.asarray(db_rs[at[match]], dtype="int64")
        return out


class Frequencies:
    """1000 Genomes GRCh38 ALT-allele frequencies for one super-population (biallelic-split SNVs)."""

    def __init__(self, root: Path, population: str):
        self.population = population
        prefix = Path(root) / population
        self.key = _memmap(Path(f"{prefix}.key.i64"), "<i8")
        self.ref = _memmap(Path(f"{prefix}.ref.u1"), "u1")
        self.alt = _memmap(Path(f"{prefix}.alt.u1"), "u1")
        self.af = _memmap(Path(f"{prefix}.af.f4"), "<f4")

    @classmethod
    def find(cls, reference_root: str | Path, population: str) -> "Frequencies | None":
        root = Path(reference_root) / "1000genomes" / "GRCh38"
        for pop in (population, "ALL"):
            if Path(root / f"{pop}.key.i64").exists():
                return cls(root, pop)
        return None

    def annotate(self, chrom: Any, pos: Any) -> pd.DataFrame:
        """Up to 3 records at each row's position, as columns A1_k (ALT base), A2_k (REF base) and AF1_k (ALT frequency)
        for k = 0, 1, 2; the caller picks the record whose alleles match the row."""
        keys = position_keys(chrom, pos)
        n = len(keys)
        lo = np.searchsorted(self.key, keys, side="left")
        hi = np.searchsorted(self.key, keys, side="right")
        out = {}
        letters = np.array(["", "A", "C", "G", "T"], dtype=object)
        for k in range(3):
            has = (keys >= 0) & (hi - lo > k)
            at = np.where(has, lo + k, 0)
            alt = np.where(has, letters[np.asarray(self.alt[at]) if len(self.alt) else np.zeros(n, dtype=int)], "")
            ref = np.where(has, letters[np.asarray(self.ref[at]) if len(self.ref) else np.zeros(n, dtype=int)], "")
            af = np.where(has, np.asarray(self.af[at], dtype="float64") if len(self.af) else np.nan, np.nan)
            out[f"A1_{k}"], out[f"A2_{k}"], out[f"AF1_{k}"] = alt.astype(object), ref.astype(object), af
        return pd.DataFrame(out)


def paired_build_sample(reference_root: str | Path) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] | None:
    """The same common SNPs in both builds, for build detection: {build: (sorted keys, ref code, alt code)}."""
    root = Path(reference_root) / "build_detection"
    if not (root / "grch38.key.i64").exists():
        return None
    out = {}
    for build, name in (("hg19", "grch37"), ("hg38", "grch38")):
        key = np.fromfile(root / f"{name}.key.i64", dtype="<i8")
        ref = np.fromfile(root / f"{name}.ref.u1", dtype="u1")
        alt = np.fromfile(root / f"{name}.alt.u1", dtype="u1")
        out[build] = (key, ref, alt)
    return out


def build_hits(sample: tuple[np.ndarray, np.ndarray, np.ndarray], chrom: Any, pos: Any, a: Any, b: Any) -> int:
    """How many rows match a sampled SNP at the same position with the same allele pair (either order, either strand)."""
    key, ref, alt = sample
    q = position_keys(chrom, pos)
    idx = np.minimum(np.searchsorted(key, q), max(len(key) - 1, 0))
    at = (q >= 0) & (len(key) > 0) & (key[idx] == q)
    if not at.any():
        return 0
    comp = np.array([0, 4, 3, 2, 1], dtype="uint8")  # A<->T, C<->G
    ca, cb = base_codes(a)[at], base_codes(b)[at]
    r, s = ref[idx[at]], alt[idx[at]]
    same = ((ca == r) & (cb == s)) | ((ca == s) & (cb == r))
    flipped = ((comp[ca] == r) & (comp[cb] == s)) | ((comp[ca] == s) & (comp[cb] == r))
    return int((same | flipped).sum())
