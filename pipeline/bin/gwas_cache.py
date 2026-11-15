#!/usr/bin/env python3
"""A cache of what the pipeline works out for a variant from its alleles and position alone, shared by every GWAS.

Studies that test the same variants (the 4,907 deCODE proteins, or the Olink and SomaScan proteins of one study) repeat
the same per-variant work in every file: the 1000 Genomes annotation of (CHR, POS38, EA, OA), the alignment of the
alleles to the GRCh38 sequence, and the left-alignment of indels. These results are pure functions of the variant, so
the first GWAS that meets a variant stores them and every later GWAS reads them back. Nothing that depends on a
study's own statistics (frequencies, betas, P values, palindromic strand decisions) is cached.

Layout, under <cache root>/<kind>-<fingerprint>/<bucket 00-63>/:
  seg-<time>-<pid>-<n>.arrow   rows sorted by h1: h1 and h2 (2 independent 64-bit hashes of the variant key, 128 bits in
                               all) and the cached values; written whole under a temporary name and renamed, so a
                               reader never sees a partial file
A lookup reads every segment of the buckets its keys fall in (memory-mapped) and takes the first match. Anything that
goes wrong (a missing file, a bucket being compacted, a full disk) is a cache miss: the value is worked out again, so
the output never depends on the cache's contents. The fingerprint holds the reference files' names, sizes and times, so
a rebuilt reference starts a new cache. `compact` merges a bucket's segments once it holds more than a few; it runs at
the end of a GWAS and skips a bucket another job is already compacting. The segments it merges are renamed to
`.old` and deleted a day later, because a job on another node may still have one memory-mapped and a deleted file
kills that job on a shared file system.

GWAS_VARIANT_CACHE sets the cache root (default <reference root>/variant_cache); `off` switches the cache off.
"""

from __future__ import annotations

import fcntl
import hashlib
import itertools
import os
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gwas_lookup import hash_keys

try:
    import pyarrow as pa
except ImportError:  # the cache is then off
    pa = None

BUCKETS = 64
SHIFT = np.uint64(58)  # the top 6 bits of h1 pick the bucket
COMPACT_AT = 8  # segments in a bucket before it is merged into 1
OLD_AFTER = 24 * 3600  # seconds a merged-away segment is kept before it is deleted
_counter = itertools.count()
_OPEN: dict[tuple[str, str], "VariantCache | None"] = {}


def cache_root(reference_root: str | Path) -> Path | None:
    setting = os.environ.get("GWAS_VARIANT_CACHE", "")
    if setting.lower() in ("off", "0", "false", "no"):
        return None
    return Path(setting) if setting else Path(reference_root) / "variant_cache"


def file_signature(*paths: Any) -> str:
    """Names, sizes and times of the reference files a cached value was worked out from."""
    parts = []
    for path in paths:
        try:
            stat = os.stat(path)
            parts.append(f"{Path(path).name}:{stat.st_size}:{int(stat.st_mtime)}")
        except OSError:
            parts.append(f"{Path(path).name}:missing")
    return "|".join(parts)


def variant_keys(*columns: Any) -> tuple[np.ndarray, np.ndarray]:
    """The 2 hashes of a row's key: its columns as text joined with `|` (so ("1|2", "3") and ("1", "2|3") differ only
    through the second hash and the check below; chromosome, position and alleles contain no `|`)."""
    text = pd.Series(np.asarray(columns[0], dtype=object)).astype(str)
    for column in columns[1:]:
        text = text + "|" + pd.Series(np.asarray(column, dtype=object)).astype(str)
    values = text.to_numpy(dtype=object)
    return hash_keys(values), hash_keys(values, second=True)


class VariantCache:
    def __init__(self, directory: Path, columns: dict[str, Any]):
        self.directory = directory
        self.columns = columns  # value name -> pyarrow type
        self.hits = 0
        self.misses = 0

    @classmethod
    def open(cls, root: Path | None, kind: str, fingerprint: str, columns: dict[str, Any]) -> "VariantCache | None":
        if root is None or pa is None:
            return None
        key = (str(root), kind + fingerprint)
        if key not in _OPEN:
            name = f"{kind}-{hashlib.sha1(fingerprint.encode()).hexdigest()[:12]}"
            directory = root / name
            try:
                directory.mkdir(parents=True, exist_ok=True)
                _OPEN[key] = cls(directory, columns)
            except OSError:
                _OPEN[key] = None
        return _OPEN[key]

    def _segments(self, bucket: int) -> list[Path]:
        try:
            return sorted(p for p in (self.directory / f"{bucket:02d}").iterdir() if p.name.startswith("seg-") and p.suffix == ".arrow")
        except OSError:
            return []

    def lookup(self, h1: np.ndarray, h2: np.ndarray) -> tuple[np.ndarray, dict[str, np.ndarray]]:
        """(found, values): for each key whether it is cached, and its cached values (object arrays; a float column
        comes back as float64)."""
        n = len(h1)
        found = np.zeros(n, dtype=bool)
        values: dict[str, np.ndarray] = {c: (np.full(n, np.nan) if pa.types.is_floating(t) else np.empty(n, dtype=object))
                                          for c, t in self.columns.items()}
        if n:
            bucket = (h1 >> SHIFT).astype(np.int64)
            for b in np.unique(bucket):
                idx = np.flatnonzero(bucket == b)
                for segment in self._segments(int(b)):
                    open_rows = idx[~found[idx]]
                    if not len(open_rows):
                        break
                    try:
                        table = pa.ipc.open_file(pa.memory_map(str(segment), "r")).read_all().combine_chunks()
                        k1 = table.column("h1").to_numpy()
                        k2 = table.column("h2").to_numpy()
                        if not len(k1):
                            continue
                        position = np.minimum(np.searchsorted(k1, h1[open_rows]), len(k1) - 1)
                        hit = (k1[position] == h1[open_rows]) & (k2[position] == h2[open_rows])
                        if not hit.any():
                            continue
                        rows, source = open_rows[hit], position[hit]
                        taken = {name: table.column(name).take(pa.array(source)).to_numpy(zero_copy_only=False)
                                 for name in self.columns}
                    except (OSError, ValueError, pa.ArrowException):
                        continue  # compacted away, stale or unreadable: a miss
                    for name in self.columns:
                        values[name][rows] = taken[name]
                    found[rows] = True
        self.hits += int(found.sum())
        self.misses += int((~found).sum())
        return found, values

    def store(self, h1: np.ndarray, h2: np.ndarray, values: dict[str, np.ndarray]) -> None:
        """Add keys (and their values) as new segments, 1 per bucket; a failure to write is ignored."""
        if not len(h1):
            return
        bucket = (h1 >> SHIFT).astype(np.int64)
        for b in np.unique(bucket):
            idx = np.flatnonzero(bucket == b)
            order = idx[np.argsort(h1[idx], kind="stable")]
            _, first = np.unique(h1[order], return_index=True)
            order = order[first]
            arrays = [pa.array(h1[order], type=pa.uint64()), pa.array(h2[order], type=pa.uint64())]
            names = ["h1", "h2"]
            for name, kind in self.columns.items():
                arrays.append(pa.array(values[name][order].tolist() if not pa.types.is_floating(kind) else values[name][order], type=kind))
                names.append(name)
            self._write(int(b), pa.Table.from_arrays(arrays, names=names))

    def _write(self, bucket: int, table: Any) -> Path | None:
        directory = self.directory / f"{bucket:02d}"
        name = f"seg-{time.time_ns()}-{os.getpid()}-{next(_counter)}"
        tmp = directory / f".{name}.tmp"
        try:
            directory.mkdir(exist_ok=True)
            with pa.OSFile(str(tmp), "wb") as sink:
                with pa.ipc.new_file(sink, table.schema) as writer:
                    writer.write_table(table)
            final = directory / f"{name}.arrow"
            os.replace(tmp, final)
            return final
        except OSError:
            try:
                tmp.unlink()
            except OSError:
                pass
            return None

    def compact(self, at: int = COMPACT_AT) -> int:
        """Merge each bucket that holds more than `at` segments into 1; returns how many buckets were merged. A bucket
        another job is merging is skipped."""
        merged = 0
        for b in range(BUCKETS):
            segments = self._segments(b)
            if len(segments) <= at:
                continue
            lock_path = self.directory / f"{b:02d}.lock"
            try:
                handle = open(lock_path, "w")
            except OSError:
                continue
            try:
                try:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except OSError:
                    continue
                segments = self._segments(b)
                tables = []
                for segment in segments:
                    try:
                        tables.append(pa.ipc.open_file(pa.memory_map(str(segment), "r")).read_all())
                    except (OSError, pa.ArrowException):
                        pass
                if len(tables) <= 1:
                    continue
                table = pa.concat_tables(tables).combine_chunks()
                h1 = table.column("h1").to_numpy()
                order = np.argsort(h1, kind="stable")
                _, first = np.unique(h1[order], return_index=True)
                keep = pa.array(order[first])
                final = self._write(b, table.take(keep))
                if final is None:
                    continue
                for segment in segments:  # renamed, not deleted: a job that has the file mapped keeps a valid handle
                    try:
                        os.replace(segment, segment.with_name(segment.name + ".old"))
                    except OSError:
                        pass
                merged += 1
            finally:
                handle.close()
        self._remove_old()
        return merged

    def _remove_old(self, age: float = OLD_AFTER) -> None:
        """Delete the segments that compaction set aside once no running job can still be using them."""
        cutoff = time.time() - age
        for bucket in self.directory.glob("[0-9][0-9]"):
            try:
                for path in bucket.glob("seg-*.arrow.old"):
                    try:
                        if path.stat().st_mtime < cutoff:
                            path.unlink()
                    except OSError:
                        pass
            except OSError:
                pass
