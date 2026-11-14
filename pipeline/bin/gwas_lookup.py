#!/usr/bin/env python3
"""An indexed copy of a large table that a study YAML joins or excludes by, so each piece of a GWAS reads only the
rows it needs.

A join or exclude item with `index: true` is read once, start to finish, into <file>.index/ beside the file:
  hash.u64    the 64-bit hash of each key, sorted
  check.u64   a second, independent 64-bit hash of the same keys, in the same order
  row.i64     for each sorted hash, the row of the table it came from (the first row when a key repeats)
  data.arrow  every column as text, in file order, memory-mapped when read
  meta.json   the source file's size and time, the key column, the columns and the row count (written last)
A piece looks up its keys by binary search on the hashes, requires the second hash to agree (128 bits in all, so a
different key matching by chance is not a practical risk), and reads only the matching rows, and only the columns it is
given. The index is rebuilt when the file changes, and built once when several GWAS start together (a lock file).
Built for the deCODE annotation file, 27 million rows, which would otherwise be read in full for every piece of every
protein. Keys are compared as text.
"""

from __future__ import annotations

import fcntl
import json
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gwas_pipeline_common import _header_line_and_metadata_count, config_error, infer_delimiter

try:
    import pyarrow as pa
except ImportError:  # reported where an index is asked for
    pa = None

VERSION = 2
CHUNK_ROWS = 1_000_000
# What read_table reads as missing, for columns that come back from the index as text
MISSING = {"", ".", "NA", "NaN", "nan", "N/A", "n/a", "NULL", "null", "#N/A", "<NA>", "-nan", "-NaN"}
_OPEN: dict[str, "Index"] = {}  # open indexes by directory, kept for the pieces that follow


def index_dir(path: Path) -> Path:
    return path.with_name(path.name + ".index")


def hash_keys(values: Any, second: bool = False) -> np.ndarray:
    key = "gwasprocess-0002" if second else "0123456789123456"
    return np.asarray(pd.util.hash_array(np.asarray(values, dtype=object), hash_key=key), dtype="<u8")


def _require_pyarrow() -> None:
    if pa is None:
        raise config_error("`index: true` in a study YAML needs pyarrow: run `bash src/source.sh` to install it")


def _meta(directory: Path) -> dict[str, Any] | None:
    try:
        return json.loads((directory / "meta.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _current(directory: Path, path: Path, key: str) -> dict[str, Any] | None:
    meta = _meta(directory)
    stat = path.stat()
    if meta and meta.get("version") == VERSION and meta.get("key") == key and meta.get("size") == stat.st_size \
            and meta.get("mtime") == int(stat.st_mtime):
        return meta
    return None


def _build(path: Path, key: str, delimiter: str | None, out: Path) -> dict[str, Any] | None:
    """Write the index into `out`; None when the file has no column `key`."""
    sep = delimiter or infer_delimiter(path)
    _, skip = _header_line_and_metadata_count(path)
    reader = pd.read_csv(path, sep=sep, compression="gzip" if str(path).endswith(".gz") else None, skiprows=skip,
                         dtype=str, keep_default_na=False, na_filter=False, chunksize=CHUNK_ROWS)
    out.mkdir(parents=True)
    hashes: list[np.ndarray] = []
    checks: list[np.ndarray] = []
    columns: list[str] = []
    writer = None
    sink = pa.OSFile(str(out / "data.arrow"), "wb")
    try:
        for chunk in reader:
            if writer is None:
                columns = [str(c) for c in chunk.columns]
                if key not in columns:
                    return None
            table = pa.Table.from_pandas(chunk.astype(object), preserve_index=False,
                                         schema=pa.schema([(c, pa.string()) for c in chunk.columns]))
            if writer is None:
                writer = pa.ipc.new_file(sink, table.schema)
            writer.write_table(table)
            keys = chunk[key].to_numpy(dtype=object)
            hashes.append(hash_keys(keys))
            checks.append(hash_keys(keys, second=True))
        if writer is None:  # a header and no rows
            columns = [str(c) for c in pd.read_csv(path, sep=sep, skiprows=skip, nrows=0).columns]
            if key not in columns:
                return None
            writer = pa.ipc.new_file(sink, pa.schema([(c, pa.string()) for c in columns]))
        writer.close()
    finally:
        sink.close()
    values = np.concatenate(hashes) if hashes else np.empty(0, dtype="<u8")
    order = np.argsort(values, kind="stable")  # stable: the first row of a repeated key comes first
    values[order].tofile(out / "hash.u64")
    (np.concatenate(checks) if checks else np.empty(0, dtype="<u8"))[order].tofile(out / "check.u64")
    order.astype("<i8").tofile(out / "row.i64")
    duplicates = int((values[order][1:] == values[order][:-1]).sum()) if len(values) > 1 else 0
    stat = path.stat()
    meta = {"version": VERSION, "source": str(path), "size": stat.st_size, "mtime": int(stat.st_mtime), "key": key,
            "columns": columns, "rows": int(len(values)), "duplicate_keys": duplicates}
    (out / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    return meta


def ensure_index(path: Path, key: str, delimiter: str | None = None) -> Path | None:
    """The index directory of `path` keyed on `key`, built (once, under a lock) when missing or out of date; None when
    the file has no column `key`."""
    _require_pyarrow()
    out = index_dir(path)
    if _current(out, path, key):
        return out
    lock = open(path.with_name(path.name + ".index.lock"), "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if _current(out, path, key):  # another GWAS built it while this one waited
            return out
        print(f"Building the lookup index of {path} (first use; read once, start to finish)", flush=True)
        tmp = out.with_name(out.name + f".tmp{os.getpid()}")
        shutil.rmtree(tmp, ignore_errors=True)
        try:
            meta = _build(path, key, delimiter, tmp)
            if meta is None:
                return None
            shutil.rmtree(out, ignore_errors=True)
            os.replace(tmp, out)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        return out
    finally:
        fcntl.flock(lock, fcntl.LOCK_UN)
        lock.close()


class Index:
    def __init__(self, directory: Path):
        self.meta = _meta(directory) or {}
        self.key = self.meta["key"]
        full = bool(self.meta["rows"])
        self.hash = np.memmap(directory / "hash.u64", dtype="<u8", mode="r") if full else np.empty(0, dtype="<u8")
        self.check = np.memmap(directory / "check.u64", dtype="<u8", mode="r") if full else np.empty(0, dtype="<u8")
        self.row = np.memmap(directory / "row.i64", dtype="<i8", mode="r") if full else np.empty(0, dtype="<i8")
        self.table = pa.ipc.open_file(pa.memory_map(str(directory / "data.arrow"), "r")).read_all()

    def rows(self, text: np.ndarray) -> np.ndarray:
        """For each key text, the row of the indexed table holding it (the first when it repeats), -1 when absent."""
        found = np.full(len(text), -1, dtype="int64")
        if not len(self.hash) or not len(text):
            return found
        wanted = hash_keys(text)
        position = np.minimum(np.searchsorted(self.hash, wanted), len(self.hash) - 1)
        hit = np.flatnonzero((self.hash[position] == wanted) & (self.check[position] == hash_keys(text, second=True)))
        found[hit] = self.row[position[hit]]
        return found

    def contains(self, keys: pd.Series) -> np.ndarray:
        return self.rows(keys.astype(str).to_numpy(dtype=object)) >= 0

    def fetch(self, keys: pd.Series, columns: list[str] | None = None) -> pd.DataFrame:
        """One row per distinct key of `keys` found in the table: the key column (as text) and `columns` (default every
        column). A column is numeric where every value is a number, else text; missing values are NaN."""
        text = keys.astype(str).to_numpy(dtype=object)
        found = self.rows(text)
        present = np.flatnonzero(found >= 0)
        rows, first = np.unique(found[present], return_index=True)
        names = [c for c in (columns if columns is not None else self.meta["columns"]) if c != self.key and c in self.meta["columns"]]
        if len(rows):
            frame = self.table.select(names).take(pa.array(rows)).to_pandas()
        else:
            frame = pd.DataFrame({c: pd.Series(dtype=object) for c in names})
        frame.insert(0, self.key, pd.Series(text[present][first], dtype=object).to_numpy())
        for column in names:
            values = frame[column]
            missing = values.isin(MISSING)
            numeric = pd.to_numeric(values.where(~missing), errors="coerce")
            frame[column] = numeric if numeric.notna().sum() == (~missing).sum() and len(values) else values.where(~missing, np.nan)
        return frame


def get_index(path: Path, key: str, delimiter: str | None = None) -> Index | None:
    """The opened index of `path` keyed on `key` (built if needed, and kept open for the pieces that follow)."""
    directory = ensure_index(path, key, delimiter)
    if directory is None:
        return None
    meta = _meta(directory) or {}
    signature = f"{meta.get('size')}:{meta.get('mtime')}:{key}"
    opened = _OPEN.get(str(directory))
    if opened is None or opened.signature != signature:
        opened = _OPEN[str(directory)] = Index(directory)
        opened.signature = signature
    return opened
