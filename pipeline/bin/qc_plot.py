#!/usr/bin/env python3
"""Raw vs standardised QQ plots (compressed TIFF) with genomic inflation factors.

Usage:
    qc_plot.py --standardise-record standardise_record.json [--force false] [--dry-run false] --output row_record.json

Works in -log10 P (LP) space so P-values beyond double precision are plotted. Every point with LP >= 2
is drawn; the dense null region below is thinned to plotting resolution, which is visually lossless.
The raw panel reads the columns and delimiter the standardiser used. An existing plot is reused only
while it is newer than the standardised VCF.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import gwas_stats  # noqa: E402
from gwas_pipeline_common import (  # noqa: E402
    BETA_ALIASES,
    LOGP_ALIASES,
    P_ALIASES,
    SE_ALIASES,
    count_metadata_lines,
    infer_delimiter,
    is_gwas_vcf,
    numeric_series,
    parse_bool,
    pick_existing_column,
    read_json,
    write_json,
)

TABLE_CHUNK = 1_000_000
THIN_BELOW = 2.0
THIN_RESOLUTION = 0.005


def lp_from_vcf(path: Path) -> np.ndarray:
    opener = gzip.open if str(path).endswith(".gz") else open
    values: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            keys = parts[8].split(":")
            if "LP" in keys:
                sample = parts[9].split(":")
                index = keys.index("LP")
                values.append(sample[index] if index < len(sample) else ".")
    return numeric_series(pd.Series(values, dtype=object)).to_numpy(dtype="float64")


def lp_from_table(path: Path, source_columns: dict[str, str] | None = None, delimiter: str | None = None) -> np.ndarray:
    sep = delimiter or infer_delimiter(path)
    skip = count_metadata_lines(path)
    compression = "gzip" if str(path).endswith(".gz") else None
    header = list(pd.read_csv(path, sep=sep, compression=compression, nrows=0, skiprows=skip).columns)
    stripped = {str(c).lstrip("#"): c for c in header}
    source_columns = source_columns or {}

    def column(field: str, aliases: list[str]) -> str | None:
        chosen = source_columns.get(field)
        if chosen in stripped:
            return stripped[chosen]
        if source_columns:
            return None
        found = pick_existing_column(list(stripped), aliases)
        return stripped[found] if found else None

    cols = {field: column(field, aliases) for field, aliases in
            (("P", P_ALIASES), ("LOG10P", LOGP_ALIASES), ("BETA", BETA_ALIASES), ("SE", SE_ALIASES))}
    usecols = sorted({c for c in cols.values() if c})
    if not (cols["P"] or cols["LOG10P"]):
        return np.array([], dtype="float64")
    chunks = []
    reader = pd.read_csv(path, sep=sep, compression=compression, usecols=usecols, chunksize=TABLE_CHUNK, low_memory=False,
                         skiprows=skip, na_values=["."])
    for chunk in reader:
        get = lambda field: numeric_series(chunk[cols[field]]) if cols[field] else None  # noqa: E731
        lp, _ = gwas_stats.compute_lp(p=get("P"), log10p=get("LOG10P"), beta=get("BETA"), se=get("SE"))
        chunks.append(lp)
    return np.concatenate(chunks) if chunks else np.array([], dtype="float64")


def extract_lp(path: str | Path, source_columns: dict[str, str] | None = None, delimiter: str | None = None) -> np.ndarray:
    path = Path(path)
    if not path.exists():
        return np.array([], dtype="float64")
    if is_gwas_vcf(path):
        return lp_from_vcf(path)
    return lp_from_table(path, source_columns, delimiter)


def qq_points(lp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(expected, observed) LP pairs, largest first, with the LP < 2 bulk thinned to plotting resolution."""
    values = np.asarray(lp, dtype="float64")
    observed = np.sort(values[np.isfinite(values) & (values >= 0)])[::-1]
    n = observed.size
    if n == 0:
        return np.array([]), np.array([])
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    bulk = (observed < THIN_BELOW) & (expected < THIN_BELOW)
    cells = np.stack([np.round(expected[bulk] / THIN_RESOLUTION), np.round(observed[bulk] / THIN_RESOLUTION)], axis=1)
    _, first = np.unique(cells, axis=0, return_index=True)
    keep = np.concatenate([np.where(~bulk)[0], np.where(bulk)[0][np.sort(first)]])
    keep.sort()
    return expected[keep], observed[keep]


def draw_qq(ax: Any, lp: np.ndarray, title: str) -> float | None:
    expected, observed = qq_points(lp)
    if observed.size == 0:
        ax.text(0.5, 0.5, "No valid P values", ha="center", va="center")
        ax.set_title(title)
        return None
    lam = gwas_stats.lambda_gc(lp)
    n = int(np.isfinite(lp).sum())
    ax.scatter(expected, observed, s=5, alpha=0.6, linewidths=0, rasterized=True)
    limit = float(max(expected.max(), observed.max(), 1.0))
    ax.plot([0, limit], [0, limit], color="black", linewidth=0.8)
    ax.set_xlim(0, limit * 1.02)
    ax.set_ylim(0, limit * 1.02)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(f"{title} (n={n:,}, lambda={lam:.3f})" if lam is not None else f"{title} (n={n:,})")
    return lam


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--standardise-record", required=True)
    parser.add_argument("--force", default="false")
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    record = read_json(args.standardise_record)
    plot_path = Path(record["plot_tiff"])
    vcf_path = Path(record["standard_vcf"])
    if parse_bool(args.dry_run):
        record["plot_created"] = False
        write_json(args.output, record)
        return
    fresh = plot_path.exists() and plot_path.stat().st_size > 0 and vcf_path.exists() \
        and plot_path.stat().st_mtime >= vcf_path.stat().st_mtime
    if fresh and not parse_bool(args.force):
        record["plot_created"] = False
        write_json(args.output, record)
        return

    qc = record.get("qc") or {}
    raw_lp = extract_lp(record["raw_path"], qc.get("source_columns"), qc.get("input_delimiter"))
    std_lp = extract_lp(vcf_path)
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    record["lambda_gc_raw"] = draw_qq(axes[0], raw_lp, "Raw QQ")
    record["lambda_gc_standardised"] = draw_qq(axes[1], std_lp, "Standardised QQ")
    fig.suptitle(f"{record.get('study_id', '')} | {record.get('phenotype', '')} | {record.get('ancestry', '')}", fontsize=11)
    tmp = plot_path.with_name(f".{plot_path.name}.tmp")
    fig.savefig(tmp, format="tiff", dpi=150, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    tmp.replace(plot_path)
    record["plot_created"] = True
    write_json(args.output, record)


if __name__ == "__main__":
    main()
