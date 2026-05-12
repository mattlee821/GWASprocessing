#!/usr/bin/env python3
"""Create raw/standardised QQ and Manhattan QC plots."""

from __future__ import annotations

import argparse
import gzip
import math
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gwas_pipeline_common import (
    CHR_ALIASES,
    LOGP_ALIASES,
    P_ALIASES,
    POS38_ALIASES,
    POS_ALIASES,
    RSID_ALIASES,
    SNPID_ALIASES,
    ancestry_to_population,
    infer_delimiter,
    is_gwas_vcf,
    logp_to_p,
    normalize_chromosome,
    numeric_series,
    parse_bool,
    pick_existing_column,
    read_json,
    read_reference_lookup,
    write_json,
)


MAX_RAW_REFERENCE_ROWS = 1_000_000
MAX_QQ_ROWS = 200_000
MAX_MANHATTAN_ROWS = 300_000
TOP_ASSOCIATION_ROWS = 50_000
PLOT_CHUNK_SIZE = 500_000


def normalize_plot_identifiers(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    snp = out["SNP"].astype("string").str.strip()
    snpid = out["SNPID"].astype("string").str.strip()
    empty_snp = snp.isna() | snp.str.lower().isin(["", "nan", "none", "na", "."])
    empty_snpid = snpid.isna() | snpid.str.lower().isin(["", "nan", "none", "na", "."])
    original_snp = snp.mask(empty_snp, pd.NA)
    out.loc[empty_snpid & original_snp.notna(), "SNPID"] = original_snp[empty_snpid & original_snp.notna()]
    rs_from_snp = original_snp.str.extract(r"\b(rs\d+)\b", flags=re.I, expand=False).str.lower()
    rs_from_snpid = snpid.mask(empty_snpid, pd.NA).str.extract(r"\b(rs\d+)\b", flags=re.I, expand=False).str.lower()
    extracted = rs_from_snp.fillna(rs_from_snpid)
    out.loc[extracted.notna(), "SNP"] = extracted[extracted.notna()]
    return out


def select_plot_subset(df: pd.DataFrame, max_rows: int, top_rows: int = TOP_ASSOCIATION_ROWS) -> pd.DataFrame:
    if len(df) <= max_rows:
        return df
    out = df.copy()
    out["P"] = numeric_series(out["P"])
    valid = out[out["P"].notna()]
    if valid.empty:
        return out.sample(n=max_rows, random_state=1)
    top_n = min(top_rows, max_rows // 2, len(valid))
    top = valid.nsmallest(top_n, "P")
    remaining = out.drop(index=top.index, errors="ignore")
    background_n = min(max_rows - len(top), len(remaining))
    background = remaining.sample(n=background_n, random_state=1) if background_n > 0 else remaining.iloc[0:0]
    return pd.concat([top, background], ignore_index=True)


def plot_frame_from_table(df: pd.DataFrame, source_cols: dict[str, str | None]) -> pd.DataFrame:
    p_col = source_cols.get("P")
    logp_col = source_cols.get("LOG10P")
    if p_col:
        p = numeric_series(df[p_col])
    elif logp_col:
        p = df[logp_col].map(logp_to_p)
    else:
        p = pd.Series(np.nan, index=df.index)
    return pd.DataFrame(
        {
            "CHR": df[source_cols["CHR"]] if source_cols.get("CHR") else pd.Series(np.nan, index=df.index),
            "POS": numeric_series(df[source_cols["POS"]]) if source_cols.get("POS") else pd.Series(np.nan, index=df.index),
            "P": numeric_series(p),
            "SNP": df[source_cols["SNP"]].astype("string")
            if source_cols.get("SNP")
            else pd.Series(pd.NA, index=df.index, dtype="string"),
            "SNPID": df[source_cols["SNPID"]].astype("string")
            if source_cols.get("SNPID")
            else pd.Series(pd.NA, index=df.index, dtype="string"),
        }
    )


def extract_table_plot_df(path: Path, record: dict[str, Any] | None = None) -> pd.DataFrame:
    sep = infer_delimiter(path)
    compression = "gzip" if str(path).endswith(".gz") else None
    header = pd.read_csv(path, sep=sep, compression=compression, nrows=0).columns
    source_cols = {
        "CHR": pick_existing_column(header, CHR_ALIASES),
        "POS": pick_existing_column(header, POS38_ALIASES + POS_ALIASES),
        "SNP": pick_existing_column(header, RSID_ALIASES),
        "SNPID": pick_existing_column(header, SNPID_ALIASES),
        "P": pick_existing_column(header, P_ALIASES),
        "LOG10P": pick_existing_column(header, LOGP_ALIASES),
    }
    usecols = [col for col in dict.fromkeys(source_cols.values()) if col]
    if not usecols:
        return pd.DataFrame(columns=["CHR", "POS", "P", "SNP", "SNPID"])

    sampled: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        path,
        sep=sep,
        compression=compression,
        usecols=usecols,
        chunksize=PLOT_CHUNK_SIZE,
        low_memory=False,
    ):
        frame = normalize_plot_identifiers(plot_frame_from_table(chunk, source_cols))
        sampled.append(select_plot_subset(frame, MAX_RAW_REFERENCE_ROWS))
        if sum(len(item) for item in sampled) > MAX_RAW_REFERENCE_ROWS * 2:
            sampled = [select_plot_subset(pd.concat(sampled, ignore_index=True), MAX_RAW_REFERENCE_ROWS)]

    if not sampled:
        return pd.DataFrame(columns=["CHR", "POS", "P", "SNP", "SNPID"])
    out = select_plot_subset(pd.concat(sampled, ignore_index=True), MAX_RAW_REFERENCE_ROWS)
    return fill_plot_positions_from_reference(out, record or {})


def extract_vcf_plot_df(path: Path) -> pd.DataFrame:
    sampled: list[pd.DataFrame] = []
    rows: list[dict[str, Any]] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            info = {}
            for item in parts[7].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info[key] = value
            lp = np.nan
            if len(parts) >= 10:
                keys = parts[8].split(":")
                values = parts[9].split(":")
                if "LP" in keys:
                    index = keys.index("LP")
                    if index < len(values):
                        lp = values[index]
            rows.append(
                {
                    "CHR": parts[0],
                    "POS": parts[1],
                    "P": logp_to_p(lp),
                    "SNP": parts[2],
                    "SNPID": info.get("SNPID", pd.NA),
                }
            )
            if len(rows) >= PLOT_CHUNK_SIZE:
                sampled.append(select_plot_subset(normalize_plot_identifiers(pd.DataFrame(rows)), MAX_RAW_REFERENCE_ROWS))
                rows = []
                if sum(len(item) for item in sampled) > MAX_RAW_REFERENCE_ROWS * 2:
                    sampled = [select_plot_subset(pd.concat(sampled, ignore_index=True), MAX_RAW_REFERENCE_ROWS)]
    if rows:
        sampled.append(select_plot_subset(normalize_plot_identifiers(pd.DataFrame(rows)), MAX_RAW_REFERENCE_ROWS))
    if not sampled:
        return pd.DataFrame(columns=["CHR", "POS", "P", "SNP", "SNPID"])
    return select_plot_subset(pd.concat(sampled, ignore_index=True), MAX_RAW_REFERENCE_ROWS)


def extract_plot_df(path: str | Path, record: dict[str, Any] | None = None) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=["CHR", "POS", "P", "SNP", "SNPID"])
    if is_gwas_vcf(path):
        return extract_vcf_plot_df(path)

    return extract_table_plot_df(path, record)


def fill_plot_positions_from_reference(df: pd.DataFrame, record: dict[str, Any]) -> pd.DataFrame:
    if df.empty or not record:
        return df
    missing_chr = df["CHR"].isna() | df["CHR"].astype("string").str.lower().isin(["", "nan", "none", "na"])
    missing_pos = df["POS"].isna()
    if not (missing_chr.any() or missing_pos.any()):
        return df
    reference_root = record.get("reference_root")
    if not reference_root:
        return df
    population = ancestry_to_population(record.get("ancestry", ""))
    ref = read_reference_lookup(reference_root, population)
    if ref.empty or "SNP" not in ref.columns:
        return df
    ref = ref.copy()
    if "CHR" in ref.columns:
        ref["REF_CHR"] = ref["CHR"]
    pos_col = next((c for c in ["POS38", "POS", "REF_POS38", "POS19"] if c in ref.columns), None)
    keep = ["SNP"]
    if "REF_CHR" in ref.columns:
        keep.append("REF_CHR")
    if pos_col:
        ref["REF_POS"] = ref[pos_col]
        keep.append("REF_POS")
    if "SNPID" in ref.columns:
        ref["REF_SNPID"] = ref["SNPID"].astype("string")
        keep.append("REF_SNPID")
    if len(keep) == 1:
        return df
    out = df.copy()
    out["SNP"] = out["SNP"].astype("string")
    out["SNPID"] = out["SNPID"].astype("string")
    ref["SNP"] = ref["SNP"].astype("string")
    merged = out.merge(ref[keep].drop_duplicates(), how="left", on="SNP")
    if "REF_SNPID" in ref.columns:
        needs_snpid = merged.get("REF_POS", pd.Series(pd.NA, index=merged.index)).isna() & out["SNPID"].notna()
        if needs_snpid.any():
            by_snpid = out.loc[needs_snpid].merge(
                ref[keep].drop_duplicates(),
                how="left",
                left_on="SNPID",
                right_on="REF_SNPID",
                suffixes=("", "_by_snpid"),
            )
            fill_index = merged.index[needs_snpid]
            for col in ["REF_CHR", "REF_POS"]:
                by_col = f"{col}_by_snpid"
                source_col = by_col if by_col in by_snpid.columns else col
                if col in merged.columns and source_col in by_snpid.columns:
                    merged.loc[fill_index, col] = by_snpid[source_col].values
    if "REF_CHR" in merged.columns:
        merged["CHR"] = merged["CHR"].astype("object")
        missing_chr = merged["CHR"].isna() | merged["CHR"].astype("string").str.lower().isin(["", "nan", "none", "na"])
        merged.loc[missing_chr, "CHR"] = merged.loc[missing_chr, "REF_CHR"]
    if "REF_POS" in merged.columns:
        missing_pos = merged["POS"].isna()
        merged.loc[missing_pos, "POS"] = merged.loc[missing_pos, "REF_POS"]
    return merged[[c for c in ["CHR", "POS", "P", "SNP", "SNPID"] if c in merged.columns]]


def clean_qq_df(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["P"] = numeric_series(out["P"])
    out = out[(out["P"] > 0) & (out["P"] <= 1)]
    out["LOGP"] = -np.log10(out["P"].astype(float))
    return out


def clean_manhattan_df(df: pd.DataFrame) -> pd.DataFrame:
    out = clean_qq_df(df)
    out["POS"] = numeric_series(out["POS"])
    out["CHR_CLEAN"] = out["CHR"].map(normalize_chromosome)
    out = out[out["POS"].notna() & out["CHR_CLEAN"].notna()]
    return out


def qqplot(ax, df: pd.DataFrame, title: str) -> None:
    df = clean_qq_df(df)
    if df.empty:
        ax.text(0.5, 0.5, "No valid P values", ha="center", va="center")
        ax.set_title(title)
        return
    df = select_plot_subset(df, MAX_QQ_ROWS)
    observed = np.sort(df["LOGP"].to_numpy())
    expected = -np.log10((np.arange(1, len(observed) + 1) - 0.5) / len(observed))
    ax.scatter(expected, observed, s=6, alpha=0.65, linewidths=0)
    x_limit = max(float(np.nanmax(expected)), 1)
    y_limit = max(float(np.nanmax(observed)), x_limit)
    ax.plot([0, x_limit], [0, x_limit], color="black", linewidth=0.8)
    ax.set_xlim(0, x_limit * 1.02)
    ax.set_ylim(0, y_limit * 1.02)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(title)


def chrom_key(value: Any) -> tuple[int, str]:
    text = normalize_chromosome(value)
    if pd.isna(text):
        return 999, ""
    text = str(text)
    mapping = {"X": 23, "Y": 24, "MT": 25, "M": 25}
    if text in mapping:
        return mapping[text], text
    try:
        return int(float(text)), text
    except ValueError:
        return 999, text


def manhattan(ax, df: pd.DataFrame, title: str) -> None:
    df = clean_manhattan_df(df)
    if df.empty:
        ax.text(0.5, 0.5, "No valid P/position values", ha="center", va="center")
        ax.set_title(title)
        return
    df = select_plot_subset(df, MAX_MANHATTAN_ROWS)
    df["CHR_ORDER"] = df["CHR_CLEAN"].map(chrom_key)
    df = df.sort_values(["CHR_ORDER", "POS"])
    current = 0
    ticks = []
    labels = []
    colors = ["#111111", "#8a8a8a"]
    for index, (chrom, group) in enumerate(df.groupby("CHR_CLEAN", sort=False)):
        xpos = group["POS"] + current
        ax.scatter(xpos, group["LOGP"], c=colors[index % 2], s=3, alpha=0.7, linewidths=0)
        ticks.append(float(current + group["POS"].max() / 2))
        labels.append(str(chrom))
        current += float(group["POS"].max()) + 1_000_000
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=0, fontsize=6)
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(P)")
    ax.set_title(title)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--standardise-record", required=True)
    parser.add_argument("--force", default="false")
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    record = read_json(args.standardise_record)
    plot_path = Path(record["plot_tiff"])
    dry_run = parse_bool(args.dry_run)
    force = parse_bool(args.force)

    if dry_run:
        record["plot_created"] = False
        write_json(args.output, record)
        return
    if plot_path.exists() and plot_path.stat().st_size > 0 and not force:
        record["plot_created"] = False
        write_json(args.output, record)
        return

    raw_df = extract_plot_df(record["raw_path"], record)
    std_df = extract_plot_df(record["standard_vcf"], record)

    plot_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)
    qqplot(axes[0, 0], raw_df, "Raw QQ")
    manhattan(axes[0, 1], raw_df, "Raw Manhattan")
    qqplot(axes[1, 0], std_df, "Standardised QQ")
    manhattan(axes[1, 1], std_df, "Standardised Manhattan")
    fig.suptitle(f"{record.get('study_id', '')} | {record.get('phenotype', '')} | {record.get('ancestry', '')}", fontsize=11)
    fig.savefig(plot_path, format="tiff", dpi=120, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    record["plot_created"] = True
    write_json(args.output, record)


if __name__ == "__main__":
    main()
