#!/usr/bin/env python3
"""Standardise a staged GWAS file to a lightweight GWAS-VCF representation."""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gwas_pipeline_common import (
    ALT_ALIASES,
    BETA_ALIASES,
    CHR_ALIASES,
    EA_ALIASES,
    EAF_ALIASES,
    LOGP_ALIASES,
    N_ALIASES,
    OA_ALIASES,
    P_ALIASES,
    POS19_ALIASES,
    POS38_ALIASES,
    POS_ALIASES,
    REF_ALIASES,
    RSID_ALIASES,
    SE_ALIASES,
    SNPID_ALIASES,
    ancestry_to_population,
    command_exists,
    ensure_dir,
    is_gwas_vcf,
    load_yaml,
    logp_to_p,
    merge_dicts,
    normalize_chromosome,
    numeric_series,
    parse_bool,
    pick_existing_column,
    read_json,
    read_reference_lookup,
    read_table,
    resolve_study_config,
    utc_now,
    vcf_records,
    write_json,
)


CANONICAL_COLUMNS = [
    "CHR",
    "POS38",
    "POS19",
    "SNP",
    "SNPID",
    "EA",
    "OA",
    "EAF",
    "BETA",
    "SE",
    "P",
    "N",
]

VCF_WRITE_CHUNK = 500_000


def resolve_config_path(path_value: str, config: dict[str, Any]) -> Path:
    path = Path(os.path.expandvars(path_value))
    if path.is_absolute():
        return path
    config_path = config.get("_config_path")
    if config_path:
        return Path(config_path).parent / path
    return path


def apply_configured_exclusions(df: pd.DataFrame, config: dict[str, Any], qc: dict[str, Any]) -> pd.DataFrame:
    exclusion = config.get("exclude") or config.get("exclusions")
    if not exclusion:
        return df
    exclusions = exclusion if isinstance(exclusion, list) else [exclusion]
    out = df
    for item in exclusions:
        if not isinstance(item, dict) or not item.get("file"):
            continue
        path = resolve_config_path(str(item["file"]), config)
        key = item.get("key", "SNP")
        if not path.exists() or key not in out.columns:
            qc.setdefault("warnings", []).append(f"exclusion skipped: {path} key={key}")
            continue
        ex = read_table(path, delimiter=item.get("delimiter", "\t"))
        if key not in ex.columns:
            qc.setdefault("warnings", []).append(f"exclusion key not found in {path}: {key}")
            continue
        before = len(out)
        out = out[~out[key].astype(str).isin(set(ex[key].astype(str)))]
        qc.setdefault("exclusions", []).append({"file": str(path), "key": key, "removed": before - len(out)})
    return out


def apply_configured_joins(df: pd.DataFrame, config: dict[str, Any], qc: dict[str, Any]) -> pd.DataFrame:
    joins = config.get("joins", [])
    out = df
    for item in joins:
        if not isinstance(item, dict) or not item.get("file"):
            continue
        path = resolve_config_path(str(item["file"]), config)
        on = item.get("on")
        left_on = item.get("left_on", on)
        right_on = item.get("right_on", on)
        if not path.exists() or not left_on or left_on not in out.columns:
            qc.setdefault("warnings", []).append(f"join skipped: {path}")
            continue
        join_df = read_table(path, delimiter=item.get("delimiter", "\t"))
        if right_on not in join_df.columns:
            qc.setdefault("warnings", []).append(f"join key not found in {path}: {right_on}")
            continue
        before_cols = set(out.columns)
        out = pd.merge(out, join_df, how=item.get("how", "left"), left_on=left_on, right_on=right_on, suffixes=("", "_join"))
        qc.setdefault("joins", []).append(
            {"file": str(path), "left_on": left_on, "right_on": right_on, "added_columns": sorted(set(out.columns) - before_cols)}
        )
    return out


def column_has_values(series: pd.Series) -> bool:
    values = series.astype("string").str.strip()
    missing = values.isna() | values.str.lower().isin(["", "nan", "none", "na", "."])
    return bool((~missing).any())


def select_column(df: pd.DataFrame, config: dict[str, Any], canonical: str, aliases: list[str]) -> str | None:
    configured = (config.get("columns") or {}).get(canonical)
    if configured and configured in df.columns:
        return configured
    for alias in aliases:
        if alias in df.columns and column_has_values(df[alias]):
            return alias
    return pick_existing_column(df.columns, aliases)


def normalise_sample_size(series: pd.Series, qc: dict[str, Any]) -> pd.Series:
    values = numeric_series(series)
    present = values.notna()
    fractional = present & ((values - np.floor(values)) != 0)
    if fractional.any():
        qc["sample_size_values_rounded"] = int(fractional.sum())
    rounded = values.copy()
    rounded.loc[present] = np.floor(values.loc[present] + 0.5)
    return rounded


def normalize_delimiter_value(value: Any) -> str | None:
    text = str(value or "").strip()
    if not text:
        return None
    lower = text.lower()
    if lower in {"tab", "\\t", "t"}:
        return "\t"
    if lower in {"space", "spaces", "whitespace", "white_space", "\\s+", "s+"}:
        return r"\s+"
    if lower in {"comma", "csv"}:
        return ","
    return text


def normalize_build_value(value: Any) -> str:
    text = str(value or "").strip().lower()
    if not text:
        return ""
    text = text.replace("_", "").replace("-", "")
    if text in {"19", "hg19", "grch37", "37", "b37", "build37"}:
        return "hg19"
    if text in {"38", "hg38", "grch38", "38", "b38", "build38"}:
        return "hg38"
    return ""


def recursive_metadata_value(payload: Any, keys: set[str]) -> Any:
    if isinstance(payload, dict):
        for key, value in payload.items():
            if str(key).lower() in keys and value not in [None, ""]:
                return value
        for value in payload.values():
            found = recursive_metadata_value(value, keys)
            if found not in [None, ""]:
                return found
    if isinstance(payload, list):
        for item in payload:
            found = recursive_metadata_value(item, keys)
            if found not in [None, ""]:
                return found
    return ""


def configured_or_metadata_build(record: dict[str, Any], config: dict[str, Any]) -> tuple[str, str]:
    configured = normalize_build_value(config.get("input_build") or record.get("input_build"))
    if configured:
        return configured, "configured"
    metadata_build = normalize_build_value(
        recursive_metadata_value(record.get("metadata", {}), {"build", "genome_build", "assembly", "coordinate_build"})
    )
    if metadata_build:
        return metadata_build, "metadata"
    return "", ""


def prepare_reference_for_merge(ref: pd.DataFrame, population: str) -> pd.DataFrame:
    ref = ref.copy()
    rename = {}
    for src, dst in [
        ("rsid", "SNP"),
        ("rsID", "SNP"),
        ("CHR", "REF_CHR"),
        ("POS19", "REF_POS19"),
        ("POS38", "REF_POS38"),
        ("SNPID", "REF_SNPID"),
    ]:
        if src in ref.columns and dst not in ref.columns:
            rename[src] = dst
    ref = ref.rename(columns=rename)
    eaf_col = next((c for c in [f"EAF_{population}", "EAF", "AF"] if c in ref.columns), None)
    if eaf_col and eaf_col != "REF_EAF":
        ref = ref.rename(columns={eaf_col: "REF_EAF"})
    return ref


def infer_build_from_reference(out: pd.DataFrame, pos: pd.Series, ref: pd.DataFrame, qc: dict[str, Any]) -> str:
    if ref.empty or "SNP" not in ref.columns:
        return ""
    candidates = pd.DataFrame({"SNP": out["SNP"], "SNPID": out["SNPID"], "POS": numeric_series(pos)})
    candidates = candidates[candidates["POS"].notna()].copy()
    if candidates.empty:
        return ""
    if len(candidates) > 20000:
        candidates = candidates.sample(n=20000, random_state=1)
    candidates["SNP"] = candidates["SNP"].astype("string")
    candidates["SNPID"] = candidates["SNPID"].astype("string")
    ref = ref.copy()
    ref["SNP"] = ref["SNP"].astype("string")
    if "REF_SNPID" in ref.columns:
        ref["REF_SNPID"] = ref["REF_SNPID"].astype("string")
    keep = [c for c in ["SNP", "REF_SNPID", "REF_POS19", "REF_POS38"] if c in ref.columns]
    ref_by_snp = ref[keep].dropna(subset=["SNP"]).drop_duplicates(subset=["SNP"])
    merged = candidates.merge(ref_by_snp, how="left", on="SNP")
    if "REF_SNPID" in ref.columns:
        missing_reference = merged.get("REF_POS19", pd.Series(pd.NA, index=merged.index)).isna() & merged["SNPID"].notna()
        if missing_reference.any():
            ref_by_snpid = ref[keep].dropna(subset=["REF_SNPID"]).drop_duplicates(subset=["REF_SNPID"])
            by_snpid = merged.loc[missing_reference, ["SNPID"]].merge(
                ref_by_snpid,
                how="left",
                left_on="SNPID",
                right_on="REF_SNPID",
                suffixes=("", "_by_snpid"),
            )
            fill_index = merged.index[missing_reference]
            for col in ["REF_POS19", "REF_POS38"]:
                source_col = f"{col}_by_snpid" if f"{col}_by_snpid" in by_snpid.columns else col
                if col in merged.columns and source_col in by_snpid.columns:
                    merged.loc[fill_index, col] = by_snpid[source_col].values
    matches: dict[str, int] = {}
    for build, ref_col in [("hg19", "REF_POS19"), ("hg38", "REF_POS38")]:
        if ref_col in merged.columns:
            matches[build] = int((numeric_series(merged["POS"]) == numeric_series(merged[ref_col])).sum())
        else:
            matches[build] = 0
    min_matches = min(10, max(1, int(len(candidates) * 0.01)))
    qc["build_detection"] = {"method": "reference_match", "matches": matches, "sampled_rows": int(len(candidates))}
    if matches["hg19"] >= min_matches and matches["hg19"] > matches["hg38"]:
        qc["build_detection"]["build"] = "hg19"
        return "hg19"
    if matches["hg38"] >= min_matches and matches["hg38"] > matches["hg19"]:
        qc["build_detection"]["build"] = "hg38"
        return "hg38"
    return ""


def manifest_config_overrides(record: dict[str, Any]) -> dict[str, Any]:
    overrides: dict[str, Any] = {}
    if str(record.get("input_build", "")).strip():
        overrides["input_build"] = record["input_build"]
    delimiter = normalize_delimiter_value(record.get("delimiter"))
    if delimiter:
        overrides["delimiter"] = delimiter
    if str(record.get("population", "")).strip():
        overrides["reference_population"] = str(record["population"]).strip().upper()
    if str(record.get("format", "")).strip():
        overrides["format"] = str(record["format"]).strip()
    return overrides


def normalize_variant_identifiers(out: pd.DataFrame, qc: dict[str, Any]) -> pd.DataFrame:
    out = out.copy()
    snp = out["SNP"].astype("string").str.strip()
    snpid = out["SNPID"].astype("string").str.strip()
    empty_snp = snp.isna() | snp.str.lower().isin(["", "nan", "none", "na", "."])
    empty_snpid = snpid.isna() | snpid.str.lower().isin(["", "nan", "none", "na", "."])

    original_snp = snp.mask(empty_snp, pd.NA)
    out.loc[empty_snpid & original_snp.notna(), "SNPID"] = original_snp[empty_snpid & original_snp.notna()]

    rs_from_snp = original_snp.str.extract(r"\b(rs\d+)\b", flags=re.I, expand=False).str.lower()
    rs_from_snpid = snpid.mask(empty_snpid, pd.NA).str.extract(r"\b(rs\d+)\b", flags=re.I, expand=False).str.lower()
    extracted = rs_from_snp.fillna(rs_from_snpid)
    has_extracted = extracted.notna()
    changed = has_extracted & (snp.str.lower() != extracted)
    out.loc[has_extracted, "SNP"] = extracted[has_extracted]

    final_snp = out["SNP"].astype("string").str.strip()
    invalid = final_snp.isna() | ~final_snp.str.lower().str.startswith("rs", na=False)
    out.loc[invalid, "SNP"] = pd.NA

    qc["variant_id_normalisation"] = {
        "rsid_extracted": int(has_extracted.sum()),
        "compound_snpid_preserved": int((empty_snpid & original_snp.notna()).sum()),
        "compound_snp_values_changed": int(changed.sum()),
    }
    return out


def canonicalise_table(df: pd.DataFrame, record: dict[str, Any], config: dict[str, Any], reference_root: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    qc: dict[str, Any] = {"warnings": [], "reference_fills": {}, "input_rows": int(len(df))}

    transforms = config.get("transforms", {})
    if transforms.get("strip_hash_header", True) and len(df.columns):
        df.columns = [str(c).lstrip("#") for c in df.columns]
        first = df.columns[0]
        if first in df.columns:
            df[first] = df[first].astype(str).str.lstrip("#")

    df = apply_configured_exclusions(df, config, qc)
    df = apply_configured_joins(df, config, qc)

    source_cols = {
        "CHR": select_column(df, config, "CHR", CHR_ALIASES),
        "POS": select_column(df, config, "POS", POS_ALIASES),
        "POS19": select_column(df, config, "POS19", POS19_ALIASES),
        "POS38": select_column(df, config, "POS38", POS38_ALIASES),
        "SNP": select_column(df, config, "SNP", RSID_ALIASES),
        "SNPID": select_column(df, config, "SNPID", SNPID_ALIASES),
        "EA": select_column(df, config, "EA", EA_ALIASES),
        "OA": select_column(df, config, "OA", OA_ALIASES),
        "EAF": select_column(df, config, "EAF", EAF_ALIASES),
        "BETA": select_column(df, config, "BETA", BETA_ALIASES),
        "SE": select_column(df, config, "SE", SE_ALIASES),
        "P": select_column(df, config, "P", P_ALIASES),
        "LOG10P": select_column(df, config, "LOG10P", LOGP_ALIASES),
        "N": select_column(df, config, "N", N_ALIASES),
        "REF": select_column(df, config, "REF", REF_ALIASES),
        "ALT": select_column(df, config, "ALT", ALT_ALIASES),
    }
    qc["source_columns"] = {k: v for k, v in source_cols.items() if v}

    out = pd.DataFrame(index=df.index)
    for canonical in CANONICAL_COLUMNS:
        source = source_cols.get(canonical)
        if source:
            out[canonical] = df[source]
        else:
            out[canonical] = np.nan
    for col in ["CHR", "SNP", "SNPID", "EA", "OA"]:
        out[col] = out[col].astype("string")
    for col in ["POS38", "POS19", "EAF", "BETA", "SE", "P", "N"]:
        out[col] = numeric_series(out[col])
    out["N"] = normalise_sample_size(out["N"], qc)
    out = normalize_variant_identifiers(out, qc)

    population = config.get("reference_population") or ancestry_to_population(record.get("ancestry", ""))
    ref = read_reference_lookup(reference_root, population)
    qc["reference_population"] = population
    if not ref.empty:
        ref = prepare_reference_for_merge(ref, population)
    else:
        qc["warnings"].append(f"reference lookup not found for population {population}")

    if source_cols["POS38"]:
        out["POS38"] = numeric_series(df[source_cols["POS38"]])
    if source_cols["POS19"]:
        out["POS19"] = numeric_series(df[source_cols["POS19"]])
    if out["POS38"].isna().all() and source_cols["POS"]:
        input_build, build_method = configured_or_metadata_build(record, config)
        if not input_build:
            input_build = infer_build_from_reference(out, df[source_cols["POS"]], ref, qc)
            build_method = "reference_match" if input_build else "default"
        if input_build == "hg19":
            out["POS19"] = numeric_series(df[source_cols["POS"]])
        else:
            out["POS38"] = numeric_series(df[source_cols["POS"]])
        qc["build_detection"] = {**qc.get("build_detection", {}), "method": build_method, "build": input_build or "hg38"}

    if out["P"].isna().all() and source_cols["LOG10P"]:
        out["P"] = df[source_cols["LOG10P"]].map(logp_to_p)
        qc["derived_p_from_log10p"] = True

    if transforms.get("derive_oa_from_ref_alt", False) or out["OA"].isna().all():
        ref_col = source_cols.get("REF")
        alt_col = source_cols.get("ALT")
        ea_col = source_cols.get("EA")
        if ref_col and alt_col and ea_col:
            ea_allele = df[ea_col].astype(str).str.upper()
            ref_allele = df[ref_col].astype(str).str.upper()
            alt_allele = df[alt_col].astype(str).str.upper()
            out["OA"] = np.where(
                ea_allele == ref_allele,
                alt_allele,
                np.where(ea_allele == alt_allele, ref_allele, out["OA"]),
            )
            qc["derived_oa_from_ref_alt"] = True

    if transforms.get("uppercase_alleles", True):
        for col in ["EA", "OA"]:
            out[col] = out[col].astype("string").str.upper()

    if not ref.empty:
        out = fill_from_reference(out, ref, population, qc)

    if out["POS38"].isna().all() and out["POS19"].notna().any():
        out["POS38"] = out["POS19"]
        qc["warnings"].append("liftover unavailable; POS38 copied from POS19")
    if out["POS19"].isna().all() and out["POS38"].notna().any():
        out["POS19"] = out["POS38"]
        qc["warnings"].append("POS19 unavailable; copied from POS38 for traceability")

    out["CHR"] = out["CHR"].map(normalize_chromosome)

    extras = config.get("extras", [])
    for extra in extras:
        if extra in df.columns and extra not in out.columns:
            out[extra] = df[extra]

    out["phenotype"] = record.get("phenotype", "")
    out["ancestry"] = record.get("ancestry", "")
    out["study_id"] = record.get("study_id", "")

    before = len(out)
    out = out.dropna(subset=["CHR", "POS38", "EA", "OA", "BETA", "SE", "P"], how="any")
    qc["rows_removed_missing_required"] = int(before - len(out))
    qc["output_rows"] = int(len(out))
    return out, qc


def fill_from_reference(out: pd.DataFrame, ref: pd.DataFrame, population: str, qc: dict[str, Any]) -> pd.DataFrame:
    # ref is already normalised by prepare_reference_for_merge in the caller.
    out = out.copy()
    eaf_col = "REF_EAF" if "REF_EAF" in ref.columns else None
    keep = [c for c in ["SNP", "REF_SNPID", "REF_CHR", "REF_POS19", "REF_POS38", eaf_col] if c]
    if "SNP" not in ref.columns or not keep:
        qc["warnings"].append("reference lookup lacks SNP/position columns")
        return out
    out["SNP"] = out["SNP"].astype("string")
    out["SNPID"] = out["SNPID"].astype("string")
    ref["SNP"] = ref["SNP"].astype("string")
    if "REF_SNPID" in ref.columns:
        ref["REF_SNPID"] = ref["REF_SNPID"].astype("string")

    ref_by_snp = ref[keep].dropna(subset=["SNP"]).drop_duplicates(subset=["SNP"])
    merged = out.merge(ref_by_snp, how="left", left_on="SNP", right_on="SNP", suffixes=("", "_ref"))
    if "REF_SNPID" in ref.columns:
        missing = merged["SNP"].isna() | ~merged["SNP"].astype(str).str.startswith("rs", na=False)
        if missing.any():
            ref_by_snpid = ref[keep].dropna(subset=["REF_SNPID"]).drop_duplicates(subset=["REF_SNPID"])
            by_snpid = out.loc[missing].merge(
                ref_by_snpid,
                how="left",
                left_on="SNPID",
                right_on="REF_SNPID",
                suffixes=("", "_ref2"),
            )
            fill_index = merged.index[missing]
            if "SNP_ref2" in by_snpid.columns:
                merged.loc[fill_index, "SNP"] = by_snpid["SNP_ref2"].values
                qc["reference_fills"]["rsid_from_snpid"] = int(pd.Series(by_snpid["SNP_ref2"]).notna().sum())
            for ref_col in ["REF_POS19", "REF_POS38", "REF_EAF"]:
                ref2 = f"{ref_col}_ref2"
                source_col = ref2 if ref2 in by_snpid.columns else ref_col
                if source_col in by_snpid.columns and ref_col in merged.columns:
                    merged.loc[fill_index, ref_col] = by_snpid[source_col].values

    if "REF_CHR" in merged.columns:
        merged["CHR"] = merged["CHR"].astype("object")
        missing = merged["CHR"].isna() | merged["CHR"].astype("string").str.lower().isin(["", "nan", "none", "na"])
        merged.loc[missing, "CHR"] = merged.loc[missing, "REF_CHR"]
        qc["reference_fills"]["CHR"] = int(missing.sum())
    for target, ref_col in [("POS19", "REF_POS19"), ("POS38", "REF_POS38")]:
        if ref_col in merged.columns:
            missing = merged[target].isna()
            merged.loc[missing, target] = merged.loc[missing, ref_col]
            qc["reference_fills"][target] = int(missing.sum())
    if eaf_col and eaf_col in merged.columns:
        missing = merged["EAF"].isna()
        merged.loc[missing, "EAF"] = merged.loc[missing, eaf_col]
        qc["reference_fills"]["EAF"] = int(missing.sum())

    drop_cols = [c for c in merged.columns if c.startswith("REF_") or c == eaf_col]
    return merged.drop(columns=drop_cols, errors="ignore")


def chrom_sort_value(chrom: Any) -> tuple[int, str]:
    text = normalize_chromosome(chrom)
    if pd.isna(text):
        return (999, "")
    text = str(text)
    mapping = {"X": 23, "Y": 24, "MT": 25, "M": 25}
    if text in mapping:
        return (mapping[text], text)
    try:
        return (int(float(text)), text)
    except ValueError:
        return (999, text)


def _vcf_float_series(series: pd.Series, index: pd.Index) -> pd.Series:
    arr = pd.to_numeric(series, errors="coerce").to_numpy(dtype="float64")
    out = np.full(arr.shape, ".", dtype=object)
    finite = np.isfinite(arr)
    out[finite] = [f"{value:.8g}" for value in arr[finite]]
    return pd.Series(out, index=index)


def _vcf_int_series(series: pd.Series, index: pd.Index) -> pd.Series:
    arr = pd.to_numeric(series, errors="coerce").to_numpy(dtype="float64")
    out = np.full(arr.shape, ".", dtype=object)
    finite = np.isfinite(arr)
    out[finite] = np.floor(arr[finite] + 0.5).astype("int64").astype(str)
    return pd.Series(out, index=index)


def _vcf_allele_series(series: pd.Series, fallback: str = "N") -> pd.Series:
    alleles = series.astype("string").str.strip().str.upper()
    invalid = alleles.isna() | alleles.isin(["", "NAN", "NONE", "."])
    alleles = alleles.mask(invalid, fallback).str.replace(r"[^ACGTN]", "N", regex=True)
    return alleles.astype(object)


def build_vcf_lines(df: pd.DataFrame) -> pd.Series:
    """Build VCF body lines as a vectorised string Series (one entry per variant)."""
    index = df.index
    chrom = df["CHR"].map(normalize_chromosome)
    chrom = chrom.where(chrom.notna(), ".").astype(object)
    pos = numeric_series(df["POS38"]).astype("int64").astype(str).astype(object)

    snp = df["SNP"].astype("string").str.strip()
    snp = snp.mask(snp.isna() | (snp.str.lower() == "nan") | (snp == ""), ".").astype(object)

    ref = _vcf_allele_series(df["OA"])
    alt = _vcf_allele_series(df["EA"])

    pos19 = numeric_series(df["POS19"])
    pos19_info = pd.Series("POS19=.", index=index, dtype=object)
    has_pos19 = pos19.notna()
    pos19_info[has_pos19] = "POS19=" + pos19[has_pos19].astype("int64").astype(str)

    snpid = df["SNPID"]
    snpid_info = pd.Series("SNPID=.", index=index, dtype=object)
    has_snpid = snpid.notna()
    snpid_info[has_snpid] = "SNPID=" + snpid[has_snpid].astype(str)

    info = pos19_info + ";" + snpid_info + ";EAF_SOURCE=input"

    p = pd.to_numeric(df["P"], errors="coerce").to_numpy(dtype="float64")
    lp = np.full(p.shape, np.nan)
    valid_p = np.isfinite(p) & (p > 0)
    lp[valid_p] = -np.log10(p[valid_p])

    sample = (
        _vcf_float_series(df["BETA"], index) + ":"
        + _vcf_float_series(df["SE"], index) + ":"
        + _vcf_float_series(pd.Series(lp, index=index), index) + ":"
        + _vcf_float_series(df["EAF"], index) + ":"
        + _vcf_int_series(df["N"], index)
    )

    return (
        chrom + "\t" + pos + "\t" + snp + "\t" + ref + "\t" + alt
        + "\t.\tPASS\t" + info + "\tES:SE:LP:AF:SS\t" + sample
    )


def write_gwas_vcf(df: pd.DataFrame, path: Path, record: dict[str, Any], qc: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df = df.copy()
    df["CHR_SORT"] = df["CHR"].map(chrom_sort_value)
    df["POS38"] = numeric_series(df["POS38"])
    df = df.dropna(subset=["POS38"]).sort_values(["CHR_SORT", "POS38"])

    sample = record.get("phenotype") or "GWAS"
    header = [
        "##fileformat=VCFv4.2",
        "##source=000_GWAS_formatting",
        "##GWASVCF=1.0",
        f"##fileDate={utc_now()}",
        *gwas_metadata_header_lines(record),
        '##INFO=<ID=POS19,Number=1,Type=Integer,Description="GRCh37/hg19 position">',
        '##INFO=<ID=SNPID,Number=1,Type=String,Description="Source SNP identifier when available">',
        '##INFO=<ID=EAF_SOURCE,Number=1,Type=String,Description="EAF source: input or 1000G reference">',
        '##FORMAT=<ID=ES,Number=1,Type=Float,Description="Effect size">',
        '##FORMAT=<ID=SE,Number=1,Type=Float,Description="Standard error">',
        '##FORMAT=<ID=LP,Number=1,Type=Float,Description="-log10 P-value">',
        '##FORMAT=<ID=AF,Number=1,Type=Float,Description="Effect allele frequency">',
        '##FORMAT=<ID=SS,Number=1,Type=Float,Description="Sample size">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + str(sample),
    ]

    lines = build_vcf_lines(df)

    plain_path = path.with_suffix("")
    handle_context = open(plain_path, "wt", encoding="utf-8") if command_exists("bgzip") else gzip.open(path, "wt", encoding="utf-8")
    with handle_context as handle:
        handle.write("\n".join(header) + "\n")
        for start in range(0, len(lines), VCF_WRITE_CHUNK):
            handle.write("\n".join(lines.iloc[start:start + VCF_WRITE_CHUNK].tolist()))
            handle.write("\n")
    if command_exists("bgzip"):
        with open(path, "wb") as out_handle:
            subprocess.run(["bgzip", "-c", str(plain_path)], check=True, stdout=out_handle)
        plain_path.unlink(missing_ok=True)


def gwas_metadata_header_lines(record: dict[str, Any]) -> list[str]:
    return [
        f"##study_id={record.get('study_id', '')}",
        f"##phenotype={record.get('phenotype', '')}",
        f"##ancestry={record.get('ancestry', '')}",
        f"##sex={record.get('sex', '')}",
        f"##PMID={record.get('PMID', '')}",
    ]


def normalise_vcf_sample_size_line(line: str) -> tuple[str, bool]:
    if line.startswith("#") or not line.strip():
        return line, False
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 10:
        return line, False
    format_keys = fields[8].split(":")
    if "SS" not in format_keys:
        return line, False
    ss_index = format_keys.index("SS")
    changed = False
    for field_index in range(9, len(fields)):
        values = fields[field_index].split(":")
        if ss_index >= len(values):
            continue
        original = values[ss_index]
        rounded = format_integer(original)
        if rounded != "." and rounded != original:
            values[ss_index] = rounded
            fields[field_index] = ":".join(values)
            changed = True
    return "\t".join(fields) + "\n", changed


def copy_vcf_with_metadata_headers(raw_path: Path, vcf_path: Path, record: dict[str, Any], qc: dict[str, Any]) -> None:
    managed_prefixes = ("##study_id=", "##phenotype=", "##ancestry=", "##sex=", "##PMID=")
    opener = gzip.open if str(raw_path).endswith(".gz") else open
    plain_path = vcf_path.with_suffix("")
    output_context = open(plain_path, "wt", encoding="utf-8") if command_exists("bgzip") else gzip.open(vcf_path, "wt", encoding="utf-8")
    inserted = False
    rounded_sample_sizes = 0
    with opener(raw_path, "rt", encoding="utf-8", errors="replace") as in_handle, output_context as out_handle:
        for line in in_handle:
            if line.startswith(managed_prefixes):
                continue
            if not inserted and line.startswith("#CHROM"):
                out_handle.write("\n".join(gwas_metadata_header_lines(record)) + "\n")
                inserted = True
            line, rounded = normalise_vcf_sample_size_line(line)
            rounded_sample_sizes += int(rounded)
            out_handle.write(line)
        if not inserted:
            out_handle.write("\n".join(gwas_metadata_header_lines(record)) + "\n")
    if rounded_sample_sizes:
        qc["sample_size_values_rounded"] = rounded_sample_sizes
    if command_exists("bgzip"):
        with open(vcf_path, "wb") as out_handle:
            subprocess.run(["bgzip", "-c", str(plain_path)], check=True, stdout=out_handle)
        plain_path.unlink(missing_ok=True)


def format_integer(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "."
    if not math.isfinite(number):
        return "."
    return str(int(math.floor(number + 0.5)))


def create_index(vcf_path: Path, index_path: Path, qc: dict[str, Any]) -> None:
    if command_exists("tabix") and command_exists("bgzip"):
        subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_path)], check=True)
        if index_path.exists():
            qc["index_method"] = "tabix"
            return
    index_path.write_text("Placeholder index. Install bgzip/tabix to create a binary tabix index.\n", encoding="utf-8")
    qc["index_method"] = "placeholder"
    qc.setdefault("warnings", []).append("binary tabix index not created because bgzip/tabix is unavailable")


def detect_vcf_build(path: str | Path) -> str:
    """Return 'GRCh38'/'GRCh37' from a VCF header, or '' when it cannot be determined."""
    opener = gzip.open if str(path).endswith(".gz") else open
    assembly_text = ""
    contig1_length: int | None = None
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("#"):
                    break
                if line.startswith("##contig"):
                    match = re.search(r"assembly=([^,>]+)", line, flags=re.I)
                    if match and not assembly_text:
                        assembly_text = match.group(1)
                    length = re.search(r"ID=(?:chr)?1,[^>]*length=(\d+)", line, flags=re.I)
                    if length and contig1_length is None:
                        contig1_length = int(length.group(1))
                elif line.startswith(("##reference", "##genomeBuild", "##genome_build", "##build")) and not assembly_text:
                    assembly_text = line.split("=", 1)[-1].strip()
    except OSError:
        return ""
    text = assembly_text.lower()
    if any(token in text for token in ("grch38", "hg38", "ncbi38", "b38", "build38")):
        return "GRCh38"
    if any(token in text for token in ("grch37", "hg19", "ncbi37", "b37", "build37")):
        return "GRCh37"
    if contig1_length == 248956422:
        return "GRCh38"
    if contig1_length == 249250621:
        return "GRCh37"
    return ""


def standardise_existing_vcf(raw_path: Path, output_prefix: Path, record: dict[str, Any], qc: dict[str, Any]) -> tuple[Path, Path]:
    vcf_path = output_prefix.with_suffix(".gwas.vcf.gz")
    index_path = Path(str(vcf_path) + ".tbi")
    copy_vcf_with_metadata_headers(raw_path, vcf_path, record, qc)
    create_index(vcf_path, index_path, qc)
    qc["input_format"] = "gwas_vcf"
    qc["metadata_headers_normalised"] = True
    return vcf_path, index_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-record", required=True)
    parser.add_argument("--config-root", required=True)
    parser.add_argument("--study-yaml", default="")
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--force", default="false")
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    record = read_json(args.stage_record)
    output_prefix = Path(record["output_prefix"])
    ensure_dir(output_prefix.parent)
    standard_vcf = output_prefix.with_suffix(".gwas.vcf.gz")
    standard_index = Path(str(standard_vcf) + ".tbi")
    metadata_json = output_prefix.with_suffix(".metadata.json")
    qc_json = output_prefix.with_suffix(".qc.json")

    force = parse_bool(args.force)
    dry_run = parse_bool(args.dry_run)
    config = resolve_study_config(args.config_root, record["study_id"], None, record.get("study_yaml") or None)
    config = merge_dicts(config, manifest_config_overrides(record))
    if args.study_yaml:
        config = merge_dicts(config, resolve_study_config(args.config_root, record["study_id"], args.study_yaml, None))
    qc: dict[str, Any] = {"warnings": [], "started_at": utc_now(), "dry_run": dry_run, "config": config.get("_config_path", "")}

    if dry_run:
        result = build_record(record, standard_vcf, standard_index, metadata_json, qc_json, output_prefix.with_suffix(".tiff"), qc)
        write_json(args.output, result)
        return

    if standard_vcf.exists() and standard_index.exists() and metadata_json.exists() and qc_json.exists() and not force:
        result = build_record(record, standard_vcf, standard_index, metadata_json, qc_json, output_prefix.with_suffix(".tiff"), read_json(qc_json))
        write_json(args.output, result)
        return

    raw_path = Path(record["raw_path"])
    if not raw_path.exists():
        raise FileNotFoundError(f"Staged raw file not found: {raw_path}")

    if is_gwas_vcf(raw_path):
        standardise_existing_vcf(raw_path, output_prefix, record, qc)
        # Passthrough VCFs keep their source coordinates; report the real build
        # rather than assuming the GRCh38 produced by the table writer.
        coordinate_build = detect_vcf_build(raw_path)
        if not coordinate_build:
            qc["warnings"].append("could not detect genome build from source VCF header; coordinate_build set to unknown")
        qc["coordinate_build"] = coordinate_build or "unknown"
    else:
        delimiter = normalize_delimiter_value(config.get("delimiter"))
        df = read_table(raw_path, delimiter=delimiter)
        table, table_qc = canonicalise_table(df, record, config, args.reference_root)
        qc.update(table_qc)
        write_gwas_vcf(table, standard_vcf, record, qc)
        create_index(standard_vcf, standard_index, qc)
        qc["input_format"] = "table"
        coordinate_build = "GRCh38"
        qc["coordinate_build"] = coordinate_build

    metadata = {
        "study_id": record["study_id"],
        "output_id": record["output_id"],
        "phenotype": record["phenotype"],
        "ancestry": record["ancestry"],
        "sex": record.get("sex", ""),
        "author": record["author"],
        "year": record["year"],
        "PMID": record["PMID"],
        "source_type": record["source_type"],
        "GWAS_location": record["GWAS_location"],
        "standard_format": "GWAS-VCF",
        "coordinate_build": coordinate_build,
        "raw_path": record["raw_path"],
    }
    qc["completed_at"] = utc_now()
    write_json(metadata_json, metadata)
    write_json(qc_json, qc)

    write_json(args.output, build_record(record, standard_vcf, standard_index, metadata_json, qc_json, output_prefix.with_suffix(".tiff"), qc))


def build_record(record: dict[str, Any], vcf: Path, index: Path, metadata: Path, qc_json: Path, plot: Path, qc: dict[str, Any]) -> dict[str, Any]:
    out = dict(record)
    out.update(
        {
            "standard_vcf": str(vcf),
            "standard_index": str(index),
            "metadata_json": str(metadata),
            "qc_json": str(qc_json),
            "plot_tiff": str(plot),
            "qc": qc,
        }
    )
    return out


if __name__ == "__main__":
    main()
