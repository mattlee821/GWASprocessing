#!/usr/bin/env python3
"""Standardise one staged GWAS (summary-statistics table or pre-built GWAS-VCF) to an analysis-ready
GWAS-VCF: GRCh38 POS, GRCh37 INFO/POS19, trait metadata in the header, ES/AF relative to ALT.

Usage (called once per manifest row by pipeline/main.nf):
    standardise_gwas.py --stage-record stage.json --config-root manifests --reference-root references \
        [--study-yaml extra.yaml] [--analysis-qc true|false] [--force true|false] [--dry-run true|false] \
        [--workers N] --output record.json

Tables are read in pieces of 1,000,000 rows (GWAS_TABLE_CHUNK overrides); --workers N (or GWAS_WORKERS)
standardises N pieces at a time. Output is the same for any piece size or number of workers.

Exit codes: 0 success; 78 configuration error (stops the Nextflow run); 65 unusable input data.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import multiprocessing
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from array import array
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np
import pandas as pd

import gwas_alleles as alleles
import gwas_coordinates as coords
import gwas_reference as grch38
import gwas_stats
import gwas_steps
import gwas_vcf
from gwas_pipeline_common import (
    ALT_ALIASES,
    BETA_ALIASES,
    CHR_ALIASES,
    CI_LOWER_ALIASES,
    CI_UPPER_ALIASES,
    EA_ALIASES,
    EAF_ALIASES,
    IMPINFO_ALIASES,
    LOGP_ALIASES,
    N_ALIASES,
    NCASE_ALIASES,
    OA_ALIASES,
    OR_ALIASES,
    OUTPUT_SCHEMA_VERSION,
    P_ALIASES,
    POS19_ALIASES,
    POS38_ALIASES,
    POS_ALIASES,
    PRIMARY_CHROMS,
    REF_ALIASES,
    RSID_ALIASES,
    SE_ALIASES,
    SNPID_ALIASES,
    SUPER_POPULATIONS,
    ancestry_to_population,
    complement_alleles,
    config_error,
    EXIT_CONFIG,
    EXIT_DATA,
    data_error,
    download_url,
    ensure_dir,
    infer_delimiter,
    is_gwas_vcf,
    make_work_dir,
    normalize_build_value,
    normalize_chromosome_series,
    normalize_delimiter_value,
    numeric_series,
    parse_bool,
    read_json,
    read_table,
    resolve_ancestry,
    resolve_row_config,
    settings_fingerprint,
    upper_alleles,
    utc_now,
    write_json,
)

FIELD_ALIASES = {
    "CHR": CHR_ALIASES, "POS": POS_ALIASES, "POS19": POS19_ALIASES, "POS38": POS38_ALIASES, "SNP": RSID_ALIASES,
    "SNPID": SNPID_ALIASES, "EA": EA_ALIASES, "OA": OA_ALIASES, "EAF": EAF_ALIASES, "BETA": BETA_ALIASES,
    "SE": SE_ALIASES, "P": P_ALIASES, "LOG10P": LOGP_ALIASES, "N": N_ALIASES, "REF": REF_ALIASES, "ALT": ALT_ALIASES,
    "OR": OR_ALIASES, "CI_LOWER": CI_LOWER_ALIASES, "CI_UPPER": CI_UPPER_ALIASES, "NCASE": NCASE_ALIASES,
    "IMPINFO": IMPINFO_ALIASES,
}
# Allele-dependent columns must come from one family (original or GWAS Catalog harmonised hm_*).
FAMILY_FIELDS = ("EA", "OA", "BETA", "EAF", "OR", "CI_LOWER", "CI_UPPER", "CHR", "POS")
KNOWN_CONFIG_KEYS = {"description", "input_build", "delimiter", "columns", "exclude", "exclusions", "joins",
                     "transforms", "extras", "format", "reference_population", "palindromic_tolerance", "sample",
                     "study_type", "ncase", "ncontrol", "sample_size", "gwas_catalog_file", "_config_path"}
KNOWN_TRANSFORMS = {"strip_hash_header", "uppercase_alleles", "derive_oa_from_ref_alt"}
REQUIRED_FIELDS = ("CHR", "POS38", "EA", "OA", "BETA", "SE", "LP")
RETENTION_WARNING = 0.5
AF_CONCORDANCE_MIN_N = 1000
BUILD_MIN_HITS = 20  # a build detected from the data needs this many matching variants (as GwasDataImport requires)
Z_975 = NormalDist().inv_cdf(0.975)
_inv_cdf = np.vectorize(NormalDist().inv_cdf, otypes=[float])
STRAND_TAG = "GWASPROCESSING_LIFT_STRAND"
MISSING_TOKENS = {"", "nan", "none", "na", "."}


# ---------------------------------------------------------------------------------------------------
# Row metadata and configuration
# ---------------------------------------------------------------------------------------------------

def resolve_trait_metadata(record: dict[str, Any]) -> dict[str, str]:
    """Mandatory header metadata: ancestry must be a super-population code; sex/phenotype get sentinels."""
    raw_ancestry = str(record.get("ancestry", "") or "").strip()
    ancestry = resolve_ancestry(raw_ancestry) or ""
    if ancestry not in SUPER_POPULATIONS:
        raise config_error(
            f"Manifest row {record.get('study_id', '')}/{record.get('output_id', '')} has missing or "
            f"unrecognised ancestry '{raw_ancestry}'. Set ancestry to a super-population code "
            f"({', '.join(SUPER_POPULATIONS)}) so the VCF header ##ancestry is populated."
        )
    sex = str(record.get("sex", "") or "").strip() or "NA"
    phenotype = str(record.get("phenotype", "") or "").strip() or str(record.get("study_id", "") or "").strip() or "unknown"
    return {"ancestry": ancestry, "sex": sex, "phenotype": phenotype}


def validate_config_keys(config: dict[str, Any], qc: dict[str, Any]) -> None:
    unknown = sorted(set(config) - KNOWN_CONFIG_KEYS)
    if unknown:
        qc["warnings"].append(f"unknown study-YAML keys ignored: {', '.join(unknown)}")
    transforms = config.get("transforms") or {}
    unknown_transforms = sorted(set(transforms) - KNOWN_TRANSFORMS)
    if unknown_transforms:
        qc["warnings"].append(f"unknown transforms ignored: {', '.join(unknown_transforms)}")


def resolve_config_path(path_value: str, config: dict[str, Any]) -> Path:
    path = Path(os.path.expandvars(path_value))
    if path.is_absolute():
        return path
    config_path = config.get("_config_path")
    return Path(config_path).parent / path if config_path else path


def apply_configured_exclusions(df: pd.DataFrame, config: dict[str, Any], qc: dict[str, Any]) -> pd.DataFrame:
    exclusion = config.get("exclude") or config.get("exclusions")
    if not exclusion:
        return df
    out = df
    for item in exclusion if isinstance(exclusion, list) else [exclusion]:
        if not isinstance(item, dict) or not item.get("file"):
            continue
        path = resolve_config_path(str(item["file"]), config)
        key = item.get("key", "SNP")
        if not path.exists() or key not in out.columns:
            qc["warnings"].append(f"exclusion skipped: {path} key={key}")
            continue
        excluded = read_table(path, delimiter=item.get("delimiter", "\t"))
        if key not in excluded.columns:
            qc["warnings"].append(f"exclusion key not found in {path}: {key}")
            continue
        before = len(out)
        out = out[~out[key].astype(str).isin(set(excluded[key].astype(str)))]
        qc.setdefault("exclusions", []).append({"file": str(path), "key": key, "removed": before - len(out)})
    return out.reset_index(drop=True)


def apply_configured_joins(df: pd.DataFrame, config: dict[str, Any], qc: dict[str, Any]) -> pd.DataFrame:
    out = df
    for item in config.get("joins", []) or []:
        if not isinstance(item, dict) or not item.get("file"):
            continue
        path = resolve_config_path(str(item["file"]), config)
        left_on = item.get("left_on", item.get("on"))
        right_on = item.get("right_on", item.get("on"))
        if not path.exists() or not left_on or left_on not in out.columns:
            qc["warnings"].append(f"join skipped: {path}")
            continue
        join_df = read_table(path, delimiter=item.get("delimiter", "\t"))
        if right_on not in join_df.columns:
            qc["warnings"].append(f"join key not found in {path}: {right_on}")
            continue
        duplicates = int(join_df.duplicated(subset=[right_on]).sum())
        if duplicates:
            join_df = join_df.drop_duplicates(subset=[right_on], keep="first")
            qc["warnings"].append(f"join file {path} had {duplicates} duplicate {right_on} keys; kept the first of each")
        before_cols, before_rows = set(out.columns), int(len(out))
        out = pd.merge(out, join_df, how=item.get("how", "left"), left_on=left_on, right_on=right_on, suffixes=("", "_join"))
        qc.setdefault("joins", []).append({"file": str(path), "left_on": left_on, "right_on": right_on,
                                           "added_columns": sorted(set(out.columns) - before_cols),
                                           "rows_before": before_rows, "rows_after": int(len(out))})
    return out.reset_index(drop=True)


def recursive_metadata_value(payload: Any, keys: set[str]) -> Any:
    if isinstance(payload, dict):
        for key, value in payload.items():
            if str(key).lower() in keys and value not in (None, ""):
                return value
        payload = list(payload.values())
    if isinstance(payload, list):
        for item in payload:
            found = recursive_metadata_value(item, keys)
            if found not in (None, ""):
                return found
    return ""


def declared_build(record: dict[str, Any], config: dict[str, Any]) -> tuple[str, str]:
    """A build stated outside the data: manifest or study YAML, then source metadata, then a VCF header's assembly."""
    configured = normalize_build_value(config.get("input_build") or record.get("input_build"))
    if configured:
        return configured, "manifest_or_yaml"
    metadata = normalize_build_value(recursive_metadata_value(record.get("metadata", {}),
                                                              {"build", "genome_build", "assembly", "coordinate_build", "genome_assembly"}))
    if metadata:
        return metadata, "source_metadata"
    header = normalize_build_value(config.get("_vcf_header_build"))
    return (header, "vcf_header") if header else ("", "")


def build_from_vcf_header(header: list[str]) -> str:
    """hg19 or hg38 from ##contig assembly= or ##reference, else ''."""
    text = " ".join(line for line in header if line.startswith(("##contig", "##reference", "##assembly"))).lower()
    for token, build in (("grch38", "hg38"), ("hg38", "hg38"), ("b38", "hg38"), ("grch37", "hg19"), ("hg19", "hg19"),
                         ("b37", "hg19"), ("human_g1k_v37", "hg19"), ("hs37d5", "hg19")):
        if token in text:
            return build
    return ""

# ---------------------------------------------------------------------------------------------------
# Table path: column selection and canonical frame
# ---------------------------------------------------------------------------------------------------

def column_has_values(series: pd.Series) -> bool:
    """True when the column holds any non-missing value (checks a sample first, then the whole column)."""
    present = series.dropna()
    if present.empty:
        return False
    if pd.api.types.is_numeric_dtype(present):
        return True

    def any_value(values: pd.Series) -> bool:
        text = values.astype(str).str.strip().str.lower()
        return bool((~text.isin(MISSING_TOKENS)).any())

    return any_value(present.head(10_000)) or any_value(present)


def _in_family(column: str, family: str | None) -> bool:
    if family is None:
        return True
    harmonised = str(column).lower().startswith("hm_")
    return harmonised if family == "hm" else not harmonised


def _pick(df: pd.DataFrame, aliases: list[str], family: str | None) -> str | None:
    lower = {str(c).lower(): c for c in df.columns}
    candidates = [a for a in aliases if a in df.columns] + [lower[a.lower()] for a in aliases if a not in df.columns and a.lower() in lower]
    candidates = [c for c in dict.fromkeys(candidates) if _in_family(c, family)]
    with_values = [c for c in candidates if column_has_values(df[c])]
    return (with_values or candidates or [None])[0]


def select_columns(df: pd.DataFrame, config: dict[str, Any]) -> dict[str, str | None]:
    configured = {k: v for k, v in (config.get("columns") or {}).items() if v in df.columns}

    def available(field: str, family: str | None) -> bool:
        if field in configured:
            return True
        column = _pick(df, FIELD_ALIASES[field], family)
        return column is not None and column_has_values(df[column])

    def complete(family: str) -> bool:
        return available("EA", family) and (available("OA", family) or available("REF", None)) and \
            (available("BETA", family) or available("OR", family))

    family = "orig" if complete("orig") else ("hm" if complete("hm") else None)
    chosen: dict[str, str | None] = {}
    for field, aliases in FIELD_ALIASES.items():
        if field in configured:
            chosen[field] = configured[field]
            continue
        fam = family if field in FAMILY_FIELDS else None
        column = _pick(df, aliases, fam)
        if column is None and field in ("CHR", "POS"):
            column = _pick(df, aliases, None)
        chosen[field] = column
    if chosen.get("SNPID") == chosen.get("SNP") and chosen.get("SNP") is not None and "SNPID" not in configured:
        chosen["SNPID"] = chosen["SNP"]
    return chosen


def normalise_sample_size(series: pd.Series, qc: dict[str, Any]) -> pd.Series:
    values = numeric_series(series).astype("float64")
    present = values.notna()
    fractional = present & ((values - np.floor(values)) != 0)
    if fractional.any():
        qc["sample_size_values_rounded"] = int(fractional.sum())
    values.loc[present] = np.floor(values.loc[present] + 0.5)
    return values


def normalize_variant_identifiers(out: pd.DataFrame, qc: dict[str, Any]) -> pd.DataFrame:
    """Keep compound IDs in SNPID, extract rsIDs into SNP, and blank SNP values that are not rsIDs."""
    out = out.copy()
    snp = out["SNP"].astype("string").str.strip()
    snpid = out["SNPID"].astype("string").str.strip()
    empty_snp = snp.isna() | snp.str.lower().isin(MISSING_TOKENS)
    empty_snpid = snpid.isna() | snpid.str.lower().isin(MISSING_TOKENS)
    original_snp = snp.mask(empty_snp, pd.NA)
    preserve = (empty_snpid & original_snp.notna()).to_numpy()
    snpid = snpid.mask(empty_snpid, pd.NA)
    snpid[preserve] = original_snp[preserve]
    out["SNPID"] = snpid.astype(object)
    rs_from_snp = original_snp.str.extract(r"\b(rs\d+)\b", flags=re.I, expand=False).str.lower()
    rs_from_snpid = snpid.str.extract(r"\b(rs\d+)\b", flags=re.I, expand=False).str.lower()
    extracted = rs_from_snp.fillna(rs_from_snpid)
    qc["variant_id_normalisation"] = {
        "rsid_extracted": int(extracted.notna().sum()),
        "compound_snpid_preserved": int(preserve.sum()),
        "compound_snp_values_changed": int((extracted.notna() & (snp.str.lower() != extracted)).sum()),
    }
    out["SNP"] = extracted.astype(object).where(extracted.notna(), None)
    return out


def build_canonical(df: pd.DataFrame, cols: dict[str, str | None], config: dict[str, Any], qc: dict[str, Any],
                    negate_log10p: bool | None = None, derive_oa: bool | None = None) -> pd.DataFrame:
    """Canonical columns from the source columns. `negate_log10p` and `derive_oa` fix decisions that otherwise
    depend on the whole column (pieces of one table share the first piece's decisions); None decides here."""
    transforms = config.get("transforms", {}) or {}

    def take(field: str) -> pd.Series:
        column = cols.get(field)
        return df[column] if column else pd.Series(np.nan, index=df.index)

    out = pd.DataFrame(index=df.index)
    out["CHR"] = normalize_chromosome_series(take("CHR"))
    for field in ("POS", "POS19", "POS38", "EAF", "BETA", "SE", "P", "LOG10P", "OR", "CI_LOWER", "CI_UPPER", "NCASE", "EZ", "IMPINFO"):
        out[field] = numeric_series(take(field)).astype("float64")
    out["N"] = normalise_sample_size(take("N"), qc)
    out["SNP"] = take("SNP").astype("string").astype(object)
    out["SNPID"] = take("SNPID").astype("string").astype(object)
    uppercase = transforms.get("uppercase_alleles", True)
    for field in ("EA", "OA"):
        values = take(field).astype("string").str.strip()
        out[field] = (values.str.upper() if uppercase else values).astype(object)
    out = normalize_variant_identifiers(out, qc)

    if derive_oa is None:
        derive_oa = bool(transforms.get("derive_oa_from_ref_alt", False) or (upper_alleles(out["OA"]) == "").all())
    if derive_oa and cols.get("REF") and cols.get("ALT"):
        ea, ref, alt = upper_alleles(out["EA"]), upper_alleles(df[cols["REF"]]), upper_alleles(df[cols["ALT"]])
        derived = np.where((ea != "") & (ea == ref), alt, np.where((ea != "") & (ea == alt), ref, ""))
        use = derived != ""
        out.loc[use, "OA"] = derived[use]
        qc["derived_oa_from_ref_alt"] = True

    derivations = {"beta_from_or": 0, "se_from_ci": 0, "se_from_p": 0}
    from_or = ~np.isfinite(out["BETA"]) & np.isfinite(out["OR"]) & (out["OR"] > 0)
    out.loc[from_or, "BETA"] = np.log(out.loc[from_or, "OR"])
    derivations["beta_from_or"] = int(from_or.sum())
    lower, upper = out["CI_LOWER"], out["CI_UPPER"]
    need_se = ~np.isfinite(out["SE"]) & np.isfinite(out["BETA"])
    ci_or = need_se & from_or & np.isfinite(lower) & np.isfinite(upper) & (lower > 0) & (upper > lower)
    out.loc[ci_or, "SE"] = (np.log(upper[ci_or]) - np.log(lower[ci_or])) / (2 * Z_975)
    ci_beta = need_se & ~from_or & np.isfinite(lower) & np.isfinite(upper) & (upper > lower)
    out.loc[ci_beta, "SE"] = (upper[ci_beta] - lower[ci_beta]) / (2 * Z_975)
    derivations["se_from_ci"] = int(ci_or.sum() + ci_beta.sum())
    p = out["P"]
    se_from_p = ~np.isfinite(out["SE"]) & np.isfinite(out["BETA"]) & np.isfinite(p) & (p > 0) & (p < 1)
    if se_from_p.any():
        z = -_inv_cdf(p[se_from_p].to_numpy() / 2.0)
        out.loc[se_from_p, "SE"] = np.abs(out.loc[se_from_p, "BETA"].to_numpy()) / z
    derivations["se_from_p"] = int(se_from_p.sum())
    qc["derivations"] = derivations

    lp, lp_info = gwas_stats.compute_lp(p=out["P"], log10p=out["LOG10P"], beta=out["BETA"], se=out["SE"],
                                        negate_log10p=negate_log10p)
    out["LP"] = lp
    qc["p_value_source"] = lp_info
    out = out.drop(columns=["P", "LOG10P", "OR", "CI_LOWER", "CI_UPPER"])
    out["EAF_SOURCE"] = np.where(np.isfinite(out["EAF"]), "input", "").astype(object)

    for extra in config.get("extras", []) or []:
        if extra in df.columns and extra not in out.columns:
            out[extra] = df[extra]
    return out


# ---------------------------------------------------------------------------------------------------
# Table path: coordinates
# ---------------------------------------------------------------------------------------------------

def fill_positions_from_ids(out: pd.DataFrame, qc: dict[str, Any]) -> pd.DataFrame:
    """Use chr:pos[:a1:a2] identifiers for rows that have no chromosome/position columns."""
    no_position = (out["CHR"].isna() | ~np.isfinite(out["POS"])).to_numpy() & ~np.isfinite(out["POS19"]).to_numpy() \
        & ~np.isfinite(out["POS38"]).to_numpy()
    if not no_position.any():
        return out
    parsed = coords.parse_variant_ids(out.loc[no_position, "SNPID"])
    usable = parsed["CHR"].notna() & np.isfinite(parsed["POS"])
    rows = parsed.index[usable]
    out.loc[rows, "CHR"] = parsed.loc[rows, "CHR"]
    out.loc[rows, "POS"] = parsed.loc[rows, "POS"]
    qc["positions_from_variant_ids"] = int(len(rows))
    return out


def detect_build_paired(rows: pd.DataFrame, paired: dict[str, Any]) -> dict[str, Any]:
    """Count a fixed-seed sample of rows matching the same common SNPs placed in GRCh37 and in GRCh38; a build is
    called when its hits are at least 10 times the other's (the minimum count is applied by the caller)."""
    sample = rows if len(rows) <= coords.DETECTION_SAMPLE else rows.sample(n=coords.DETECTION_SAMPLE, random_state=coords.DETECTION_SEED)
    hits = {build: grch38.build_hits(paired[build], sample["CHR"], sample["POS"], sample["EA"], sample["OA"]) for build in ("hg19", "hg38")}
    result = {"build": "", "hits_hg19": hits["hg19"], "hits_hg38": hits["hg38"], "sampled": int(len(sample)),
              "method_detail": "paired GRCh37/GRCh38 sample of common 1000 Genomes SNPs"}
    if hits["hg19"] > 0 and hits["hg19"] >= 10 * hits["hg38"]:
        result["build"] = "hg19"
    elif hits["hg38"] > 0 and hits["hg38"] >= 10 * hits["hg19"]:
        result["build"] = "hg38"
    return result


def resolve_input_build(out: pd.DataFrame, record: dict[str, Any], config: dict[str, Any], lookup: pd.DataFrame,
                        reference_root: str, qc: dict[str, Any]) -> str:
    """The input build, detected from the data in the first instance: a fixed-seed sample of positions and alleles
    is matched to the reference under each build (gwas_coordinates.detect_build). A declared build (manifest, YAML,
    source metadata, VCF header) is used only when detection makes no call; a declared build that disagrees with a
    detected one is reported as a warning, not used."""
    declared, declared_source = declared_build(record, config)
    rows = out[np.isfinite(out["POS"]) & ~np.isfinite(out["POS19"]) & ~np.isfinite(out["POS38"])]
    detection: dict[str, Any] = {"build": "", "reason": "reference lookup without alleles"}
    paired = grch38.paired_build_sample(reference_root)
    if paired is not None:
        detection = detect_build_paired(rows, paired)
    elif alleles.lookup_has_alleles(lookup):
        lifter, _ = coords.load_liftover(reference_root, "hg38", "hg19")
        detection = coords.detect_build(rows["CHR"], rows["POS"], rows["EA"], rows["OA"], lookup, lifter)
    qc["build_detection"] = {**detection, "declared": declared, "declared_source": declared_source}
    names = {"hg19": "GRCh37", "hg38": "GRCh38"}
    detected = detection.get("build", "")
    if detected and int(detection.get(f"hits_{detected}", 0)) < BUILD_MIN_HITS:
        qc["build_detection"]["below_minimum_hits"] = BUILD_MIN_HITS
        detected = ""
    if detected:
        if declared and declared != detected:
            qc["warnings"].append(f"declared build {names.get(declared, declared)} ({declared_source}) disagrees with the "
                                  f"build detected from the data, {names[detected]}; the detected build was used")
        qc["build_detection"].update({"build": detected, "method": "detected"})
        return detected
    if not declared:
        raise config_error(
            f"Row {record.get('study_id')}/{record.get('output_id')}: could not detect the genome build from the data "
            f"(GRCh37 hits {detection.get('hits_hg19', 0)}, GRCh38 hits {detection.get('hits_hg38', 0)} of "
            f"{detection.get('sampled', 0)} sampled positions) and none is declared. Set input_build (hg19 or hg38) in the "
            "manifest or study YAML."
        )
    qc["build_detection"].update({"build": declared, "method": "declared"})
    return declared


def fill_coordinates_from_dbsnp(out: pd.DataFrame, dbsnp: Any, qc: dict[str, Any]) -> pd.DataFrame:
    """GRCh38 chromosome and position straight from dbSNP for rows with an rsID and no position (no liftover)."""
    fills = qc.setdefault("reference_fills", {})
    fills.setdefault("POS38", 0)
    snp = out["SNP"].astype("string")
    need = (~np.isfinite(out["POS19"]) & ~np.isfinite(out["POS38"])).to_numpy() & snp.str.match(r"^rs\d+$").fillna(False).to_numpy(dtype=bool)
    qc["coordinate_source"] = qc.get("coordinate_source") or ("dbsnp" if need.any() else "input")
    if not need.any():
        return out
    numbers = np.asarray(coords.rs_number(out.loc[need, "SNP"]), dtype="int64")
    code, pos = dbsnp.positions(numbers)
    fills.setdefault("POS38_merged_rsid", 0)
    via_merge = np.zeros(len(numbers), dtype=bool)
    if dbsnp.has_merges and (code == 0).any():
        # rsIDs retired by a dbSNP merge (rs12016871 into rs9581854): position of the current rsID. SNPID keeps the
        # study's rsID; the ID column is matched on position and alleles later, so it carries the current one.
        absent = np.where(code == 0)[0]
        current = dbsnp.current_rsids(numbers[absent])
        merged = absent[current > 0]
        if len(merged):
            code_m, pos_m = dbsnp.positions(current[current > 0])
            code[merged], pos[merged] = code_m, pos_m
            via_merge[merged] = code_m > 0
    elif (code == 0).any():
        qc["warnings"].append("dbSNP merge map not built (src/reference_data/build_dbsnp_merges.py): rsIDs that dbSNP has merged "
                              "into another rsID were not given positions")
    found_chr = np.array([grch38.CHROMS[c - 1] if c > 0 else "" for c in code], dtype=object)
    current = out.loc[need, "CHR"]
    chr_ok = current.isna().to_numpy() | (current.astype(str).to_numpy() == found_chr)
    use = (code > 0) & chr_ok
    rows = out.index[need][use]
    out.loc[rows, "CHR"] = found_chr[use]
    out.loc[rows, "POS38"] = pos[use].astype("float64")
    fills["POS38"] += int(use.sum())
    fills["POS38_merged_rsid"] += int((use & via_merge).sum())
    if (use & via_merge).any():
        out["RS_MERGED"] = out.get("RS_MERGED", False)
        out.loc[out.index[need][use & via_merge], "RS_MERGED"] = True
    qc["rsid_chromosome_conflicts"] = qc.get("rsid_chromosome_conflicts", 0) + int(((code > 0) & ~chr_ok).sum())
    return out


def fill_coordinates_by_rsid(out: pd.DataFrame, reference_root: str, lookup: pd.DataFrame, qc: dict[str, Any]) -> pd.DataFrame:
    fills = qc.setdefault("reference_fills", {})
    fills.setdefault("CHR", 0)
    fills.setdefault("POS19", 0)
    snp = out["SNP"].astype("string")
    need = (~np.isfinite(out["POS19"]) & ~np.isfinite(out["POS38"])).to_numpy() & snp.str.match(r"^rs\d+$").fillna(False).to_numpy(dtype=bool)
    if not need.any():
        qc.setdefault("coordinate_source", "input")
        return out
    index = coords.CoordinateIndex.load(coords.coordinate_index_dir(reference_root))
    if index is not None:
        found = index.lookup(out.loc[need, "SNP"])
        qc["coordinate_source"] = "coordinate_index"
    else:
        ref = alleles.prepare_lookup(lookup)
        ref = ref[ref["RS"].to_numpy() >= 0].drop_duplicates("RS").set_index("RS")
        numbers = pd.Series(coords.rs_number(out.loc[need, "SNP"]), index=out.index[need])
        found = pd.DataFrame({"CHR": numbers.map(ref["CHR"].astype(object)), "POS19": numbers.map(ref["POS19"])}, index=out.index[need])
        qc["coordinate_source"] = "population_lookup"
        qc["warnings"].append("rsID coordinate index not found; positions for rows without coordinates come from the population "
                              "lookup (common biallelic SNPs only). Build it with src/reference_data/make_1000genomes_lookup_tables.py")
    found_pos = numeric_series(found["POS19"]).to_numpy(dtype="float64")
    current_chr = out.loc[need, "CHR"]
    chr_ok = current_chr.isna().to_numpy() | (current_chr.astype(str).to_numpy() == found["CHR"].astype(str).to_numpy())
    use = np.isfinite(found_pos) & chr_ok
    rows = out.index[need][use]
    fills["CHR"] += int(current_chr.isna().to_numpy()[use].sum())
    out.loc[rows, "CHR"] = found.loc[rows, "CHR"].to_numpy()
    out.loc[rows, "POS19"] = found_pos[use]
    fills["POS19"] += int(use.sum())
    # Counted here and reported once per table (count_warnings), so pieces of a table add up.
    qc["rsid_chromosome_conflicts"] = qc.get("rsid_chromosome_conflicts", 0) + int((np.isfinite(found_pos) & ~chr_ok).sum())
    return out


def _snv_mask(ea: pd.Series, oa: pd.Series) -> np.ndarray:
    ea, oa = upper_alleles(ea), upper_alleles(oa)
    return (ea.isin(list("ACGT")) & oa.isin(list("ACGT"))).to_numpy()


_BCFTOOLS_LIFTOVER: dict[str, Any] = {}


def bcftools_liftover() -> str:
    """Path of a bcftools with the +liftover plugin (GWAS_BCFTOOLS, else bcftools on PATH; BCFTOOLS_PLUGINS must point at
    the plugin), or '' when unavailable."""
    if "path" not in _BCFTOOLS_LIFTOVER:
        candidate = os.environ.get("GWAS_BCFTOOLS", "") or shutil.which("bcftools") or ""
        ok = False
        if candidate:
            try:
                listed = subprocess.run([candidate, "plugin", "-l"], capture_output=True, text=True, timeout=60).stdout
                ok = any(line.strip() == "liftover" for line in listed.splitlines())
            except (OSError, subprocess.SubprocessError):
                ok = False
        _BCFTOOLS_LIFTOVER["path"] = candidate if ok else ""
    return _BCFTOOLS_LIFTOVER["path"]


def lift_with_bcftools(out: pd.DataFrame, rows: pd.Index, reference_root: str, report: dict[str, Any]) -> pd.DataFrame:
    """Lift GRCh37 rows to GRCh38 with BCFtools/liftover (Genovese et al. 2024): alleles are first oriented to the GRCh37
    reference (REF must match the source genome), then lifted with both FASTAs and the UCSC chain; the plugin swaps
    alleles where the reference allele changed between builds and complements them on minus-strand blocks. The effect
    allele is followed through, so ES and EAF need no change here (GRCh38 alignment orients them afterwards)."""
    fasta37 = alleles.Fasta.find(reference_root, "GRCh37")
    fasta38 = alleles.Fasta.find(reference_root)
    chain = Path(reference_root) / "genome" / "hg19ToHg38.over.chain.gz"
    sub = out.loc[rows]
    action = alleles.reference_alignment(sub["CHR"], sub["POS19"], sub["EA"], sub["OA"], fasta37)["action"].to_numpy()
    ea, oa = upper_alleles(sub["EA"]).to_numpy(), upper_alleles(sub["OA"]).to_numpy()
    cea, coa = complement_alleles(sub["EA"]).to_numpy(), complement_alleles(sub["OA"]).to_numpy()
    ref = np.select([action == "keep", action == "swap", action == "complement", action == "complement_swap"], [oa, ea, coa, cea], "")
    alt = np.select([action == "keep", action == "swap", action == "complement", action == "complement_swap"], [ea, oa, cea, coa], "")
    ea_is_alt = np.isin(action, ["keep", "complement"])
    usable = ref != ""
    valid = pd.Series(ea).str.fullmatch(r"[ACGT]+").fillna(False).to_numpy(bool) & pd.Series(oa).str.fullmatch(r"[ACGT]+").fillna(False).to_numpy(bool) & (ea != oa)
    report["grch37_allele_mismatch"] = int((~usable & valid).sum())  # rows with invalid alleles are counted as such instead
    out.loc[rows, "POS38"] = np.nan
    work = make_work_dir()
    try:
        vcf = work / "grch37.vcf"
        contigs = [f"##contig=<ID={name},length={length}>" for name, (length, *_rest) in fasta37.index.items()]
        lines = ["##fileformat=VCFv4.2", *contigs, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
        chrom = normalize_chromosome_series(sub["CHR"]).astype(str).to_numpy()
        pos = sub["POS19"].to_numpy(dtype="float64")
        idx = np.where(usable & np.isfinite(pos))[0]
        order = idx[np.lexsort((pos[idx], gwas_vcf.chrom_ranks(pd.Series(chrom[idx]))))]
        body = [f"{chrom[i]}\t{int(pos[i])}\t{i}\t{ref[i]}\t{alt[i]}\t.\t.\t." for i in order]
        vcf.write_text("\n".join(lines + body) + "\n", encoding="utf-8")
        command = [bcftools_liftover(), "+liftover", "--no-version", str(vcf), "--", "-s", str(fasta37.path), "-f",
                   str(fasta38.path), "-c", str(chain), "--reject", str(work / "reject.vcf"), "--write-reject"]
        lifted = subprocess.run(command, capture_output=True, text=True, check=True).stdout
        new_chr, new_pos = np.full(len(sub), None, dtype=object), np.full(len(sub), np.nan)
        new_ea, new_oa = ea.copy(), oa.copy()
        swapped = flipped = added = 0
        for line in lifted.splitlines():
            if line.startswith("#"):
                continue
            f = line.split("\t")
            i, info = int(f[2]), f[7]
            if "," in f[4] or "SWAP=-1" in info:  # a new reference allele was added: not a simple biallelic variant any more
                added += 1
                continue
            swap = "SWAP=1" in info
            swapped += swap
            flipped += "FLIP" in info.split(";")
            # the effect allele was ALT before liftover iff ea_is_alt; a swap exchanges REF and ALT
            effect_is_alt = bool(ea_is_alt[i]) != swap
            effect, other = (f[4], f[3]) if effect_is_alt else (f[3], f[4])
            new_chr[i], new_pos[i], new_ea[i], new_oa[i] = f[0], float(f[1]), effect, other
        rejected: dict[str, int] = {}
        reject = work / "reject.vcf"
        if reject.exists():
            for line in reject.read_text(encoding="utf-8").splitlines():
                if line.startswith("#"):
                    continue
                fields = line.split("\t")  # --write-reject puts the reason in FILTER (e.g. UnmappedAnchors)
                reason = fields[6] if len(fields) > 6 and fields[6] not in ("", ".", "PASS") else "rejected"
                rejected[reason] = rejected.get(reason, 0) + 1
    finally:
        shutil.rmtree(work, ignore_errors=True)
    lifted_ok = np.isfinite(new_pos)
    same_chrom = lifted_ok & (normalize_chromosome_series(pd.Series(new_chr)).astype(str).to_numpy() == chrom)
    out.loc[rows[same_chrom], "POS38"] = new_pos[same_chrom]
    out.loc[rows[same_chrom], "EA"] = new_ea[same_chrom]
    out.loc[rows[same_chrom], "OA"] = new_oa[same_chrom]
    report.update({"tool": "bcftools +liftover", "ok": int(same_chrom.sum()), "other_chrom": int((lifted_ok & ~same_chrom).sum()),
                   "reference_allele_added": added, "swapped_between_builds": swapped, "minus_strand": flipped,
                   "rejected": rejected})
    return out


def establish_dual_build(out: pd.DataFrame, reference_root: str, qc: dict[str, Any], need_pos19: bool = True) -> pd.DataFrame:
    """Give every row a lifted GRCh38 POS38 (dropping failures) and, where needed (the GRCh37 lookup), a GRCh37 POS19."""
    out = out.copy()
    out["LIFT_STRAND"] = "+"
    report = qc.setdefault("liftover", {})
    need38 = (~np.isfinite(out["POS38"]) & np.isfinite(out["POS19"])).to_numpy()
    if need38.any() and bcftools_liftover() and alleles.Fasta.find(reference_root, "GRCh37") is not None \
            and alleles.Fasta.find(reference_root) is not None and (Path(reference_root) / "genome" / "hg19ToHg38.over.chain.gz").exists():
        report["hg19_to_hg38"] = {}
        out = lift_with_bcftools(out, out.index[need38], reference_root, report["hg19_to_hg38"])
    elif need38.any():
        lifter, reason = coords.load_liftover(reference_root, "hg19", "hg38")
        if lifter is None:
            raise config_error(
                f"{int(need38.sum())} variants have only a GRCh37 position and need hg19->hg38 liftover, but {reason}. "
                "Add hg19ToHg38.over.chain.gz (see src/reference_data/download_fasta_and_chains.sh)."
            )
        result = coords.lift(out.loc[need38, "CHR"], out.loc[need38, "POS19"], lifter)
        rows = out.index[need38]
        minus = result.ok & (result.strand == "-")
        snv = _snv_mask(out.loc[rows, "EA"], out.loc[rows, "OA"])
        positions = result.pos.copy()
        positions[minus & ~snv] = np.nan
        out.loc[rows, "POS38"] = positions
        flip_rows = rows[minus & snv]
        out.loc[flip_rows, "EA"] = complement_alleles(out.loc[flip_rows, "EA"]).to_numpy()
        out.loc[flip_rows, "OA"] = complement_alleles(out.loc[flip_rows, "OA"]).to_numpy()
        out.loc[rows[minus], "LIFT_STRAND"] = "-"
        report["hg19_to_hg38"] = {**result.counts(), "minus_strand": int(minus.sum()), "minus_strand_non_snv_dropped": int((minus & ~snv).sum())}
    need19 = (~np.isfinite(out["POS19"]) & np.isfinite(out["POS38"])).to_numpy()
    if need_pos19 and need19.any():
        lifter, reason = coords.load_liftover(reference_root, "hg38", "hg19")
        if lifter is None:
            qc["pos19_missing_reason"] = reason
            qc["pos19_missing_variants"] = qc.get("pos19_missing_variants", 0) + int(need19.sum())
        else:
            result = coords.lift(out.loc[need19, "CHR"], out.loc[need19, "POS38"], lifter)
            rows = out.index[need19]
            out.loc[rows, "POS19"] = result.pos
            out.loc[rows[result.ok & (result.strand == "-")], "LIFT_STRAND"] = "-"
            report["hg38_to_hg19"] = {**result.counts(), "minus_strand": int((result.ok & (result.strand == "-")).sum())}
    return out


# Why a row without its required fields was dropped: the first cause that applies, in this order, so the causes
# add up to the rows removed. Coordinate causes come first because a row without a position cannot be used whatever
# its statistics.
INCOMPLETE_CAUSES = ("rsid_not_in_coordinate_index", "no_position_or_rsid", "missing_allele", "invalid_allele",
                     "liftover_to_grch38_failed", "chromosome_not_primary", "invalid_position", "missing_beta",
                     "missing_or_invalid_se", "missing_or_invalid_p")


def incomplete_causes(out: pd.DataFrame, checks: dict[str, np.ndarray]) -> np.ndarray:
    """The first cause (INCOMPLETE_CAUSES) for each row failing a check in `checks`; '' for complete rows."""
    pos19, pos38 = np.isfinite(out["POS19"]).to_numpy(), np.isfinite(out["POS38"]).to_numpy()
    rsid = out["SNP"].astype("string").str.match(r"^rs\d+$").fillna(False).to_numpy(dtype=bool)
    no_position = ~pos19 & ~pos38
    conditions = [
        no_position & rsid,
        no_position & ~rsid,
        ~checks["EA"] | ~checks["OA"],
        ~checks["ALLELES"],  # before liftover: I/D codes or '-' cannot be lifted, and that is not a liftover failure
        pos19 & ~pos38,
        ~checks["CHR"],
        ~checks["POS38"],
        ~checks["BETA"],
        ~checks["SE"],
        ~checks["LP"],
    ]
    return np.select(conditions, list(INCOMPLETE_CAUSES), default="").astype(object)


def drop_incomplete_rows(out: pd.DataFrame, qc: dict[str, Any], maf_threshold: float = 0.01) -> pd.DataFrame:
    ea, oa = upper_alleles(out["EA"]), upper_alleles(out["OA"])
    checks = {
        "CHR": out["CHR"].astype(str).isin(PRIMARY_CHROMS).to_numpy(),
        "POS38": (np.isfinite(out["POS38"]) & (out["POS38"] > 0)).to_numpy(),
        "EA": (upper_alleles(out["EA"]) != "").to_numpy(),
        "OA": (upper_alleles(out["OA"]) != "").to_numpy(),
        # VCF alleles must be bases: I/D codes, "-" and identical alleles cannot be written as REF/ALT
        "ALLELES": (ea.str.fullmatch(r"[ACGT]+").fillna(False) & oa.str.fullmatch(r"[ACGT]+").fillna(False) & (ea != oa)).to_numpy(dtype=bool),
        "BETA": np.isfinite(out["BETA"]).to_numpy(),
        "SE": (np.isfinite(out["SE"]) & (out["SE"] > 0)).to_numpy(),
        "LP": (np.isfinite(out["LP"]) & (out["LP"] >= 0)).to_numpy(),
    }
    keep = np.logical_and.reduce(list(checks.values())) if len(out) else np.zeros(0, dtype=bool)
    qc["missing_required_by_field"] = {field: int((~mask).sum()) for field, mask in checks.items()}
    qc["rows_removed_missing_required"] = int((~keep).sum())
    if len(out) and not keep.all():
        cause = incomplete_causes(out, checks)
        eaf = numeric_series(out["EAF"]).to_numpy(dtype="float64")
        rare = np.isfinite(eaf) & (np.minimum(eaf, 1.0 - eaf) < maf_threshold)
        no_eaf = ~np.isfinite(eaf)
        detail = qc.setdefault("incomplete_by_cause", {})
        for name in INCOMPLETE_CAUSES:
            hit = cause == name
            if hit.any():
                detail[name] = {"rows": int(hit.sum()), "study_maf_below_threshold": int((hit & rare).sum()),
                                "no_study_eaf": int((hit & no_eaf).sum())}
    return out[keep]


def annotate_with_reference(out: pd.DataFrame, lookup: pd.DataFrame, qc: dict[str, Any], freqs: Any = None) -> tuple[pd.DataFrame, dict[str, float]]:
    """Reference annotation and EAF fills (GRCh38 1000 Genomes when installed, else the GRCh37 lookup, which also
    fills rsIDs); returns the frame and its allele-frequency concordance moments."""
    if freqs is not None:
        out, rsid_fills = alleles.annotate_with_frequencies(out, freqs), 0
    else:
        out = alleles.annotate_with_reference(out, lookup)
        out, rsid_fills = alleles.fill_missing_rsids(out)
    out, eaf_fills = alleles.fill_eaf_from_reference(out)
    fills = qc.setdefault("reference_fills", {})
    fills["rsid"], fills["EAF"] = rsid_fills, eaf_fills
    return out, alleles.af_moments(out)


def check_af_concordance(moments: dict[str, float], qc: dict[str, Any]) -> None:
    """Fail when study and reference frequencies anti-correlate (swapped allele columns); warn when they barely agree."""
    concordance = alleles.af_concordance_from_moments(moments)
    qc["af_concordance"] = concordance
    r = concordance.get("pearson_r")
    if r is not None and concordance["n"] >= AF_CONCORDANCE_MIN_N and r < 0:
        raise data_error(
            f"Effect-allele frequencies anti-correlate with the 1000 Genomes reference (r={r:.2f} over {concordance['n']} SNPs). "
            "The effect/other allele columns are probably swapped (e.g. SAIGE Allele2 is the effect allele); "
            "set columns: EA/OA in the study YAML."
        )
    if r is not None and concordance["n"] >= 100 and r < 0.5:
        qc["warnings"].append(f"weak agreement between study EAF and the reference (r={r:.2f}); check allele columns and ancestry")


def check_output_rows(written: int, record: dict[str, Any], qc: dict[str, Any]) -> None:
    if written == 0:
        missing = {k: v for k, v in qc.get("missing_required_by_field", {}).items() if v}
        raise data_error(
            f"No variants remain for {record.get('study_id')}/{record.get('output_id')} after standardisation "
            f"(rows missing required fields: {missing or 'none'}; source columns: {qc.get('source_columns')}). "
            "Map the missing columns (e.g. SE, P) in the study YAML `columns` block."
        )
    input_rows = qc.get("input_rows", 0)
    retained = written / input_rows if input_rows else 1.0
    if retained < RETENTION_WARNING:
        qc["warnings"].append(f"retained {written} of {input_rows} input rows ({retained:.0%}); see the step log (qc.json steps, "
                              "incomplete_by_cause and liftover)")
    qc["output_rows"] = int(written)


def candidate_columns(raw_path: Path, sep: str, config: dict[str, Any]) -> list[str] | None:
    """Columns worth reading: any alias match, configured column, extra, or exclusion key (None = read all)."""
    if config.get("joins"):
        return None  # joined files may supply any column, so keep everything
    header = list(read_table(raw_path, delimiter=sep, nrows=0).columns)
    wanted = {alias.lower() for aliases in FIELD_ALIASES.values() for alias in aliases}
    wanted |= {str(v).lower() for v in (config.get("columns") or {}).values()}
    wanted |= {str(e).lower() for e in (config.get("extras") or [])}
    exclusion = config.get("exclude") or config.get("exclusions") or []
    for item in exclusion if isinstance(exclusion, list) else [exclusion]:
        if isinstance(item, dict):
            wanted.add(str(item.get("key", "SNP")).lower())
    keep = [c for c in header if str(c).lstrip("#").lower() in wanted]
    return keep or None


# Tables are standardised in pieces of TABLE_CHUNK rows so memory stays bounded whatever the table's size,
# and pieces can run in parallel (--workers). Each piece goes through every per-row step and writes sortable
# lines; the steps that need the whole table run once over the merged pieces: the sort by position,
# deduplication (first record in input order, as for a whole table), the alignment summary, lambda GC and
# the allele-frequency concordance check (pooled moments). Decisions that read a whole column (which columns
# to use, LOG10P sign, deriving OA from REF/ALT, the genome build) are taken from the first piece.

TABLE_CHUNK = 1_000_000  # table rows per piece (GWAS_TABLE_CHUNK overrides)
BUILD_HEAD_ROWS = 200_000  # rows read from the start of the file to detect the build, whatever the piece size
ROW_COLUMN = "_gwasprocessing_row"  # input row number, carried through each piece to restore input order
MERGED_RSID_ORDER_OFFSET = 10 ** 12  # larger than any input row number

FIRST_PIECE_ONLY = ("build_detection", "source_columns")


def table_chunk_rows(config: dict[str, Any]) -> int | None:
    """Rows per piece; None reads the table whole, because a right or outer join needs every row at once."""
    for item in config.get("joins", []) or []:
        if isinstance(item, dict) and str(item.get("how", "left")).lower() not in ("left", "inner"):
            return None
    return max(1, int(os.environ.get("GWAS_TABLE_CHUNK", TABLE_CHUNK)))


def read_pieces(raw_path: Path, sep: str, usecols: list[str] | None, chunk_rows: int | None, config: dict[str, Any]):
    """Yield the table in pieces with ROW_COLUMN numbering the input rows; an empty table yields one empty piece."""
    strip = (config.get("transforms", {}) or {}).get("strip_hash_header", True)
    offset = 0

    def numbered(df: pd.DataFrame) -> pd.DataFrame:
        df = df.reset_index(drop=True)
        if strip and len(df.columns):
            df.columns = [str(c).lstrip("#") for c in df.columns]
        df[ROW_COLUMN] = np.arange(offset, offset + len(df), dtype="int64")
        return df

    if chunk_rows is None:
        yield numbered(read_table(raw_path, delimiter=sep, usecols=usecols))
        return
    yielded = False
    with read_table(raw_path, delimiter=sep, usecols=usecols, chunksize=chunk_rows) as reader:
        for df in reader:
            piece = numbered(df)
            offset += len(piece)
            yielded = True
            yield piece
    if not yielded:
        yield numbered(read_table(raw_path, delimiter=sep, usecols=usecols, nrows=0))


def standardise_piece(df: pd.DataFrame, settings: dict[str, Any], lookup: pd.DataFrame, fasta: Any, part: Path) -> dict[str, Any]:
    """Run every per-row step on one piece and write its lines to `part` (see write_piece_lines).

    Returns the piece's QC counts, frequency moments, the decisions later pieces must share, and seconds per step."""
    config, record, reference_root = settings["config"], settings["record"], settings["reference_root"]
    qc: dict[str, Any] = {"warnings": [], "reference_fills": {}, "input_rows": int(len(df))}
    timings: dict[str, float] = {}
    clock = [time.perf_counter()]

    def lap(step: str) -> None:
        now = time.perf_counter()
        timings[step] = timings.get(step, 0.0) + now - clock[0]
        clock[0] = now

    df = apply_configured_exclusions(df, config, qc)
    df = apply_configured_joins(df, config, qc)
    rows = df[ROW_COLUMN].to_numpy(dtype="int64")
    cols = settings.get("cols")
    if cols is None:
        cols = select_columns(df.drop(columns=[ROW_COLUMN]), config)
        qc["source_columns"] = {k: v for k, v in cols.items() if v}
    out = build_canonical(df, cols, config, qc, negate_log10p=settings.get("negate_log10p"), derive_oa=settings.get("derive_oa"))
    extras = [e for e in (config.get("extras", []) or []) if e in out.columns]
    del df
    lap("columns_and_statistics")
    out = fill_positions_from_ids(out, qc)
    build = settings.get("build")
    if np.isfinite(out["POS"]).any() and (~np.isfinite(out["POS19"]) & ~np.isfinite(out["POS38"]) & np.isfinite(out["POS"])).any():
        build = build or resolve_input_build(out, record, config, lookup, reference_root, qc)
        target = "POS19" if build == "hg19" else "POS38"
        assign = np.isfinite(out["POS"]) & ~np.isfinite(out[target])
        out.loc[assign, target] = out.loc[assign, "POS"]
    out = out.drop(columns=["POS"])
    lap("positions_and_build")
    dbsnp, freqs, _ = grch38_references(reference_root, settings["population"])
    if dbsnp is not None:
        out = fill_coordinates_from_dbsnp(out, dbsnp, qc)
    else:
        out = fill_coordinates_by_rsid(out, reference_root, lookup, qc)
    lap("rsid_coordinates")
    out = establish_dual_build(out, reference_root, qc, need_pos19=freqs is None)
    lap("liftover")
    out = drop_incomplete_rows(out, qc)
    out, moments = annotate_with_reference(out, lookup, qc, freqs)
    lap("reference_annotation")
    out = alleles.resolve_palindromic_strand(out, qc, settings["palindromic_tolerance"])
    lap("palindromic_strand")
    if fasta is not None:
        out, action = alleles.orient_to_reference(out, fasta)
        out, moved = alleles.normalise_indels(out, fasta, action)
        qc["indels_left_aligned"] = moved
    else:
        action = np.full(len(out), "unalignable", dtype=object)
    low, high = alleles.dedup_alleles(out)  # unordered allele pair: the same variant reported either way round matches
    out["HM"] = harmonisation_action(out["PAL"].to_numpy(dtype=object), action)
    lap("grch38_alignment")
    # optional FORMAT fields (SS, NC, EZ, SI) are written when the first piece has values for them (a VCF header may
    # declare fields its records never use)
    keys = settings.get("format_keys") or gwas_vcf.format_keys(
        c for c in gwas_vcf.OPTIONAL_FORMAT if cols.get(c) and c in out.columns and np.isfinite(numeric_series(out[c])).any())
    # each piece is written in genome order (then input order), so the merge only interleaves pre-sorted pieces. A row
    # placed through a retired (merged) rsID sorts after every row with a current rsID, so when the file also reports the
    # variant under its current rsID, the duplicate step keeps that row
    row_ids = rows[out.index.to_numpy()]
    if "RS_MERGED" in out.columns:
        row_ids = row_ids + MERGED_RSID_ORDER_OFFSET * out["RS_MERGED"].fillna(False).to_numpy(dtype=bool)
    order = np.lexsort((row_ids, out["POS38"].to_numpy(dtype="float64"), gwas_vcf.chrom_ranks(out["CHR"])))
    write_piece_lines(part, out.iloc[order], row_ids[order], low[order], high[order], action[order], extras, keys,
                      settings["afs_label"])
    lap("write_piece")
    decisions = {"cols": cols, "extras": extras, "build": build, "format_keys": keys,
                 "negate_log10p": bool(qc.get("p_value_source", {}).get("log10p_negated", False)),
                 "derive_oa": bool(qc.get("derived_oa_from_ref_alt", False))}
    return {"qc": qc, "moments": moments, "decisions": decisions, "part": str(part), "timings": timings}


def harmonisation_action(pal: np.ndarray, action: np.ndarray) -> np.ndarray:
    """INFO/HM per row from the palindromic decision (PAL) and the GRCh38 alignment action."""
    base = {"keep": "forward", "swap": "swapped", "complement": "complemented", "complement_swap": "complemented_swapped",
            "indel_unmatched": "unverified", "unalignable": "unverified", "drop": "dropped"}
    words = np.array([base.get(str(a), "unverified") for a in action], dtype=object)
    swapped = np.isin(action, ["swap", "complement_swap"])
    pal_words = np.array(["palindrome_" + str(p) + ("_swapped" if s else "") for p, s in zip(pal, swapped)], dtype=object)
    return np.where(pal != "", pal_words, words)


def write_piece_lines(part: Path, out: pd.DataFrame, rows: np.ndarray, low: np.ndarray, high: np.ndarray,
                      action: np.ndarray, extras: list[str], keys: list[str], afs_label: str) -> None:
    """One line per row: chromosome rank, POS38, input row, dedup allele pair, alignment action, LP, AF-filled flag,
    AF-missing flag, ES, SS, VCF record.

    The leading fields let `sort` restore genome order with input order as the tie-break; LP is kept at full
    precision for lambda GC. Rows that alignment will drop are kept so they still claim their dedup key."""
    with open(part, "w", encoding="utf-8") as handle:
        step = 500_000
        for start in range(0, len(out), step):
            piece = out.iloc[start:start + step]
            sl = slice(start, start + step)
            lines = np.array(gwas_vcf.table_body_lines(piece, extras, keys, afs_label), dtype=object)
            rank = gwas_vcf.chrom_ranks(piece["CHR"]).astype(str).astype(object)
            pos = piece["POS38"].to_numpy(dtype="float64").astype("int64").astype(str).astype(object)
            lp = np.char.mod("%.17g", piece["LP"].to_numpy(dtype="float64")).astype(object)
            eaf = piece["EAF"].to_numpy(dtype="float64")
            filled = ((piece["EAF_SOURCE"].astype(str).to_numpy() == "reference") & np.isfinite(eaf)).astype(int).astype(str).astype(object)
            no_af = (~np.isfinite(eaf)).astype(int).astype(str).astype(object)
            es = np.char.mod("%.17g", piece["BETA"].to_numpy(dtype="float64")).astype(object)
            ss = np.char.mod("%.17g", numeric_series(piece["N"]).to_numpy(dtype="float64") if "N" in piece.columns
                             else np.full(len(piece), np.nan)).astype(object)
            keys_ = (rank + "\t" + pos + "\t" + rows[sl].astype(str).astype(object) + "\t" + low[sl].astype(object)
                     + "\t" + high[sl].astype(object) + "\t" + action[sl].astype(object) + "\t" + lp + "\t" + filled
                     + "\t" + no_af + "\t" + es + "\t" + ss + "\t" + lines)
            handle.write("\n".join(keys_.tolist()) + "\n")


def merge_pieces(parts: list[Path], work: Path, vcf_path: Path, record: dict[str, Any], config: dict[str, Any],
                 decisions: dict[str, Any], fasta: Any, qc: dict[str, Any], timings: dict[str, float] | None = None,
                 reference_root: str = "", population: str = "") -> None:
    """Sort all pieces by position, remove duplicates, drop alignment failures, then write the header (whose
    ##SAMPLE counts are final only now) in front of the body, bgzip and index."""
    timings = {} if timings is None else timings
    started = time.perf_counter()
    body = work / "body.vcf"
    actions: dict[str, int] = {}
    lp_values = array("d")
    duplicates = written = filled_written = 0
    stats = {"af_missing": 0, "effect_sum": 0.0, "effect_abs_sum": 0.0, "effect_n": 0, "max_ss": -math.inf}
    current: tuple[str, str] | None = None
    seen: set[tuple[str, str]] = set()
    command = ["sort", "-m", "-t", "\t", "-k1,1n", "-k2,2n", "-k3,3n", "-T", str(work), *map(str, parts)]
    dbsnp = grch38_references(reference_root, population)[0] if reference_root else None
    batch: list[str] = []
    matched = 0

    def flush(handle) -> None:
        """Write a batch of records, with ID and FORMAT/ID set to the dbSNP rsID matching (CHROM, POS, REF, ALT)."""
        nonlocal matched
        if dbsnp is not None and batch:
            fields = [line.split("\t", 5) for line in batch]
            found = dbsnp.rsids([f[0] for f in fields], [int(f[1]) for f in fields], [f[3] for f in fields], [f[4] for f in fields])
            matched += int((found > 0).sum())
            for i, f in enumerate(fields):
                rsid = f"rs{found[i]}" if found[i] > 0 else "."
                rest = f[5]
                head, sample = rest.rsplit(":", 1)  # FORMAT/ID is the last value of the sample column
                batch[i] = "\t".join((f[0], f[1], rsid, f[3], f[4], f"{head}:{rsid}\n"))
        handle.writelines(batch)
        batch.clear()
    env = {**os.environ, "LC_ALL": "C"}
    with open(body, "w", encoding="utf-8") as out, \
            subprocess.Popen(command, stdout=subprocess.PIPE, text=True, encoding="utf-8", env=env) as proc:
        for line in proc.stdout:
            rank, pos, _row, low, high, action, lp, filled, no_af, es, ss, record_line = line.split("\t", 11)
            if (rank, pos) != current:
                current, seen = (rank, pos), set()
            if (low, high) in seen:
                duplicates += 1
                continue
            seen.add((low, high))
            actions[action] = actions.get(action, 0) + 1
            if action in alleles.REMOVED_ACTIONS:
                continue
            batch.append(record_line)
            if len(batch) >= 500_000:
                flush(out)
            lp_values.append(float(lp))
            written += 1
            filled_written += filled == "1"
            stats["af_missing"] += no_af == "1"
            effect, size = float(es), float(ss)
            if math.isfinite(effect):
                stats["effect_sum"] += effect
                stats["effect_abs_sum"] += abs(effect)
                stats["effect_n"] += 1
            if math.isfinite(size):
                stats["max_ss"] = max(stats["max_ss"], size)
        flush(out)
    if proc.returncode != 0:
        raise RuntimeError(f"sort exited {proc.returncode} while merging {len(parts)} table pieces")
    if dbsnp is not None:
        qc["rsid_dbsnp"] = {"matched": matched, "unmatched": written - matched, "source": dbsnp.label}
    qc["duplicates_removed"] = duplicates
    qc["af_filled_written"] = filled_written
    if fasta is not None:
        qc["reference_alignment"] = alleles.alignment_summary_from_counts(actions, fasta)
    else:
        qc["reference_alignment"] = {"performed": False, "reason": "no GRCh38 FASTA under <reference_root>/genome"}
        qc["warnings"].append("REF/ALT not aligned to the GRCh38 reference: no FASTA found (run src/reference_data/download_fasta_and_chains.sh)")
    check_output_rows(written, record, qc)
    lp_array = np.frombuffer(lp_values, dtype="float64") if written else np.array([])
    qc["lambda_gc"] = gwas_stats.lambda_gc(lp_array)
    qc["qc_metrics"] = {
        "variants": written, "n_p_below_5e-8": int((lp_array > -math.log10(5e-8)).sum()),
        "mean_effect": stats["effect_sum"] / stats["effect_n"] if stats["effect_n"] else None,
        "mean_abs_effect": stats["effect_abs_sum"] / stats["effect_n"] if stats["effect_n"] else None,
        "af_missing_fraction": stats["af_missing"] / written if written else None,
        "max_sample_size": stats["max_ss"] if math.isfinite(stats["max_ss"]) else None,
    }
    header = work / "header.vcf"
    fields = sample_fields(record, config, qc, decisions["format_keys"])
    qc["sample_header"] = fields
    qc["qc_flags"] = qc_flags(qc, fields)
    reference = fasta.path.name if fasta is not None else "GRCh38"
    header.write_text("\n".join(gwas_vcf.table_header(record, decisions["extras"], decisions["format_keys"], fields,
                                                       reference, qc.get("settings_hash", ""))) + "\n", encoding="utf-8")
    timings["sort_merge_deduplicate"] = time.perf_counter() - started
    started = time.perf_counter()
    gwas_vcf.bgzip_and_index([header, body], vcf_path)
    timings["compress_and_index"] = time.perf_counter() - started


# OpenGWAS QC flags (opengwas-reports funcs/flags.R): flags for a reviewer, never a reason to stop.
FLAG_RULES = {
    "af_correlation": "absolute correlation of study and 1000 Genomes allele frequencies below 0.7",
    "inflation_factor": "lambda GC above 1.2",
    "n": "maximum sample size below 10,000",
    "mean_EFFECT_05": "absolute mean effect above 0.5",
    "mean_EFFECT_01": "absolute mean effect above 0.1",
    "n_p_sig": "more than 1,000 variants at P < 5e-8",
}


def qc_flags(qc: dict[str, Any], fields: dict[str, Any]) -> dict[str, Any]:
    """{flag: description} for each OpenGWAS QC flag raised (LDSC and clumping flags are not computed)."""
    metrics = qc.get("qc_metrics", {}) or {}
    r = (qc.get("af_concordance") or {}).get("pearson_r")
    lam = qc.get("lambda_gc")
    max_n = metrics.get("max_sample_size")
    if max_n is None:
        totals = [fields.get("TotalControls"), fields.get("TotalCases")]
        max_n = sum(float(t) for t in totals if t is not None) or None
    mean_effect = metrics.get("mean_effect")
    raised = {
        "af_correlation": r is not None and abs(r) < 0.7,
        "inflation_factor": lam is not None and lam > 1.2,
        "n": max_n is not None and max_n < 10_000,
        "mean_EFFECT_05": mean_effect is not None and abs(mean_effect) > 0.5,
        "mean_EFFECT_01": mean_effect is not None and abs(mean_effect) > 0.1,
        "n_p_sig": int(metrics.get("n_p_below_5e-8", 0)) > 1000,
    }
    return {flag: FLAG_RULES[flag] for flag, on in raised.items() if on}


def _number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0 else None


def study_totals(record: dict[str, Any], config: dict[str, Any], has_cases: bool) -> tuple[str, float | None, float | None]:
    """(StudyType, TotalCases, TotalControls) as gwas2vcf writes them: case-control studies give both totals,
    continuous traits put the total sample size in TotalControls. Study YAML values win over source metadata."""
    metadata = record.get("metadata", {}) or {}

    def find(*keys: str) -> float | None:
        for key in keys:
            if config.get(key) not in (None, ""):
                return _number(config.get(key))
        return _number(recursive_metadata_value(metadata, set(keys)))

    cases, controls = find("ncase", "n_cases"), find("ncontrol", "n_controls")
    total = find("sample_size")
    declared = str(config.get("study_type", "") or "").strip().lower().replace("-", "").replace("_", "")
    case_control = declared == "casecontrol" or (declared != "continuous" and (cases is not None or has_cases))
    if case_control:
        return "CaseControl", cases, controls
    return "Continuous", None, total


def sample_fields(record: dict[str, Any], config: dict[str, Any], qc: dict[str, Any], keys: list[str]) -> dict[str, Any]:
    """The ##SAMPLE values: harmonisation counts that add up to TotalVariants, study type and totals, trait metadata."""
    total, written = int(qc.get("input_rows", 0)), int(qc.get("output_rows", 0))
    alignment = qc.get("reference_alignment", {}) or {}
    not_harmonised = int(alignment.get("dropped_mismatch", 0) or 0) + int(alignment.get("indels_unmatched_dropped", 0) or 0) \
        + int(qc.get("duplicates_removed", 0))
    fields: dict[str, Any] = {
        "TotalVariants": total, "VariantsNotRead": total - written - not_harmonised, "HarmonisedVariants": written,
        "VariantsNotHarmonised": not_harmonised, "SwitchedAlleles": int(alignment.get("swapped", 0) or 0),
        "AFFilled": int(qc.get("af_filled_written", 0)),
    }
    study_type, cases, controls = study_totals(record, config, "NC" in keys)
    if controls is not None:
        fields["TotalControls"] = int(controls) if float(controls).is_integer() else controls
    if cases is not None:
        fields["TotalCases"] = int(cases) if float(cases).is_integer() else cases
    fields["StudyType"] = study_type
    fields.update({"StudyID": record.get("study_id", ""), "Phenotype": record.get("phenotype", ""),
                   "Ancestry": record.get("ancestry", ""), "Sex": record.get("sex", ""), "PMID": record.get("PMID", "")})
    return fields


def _merge_counts(a: Any, b: Any) -> Any:
    """Pool a piece's QC value into the running total: counts add, dicts merge, lists union; anything else keeps the first."""
    if a is None:
        return b
    if isinstance(a, bool) or isinstance(b, bool):
        return a
    if isinstance(a, (int, np.integer)) and isinstance(b, (int, np.integer)):
        return int(a) + int(b)
    if isinstance(a, dict) and isinstance(b, dict):
        return {k: _merge_counts(a.get(k), b.get(k)) for k in dict.fromkeys([*a, *b])}
    if isinstance(a, list) and isinstance(b, list):
        if any(isinstance(x, dict) for x in a + b):  # e.g. one entry per configured exclusion file, in config order
            return [_merge_counts(x, y) for x, y in zip(a, b)] if len(a) == len(b) else a
        return list(dict.fromkeys([*a, *b]))
    return a


def merge_piece_qc(qc: dict[str, Any], results: list[dict[str, Any]], record: dict[str, Any]) -> dict[str, float]:
    """Pool the pieces' QC into `qc`; returns the pooled frequency moments. Fails if pieces disagree on the build."""
    builds = {r["decisions"]["build"] for r in results if r["decisions"].get("build")}
    if len(builds) > 1:
        raise config_error(f"Row {record.get('study_id')}/{record.get('output_id')}: pieces of the table match different "
                           f"genome builds ({', '.join(sorted(builds))}). Set input_build in the manifest or study YAML.")
    moments = {"n": 0, "mean_x": 0.0, "mean_y": 0.0, "m2_x": 0.0, "m2_y": 0.0, "c_xy": 0.0, "far": 0}
    for result in results:
        moments = alleles.combine_af_moments(moments, result["moments"])
        for key, value in result["qc"].items():
            if key in FIRST_PIECE_ONLY:
                qc.setdefault(key, value)
            elif key == "coordinate_source":
                qc[key] = value if qc.get(key) in (None, "input") else qc[key]
            else:
                qc[key] = _merge_counts(qc.get(key), value)
    conflicts = qc.pop("rsid_chromosome_conflicts", 0)
    if conflicts:
        qc["warnings"].append(f"{conflicts} rsIDs map to a different chromosome than the input reports; their coordinates were not filled")
    unlifted = qc.pop("pos19_missing_variants", 0)
    reason = qc.pop("pos19_missing_reason", "")
    if unlifted:
        qc["warnings"].append(f"{reason}; POS19 left missing for {unlifted} variants")
    return moments


_WORKER: dict[str, Any] = {}
_GRCH38: dict[tuple[str, str], tuple[Any, Any, Any]] = {}


def grch38_references(reference_root: str, population: str) -> tuple[Any, Any, Any]:
    """(dbSNP index, 1000 Genomes frequencies, paired build-detection sample), each None when not installed.
    Memory-mapped, opened once per process."""
    key = (str(reference_root), str(population))
    if key not in _GRCH38:
        _GRCH38[key] = (grch38.Dbsnp.find(reference_root), grch38.Frequencies.find(reference_root, population),
                        grch38.paired_build_sample(reference_root))
    return _GRCH38[key]


def _init_worker(reference_root: str, population: str, fasta_path: str) -> None:
    if grch38_references(reference_root, population)[1] is not None:
        _WORKER["lookup"] = pd.DataFrame(columns=alleles.LOOKUP_COLUMNS)
    else:
        _WORKER["lookup"], _ = alleles.read_reference_lookup(reference_root, population)
    _WORKER["fasta"] = alleles.Fasta(fasta_path) if fasta_path else None


def _piece_in_worker(df: pd.DataFrame, settings: dict[str, Any], part: str) -> dict[str, Any]:
    return standardise_piece(df, settings, _WORKER["lookup"], _WORKER["fasta"], Path(part))


def run_pieces(pieces: Any, settings: dict[str, Any], lookup: pd.DataFrame, fasta: Any, work: Path, workers: int) -> list[dict[str, Any]]:
    """Standardise the first piece here (it fixes the shared decisions), then the rest here or in `workers` processes."""
    part = lambda i: work / f"piece-{i:06d}.txt"  # noqa: E731
    pieces = iter(pieces)
    results = [standardise_piece(next(pieces), settings, lookup, fasta, part(0))]
    settings = {**settings, **results[0]["decisions"]}
    if workers <= 1:
        for i, df in enumerate(pieces, start=1):
            results.append(standardise_piece(df, settings, lookup, fasta, part(i)))
            settings["build"] = settings["build"] or results[-1]["decisions"]["build"]
        return results
    context = multiprocessing.get_context("spawn")  # fork is unsafe on macOS; workers load the cached lookup themselves
    initargs = (settings["reference_root"], settings["population"], str(fasta.path) if fasta is not None else "")
    with ProcessPoolExecutor(max_workers=workers, mp_context=context, initializer=_init_worker, initargs=initargs) as pool:
        pending: set[Any] = set()
        for i, df in enumerate(pieces, start=1):
            if len(pending) >= workers:  # at most `workers` pieces in flight, so memory stays bounded
                done, pending = wait(pending, return_when=FIRST_COMPLETED)
                results.extend(future.result() for future in done)
            pending.add(pool.submit(_piece_in_worker, df, settings, str(part(i))))
            del df
        results.extend(future.result() for future in wait(pending).done)
    return results


def detect_build_from_head(head: pd.DataFrame, settings: dict[str, Any], lookup: pd.DataFrame, qc: dict[str, Any]) -> str:
    """The input build, decided once from the first rows of the file (so it does not depend on the piece size)."""
    config, record = settings["config"], settings["record"]
    scratch: dict[str, Any] = {"warnings": [], "reference_fills": {}}
    cols = settings.get("cols") or select_columns(head.drop(columns=[ROW_COLUMN], errors="ignore"), config)
    out = fill_positions_from_ids(build_canonical(head, cols, config, scratch), scratch)
    needs = np.isfinite(out["POS"]) & ~np.isfinite(out["POS19"]) & ~np.isfinite(out["POS38"])
    if not needs.any():
        return ""
    build = resolve_input_build(out, record, config, lookup, settings["reference_root"], scratch)
    qc["build_detection"] = scratch["build_detection"]
    qc["warnings"].extend(scratch["warnings"])
    return build


def standardise_table(raw_path: Path, vcf_path: Path, record: dict[str, Any], config: dict[str, Any], reference_root: str,
                      lookup: pd.DataFrame, fasta: Any, qc: dict[str, Any], population: str = "ALL", population_used: str = "",
                      workers: int = 1) -> None:
    """Standardise a table or a GWAS-VCF (read into the same canonical columns) to GRCh38 GWAS-VCF."""
    vcf_input = use_vcf_reader(raw_path, config)
    chunk_rows = table_chunk_rows(config)
    tolerance = float(config.get("palindromic_tolerance", alleles.PALINDROMIC_TOLERANCE))
    settings = {"config": config, "record": record, "reference_root": reference_root, "population": population,
                "palindromic_tolerance": tolerance, "afs_label": f"1000G_{population_used or population}"}
    if vcf_input:
        header = read_vcf_header(raw_path)
        config = {**config, "_vcf_header_build": build_from_vcf_header(header)}
        settings["config"] = config
        declared = {line.split("ID=", 1)[1].split(",", 1)[0] for line in header if line.startswith("##FORMAT=<ID=")}
        optional = {column for key, column in VCF_FORMAT_COLUMNS.items() if column in gwas_vcf.OPTIONAL_FORMAT}
        settings["cols"] = {field: column for field, column in VCF_COLUMNS.items()
                            if column not in optional or any(VCF_FORMAT_COLUMNS[k] == column for k in declared if k in VCF_FORMAT_COLUMNS)}
        qc["source_columns"] = {key: f"FORMAT/{fmt}" for fmt, key in VCF_FORMAT_COLUMNS.items()}
        pieces = read_vcf_pieces(raw_path, chunk_rows, config)
    else:
        sep = normalize_delimiter_value(config.get("delimiter")) or infer_delimiter(raw_path)
        qc["input_delimiter"] = sep
        pieces = read_pieces(raw_path, sep, candidate_columns(raw_path, sep, config), chunk_rows, config)
    head_rows = max(BUILD_HEAD_ROWS, chunk_rows or 0) if chunk_rows else None
    head_pieces = read_vcf_pieces(raw_path, head_rows, config) if vcf_input else \
        read_pieces(raw_path, qc["input_delimiter"], candidate_columns(raw_path, qc["input_delimiter"], config), head_rows, config)
    settings["build"] = detect_build_from_head(next(iter(head_pieces)), settings, lookup, qc)
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    work = make_work_dir()
    try:
        started = time.perf_counter()
        results = sorted(run_pieces(pieces, settings, lookup, fasta, work, workers), key=lambda r: r["part"])
        wall = {"read_and_standardise_pieces": time.perf_counter() - started}
        qc["table_pieces"] = {"pieces": len(results), "rows_per_piece": chunk_rows or "all", "workers": workers}
        moments = merge_piece_qc(qc, results, record)
        check_af_concordance(moments, qc)
        merge_pieces([Path(r["part"]) for r in results], work, vcf_path, record, config, results[0]["decisions"], fasta, qc, wall,
                     reference_root, population)
        per_step: dict[str, float] = {}
        for result in results:
            for step, seconds in result["timings"].items():
                per_step[step] = per_step.get(step, 0.0) + seconds
        qc["timings_s"] = {"wall": {k: round(v, 1) for k, v in wall.items()},
                           "piece_steps_summed_over_pieces": {k: round(v, 1) for k, v in per_step.items()}}
    finally:
        shutil.rmtree(work, ignore_errors=True)
    qc["input_format"] = "gwas_vcf" if vcf_input else "table"


# ---------------------------------------------------------------------------------------------------
# GWAS-VCF inputs: read into the table pipeline
# ---------------------------------------------------------------------------------------------------
# A GWAS-VCF (e.g. from OpenGWAS) is read record by record into the canonical table columns and then goes
# through exactly the same steps as a table: build detection from the data, coordinates, alignment, palindromic
# resolution, deduplication and writing. Multi-allelic records are split into 1 row per ALT allele.

VCF_FORMAT_COLUMNS = {"ES": "BETA", "SE": "SE", "LP": "LOG10P", "AF": "EAF", "SS": "N", "NC": "NCASE", "EZ": "EZ", "SI": "IMPINFO"}
VCF_COLUMNS = {"CHR": "CHR", "POS": "POS", "SNP": "SNP", "SNPID": "SNPID", "EA": "EA", "OA": "OA", "BETA": "BETA", "SE": "SE",
               "LOG10P": "LOG10P", "EAF": "EAF", "N": "N", "NCASE": "NCASE", "EZ": "EZ", "IMPINFO": "IMPINFO"}


def _open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8", errors="replace")


def read_vcf_header(path: Path) -> list[str]:
    header: list[str] = []
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            header.append(line.rstrip("\n"))
    return header


def use_vcf_reader(raw_path: Path, config: dict[str, Any]) -> bool:
    declared = str(config.get("format", "") or "").strip().lower()
    if declared in {"vcf", "gwas_vcf", "gwas-vcf", "gwasvcf"}:
        return True
    if declared in {"table", "tsv", "csv", "txt", "text"}:
        return False
    return is_gwas_vcf(raw_path)


def vcf_records_to_table(chunk: pd.DataFrame, sample_column: int) -> pd.DataFrame:
    """Canonical columns from raw VCF records (columns 0-8 plus the chosen sample), 1 row per ALT allele."""
    chunk = chunk.reset_index(drop=True)
    for column in range(9):
        if column not in chunk.columns:
            chunk[column] = ""
    chunk = chunk.fillna("")
    out = pd.DataFrame({"CHR": chunk[0], "POS": chunk[1], "SNPID": chunk[2], "OA": chunk[3], "EA": chunk[4]})
    keys = chunk[8].astype(str).to_numpy()
    values = chunk[sample_column].astype(str).to_numpy() if sample_column in chunk.columns else np.full(len(chunk), "", dtype=object)
    parsed: dict[str, list[str]] = {name: [""] * len(chunk) for name in VCF_FORMAT_COLUMNS.values()}
    parsed["FORMAT_ID"] = [""] * len(chunk)
    for i, (k, v) in enumerate(zip(keys, values)):
        for key, value in zip(k.split(":"), v.split(":")):
            if key in VCF_FORMAT_COLUMNS:
                parsed[VCF_FORMAT_COLUMNS[key]][i] = value
            elif key == "ID":
                parsed["FORMAT_ID"][i] = value
    for name, column in parsed.items():
        out[name] = column
    ids = out["SNPID"].astype(str)
    out["SNP"] = np.where(ids.str.match(r"(?i)^rs\d+$"), ids, out["FORMAT_ID"])
    out = out.drop(columns=["FORMAT_ID"])
    multi = out["EA"].astype(str).str.contains(",", regex=False)
    if multi.any():
        split = ["EA", *VCF_FORMAT_COLUMNS.values()]
        for column in split:
            out[column] = out[column].astype(str).str.split(",")
        exploded = out[multi].explode(split)
        out = pd.concat([out[~multi].assign(**{c: out.loc[~multi, c].str[0] for c in split}), exploded])
        out = out.sort_index(kind="stable").reset_index(drop=True)
    out = out.replace({".": "", "NA": ""})
    # OpenGWAS writes LP = 999 for P = 0 (pygwasvcf reads it back as P = 0); treat it as missing so LP is recovered
    # from ES/SE like a P of 0 in a table
    out.loc[pd.to_numeric(out["LOG10P"], errors="coerce") == 999, "LOG10P"] = ""
    return out


def read_vcf_pieces(raw_path: Path, chunk_rows: int | None, config: dict[str, Any]):
    """Yield a GWAS-VCF in pieces of canonical columns, with ROW_COLUMN numbering the input records."""
    header = read_vcf_header(raw_path)
    if not header or not header[-1].startswith("#CHROM"):
        raise data_error(f"{raw_path} has no #CHROM header line; it is not a VCF")
    samples = header[-1].split("\t")[9:]
    wanted = str(config.get("sample", "") or "")
    sample_column = 9 + (samples.index(wanted) if wanted in samples else 0)
    offset = 0
    reader = pd.read_csv(raw_path, sep="\t", header=None, skiprows=len(header), dtype=str, na_filter=False,
                         quoting=csv.QUOTE_NONE, chunksize=chunk_rows or TABLE_CHUNK, comment=None)
    yielded = False
    for chunk in reader:
        piece = vcf_records_to_table(chunk, sample_column)
        piece[ROW_COLUMN] = np.arange(offset, offset + len(piece), dtype="int64")
        offset += len(piece)
        yielded = True
        yield piece
    if not yielded:
        empty = pd.DataFrame(columns=[*VCF_COLUMNS.values(), ROW_COLUMN])
        yield empty


# ---------------------------------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------------------------------

def provenance() -> dict[str, Any]:
    info = {"schema_version": OUTPUT_SCHEMA_VERSION, "python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "git_commit": ""}
    try:
        info["git_commit"] = subprocess.run(["git", "-C", str(Path(__file__).resolve().parent), "rev-parse", "HEAD"],
                                            capture_output=True, text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return info


def standardise_with_fallback(raw_path: Path, vcf_path: Path, record: dict[str, Any], config: dict[str, Any],
                              reference_root: str, lookup: pd.DataFrame, fasta: Any, qc: dict[str, Any], **options: Any) -> Path:
    """Standardise the staged file. A raw GWAS Catalog file that cannot be standardised (data or configuration error)
    falls back to the Catalog's harmonised file, from a fresh QC record; the fallback is recorded. Returns the file used."""
    base = {key: (list(value) if isinstance(value, list) else dict(value) if isinstance(value, dict) else value)
            for key, value in qc.items()}
    try:
        standardise_table(raw_path, vcf_path, record, config, reference_root, lookup, fasta, qc, **options)
        return raw_path
    except SystemExit as exc:
        metadata = record.get("metadata", {}) or {}
        url = metadata.get("harmonised_url", "")
        if exc.code not in (EXIT_DATA, EXIT_CONFIG) or metadata.get("source_type") != "gwas_catalog" \
                or metadata.get("file_kind") != "raw" or not url:
            raise
        harmonised = download_url(url, raw_path.parent)
        qc.clear()
        qc.update(base)
        qc["gwas_catalog_fallback"] = {"raw_file": raw_path.name, "harmonised_file": harmonised.name, "exit_code": exc.code}
        qc["warnings"].append(f"the raw GWAS Catalog file {raw_path.name} could not be standardised (exit {exc.code}); "
                              f"used the harmonised file {harmonised.name} instead")
        standardise_table(harmonised, vcf_path, record, config, reference_root, lookup, fasta, qc, **options)
        return harmonised


def build_record(record: dict[str, Any], vcf: Path, index: Path, metadata: Path, qc_json: Path, plot: Path, qc: dict[str, Any]) -> dict[str, Any]:
    out = dict(record)
    out.update({"standard_vcf": str(vcf), "standard_index": str(index), "metadata_json": str(metadata), "qc_json": str(qc_json),
                "report_html": str(plot), "qc": qc})
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-record", required=True)
    parser.add_argument("--config-root", required=True)
    parser.add_argument("--study-yaml", default="")
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--force", default="false")
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--workers", type=int, default=int(os.environ.get("GWAS_WORKERS", "1") or 1),
                        help="processes standardising pieces of a table in parallel (each holds its own reference lookup)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    record = read_json(args.stage_record)
    manifest_values = dict(record)  # the fingerprint uses the values exactly as prepare_manifest saw them
    record.update(resolve_trait_metadata(record))
    output_prefix = Path(record["output_prefix"])
    ensure_dir(output_prefix.parent)
    standard_vcf = output_prefix.with_suffix(".gwas.vcf.gz")
    standard_index = Path(str(standard_vcf) + ".tbi")
    metadata_json = output_prefix.with_suffix(".metadata.json")
    qc_json = output_prefix.with_suffix(".qc.json")
    plot = output_prefix.with_suffix(".report.html")
    force, dry_run = parse_bool(args.force), parse_bool(args.dry_run)
    config = resolve_row_config(args.config_root, record, args.study_yaml)
    fingerprint = settings_fingerprint(manifest_values, config, args.reference_root)

    if dry_run:
        write_json(args.output, build_record(record, standard_vcf, standard_index, metadata_json, qc_json, plot,
                                             {"dry_run": True, "settings_hash": fingerprint}))
        return
    if not force and all(p.exists() for p in (standard_vcf, standard_index, metadata_json, qc_json)):
        existing = read_json(qc_json)
        if existing.get("settings_hash") == fingerprint:
            result = build_record(record, standard_vcf, standard_index, metadata_json, qc_json, plot, existing)
            result["reused_existing_outputs"] = True
            write_json(args.output, result)
            return
    qc_json.unlink(missing_ok=True)  # a failed regeneration must not look complete

    gwas_vcf.require_htslib()
    raw_path = Path(record["raw_path"])
    if not raw_path.exists():
        raise FileNotFoundError(f"Staged raw file not found: {raw_path}")
    qc: dict[str, Any] = {"warnings": [], "started_at": utc_now(), "dry_run": False, "config": config.get("_config_path", ""),
                          "settings_hash": fingerprint, "provenance": provenance(), "reference_fills": {}}
    validate_config_keys(config, qc)
    population = config.get("reference_population") or ancestry_to_population(record.get("ancestry", ""))
    _, freqs38, _ = grch38_references(args.reference_root, population)
    if freqs38 is not None:
        lookup, used = pd.DataFrame(columns=alleles.LOOKUP_COLUMNS), freqs38.population
        qc["reference_frequencies"] = "1000 Genomes GRCh38 (2,504 unrelated phase 3 samples)"
    else:
        lookup, used = alleles.read_reference_lookup(args.reference_root, population)
    qc["reference_population"], qc["reference_population_used"] = population, used
    if not used:
        qc["warnings"].append(f"no 1000 Genomes lookup found under {args.reference_root}; no rsID/EAF fills or strand checks")
    elif used != population:
        qc["warnings"].append(f"lookup for population {population} not found; used {used}")
    if used and lookup.empty and freqs38 is None:
        qc["warnings"].append(f"the {used} reference lookup is empty")
    elif used and freqs38 is None and not alleles.lookup_has_alleles(lookup):
        qc["warnings"].append("reference lookup lacks A1/A2/AF1; rebuild it with make_1000genomes_lookup_tables.py")
    fasta = alleles.Fasta.find(args.reference_root)

    used_path = standardise_with_fallback(raw_path, standard_vcf, record, config, args.reference_root, lookup, fasta, qc,
                                          population=population, population_used=used or population, workers=max(1, args.workers))
    record["raw_path"] = str(used_path)
    qc["coordinate_build"] = "GRCh38"
    qc["completed_at"] = utc_now()
    qc["steps"] = gwas_steps.build_steps(qc)
    print(gwas_steps.format_steps(qc, record), file=sys.stderr)
    metadata = {
        "study_id": record["study_id"], "output_id": record["output_id"], "phenotype": record["phenotype"],
        "ancestry": record["ancestry"], "sex": record.get("sex", ""), "author": record["author"], "year": record["year"],
        "PMID": record["PMID"], "source_type": record["source_type"], "GWAS_location": record["GWAS_location"],
        "standard_format": "GWAS-VCF", "coordinate_build": "GRCh38", "raw_path": record["raw_path"],
        "pipeline_version": OUTPUT_SCHEMA_VERSION, "settings_hash": fingerprint,
    }
    write_json(metadata_json, metadata)
    write_json(qc_json, qc)
    write_json(args.output, build_record(record, standard_vcf, standard_index, metadata_json, qc_json, plot, qc))


if __name__ == "__main__":
    main()
