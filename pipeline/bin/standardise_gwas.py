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
import gwas_stats
import gwas_steps
import gwas_vcf
from gwas_pipeline_common import (
    ALT_ALIASES,
    BETA_ALIASES,
    CHR_ALIASES,
    CI_LOWER_ALIASES,
    CI_UPPER_ALIASES,
    COMPLEMENT,
    EA_ALIASES,
    EAF_ALIASES,
    LOGP_ALIASES,
    N_ALIASES,
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
    data_error,
    ensure_dir,
    infer_delimiter,
    is_gwas_vcf,
    make_work_dir,
    normalize_build_value,
    normalize_chromosome,
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
    "OR": OR_ALIASES, "CI_LOWER": CI_LOWER_ALIASES, "CI_UPPER": CI_UPPER_ALIASES,
}
# Allele-dependent columns must come from one family (original or GWAS Catalog harmonised hm_*).
FAMILY_FIELDS = ("EA", "OA", "BETA", "EAF", "OR", "CI_LOWER", "CI_UPPER", "CHR", "POS")
KNOWN_CONFIG_KEYS = {"description", "input_build", "delimiter", "columns", "exclude", "exclusions", "joins",
                     "transforms", "extras", "analysis_qc", "format", "reference_population", "_config_path"}
KNOWN_TRANSFORMS = {"strip_hash_header", "uppercase_alleles", "derive_oa_from_ref_alt"}
REQUIRED_FIELDS = ("CHR", "POS38", "EA", "OA", "BETA", "SE", "LP")
RETENTION_WARNING = 0.5
CONTRADICTION_MIN_HITS = 10  # a declared build is only overruled by this many matching positions
AF_CONCORDANCE_MIN_N = 1000
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
    configured = normalize_build_value(config.get("input_build") or record.get("input_build"))
    if configured:
        return configured, "manifest_or_yaml"
    metadata = normalize_build_value(recursive_metadata_value(record.get("metadata", {}),
                                                              {"build", "genome_build", "assembly", "coordinate_build", "genome_assembly"}))
    return (metadata, "source_metadata") if metadata else ("", "")


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
    for field in ("POS", "POS19", "POS38", "EAF", "BETA", "SE", "P", "LOG10P", "OR", "CI_LOWER", "CI_UPPER"):
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


def contradicts(declared: str, detection: dict[str, Any]) -> bool:
    detected = detection.get("build", "")
    if not declared or not detected or detected == declared:
        return False
    return int(detection.get(f"hits_{detected}", 0)) >= CONTRADICTION_MIN_HITS


def resolve_input_build(out: pd.DataFrame, record: dict[str, Any], config: dict[str, Any], lookup: pd.DataFrame,
                        reference_root: str, qc: dict[str, Any]) -> str:
    declared, declared_source = declared_build(record, config)
    rows = out[np.isfinite(out["POS"]) & ~np.isfinite(out["POS19"]) & ~np.isfinite(out["POS38"])]
    detection: dict[str, Any] = {"build": "", "reason": "reference lookup without alleles"}
    if alleles.lookup_has_alleles(lookup):
        lifter, _ = coords.load_liftover(reference_root, "hg38", "hg19")
        detection = coords.detect_build(rows["CHR"], rows["POS"], rows["EA"], rows["OA"], lookup, lifter)
    qc["build_detection"] = {**detection, "declared": declared, "declared_source": declared_source}
    names = {"hg19": "GRCh37", "hg38": "GRCh38"}
    if contradicts(declared, detection):
        raise config_error(
            f"Row {record.get('study_id')}/{record.get('output_id')}: input_build is declared as {names[declared]} "
            f"({declared_source}) but the positions match {names[detection['build']]} "
            f"(GRCh37 hits {detection['hits_hg19']}, GRCh38 hits {detection['hits_hg38']} of {detection['sampled']} sampled). "
            "Correct input_build in the manifest or study YAML."
        )
    build = declared or detection.get("build", "")
    if not build:
        raise config_error(
            f"Row {record.get('study_id')}/{record.get('output_id')}: could not determine the genome build "
            f"(GRCh37 hits {detection.get('hits_hg19', 0)}, GRCh38 hits {detection.get('hits_hg38', 0)} of "
            f"{detection.get('sampled', 0)} sampled positions). Set input_build (hg19 or hg38) in the manifest or study YAML."
        )
    qc["build_detection"].update({"build": build, "method": "declared" if declared else "detected"})
    return build


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


def establish_dual_build(out: pd.DataFrame, reference_root: str, qc: dict[str, Any]) -> pd.DataFrame:
    """Give every row a lifted GRCh38 POS38 (dropping failures) and, where possible, a GRCh37 POS19."""
    out = out.copy()
    out["LIFT_STRAND"] = "+"
    report = qc.setdefault("liftover", {})
    need38 = (~np.isfinite(out["POS38"]) & np.isfinite(out["POS19"])).to_numpy()
    if need38.any():
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
    if need19.any():
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
INCOMPLETE_CAUSES = ("rsid_not_in_coordinate_index", "no_position_or_rsid", "liftover_to_grch38_failed", "chromosome_not_primary",
                     "invalid_position", "missing_allele", "missing_beta", "missing_or_invalid_se", "missing_or_invalid_p")


def incomplete_causes(out: pd.DataFrame, checks: dict[str, np.ndarray]) -> np.ndarray:
    """The first cause (INCOMPLETE_CAUSES) for each row failing a check in `checks`; '' for complete rows."""
    pos19, pos38 = np.isfinite(out["POS19"]).to_numpy(), np.isfinite(out["POS38"]).to_numpy()
    rsid = out["SNP"].astype("string").str.match(r"^rs\d+$").fillna(False).to_numpy(dtype=bool)
    no_position = ~pos19 & ~pos38
    conditions = [
        no_position & rsid,
        no_position & ~rsid,
        pos19 & ~pos38,
        ~checks["CHR"],
        ~checks["POS38"],
        ~checks["EA"] | ~checks["OA"],
        ~checks["BETA"],
        ~checks["SE"],
        ~checks["LP"],
    ]
    return np.select(conditions, list(INCOMPLETE_CAUSES), default="").astype(object)


def drop_incomplete_rows(out: pd.DataFrame, qc: dict[str, Any], maf_threshold: float = 0.01) -> pd.DataFrame:
    checks = {
        "CHR": out["CHR"].astype(str).isin(PRIMARY_CHROMS).to_numpy(),
        "POS38": (np.isfinite(out["POS38"]) & (out["POS38"] > 0)).to_numpy(),
        "EA": (upper_alleles(out["EA"]) != "").to_numpy(),
        "OA": (upper_alleles(out["OA"]) != "").to_numpy(),
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


def annotate_with_reference(out: pd.DataFrame, lookup: pd.DataFrame, qc: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, float]]:
    """Reference annotation and rsID/EAF fills; returns the frame and its allele-frequency concordance moments."""
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
        qc["warnings"].append(f"retained {written} of {input_rows} input rows ({retained:.0%}); see missing_required_by_field, "
                              "liftover and qc_removed")
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
ROW_COLUMN = "_gwasprocessing_row"  # input row number, carried through each piece to restore input order
FIRST_PIECE_ONLY = ("build_detection", "analysis_qc_settings", "source_columns")


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
    out = fill_coordinates_by_rsid(out, reference_root, lookup, qc)
    lap("rsid_coordinates")
    out = establish_dual_build(out, reference_root, qc)
    lap("liftover")
    maf_threshold = float({**alleles.ANALYSIS_QC_DEFAULTS, **(config.get("analysis_qc") or {})}.get("maf_threshold", 0.01))
    out = drop_incomplete_rows(out, qc, maf_threshold)
    out, moments = annotate_with_reference(out, lookup, qc)
    lap("reference_annotation")
    if settings["analysis_qc"]:
        out = alleles.apply_analysis_qc(out, config.get("analysis_qc"), qc, deduplicate_rows=False)
    low, high = alleles.dedup_alleles(out)  # the key deduplication uses, before alignment can complement alleles
    lap("analysis_qc")
    if fasta is not None:
        out, action = alleles.orient_to_reference(out, fasta)
    else:
        action = np.full(len(out), "keep", dtype=object)
    lap("grch38_alignment")
    write_piece_lines(part, out, rows[out.index.to_numpy()], low, high, action, extras)
    lap("write_piece")
    decisions = {"cols": cols, "extras": extras, "build": build,
                 "negate_log10p": bool(qc.get("p_value_source", {}).get("log10p_negated", False)),
                 "derive_oa": bool(qc.get("derived_oa_from_ref_alt", False))}
    return {"qc": qc, "moments": moments, "decisions": decisions, "part": str(part), "timings": timings}


def write_piece_lines(part: Path, out: pd.DataFrame, rows: np.ndarray, low: np.ndarray, high: np.ndarray,
                      action: np.ndarray, extras: list[str]) -> None:
    """One line per row: chromosome rank, POS38, input row, dedup allele pair, alignment action, LP, VCF record.

    The leading fields let `sort` restore genome order with input order as the tie-break; LP is kept at full
    precision for lambda GC. Rows that alignment will drop are kept so they still claim their dedup key."""
    with open(part, "w", encoding="utf-8") as handle:
        step = 500_000
        for start in range(0, len(out), step):
            piece = out.iloc[start:start + step]
            sl = slice(start, start + step)
            lines = np.array(gwas_vcf.table_body_lines(piece, extras), dtype=object)
            rank = gwas_vcf.chrom_ranks(piece["CHR"]).astype(str).astype(object)
            pos = piece["POS38"].to_numpy(dtype="float64").astype("int64").astype(str).astype(object)
            lp = np.char.mod("%.17g", piece["LP"].to_numpy(dtype="float64")).astype(object)
            keys = (rank + "\t" + pos + "\t" + rows[sl].astype(str).astype(object) + "\t" + low[sl].astype(object)
                    + "\t" + high[sl].astype(object) + "\t" + action[sl].astype(object) + "\t" + lp + "\t" + lines)
            handle.write("\n".join(keys.tolist()) + "\n")


def merge_pieces(parts: list[Path], work: Path, vcf_path: Path, record: dict[str, Any], extras: list[str],
                 deduplicate: bool, fasta: Any, qc: dict[str, Any], timings: dict[str, float] | None = None) -> None:
    """Sort all pieces by position, deduplicate, drop alignment failures, write and index the VCF."""
    timings = {} if timings is None else timings
    started = time.perf_counter()
    plain = work / "table.vcf"
    actions: dict[str, int] = {}
    lp_values = array("d")
    duplicates = written = 0
    current: tuple[str, str] | None = None
    seen: set[tuple[str, str]] = set()
    command = ["sort", "-t", "\t", "-k1,1n", "-k2,2n", "-k3,3n", "-T", str(work), *map(str, parts)]
    env = {**os.environ, "LC_ALL": "C"}
    with open(plain, "w", encoding="utf-8") as out, \
            subprocess.Popen(command, stdout=subprocess.PIPE, text=True, encoding="utf-8", env=env) as proc:
        out.write("\n".join(gwas_vcf.table_header(record, extras)) + "\n")
        for line in proc.stdout:
            rank, pos, _row, low, high, action, lp, body = line.split("\t", 7)
            if deduplicate:
                if (rank, pos) != current:
                    current, seen = (rank, pos), set()
                if (low, high) in seen:
                    duplicates += 1
                    continue
                seen.add((low, high))
            actions[action] = actions.get(action, 0) + 1
            if action == "drop":
                continue
            out.write(body)
            lp_values.append(float(lp))
            written += 1
    if proc.returncode != 0:
        raise RuntimeError(f"sort exited {proc.returncode} while merging {len(parts)} table pieces")
    if deduplicate:
        qc["qc_removed"]["duplicate"] = qc["qc_removed"].get("duplicate", 0) + duplicates
    if fasta is not None:
        qc["reference_alignment"] = alleles.alignment_summary_from_counts(actions, fasta)
    else:
        qc["reference_alignment"] = {"performed": False, "reason": "no GRCh38 FASTA under <reference_root>/genome"}
        qc["warnings"].append("REF/ALT not aligned to the GRCh38 reference: no FASTA found (run src/reference_data/download_fasta_and_chains.sh)")
    check_output_rows(written, record, qc)
    qc["lambda_gc"] = gwas_stats.lambda_gc(np.frombuffer(lp_values, dtype="float64") if written else np.array([]))
    timings["sort_merge_deduplicate"] = time.perf_counter() - started
    started = time.perf_counter()
    gwas_vcf.bgzip_and_index(plain, vcf_path)
    timings["compress_and_index"] = time.perf_counter() - started


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


def _init_worker(reference_root: str, population: str, fasta_path: str) -> None:
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


def standardise_table(raw_path: Path, vcf_path: Path, record: dict[str, Any], config: dict[str, Any], reference_root: str,
                      analysis_qc: bool, lookup: pd.DataFrame, fasta: Any, qc: dict[str, Any], population: str = "ALL",
                      workers: int = 1) -> None:
    sep = normalize_delimiter_value(config.get("delimiter")) or infer_delimiter(raw_path)
    usecols = candidate_columns(raw_path, sep, config)
    qc["input_delimiter"] = sep
    chunk_rows = table_chunk_rows(config)
    settings = {"config": config, "record": record, "reference_root": reference_root, "analysis_qc": bool(analysis_qc),
                "population": population}
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    work = make_work_dir()
    try:
        started = time.perf_counter()
        pieces = read_pieces(raw_path, sep, usecols, chunk_rows, config)
        results = sorted(run_pieces(pieces, settings, lookup, fasta, work, workers), key=lambda r: r["part"])
        wall = {"read_and_standardise_pieces": time.perf_counter() - started}
        qc["table_pieces"] = {"pieces": len(results), "rows_per_piece": chunk_rows or "all", "workers": workers}
        moments = merge_piece_qc(qc, results, record)
        qc["analysis_qc"] = bool(analysis_qc)
        check_af_concordance(moments, qc)
        dedup = bool(analysis_qc) and bool(qc.get("analysis_qc_settings", {}).get("deduplicate", True))
        merge_pieces([Path(r["part"]) for r in results], work, vcf_path, record, results[0]["decisions"]["extras"],
                     dedup, fasta, qc, wall)
        per_step: dict[str, float] = {}
        for result in results:
            for step, seconds in result["timings"].items():
                per_step[step] = per_step.get(step, 0.0) + seconds
        qc["timings_s"] = {"wall": {k: round(v, 1) for k, v in wall.items()},
                           "piece_steps_summed_over_pieces": {k: round(v, 1) for k, v in per_step.items()}}
    finally:
        shutil.rmtree(work, ignore_errors=True)
    qc["input_format"] = "table"


# ---------------------------------------------------------------------------------------------------
# Passthrough path (pre-built GWAS-VCF)
# ---------------------------------------------------------------------------------------------------

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


def build_from_header(header: list[str]) -> str:
    assembly_text = ""
    contig1_length = None
    for line in header:
        if line.startswith("##contig"):
            match = re.search(r"assembly=([^,>]+)", line, flags=re.I)
            if match and not assembly_text:
                assembly_text = match.group(1)
            length = re.search(r"ID=(?:chr)?1,[^>]*length=(\d+)", line, flags=re.I)
            if length and contig1_length is None:
                contig1_length = int(length.group(1))
        elif line.startswith(("##reference", "##genomeBuild", "##genome_build", "##build", "##assembly")) and not assembly_text:
            assembly_text = line.split("=", 1)[-1].strip()
    build = normalize_build_value(assembly_text)
    if not build:
        text = assembly_text.lower()
        build = "hg38" if ("grch38" in text or "hg38" in text or "/b38/" in text) else ("hg19" if ("grch37" in text or "hg19" in text or "b37" in text) else "")
    if not build and contig1_length:
        build = {248956422: "hg38", 249250621: "hg19"}.get(contig1_length, "")
    return {"hg19": "GRCh37", "hg38": "GRCh38"}.get(build, "")


def detect_vcf_build(path: str | Path) -> str:
    """'GRCh38' / 'GRCh37' from a VCF header, or '' when the header does not say."""
    try:
        return build_from_header(read_vcf_header(Path(path)))
    except OSError:
        return ""


def sample_vcf_records(path: Path, limit: int) -> pd.DataFrame:
    rows = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t", 5)
            if len(fields) >= 5:
                rows.append((normalize_chromosome(fields[0]), fields[1], fields[3], fields[4]))
            if len(rows) >= limit:
                break
    return pd.DataFrame(rows, columns=["CHR", "POS", "REF", "ALT"])


def resolve_passthrough_build(raw_path: Path, header: list[str], record: dict[str, Any], config: dict[str, Any],
                              lookup: pd.DataFrame, reference_root: str, qc: dict[str, Any]) -> str:
    manifest = normalize_build_value(config.get("input_build") or record.get("input_build"))
    header_build = {"GRCh37": "hg19", "GRCh38": "hg38"}.get(build_from_header(header), "")
    metadata = declared_build({"metadata": record.get("metadata", {})}, {})[0]
    declared, source = next(((b, s) for b, s in ((manifest, "manifest_or_yaml"), (header_build, "vcf_header"),
                                                   (metadata, "source_metadata")) if b), ("", ""))
    detection: dict[str, Any] = {"build": "", "reason": "reference lookup without alleles"}
    if alleles.lookup_has_alleles(lookup):
        sample = sample_vcf_records(raw_path, coords.DETECTION_SAMPLE)
        lifter, _ = coords.load_liftover(reference_root, "hg38", "hg19")
        detection = coords.detect_build(sample["CHR"], sample["POS"], sample["ALT"], sample["REF"], lookup, lifter)
    qc["build_detection"] = {**detection, "declared": declared, "declared_source": source}
    names = {"hg19": "GRCh37", "hg38": "GRCh38"}
    if contradicts(declared, detection):
        raise config_error(
            f"Row {record.get('study_id')}/{record.get('output_id')}: the source VCF is declared {names[declared]} ({source}) "
            f"but its positions match {names[detection['build']]}. Correct input_build in the manifest or study YAML."
        )
    build = declared or detection.get("build", "")
    if not build:
        raise config_error(
            f"Row {record.get('study_id')}/{record.get('output_id')}: the source VCF does not declare its build and it could not "
            "be detected. Set input_build (hg19 or hg38) in the manifest or study YAML."
        )
    qc["build_detection"].update({"build": build, "method": "declared" if declared else "detected"})
    return names[build]


def sample_size_from_metadata(record: dict[str, Any]) -> int | None:
    value = recursive_metadata_value(record.get("metadata", {}), {"sample_size"})
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(round(number)) if math.isfinite(number) and number > 0 else None


def format_integer(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "."
    return str(int(math.floor(number + 0.5))) if math.isfinite(number) else "."


def _replace_info(info: str, updates: dict[str, str], remove: tuple[str, ...] = ()) -> str:
    parts = [p for p in info.split(";") if p and p != "." and p.split("=", 1)[0] not in set(updates) | set(remove)]
    parts.extend(f"{key}={value}" for key, value in updates.items())
    return ";".join(parts) if parts else "."


MANAGED_HEADER = ("##fileformat", "##fileDate", "##contig", "##reference", "##study_id", "##phenotype", "##ancestry",
                  "##sex", "##PMID", "##source=GWASprocessing", "##GWASprocessing_version", "##INFO=<ID=POS19,")


def passthrough_header(header: list[str], record: dict[str, Any], add_ss: bool) -> list[str]:
    kept = [line for line in header[:-1] if not line.startswith(MANAGED_HEADER)]
    extra = [gwas_vcf.POS19_INFO]
    if add_ss:
        extra.append(gwas_vcf.FORMAT_HEADER[-1])
    return ["##fileformat=VCFv4.2", f"##fileDate={utc_now()[:10].replace('-', '')}", "##source=GWASprocessing",
            f"##GWASprocessing_version={OUTPUT_SCHEMA_VERSION}", *gwas_vcf.assembly_header_lines(),
            *gwas_vcf.metadata_header_lines(record), *kept, *extra, header[-1]]


def lift_passthrough_body(raw_path: Path, plain_path: Path, header_lines: list[str], build: str, reference_root: str,
                          ss_value: int | None, tag_strand: bool, qc: dict[str, Any]) -> int:
    lifter, reason = coords.load_liftover(reference_root, "hg19", "hg38") if build == "GRCh37" else coords.load_liftover(reference_root, "hg38", "hg19")
    if build == "GRCh37" and lifter is None:
        raise config_error(f"Passthrough source is GRCh37 but {reason}. Add hg19ToHg38.over.chain.gz "
                           "(see src/reference_data/download_fasta_and_chains.sh).")
    if build == "GRCh38" and lifter is None:
        qc["warnings"].append(f"{reason}; POS19 left missing")
    counts = {"malformed_record": 0, "non_primary_chromosome": 0, "bad_position": 0, "minus_strand": 0,
              "minus_strand_non_snv_dropped": 0}
    status_counts: dict[str, int] = {}
    written = rounded = records = 0

    def flush(chunk: list[list[str]], handle) -> None:
        nonlocal written, rounded
        if not chunk:
            return
        chroms = pd.Series([f[0] for f in chunk])
        positions = pd.Series([f[1] for f in chunk])
        result = coords.lift(chroms, positions, lifter) if lifter is not None else None
        for i, fields in enumerate(chunk):
            strand = "+"
            if build == "GRCh37":
                if result.status[i] != "ok":
                    continue
                strand = result.strand[i]
                original = fields[1]
                if strand == "-":
                    if not (fields[3] in COMPLEMENT and fields[4] in COMPLEMENT):
                        counts["minus_strand_non_snv_dropped"] += 1
                        continue
                    fields[3], fields[4] = COMPLEMENT[fields[3]], COMPLEMENT[fields[4]]
                    counts["minus_strand"] += 1
                fields[1] = str(int(result.pos[i]))
                pos19 = original
            else:
                pos19 = str(int(result.pos[i])) if result is not None and result.status[i] == "ok" else "."
                if result is not None and result.status[i] == "ok" and result.strand[i] == "-":
                    strand = "-"
            updates = {"POS19": pos19}
            if tag_strand and strand == "-":
                updates[STRAND_TAG] = "-"
            fields[7] = _replace_info(fields[7], updates)
            if "SS" not in fields[8].split(":"):
                if ss_value is None:
                    counts["records_without_ss"] = counts.get("records_without_ss", 0) + 1
                else:
                    fields[8] += ":SS"
                    for j in range(9, len(fields)):
                        fields[j] += f":{ss_value}"
                    counts["sample_size_filled_records"] = counts.get("sample_size_filled_records", 0) + 1
            else:
                ss_index = fields[8].split(":").index("SS")
                for j in range(9, len(fields)):
                    values = fields[j].split(":")
                    if ss_index < len(values):
                        normalised = format_integer(values[ss_index])
                        if normalised not in (".", values[ss_index]):
                            values[ss_index] = normalised
                            fields[j] = ":".join(values)
                            rounded += 1
            handle.write("\t".join(fields) + "\n")
            written += 1
        if result is not None:
            for key, value in result.counts().items():
                status_counts[key] = status_counts.get(key, 0) + value

    with _open_text(raw_path) as src, open(plain_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(header_lines) + "\n")
        chunk: list[list[str]] = []
        for line in src:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            records += 1
            if len(fields) < 8:
                counts["malformed_record"] += 1
                continue
            chrom = normalize_chromosome(fields[0])
            if pd.isna(chrom) or str(chrom) not in PRIMARY_CHROMS:
                counts["non_primary_chromosome"] += 1
                continue
            if not re.fullmatch(r"\d+", fields[1]):
                counts["bad_position"] += 1
                continue
            fields[0] = str(chrom)
            fields[3], fields[4] = fields[3].upper(), fields[4].upper()
            chunk.append(fields)
            if len(chunk) >= 200_000:
                flush(chunk, handle)
                chunk = []
        flush(chunk, handle)
    key = "hg19_to_hg38" if build == "GRCh37" else "hg38_to_hg19"
    qc["input_rows"] = records
    qc["lifted_rows"] = written
    qc["sample_size_filled_records"] = counts.pop("sample_size_filled_records", 0)
    missing_ss = counts.pop("records_without_ss", 0)
    if missing_ss:
        qc["warnings"].append(f"{missing_ss} records have no FORMAT/SS and the source metadata gives no sample size")
    qc.setdefault("liftover", {})[key] = {**status_counts, **counts}
    if rounded:
        qc["sample_size_values_rounded"] = rounded
    return written


def _flip_number(value: str, how: str) -> str:
    try:
        number = float(value)
    except ValueError:
        return value
    if not math.isfinite(number):
        return value
    return f"{(-number if how == 'negate' else 1.0 - number):.8g}"


PASSTHROUGH_CHUNK = 1_000_000  # records per chunk when rewriting a pre-built VCF (GWAS_PASSTHROUGH_CHUNK overrides)


def _format_values(chunk: pd.DataFrame, key: str, sample_column: int = 9) -> np.ndarray:
    """FORMAT field `key` of one sample column as strings ('' where absent), grouped by FORMAT layout."""
    values = np.full(len(chunk), "", dtype=object)
    if sample_column not in chunk.columns:
        return values
    for layout in pd.unique(chunk[8]):
        keys = str(layout).split(":")
        if key not in keys:
            continue
        rows = (chunk[8] == layout).to_numpy()
        values[rows] = chunk.loc[rows, sample_column].str.split(":").str[keys.index(key)].fillna("").to_numpy(dtype=object)
    return values


def _passthrough_canon(chunk: pd.DataFrame) -> pd.DataFrame:
    info = chunk[7]
    canon = pd.DataFrame({
        "CHR": chunk[0].to_numpy(dtype=object), "POS38": numeric_series(chunk[1]).to_numpy(),
        "POS19": numeric_series(info.str.extract(r"(?:^|;)POS19=(\d+)", expand=False)).to_numpy(),
        "SNP": chunk[2].str.extract(r"^(rs\d+)$", flags=re.I, expand=False).str.lower().to_numpy(dtype=object),
        "EA": chunk[4].to_numpy(dtype=object), "OA": chunk[3].to_numpy(dtype=object),
        "EAF": numeric_series(pd.Series(_format_values(chunk, "AF"), dtype=object)).to_numpy(),
        "LIFT_STRAND": np.where(info.str.contains(f"{STRAND_TAG}=-", regex=False).to_numpy(), "-", "+"),
        "BETA": np.nan,
    })
    canon["EAF_SOURCE"] = np.where(np.isfinite(canon["EAF"]), "input", "")
    return canon


def _flip_sample(sample: str, keys: list[str]) -> str:
    values = sample.split(":")
    for k, key in enumerate(keys):
        if k < len(values) and key in ("ES", "EZ"):
            values[k] = _flip_number(values[k], "negate")
        elif k < len(values) and key == "AF":
            values[k] = _flip_number(values[k], "complement")
    return ":".join(values)


def _apply_alignment(chunk: pd.DataFrame, action: np.ndarray, ea: np.ndarray, oa: np.ndarray) -> None:
    """Rewrite REF/ALT, and for swapped records ES/EZ/AF (all samples) and INFO/AF, in place."""
    complement = np.isin(action, ["complement", "complement_swap"])
    if complement.any():
        ea[complement] = complement_alleles(pd.Series(ea[complement], dtype=object)).to_numpy(dtype=object)
        oa[complement] = complement_alleles(pd.Series(oa[complement], dtype=object)).to_numpy(dtype=object)
    swap = np.isin(action, ["swap", "complement_swap"])
    ref, alt = np.where(swap, ea, oa), np.where(swap, oa, ea)
    chunk[3], chunk[4] = ref, alt
    for i in np.where(swap)[0]:
        keys = str(chunk.iat[i, 8]).split(":")
        for column in range(9, chunk.shape[1]):
            chunk.iat[i, column] = _flip_sample(str(chunk.iat[i, column]), keys)
        chunk.iat[i, 7] = ";".join(f"AF={_flip_number(p[3:], 'complement')}" if p.startswith("AF=") else p
                                   for p in str(chunk.iat[i, 7]).split(";"))


def _has_body(plain_path: Path, header_count: int) -> str:
    """First data line after the header ('' when the file has no records)."""
    with open(plain_path, "r", encoding="utf-8") as handle:
        return next((line for i, line in enumerate(handle) if i >= header_count and line.strip()), "")


def _add_counts(total: dict[str, int], part: dict[str, Any]) -> None:
    for key, value in part.items():
        if isinstance(value, (int, np.integer)):
            total[key] = total.get(key, 0) + int(value)


def rewrite_passthrough_body(plain_path: Path, final_path: Path, header_count: int, config: dict[str, Any], analysis_qc: bool,
                             lookup: pd.DataFrame, fasta: Any, qc: dict[str, Any]) -> int:
    """Analysis QC and GRCh38 alignment over the lifted, sorted body, chunk by chunk (bounded memory).

    The body is sorted by position, so duplicates can only span a chunk boundary at the boundary position;
    the keys kept at the previous chunk's last position are carried forward to catch them.
    """
    chunk_size = int(os.environ.get("GWAS_PASSTHROUGH_CHUNK", PASSTHROUGH_CHUNK))
    removed: dict[str, int] = {}
    palindromic: dict[str, int] = {}
    maf_detail: dict[str, int] = {}
    actions: dict[str, int] = {}
    lp_parts: list[np.ndarray] = []
    carry: set[tuple[int, str, str]] = set()
    written = 0
    dedup = analysis_qc and {**alleles.ANALYSIS_QC_DEFAULTS, **(config.get("analysis_qc") or {})}.get("deduplicate", True)
    with open(plain_path, "r", encoding="utf-8") as src:
        header = [next(src).rstrip("\n") for _ in range(header_count)]
    if not _has_body(plain_path, header_count):
        return 0
    with open(final_path, "w", encoding="utf-8") as out:
        out.write("\n".join(header) + "\n")
        reader = pd.read_csv(plain_path, sep="\t", header=None, skiprows=header_count, dtype=str, na_filter=False,
                             quoting=csv.QUOTE_NONE, chunksize=chunk_size)
        for chunk in reader:
            chunk = chunk.reset_index(drop=True)
            canon = _passthrough_canon(chunk)
            keep = np.ones(len(chunk), dtype=bool)
            ea, oa = canon["EA"].to_numpy(dtype=object).copy(), canon["OA"].to_numpy(dtype=object).copy()
            if analysis_qc:
                chunk_qc: dict[str, Any] = {}
                kept = alleles.apply_analysis_qc(alleles.annotate_with_reference(canon, lookup), config.get("analysis_qc"), chunk_qc)
                keep[:] = False
                keep[kept.index.to_numpy()] = True
                ea[kept.index.to_numpy()] = kept["EA"].to_numpy(dtype=object)
                oa[kept.index.to_numpy()] = kept["OA"].to_numpy(dtype=object)
                _add_counts(removed, chunk_qc["qc_removed"])
                _add_counts(palindromic, chunk_qc.get("palindromic_resolution", {}))
                _add_counts(maf_detail, chunk_qc.get("maf_detail", {}))
                qc["analysis_qc_settings"] = chunk_qc["analysis_qc_settings"]
                if dedup:
                    key = coords.position_key(canon["CHR"], canon["POS38"])
                    low = np.where(ea <= oa, ea, oa)
                    high = np.where(ea <= oa, oa, ea)
                    crossing = keep & np.array([(int(k), str(lo), str(hi)) in carry for k, lo, hi in zip(key, low, high)], dtype=bool)
                    keep &= ~crossing
                    removed["duplicate"] = removed.get("duplicate", 0) + int(crossing.sum())
                    if keep.any():
                        last = key[keep].max()
                        carry = {(int(k), str(lo), str(hi)) for k, lo, hi, kp in zip(key, low, high, keep) if kp and k == last}
            action = np.full(len(chunk), "keep", dtype=object)
            if fasta is not None and keep.any():
                rows = np.where(keep)[0]
                action[rows] = alleles.reference_alignment(canon["CHR"].iloc[rows], canon["POS38"].iloc[rows], pd.Series(ea[rows]),
                                                           pd.Series(oa[rows]), fasta)["action"].to_numpy()
                _add_counts(actions, {str(k): int(v) for k, v in zip(*np.unique(action[rows], return_counts=True))})
            _apply_alignment(chunk, action, ea, oa)
            chunk[7] = chunk[7].str.replace(f";{STRAND_TAG}=-", "", regex=False)
            final = keep & (action != "drop")
            selected = chunk[final]
            if len(selected):
                columns = [selected[c].astype(str) for c in selected.columns]
                out.write("\n".join(columns[0].str.cat(columns[1:], sep="\t").tolist()) + "\n")
                lp_parts.append(numeric_series(pd.Series(_format_values(selected.reset_index(drop=True), "LP"), dtype=object)).to_numpy())
            written += int(final.sum())
    if analysis_qc:
        qc["qc_removed"] = removed
        qc["palindromic_resolution"] = palindromic
        if maf_detail:
            qc["maf_detail"] = maf_detail
    if fasta is not None:
        expanded = np.concatenate([np.full(v, k, dtype=object) for k, v in actions.items()]) if actions else np.array([], dtype=object)
        qc["reference_alignment"] = alleles.alignment_summary(expanded, fasta)
    else:
        qc["reference_alignment"] = {"performed": False, "reason": "no GRCh38 FASTA under <reference_root>/genome"}
        qc["warnings"].append("REF/ALT not re-checked against the GRCh38 reference: no FASTA found")
    qc["lambda_gc"] = gwas_stats.lambda_gc(np.concatenate(lp_parts) if lp_parts else np.array([]))
    return written


def scan_passthrough_body(plain_path: Path, header_count: int, qc: dict[str, Any]) -> int:
    """Count records and compute lambda GC without rewriting (nothing to change)."""
    lp_parts: list[np.ndarray] = []
    count = 0
    first = _has_body(plain_path, header_count)
    if not first:
        return 0
    usecols = [8, 9] if first.count("\t") >= 9 else [8]
    reader = pd.read_csv(plain_path, sep="\t", header=None, skiprows=header_count, dtype=str, na_filter=False,
                         quoting=csv.QUOTE_NONE, usecols=usecols, chunksize=int(os.environ.get("GWAS_PASSTHROUGH_CHUNK", PASSTHROUGH_CHUNK)))
    for chunk in reader:
        chunk = chunk.reset_index(drop=True)
        count += len(chunk)
        lp_parts.append(numeric_series(pd.Series(_format_values(chunk, "LP"), dtype=object)).to_numpy())
    qc["lambda_gc"] = gwas_stats.lambda_gc(np.concatenate(lp_parts) if lp_parts else np.array([]))
    qc["reference_alignment"] = {"performed": False, "reason": "no GRCh38 FASTA under <reference_root>/genome"}
    qc["warnings"].append("REF/ALT not re-checked against the GRCh38 reference: no FASTA found")
    return count


def standardise_existing_vcf(raw_path: Path, vcf_path: Path, record: dict[str, Any], config: dict[str, Any], reference_root: str,
                             analysis_qc: bool, lookup: pd.DataFrame, fasta: Any, qc: dict[str, Any]) -> None:
    header = read_vcf_header(raw_path)
    if not header or not header[-1].startswith("#CHROM"):
        raise data_error(f"{raw_path} has no #CHROM header line; it is not a VCF")
    build = resolve_passthrough_build(raw_path, header, record, config, lookup, reference_root, qc)
    qc["source_build"] = build
    format_line_has_ss = any(line.startswith("##FORMAT=<ID=SS,") for line in header)
    # Headers can declare SS while records omit it (OpenGWAS), so the fill is decided per record.
    ss_value = sample_size_from_metadata(record)
    new_header = passthrough_header(header, record, add_ss=ss_value is not None and not format_line_has_ss)
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    work = make_work_dir()
    try:
        plain = work / "lifted.vcf"
        lift_passthrough_body(raw_path, plain, new_header, build, reference_root, ss_value, analysis_qc, qc)
        qc["sample_size_source"] = "gwasinfo" if qc["sample_size_filled_records"] else ("source" if format_line_has_ss else "missing")
        gwas_vcf.sort_vcf_body(plain)
        qc["analysis_qc"] = bool(analysis_qc)
        if analysis_qc or fasta is not None:
            final = work / "final.vcf"
            written = rewrite_passthrough_body(plain, final, len(new_header), config, analysis_qc, lookup, fasta, qc)
        else:
            final = plain
            written = scan_passthrough_body(plain, len(new_header), qc)
        if written == 0:
            raise data_error(f"No variants remain for {record.get('study_id')}/{record.get('output_id')} after lifting the "
                             f"source VCF to GRCh38 (liftover: {qc.get('liftover')}).")
        qc["output_rows"] = written
        gwas_vcf.bgzip_and_index(final, vcf_path)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    qc["input_format"] = "gwas_vcf"
    qc["metadata_headers_normalised"] = True


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


def use_passthrough(raw_path: Path, config: dict[str, Any]) -> bool:
    declared = str(config.get("format", "") or "").strip().lower()
    if declared in {"vcf", "gwas_vcf", "gwas-vcf", "gwasvcf"}:
        return True
    if declared in {"table", "tsv", "csv", "txt", "text"}:
        return False
    return is_gwas_vcf(raw_path)


def build_record(record: dict[str, Any], vcf: Path, index: Path, metadata: Path, qc_json: Path, plot: Path, qc: dict[str, Any]) -> dict[str, Any]:
    out = dict(record)
    out.update({"standard_vcf": str(vcf), "standard_index": str(index), "metadata_json": str(metadata), "qc_json": str(qc_json),
                "plot_tiff": str(plot), "qc": qc})
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage-record", required=True)
    parser.add_argument("--config-root", required=True)
    parser.add_argument("--study-yaml", default="")
    parser.add_argument("--reference-root", required=True)
    parser.add_argument("--force", default="false")
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--analysis-qc", default="false")
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
    plot = output_prefix.with_suffix(".tiff")
    force, dry_run, analysis_qc = parse_bool(args.force), parse_bool(args.dry_run), parse_bool(args.analysis_qc)
    config = resolve_row_config(args.config_root, record, args.study_yaml)
    fingerprint = settings_fingerprint(manifest_values, config, analysis_qc, args.reference_root)

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
    output_prefix.with_suffix(".steps.txt").unlink(missing_ok=True)

    gwas_vcf.require_htslib()
    raw_path = Path(record["raw_path"])
    if not raw_path.exists():
        raise FileNotFoundError(f"Staged raw file not found: {raw_path}")
    qc: dict[str, Any] = {"warnings": [], "started_at": utc_now(), "dry_run": False, "config": config.get("_config_path", ""),
                          "settings_hash": fingerprint, "provenance": provenance(), "reference_fills": {}}
    validate_config_keys(config, qc)
    population = config.get("reference_population") or ancestry_to_population(record.get("ancestry", ""))
    lookup, used = alleles.read_reference_lookup(args.reference_root, population)
    qc["reference_population"], qc["reference_population_used"] = population, used
    if not used:
        qc["warnings"].append(f"no 1000 Genomes lookup found under {args.reference_root}; no rsID/EAF fills or strand checks")
    elif used != population:
        qc["warnings"].append(f"lookup for population {population} not found; used {used}")
    if used and lookup.empty:
        qc["warnings"].append(f"the {used} reference lookup is empty")
    elif used and not alleles.lookup_has_alleles(lookup):
        qc["warnings"].append("reference lookup lacks A1/A2/AF1; rebuild it with make_1000genomes_lookup_tables.py")
    fasta = alleles.Fasta.find(args.reference_root)

    if use_passthrough(raw_path, config):
        standardise_existing_vcf(raw_path, standard_vcf, record, config, args.reference_root, analysis_qc, lookup, fasta, qc)
    else:
        standardise_table(raw_path, standard_vcf, record, config, args.reference_root, analysis_qc, lookup, fasta, qc,
                          population=population, workers=max(1, args.workers))
    qc["coordinate_build"] = "GRCh38"
    qc["completed_at"] = utc_now()
    qc["steps"] = gwas_steps.build_steps(qc)
    steps_log = output_prefix.with_suffix(".steps.txt")
    steps_log.write_text(gwas_steps.format_steps(qc, record), encoding="utf-8")
    print(steps_log.read_text(encoding="utf-8"), file=sys.stderr)
    metadata = {
        "study_id": record["study_id"], "output_id": record["output_id"], "phenotype": record["phenotype"],
        "ancestry": record["ancestry"], "sex": record.get("sex", ""), "author": record["author"], "year": record["year"],
        "PMID": record["PMID"], "source_type": record["source_type"], "GWAS_location": record["GWAS_location"],
        "standard_format": "GWAS-VCF", "coordinate_build": "GRCh38", "raw_path": record["raw_path"],
        "pipeline_version": OUTPUT_SCHEMA_VERSION, "settings_hash": fingerprint,
    }
    write_json(metadata_json, metadata)
    write_json(qc_json, qc)
    write_json(args.output, {**build_record(record, standard_vcf, standard_index, metadata_json, qc_json, plot, qc),
                             "steps_log": str(steps_log)})


if __name__ == "__main__":
    main()
