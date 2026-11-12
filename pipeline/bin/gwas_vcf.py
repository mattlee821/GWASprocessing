#!/usr/bin/env python3
"""GWAS-VCF writing (specification 1.0.0, as gwas2vcf writes it): header, body, chromosome order, bgzip + tabix.

Every output is GRCh38 (##reference and ##contig lines with assembly=GRCh38) and has 1 sample named by the output
ID. A ##SAMPLE line carries the harmonisation counts, study type, case/control totals and trait metadata, each key
defined by a ##META line. REF is the non-effect allele and ALT the effect allele; ES and AF are relative to ALT.
The file is bgzip-compressed with a tabix index.
"""

from __future__ import annotations

import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from gwas_pipeline_common import (
    OUTPUT_SCHEMA_VERSION,
    PRIMARY_CHROMS,
    command_exists,
    config_error,
    normalize_chromosome,
    normalize_chromosome_series,
    numeric_series,
)

GRCH38_LENGTHS = {
    "1": 248956422, "2": 242193529, "3": 198295559, "4": 190214555, "5": 181538259, "6": 170805979,
    "7": 159345973, "8": 145138636, "9": 138394717, "10": 133797422, "11": 135086622, "12": 133275309,
    "13": 114364328, "14": 107043718, "15": 101991189, "16": 90338345, "17": 83257441, "18": 80373285,
    "19": 58617616, "20": 64444167, "21": 46709983, "22": 50818468, "X": 156040895, "Y": 57227415, "MT": 16569,
}
RESERVED_INFO = {"SNPID", "HM"}
_RANK = {chrom: rank for rank, chrom in enumerate(PRIMARY_CHROMS, start=1)}

FORMAT_DEFINITIONS = {
    "ES": ("A", "Float", "Effect size estimate relative to the alternative allele"),
    "SE": ("A", "Float", "Standard error of effect size estimate"),
    "LP": ("A", "Float", "-log10 p-value for effect estimate"),
    "AF": ("A", "Float", "Alternate allele frequency in the association study, or the 1000 Genomes frequency where AFS says so"),
    "AFS": ("A", "String", "Source of AF: study (reported by the study) or 1000G_<POP> (1000 Genomes frequency of that "
            "super-population, used because the study reported none)"),
    "SS": ("A", "Float", "Sample size used to estimate genetic effect"),
    "NC": ("A", "Float", "Number of cases used to estimate genetic effect"),
    "EZ": ("A", "Float", "Z-score provided if it was used to derive the EFFECT and SE fields"),
    "SI": ("A", "Float", "Accuracy score of summary data imputation"),
    "ID": ("1", "String", "Study variant identifier: the dbSNP rsID, as in the ID column"),
}
# FORMAT keys always written, then those written only when the input provides them (canonical column: key), then ID.
BASE_FORMAT = ("ES", "SE", "LP", "AF", "AFS")
OPTIONAL_FORMAT = {"N": "SS", "NCASE": "NC", "EZ": "EZ", "IMPINFO": "SI"}
HM_DESCRIPTION = ("Harmonisation action. Non-palindromic: forward (as reported), swapped (effect allele was the reference "
                  "base; ES negated and AF to 1-AF), complemented (alleles on the other strand from the "
                  "study's: reported on the reverse strand, or complemented by liftover), complemented_swapped, "
                  "unverified (no reference sequence at the position). "
                  "Palindromic (A/T or C/G): palindrome_forward or palindrome_complemented (strand inferred from study vs "
                  "1000 Genomes frequency with tolerance 0.08), or palindrome_ambiguous (strand could not be inferred; kept "
                  "as reported), palindrome_ambiguous_complemented (strand could not be inferred, and liftover complemented the "
                  "alleles on a minus-strand block); each with _swapped when the alleles were then swapped")
INFO_DEFINITIONS = [
    '##INFO=<ID=SNPID,Number=1,Type=String,Description="Variant identifier in the source file">',
    f'##INFO=<ID=HM,Number=1,Type=String,Description="{HM_DESCRIPTION}">',
]
# ##SAMPLE keys with their ##META type and description, in the order written.
SAMPLE_KEYS = {
    "TotalVariants": ("Integer", "Total number of variants in input"),
    "VariantsNotRead": ("Integer", "Number of variants that could not be read: no GRCh38 position, or missing or invalid "
                        "alleles, effect size, standard error or P value"),
    "HarmonisedVariants": ("Integer", "Total number of harmonised variants (written to this file)"),
    "VariantsNotHarmonised": ("Integer", "Total number of variants that could not be harmonised: SNV alleles matching neither "
                              "strand of the GRCh38 reference, indel alleles not matching its sequence, or duplicates of an earlier record"),
    "SwitchedAlleles": ("Integer", "Total number of variants whose effect and other alleles were swapped to match the reference"),
    "AFFilled": ("Integer", "Number of variants whose AF is the 1000 Genomes frequency because the study reported none"),
    "TotalControls": ("Float", "Total number of controls in the association study (total sample size for a continuous trait)"),
    "TotalCases": ("Float", "Total number of cases in the association study"),
    "StudyType": ("String", "Type of GWAS study [Continuous or CaseControl]"),
    "StudyID": ("String", "Study identifier"),
    "Phenotype": ("String", "Trait"),
    "Ancestry": ("String", "Ancestry of the study sample (1000 Genomes super-population code)"),
    "Sex": ("String", "Sex of the study sample"),
    "PMID": ("String", "PubMed identifier of the study"),
}


def chrom_rank(chrom: Any) -> int:
    return _RANK.get(str(normalize_chromosome(chrom)), 99)


def chrom_ranks(chrom: pd.Series) -> np.ndarray:
    """Numeric sort rank of each chromosome (1-22, X, Y, MT; 99 for anything else)."""
    return normalize_chromosome_series(chrom).map(lambda c: _RANK.get(str(c), 99)).to_numpy(dtype="int64")


def sort_frame(df: pd.DataFrame) -> pd.DataFrame:
    rank = chrom_ranks(df["CHR"])
    order = np.lexsort((numeric_series(df["POS38"]).to_numpy(), rank))
    return df.iloc[order]


def sanitize_header_value(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value if value is not None else "")).strip()


def sample_name(record: dict[str, Any]) -> str:
    """The output ID (e.g. GWAScatalog_GCST009004): unique per dataset; the sample column and the ##SAMPLE ID."""
    value = sanitize_header_value(record.get("output_id")) or sanitize_header_value(record.get("study_id")) or "GWAS"
    return re.sub(r"[\s,;=<>]", "_", value)


def info_key(name: Any) -> str:
    key = re.sub(r"[^A-Za-z0-9_.]", "_", str(name))
    if not re.match(r"[A-Za-z_]", key):
        key = f"X_{key}"
    return f"SRC_{key}" if key in RESERVED_INFO else key


def info_value(value: Any) -> str:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return "."
    text = re.sub(r"[;=,\s]", "_", str(value).strip())
    return text or "."


def header_value(value: Any) -> str:
    """A ##SAMPLE value: bare when it is a simple token, otherwise double-quoted."""
    text = sanitize_header_value(value)
    if text == "":
        return "NA"
    if re.fullmatch(r"[A-Za-z0-9_.+:/-]+", text):
        return text
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def contig_header_lines() -> list[str]:
    return [f"##contig=<ID={chrom},length={length},assembly=GRCh38>" for chrom, length in GRCH38_LENGTHS.items()]


def format_keys(columns: Iterable[str]) -> list[str]:
    """FORMAT keys for a file: the base keys, the optional keys whose source column is present, then ID."""
    columns = set(columns)
    return [*BASE_FORMAT, *[key for column, key in OPTIONAL_FORMAT.items() if column in columns], "ID"]


def sample_header_lines(sample: str, fields: dict[str, Any]) -> list[str]:
    meta = [f'##META=<ID={key},Number=1,Type={kind},Description="{description}">'
            for key, (kind, description) in SAMPLE_KEYS.items() if key in fields]
    values = "".join(f",{key}={header_value(fields[key])}" for key in SAMPLE_KEYS if key in fields)
    return [*meta, f"##SAMPLE=<ID={sample}{values}>"]


def table_header(record: dict[str, Any], extras: Iterable[str], keys: list[str], sample_fields: dict[str, Any],
                 reference: str = "GRCh38", settings_hash: str = "") -> list[str]:
    extras_info = [f'##INFO=<ID={info_key(name)},Number=1,Type=String,Description="Source column {sanitize_header_value(name)}">'
                   for name in extras]
    sample = sample_name(record)
    return [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now(timezone.utc).strftime('%Y%m%d')}",
        "##source=GWASprocessing",
        f"##GWASprocessing_version={OUTPUT_SCHEMA_VERSION}",
        "##GWASVCF=1.0",
        *([f"##settings_hash={settings_hash}"] if settings_hash else []),
        f"##reference={sanitize_header_value(reference) or 'GRCh38'}",
        *contig_header_lines(),
        '##FILTER=<ID=PASS,Description="All filters passed">',
        *INFO_DEFINITIONS,
        *extras_info,
        *[f'##FORMAT=<ID={key},Number={FORMAT_DEFINITIONS[key][0]},Type={FORMAT_DEFINITIONS[key][1]},'
          f'Description="{FORMAT_DEFINITIONS[key][2]}">' for key in keys],
        *sample_header_lines(sample, sample_fields),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample,
    ]


def _fmt_float(values: Any) -> np.ndarray:
    arr = numeric_series(pd.Series(values)).to_numpy(dtype="float64")
    out = np.full(arr.shape, ".", dtype=object)
    finite = np.isfinite(arr)
    out[finite] = [f"{v:.8g}" for v in arr[finite]]
    return out


def _fmt_int(values: Any) -> np.ndarray:
    arr = numeric_series(pd.Series(values)).to_numpy(dtype="float64")
    out = np.full(arr.shape, ".", dtype=object)
    finite = np.isfinite(arr)
    out[finite] = np.floor(arr[finite] + 0.5).astype("int64").astype(str)
    return out


def _alleles(values: Any) -> np.ndarray:
    return pd.Series(values).astype("string").str.strip().str.upper().fillna("N").to_numpy(dtype=object)


def table_body_lines(df: pd.DataFrame, extras: Iterable[str], keys: list[str], afs_label: str = "1000G") -> list[str]:
    """VCF records for oriented rows (OA = REF, EA = ALT); `keys` from format_keys(), `afs_label` e.g. 1000G_EUR."""
    snp = df["SNP"].astype("string").str.strip()
    is_rs = snp.str.match(r"(?i)^rs\d+$").fillna(False).to_numpy(dtype=bool)
    ids = np.where(is_rs, snp.str.lower().fillna(".").to_numpy(dtype=object), ".").astype(object)
    snpid = [info_value(v) for v in df["SNPID"]] if "SNPID" in df.columns else ["."] * len(df)
    hm = df["HM"].fillna(".").astype(str).to_numpy(dtype=object) if "HM" in df.columns else np.full(len(df), ".", dtype=object)
    info = np.array([f"SNPID={s};HM={h}" for s, h in zip(snpid, hm)], dtype=object)
    for name in extras:
        values = [info_value(v) for v in df[name]] if name in df.columns else ["."] * len(df)
        info = info + np.array([f";{info_key(name)}={v}" for v in values], dtype=object)
    source = df["EAF_SOURCE"].fillna("").astype(str).to_numpy() if "EAF_SOURCE" in df.columns else np.full(len(df), "")
    af = _fmt_float(df["EAF"])
    afs = np.where(af == ".", ".", np.where(source == "reference", afs_label, "study")).astype(object)
    columns = {"ES": _fmt_float(df["BETA"]), "SE": _fmt_float(df["SE"]), "LP": _fmt_float(df["LP"]), "AF": af, "AFS": afs,
               "ID": ids}
    for column, key in OPTIONAL_FORMAT.items():
        if key in keys:
            columns[key] = _fmt_float(df[column]) if column in df.columns else np.full(len(df), ".", dtype=object)
    sample = columns[keys[0]]
    for key in keys[1:]:
        sample = sample + ":" + columns[key]
    chrom = normalize_chromosome_series(df["CHR"]).astype(str).to_numpy(dtype=object)
    pos = _fmt_int(df["POS38"])
    lines = (chrom + "\t" + pos + "\t" + ids + "\t" + _alleles(df["OA"]) + "\t" + _alleles(df["EA"])
             + "\t.\tPASS\t" + info + "\t" + ":".join(keys) + "\t" + sample)
    return list(lines)


def require_htslib() -> None:
    if not (command_exists("bgzip") and command_exists("tabix")):
        raise config_error(
            "bgzip and tabix (htslib) are required to write indexed GWAS-VCF files but were not found on PATH. "
            "Install htslib (e.g. `brew install htslib`, `conda install -c bioconda htslib`, or `module load htslib`) "
            "or run `bash src/source.sh --tools-only`."
        )


def bgzip_and_index(sources: Path | list[Path], vcf_path: Path) -> None:
    """Concatenate `sources` (e.g. header then body), bgzip to `vcf_path` atomically, tabix-index it, remove the sources."""
    require_htslib()
    parts = [Path(sources)] if isinstance(sources, (str, Path)) else [Path(p) for p in sources]
    tmp = vcf_path.with_name(f".{vcf_path.name}.tmp-{os.getpid()}")
    try:
        with open(tmp, "wb") as handle:
            cat = subprocess.Popen(["cat", *map(str, parts)], stdout=subprocess.PIPE)
            subprocess.run(["bgzip", "-c"], stdin=cat.stdout, stdout=handle, check=True)
            cat.stdout.close()
            if cat.wait() != 0:
                raise RuntimeError(f"cat exited {cat.returncode} while compressing {vcf_path}")
        os.replace(tmp, vcf_path)
    finally:
        tmp.unlink(missing_ok=True)
    index = Path(str(vcf_path) + ".tbi")
    index.unlink(missing_ok=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_path)], check=True)
    if not index.exists():
        raise RuntimeError(f"tabix did not create {index}")
    for part in parts:
        part.unlink(missing_ok=True)
