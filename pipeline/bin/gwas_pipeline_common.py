#!/usr/bin/env python3
"""Shared helpers for the GWAS-VCF standardisation pipeline."""

from __future__ import annotations

import csv
import gzip
import hashlib
import html.parser
import json
import math
import os
import re
import shutil
import ssl
import subprocess
import sys
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml


MANIFEST_COLUMNS = ["GWAS_location", "phenotype", "ancestry", "author", "year", "PMID"]
MANIFEST_OPTIONAL_COLUMNS = ["sex", "source_type", "input_build", "population", "delimiter", "format", "study_yaml"]
STATE_COLUMNS = [
    "row_hash",
    "study_id",
    "output_id",
    "gwas_location",
    "standard_vcf",
    "standard_index",
    "metadata_json",
    "qc_json",
    "plot_tiff",
    "updated_at",
]

CHR_ALIASES = ["CHR", "CHROM", "#CHROM", "chrom", "chr", "chromosome", "hm_chrom", "Chrom", "Chrom_x"]
POS_ALIASES = ["POS", "BP", "base_pair_location", "hm_pos", "position", "pos", "Pos", "Pos_x", "chrom_start"]
POS19_ALIASES = ["POS19", "pos19", "POS_hg19", "START_hg19"]
POS38_ALIASES = ["POS38", "pos38", "POS_hg38", "START_hg38"]
RSID_ALIASES = ["SNP", "rsid", "rsids", "hm_rsid", "rsids_y", "ID", "MarkerName", "variant", "variant_id"]
SNPID_ALIASES = ["SNPID", "source_snp_id", "Name", "ID", "MarkerName", "variant_id"]
EA_ALIASES = ["EA", "effect_allele", "hm_effect_allele", "effectAllele", "effectAllele_y", "ALLELE1", "A1", "Tested_Allele", "Allele1", "ALT"]
OA_ALIASES = ["OA", "other_allele", "hm_other_allele", "otherAllele", "otherAllele_y", "ALLELE0", "A2", "Other_Allele", "Allele2", "REF"]
EAF_ALIASES = ["EAF", "effectAlleleFreq", "effect_allele_frequency", "hm_effect_allele_frequency", "eaf_ref", "Freq_Tested_Allele", "A1FREQ", "A1_FREQ", "Freq1", "Freq1.Hapmap", "freq", "AF"]
BETA_ALIASES = ["BETA", "Beta", "Effect", "ES", "beta", "hm_beta", "b", "effect"]
SE_ALIASES = ["SE", "StdErr", "standard_error", "se"]
P_ALIASES = ["P", "Pvalue", "Pval", "pval", "p_value", "P-value", "P.value", "p-value", "p.value", "p"]
LOGP_ALIASES = ["LOG10P", "minus_log10_pval", "Plog10", "LP"]
N_ALIASES = ["N", "TotalSampleSize", "OBS_CT", "n", "samplesize", "SS"]
REF_ALIASES = ["REF", "ref"]
ALT_ALIASES = ["ALT", "alt"]

ANCESTRY_TO_POPULATION = {
    "african": "AFR",
    "afr": "AFR",
    "african-american": "AFR",
    "aa": "AFR",
    "american": "AMR",
    "amr": "AMR",
    "east-asian": "EAS",
    "eas": "EAS",
    "asian": "EAS",
    "european": "EUR",
    "eur": "EUR",
    "white": "EUR",
    "central-south-asian": "SAS",
    "south-asian": "SAS",
    "sas": "SAS",
    "middle-east": "ALL",
    "combined": "ALL",
    "multi": "ALL",
    "all": "ALL",
    "unknown": "ALL",
    "": "ALL",
}


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.links.append(value)


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def sanitize(value: Any, default: str = "unknown") -> str:
    text = str(value if value is not None else "").strip().lower()
    text = re.sub(r"^pmid", "", text)
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    return text or default


def safe_path_component(value: Any, default: str = "unknown") -> str:
    text = str(value if value is not None else "").strip()
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("._-")
    return text or default


def short_hash(value: str, n: int = 10) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:n]


def file_sha256(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def study_id_from_row(row: dict[str, str]) -> str:
    return "_".join(
        [
            sanitize(row.get("author")),
            sanitize(row.get("year")),
            sanitize(row.get("PMID")),
        ]
    )


def infer_source_type(location: str) -> str:
    loc = str(location).strip()
    lower = loc.lower()
    if "ebi.ac.uk/gwas/studies/" in lower or re.fullmatch(r"GCST\d+", loc, flags=re.I):
        return "gwas_catalog"
    if "opengwas.io/datasets/" in lower or re.fullmatch(r"[a-z]+-[a-z]-\d+", lower):
        return "opengwas"
    if lower.startswith("syn"):
        return "synapse"
    if lower.startswith("s3://"):
        return "s3"
    if lower.startswith("gs://"):
        return "gs"
    if lower.startswith(("http://", "https://", "ftp://")):
        return "direct_url"
    return "local"


def normalize_source_type(source_type: str) -> str:
    source = sanitize(source_type).replace("-", "_")
    aliases = {
        "gwascatalog": "gwas_catalog",
        "gwas_catalogue": "gwas_catalog",
        "catalog": "gwas_catalog",
        "catalogue": "gwas_catalog",
        "opengwas": "opengwas",
        "open_gwas": "opengwas",
        "direct": "direct_url",
        "directurl": "direct_url",
        "url": "direct_url",
        "http": "direct_url",
        "https": "direct_url",
        "local_manual": "local",
        "manual": "local",
    }
    return aliases.get(source, source)


def extract_accession_or_basename(location: str, source_type: str) -> str:
    if source_type == "gwas_catalog":
        match = re.search(r"GCST\d+", location, flags=re.I)
        return match.group(0).upper() if match else "gwas_catalog"
    if source_type == "opengwas":
        match = re.search(r"([a-z]+-[a-z]-\d+)", location, flags=re.I)
        return match.group(1).lower() if match else "opengwas"
    parsed = urllib.parse.urlparse(location)
    name = Path(parsed.path).name if parsed.path else Path(location).name
    if not name:
        name = source_type
    for suffix in [".tsv.gz", ".txt.gz", ".tbl.gz", ".csv.gz", ".vcf.gz", ".tsv", ".txt", ".tbl", ".csv", ".vcf", ".gz", ".zip", ".tar"]:
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name or source_type


def output_id_from_row(row: dict[str, str], source_type: str) -> str:
    source_id = extract_accession_or_basename(row["GWAS_location"], source_type)
    prefix = {
        "gwas_catalog": "GWAScatalog",
        "opengwas": "openGWAS",
        "direct_url": "directURL",
        "local": "local",
        "synapse": "synapse",
        "s3": "s3",
        "gs": "gs",
    }.get(source_type, source_type)
    unique_id = safe_path_component(source_id)
    if source_type not in {"gwas_catalog", "opengwas", "synapse"}:
        unique_id = f"{unique_id}_{short_hash(row['GWAS_location'], 8)}"
    return f"{prefix}_{unique_id}"


def row_hash(row: dict[str, str]) -> str:
    payload = {col: str(row.get(col, "")).strip() for col in MANIFEST_COLUMNS}
    payload.update({col: str(row.get(col, "")).strip() for col in MANIFEST_OPTIONAL_COLUMNS if str(row.get(col, "")).strip()})
    return short_hash(json.dumps(payload, sort_keys=True), 16)


def read_manifest(path: str | Path) -> list[dict[str, str]]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [col for col in MANIFEST_COLUMNS if col not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"Manifest is missing required columns: {', '.join(missing)}")
        rows = []
        for i, row in enumerate(reader, start=2):
            if not any(str(row.get(col, "")).strip() for col in MANIFEST_COLUMNS):
                continue
            clean = {col: str(row.get(col, "")).strip() for col in MANIFEST_COLUMNS}
            for col in MANIFEST_OPTIONAL_COLUMNS:
                clean[col] = str(row.get(col, "")).strip() if row.get(col) is not None else ""
            if not clean["GWAS_location"]:
                raise ValueError(f"Manifest row {i} has empty GWAS_location")
            rows.append(clean)
    return rows


def read_state(path: str | Path) -> dict[str, dict[str, str]]:
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["row_hash"]: row for row in reader if row.get("row_hash")}


def write_state(path: str | Path, rows: Iterable[dict[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=STATE_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in STATE_COLUMNS})


def state_row_complete(row: dict[str, str], qcplot: bool) -> bool:
    required = ["standard_vcf", "standard_index", "metadata_json", "qc_json"]
    if qcplot:
        required.append("plot_tiff")
    for key in required:
        value = row.get(key, "")
        if not value or not Path(value).exists() or Path(value).stat().st_size == 0:
            return False
    return True


def write_json(path: str | Path, payload: dict[str, Any]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def read_json(path: str | Path) -> dict[str, Any]:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def load_yaml(path: str | Path | None) -> dict[str, Any]:
    if not path:
        return {}
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"YAML config must be a mapping: {path}")
    data["_config_path"] = str(path)
    return data


def resolve_config_yaml_path(config_root: str | Path, path_value: str) -> Path:
    path = Path(os.path.expandvars(path_value))
    if path.is_absolute():
        return path
    root = Path(config_root)
    root_candidate = root / path
    if root_candidate.exists():
        return root_candidate
    cwd_candidate = Path.cwd() / path
    if cwd_candidate.exists():
        return cwd_candidate
    return root_candidate


def resolve_study_config(
    config_root: str | Path,
    study_id: str,
    explicit_yaml: str | None = None,
    manifest_yaml: str | None = None,
) -> dict[str, Any]:
    if explicit_yaml:
        return load_yaml(resolve_config_yaml_path(config_root, explicit_yaml))
    if manifest_yaml:
        return load_yaml(resolve_config_yaml_path(config_root, manifest_yaml))
    candidate = Path(config_root) / f"{study_id}.yaml"
    return load_yaml(candidate)


def merge_dicts(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge_dicts(out[key], value)
        else:
            out[key] = value
    return out


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def urlopen(url: str, timeout: int = 120):
    """Open HTTPS URLs with certifi CA roots when the Python install needs them."""
    context = None
    try:
        import certifi

        context = ssl.create_default_context(cafile=certifi.where())
    except Exception:
        context = None
    return urllib.request.urlopen(url, timeout=timeout, context=context)


def url_links(url: str) -> list[str]:
    with urlopen(url, timeout=60) as response:
        html = response.read().decode("utf-8", errors="replace")
    parser = LinkParser()
    parser.feed(html)
    return parser.links


def guess_filename_from_url(url: str, fallback: str = "downloaded_gwas") -> str:
    parsed = urllib.parse.urlparse(url)
    name = Path(parsed.path).name
    if name:
        return urllib.parse.unquote(name)
    query = urllib.parse.parse_qs(parsed.query)
    for key in ["filename", "file", "download"]:
        if query.get(key):
            return Path(query[key][0]).name
    return fallback


def download_url(url: str, out_dir: str | Path, filename: str | None = None, retries: int = 3, timeout: int = 600) -> Path:
    out_dir = ensure_dir(out_dir)
    filename = filename or guess_filename_from_url(url)
    out = out_dir / sanitize_filename(filename)
    if out.exists() and out.stat().st_size > 0:
        return out
    part = out.with_name(out.name + ".part")
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        part.unlink(missing_ok=True)
        try:
            with urlopen(url, timeout=timeout) as response, open(part, "wb") as handle:
                expected_size = response.headers.get("Content-Length")
                shutil.copyfileobj(response, handle)
            if part.stat().st_size == 0:
                raise IOError(f"Downloaded empty file from {url}")
            if expected_size and part.stat().st_size != int(expected_size):
                raise IOError(
                    f"Incomplete download from {url}: expected {expected_size} bytes, got {part.stat().st_size}"
                )
            part.replace(out)
            return out
        except Exception as exc:
            last_error = exc
            part.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(min(30, 2**attempt))
    if last_error:
        raise last_error
    return out


def sanitize_filename(filename: str) -> str:
    filename = Path(filename).name
    filename = re.sub(r"[^A-Za-z0-9._+-]+", "_", filename)
    return filename or "gwas"


def shell_copy_or_download(location: str, out_dir: str | Path) -> Path:
    out_dir = ensure_dir(out_dir)
    source_type = infer_source_type(location)
    if source_type == "local":
        src = Path(location.replace("file://", ""))
        if not src.exists():
            raise FileNotFoundError(f"Local GWAS file not found: {src}")
        dest = out_dir / sanitize_filename(src.name)
        if not dest.exists():
            shutil.copy2(src, dest)
        return dest
    if source_type in {"direct_url", "gwas_catalog", "opengwas"}:
        return download_url(location, out_dir)
    if source_type in {"s3", "gs"}:
        cmd = ["aws", "s3", "cp", location, str(out_dir)] if source_type == "s3" else ["gsutil", "cp", location, str(out_dir)]
        subprocess.run(cmd, check=True)
        return out_dir / Path(urllib.parse.urlparse(location).path).name
    raise ValueError(f"Unsupported generic source type: {source_type}")


def pick_existing_column(columns: Iterable[str], aliases: list[str]) -> str | None:
    columns_list = list(columns)
    lower_to_real = {c.lower(): c for c in columns_list}
    for alias in aliases:
        if alias in columns_list:
            return alias
        if alias.lower() in lower_to_real:
            return lower_to_real[alias.lower()]
    return None


def infer_delimiter(path: str | Path) -> str:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                sample = line
                break
        else:
            return "\t"
    counts = {"\t": sample.count("\t"), ",": sample.count(","), " ": len(re.findall(r"\s+", sample))}
    delimiter = max(counts, key=counts.get)
    return r"\s+" if delimiter == " " and counts[delimiter] > 0 else delimiter


def read_table(path: str | Path, delimiter: str | None = None, comment: str | None = None) -> pd.DataFrame:
    path = Path(path)
    sep = delimiter or infer_delimiter(path)
    compression = "gzip" if str(path).endswith(".gz") else None
    return pd.read_csv(path, sep=sep, compression=compression, comment=comment, low_memory=False)


def read_reference_lookup(reference_root: str | Path, population: str) -> pd.DataFrame:
    root = Path(reference_root) / "1000genomes" / "phase3" / "lookup"
    for pop in [population, "ALL", "EUR"]:
        for suffix in [".tsv.gz", ".tsv"]:
            path = root / f"{pop}{suffix}"
            if path.exists():
                return read_table(path, delimiter="\t")
    return pd.DataFrame()


def ancestry_to_population(ancestry: str) -> str:
    return ANCESTRY_TO_POPULATION.get(sanitize(ancestry).replace("_", "-"), "ALL")


def normalize_chromosome(value: Any) -> Any:
    if value is None or pd.isna(value):
        return pd.NA
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na", "null", "."}:
        return pd.NA
    text = re.sub(r"^chr", "", text, flags=re.I).strip().upper()
    aliases = {
        "23": "X",
        "24": "Y",
        "25": "MT",
        "26": "MT",
        "M": "MT",
        "MITO": "MT",
        "MITOCHONDRIAL": "MT",
        "XY": "X",
        "PAR": "X",
    }
    if text in aliases:
        return aliases[text]
    try:
        number = float(text)
        if number.is_integer():
            integer = int(number)
            if 1 <= integer <= 22:
                return str(integer)
            if integer == 23:
                return "X"
            if integer == 24:
                return "Y"
            if integer in {25, 26}:
                return "MT"
    except ValueError:
        pass
    if text in {"X", "Y", "MT"}:
        return text
    return text


def numeric_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def p_to_logp(value: Any) -> float | None:
    try:
        p = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(p) or p <= 0:
        return None
    return -math.log10(p)


def logp_to_p(value: Any) -> float | None:
    try:
        lp = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(lp):
        return None
    return 10 ** (-lp)


def parse_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def is_gwas_vcf(path: str | Path) -> bool:
    name = str(path).lower()
    if name.endswith((".vcf", ".vcf.gz", ".gwas.vcf", ".gwas.vcf.gz")):
        return True
    try:
        opener = gzip.open if name.endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            first = handle.readline()
        return first.startswith("##fileformat=VCF")
    except Exception:
        return False


def vcf_records(path: str | Path) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    sample_name = "sample"
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                if len(parts) > 9:
                    sample_name = parts[9]
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            row: dict[str, Any] = {
                "CHR": parts[0],
                "POS38": parts[1],
                "SNP": parts[2],
                "OA": parts[3],
                "EA": parts[4],
            }
            info = {}
            for item in parts[7].split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    info[key] = value
            row.update(info)
            if len(parts) >= 10:
                keys = parts[8].split(":")
                values = parts[9].split(":")
                row.update({k: values[i] if i < len(values) else "" for i, k in enumerate(keys)})
            row["sample"] = sample_name
            rows.append(row)
    return pd.DataFrame(rows)


def command_exists(command: str) -> bool:
    return shutil.which(command) is not None
