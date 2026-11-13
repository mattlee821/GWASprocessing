#!/usr/bin/env python3
"""Shared helpers for the GWAS-VCF standardisation pipeline."""

from __future__ import annotations

import csv
import gzip
import hashlib
import html.parser
import io
import json
import math
import os
import re
import shutil
import ssl
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
import yaml


# Version of the output contract (VCF layout, QC JSON fields, standardisation rules). Bump it when a
# change alters what is written, so completed rows are regenerated on the next run.
OUTPUT_SCHEMA_VERSION = "3.0.0"

# sysexits.h codes. Nextflow terminates the whole run on EX_CONFIG (an operator must fix configuration)
# and records the row as failed without retrying on EX_DATAERR (the input itself is unusable).
EXIT_CONFIG = 78
EXIT_DATA = 65


def config_error(message: str) -> SystemExit:
    print(message, file=sys.stderr)
    return SystemExit(EXIT_CONFIG)


def data_error(message: str) -> SystemExit:
    print(message, file=sys.stderr)
    return SystemExit(EXIT_DATA)


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
    "report_html",
    "settings_hash",
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
# Effect-allele frequency names (never MAF, which is not tied to the effect allele); allele-tagged forms as in
# MungeSumstats' dictionary (FREQ1, A1FREQ, FRQ_A1, ...). Freq_Tested_Allele_in_HRS: Yengo et al. 2018.
EAF_ALIASES = ["EAF", "effectAlleleFreq", "effect_allele_frequency", "hm_effect_allele_frequency", "eaf_ref", "Freq_Tested_Allele",
               "Freq_Tested_Allele_in_HRS", "A1FREQ", "A1_FREQ", "Freq1", "Freq1.Hapmap", "freq", "AF", "FRQ", "FRQ_A1", "FREQ_A1",
               "A1_AF", "AF1", "AF_Allele2", "A2FREQ", "EFFECT_ALLELE_FREQ", "Freq_Effect_Allele", "effect_allele_freq", "ALT_FREQS",
               "EAF_HRC", "EAF_UKB"]
BETA_ALIASES = ["BETA", "Beta", "Effect", "ES", "beta", "hm_beta", "b", "effect"]
SE_ALIASES = ["SE", "StdErr", "standard_error", "se"]
P_ALIASES = ["P", "Pvalue", "Pval", "pval", "p_value", "P-value", "P.value", "p-value", "p.value", "p"]
LOGP_ALIASES = ["LOG10P", "minus_log10_pval", "Plog10", "LP"]
N_ALIASES = ["N", "TotalSampleSize", "OBS_CT", "n", "samplesize", "SS"]
NCASE_ALIASES = ["NCASE", "N_CASE", "N_CASES", "NCASES", "Ncases", "n_cases", "N_CAS", "cases", "NC"]
IMPINFO_ALIASES = ["IMPINFO", "INFO", "info", "Info", "IMP_INFO", "imputation_info", "INFO_SCORE", "info_score", "Rsq", "SI"]
REF_ALIASES = ["REF", "ref"]
ALT_ALIASES = ["ALT", "alt"]
OR_ALIASES = ["OR", "odds_ratio", "hm_odds_ratio", "OddsRatio"]
CI_LOWER_ALIASES = ["OR_95L", "OR_L95", "L95", "LCI", "ci_lower", "hm_ci_lower", "CI_lower", "lower_ci"]
CI_UPPER_ALIASES = ["OR_95U", "OR_U95", "U95", "UCI", "ci_upper", "hm_ci_upper", "CI_upper", "upper_ci"]

SUPER_POPULATIONS = ("ALL", "AFR", "AMR", "EAS", "EUR", "SAS")

# Unambiguous labels only. "asian", "american", "middle east", "unknown" and blanks are rejected so
# the operator chooses the 1000 Genomes super-population explicitly.
ANCESTRY_SYNONYMS = {
    "all": "ALL", "combined": "ALL", "multi": "ALL", "multi-ancestry": "ALL", "trans-ancestry": "ALL",
    "cross-ancestry": "ALL", "mixed": "ALL",
    "afr": "AFR", "african": "AFR", "african-american": "AFR", "aa": "AFR", "black": "AFR",
    "amr": "AMR", "admixed-american": "AMR", "hispanic": "AMR", "latino": "AMR", "hispanic-latino": "AMR",
    "latin-american": "AMR",
    "eas": "EAS", "east-asian": "EAS",
    "eur": "EUR", "european": "EUR", "white": "EUR", "white-british": "EUR",
    "sas": "SAS", "south-asian": "SAS", "central-south-asian": "SAS",
}

# OpenGWAS dataset IDs: <batch>-<sub-batch>-<identifier>, e.g. ieu-b-5129, ukb-d-30780_irnt, ebi-a-GCST90000123.
OPENGWAS_ID_PATTERN = r"[a-z]+-[a-z]-[A-Za-z0-9_]+"


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


def opengwas_dataset_id(location: str) -> str | None:
    """Return the OpenGWAS dataset ID in a bare ID or dataset URL, preserving its case."""
    text = str(location).strip()
    match = re.search(rf"(?:^|/)({OPENGWAS_ID_PATTERN})(?:$|[/?#])", text)
    return match.group(1) if match else None


def infer_source_type(location: str) -> str:
    loc = str(location).strip()
    lower = loc.lower()
    if "ebi.ac.uk/gwas/studies/" in lower or re.fullmatch(r"GCST\d+", loc, flags=re.I):
        return "gwas_catalog"
    if "opengwas.io/datasets/" in lower or "gwas.mrcieu.ac.uk/datasets/" in lower:
        return "opengwas"
    if lower.startswith(("http://", "https://", "ftp://")):
        return "direct_url"
    if re.fullmatch(OPENGWAS_ID_PATTERN, loc) and not Path(loc).exists():
        return "opengwas"
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
        return opengwas_dataset_id(location) or "opengwas"
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
    }.get(source_type, source_type)
    unique_id = safe_path_component(source_id)
    if source_type not in {"gwas_catalog", "opengwas"}:
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
            clean["_line"] = i
            rows.append(clean)
    return rows


def read_state(path: str | Path) -> dict[str, dict[str, str]]:
    path = Path(path)
    if not path.exists():
        return {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["row_hash"]: row for row in reader if row.get("row_hash")}


def make_work_dir(prefix: str = "gwasprocessing-") -> Path:
    """Private directory for large intermediates: $GWAS_TMPDIR if set, else the system temporary directory
    (node-local scratch under SLURM). Keeps multi-GB files out of the output folder."""
    return Path(tempfile.mkdtemp(prefix=prefix, dir=os.environ.get("GWAS_TMPDIR") or None))


def atomic_write_text(path: str | Path, text: str) -> None:
    """Write text to a sibling temporary file, then rename it over `path` (never leaves a partial file)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with open(tmp, "w", newline="", encoding="utf-8") as handle:
            handle.write(text)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def write_state(path: str | Path, rows: Iterable[dict[str, str]]) -> None:
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=STATE_COLUMNS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({col: row.get(col, "") for col in STATE_COLUMNS})
    atomic_write_text(path, buffer.getvalue())


def state_row_complete(row: dict[str, str], report: bool) -> bool:
    required = ["standard_vcf", "standard_index", "metadata_json", "qc_json"]
    if report:
        required.append("report_html")
    for key in required:
        value = row.get(key, "")
        if not value or not Path(value).exists() or Path(value).stat().st_size == 0:
            return False
    return True


def write_json(path: str | Path, payload: dict[str, Any]) -> None:
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    atomic_write_text(path, text)


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
    repo_candidate = root.parent / path  # paths written relative to the repository root, e.g. manifests/x.yaml
    if repo_candidate.exists():
        return repo_candidate
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
    for named in (explicit_yaml, manifest_yaml):
        if named:  # a study YAML that is named must exist: silently ignoring it would drop its settings
            path = resolve_config_yaml_path(config_root, named)
            if not path.exists():
                raise config_error(f"Study YAML {named} for {study_id} not found (looked in {config_root}, its parent and the "
                                   "working directory)")
            return load_yaml(path)
    candidate = Path(config_root) / f"{study_id}.yaml"
    return load_yaml(candidate)


def manifest_config_overrides(record: dict[str, Any]) -> dict[str, Any]:
    """Manifest columns that override study-YAML settings for one row."""
    overrides: dict[str, Any] = {}
    if str(record.get("input_build", "") or "").strip():
        overrides["input_build"] = str(record["input_build"]).strip()
    delimiter = normalize_delimiter_value(record.get("delimiter"))
    if delimiter:
        overrides["delimiter"] = delimiter
    if str(record.get("population", "") or "").strip():
        overrides["reference_population"] = str(record["population"]).strip().upper()
    if str(record.get("format", "") or "").strip():
        overrides["format"] = str(record["format"]).strip()
    return overrides


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


def resolve_row_config(config_root: str | Path, record: dict[str, Any], cli_study_yaml: str = "") -> dict[str, Any]:
    """Study YAML (manifest study_yaml, else <study_id>.yaml) + manifest overrides + an optional CLI YAML."""
    config = resolve_study_config(config_root, record["study_id"], None, record.get("study_yaml") or None)
    config = merge_dicts(config, manifest_config_overrides(record))
    if cli_study_yaml:
        config = merge_dicts(config, resolve_study_config(config_root, record["study_id"], cli_study_yaml, None))
    return config


def without_appledouble(paths: Iterable[Path]) -> list[Path]:
    """Drop the ._ metadata files macOS writes beside each file on exFAT and network drives: Path.glob matches them."""
    return [path for path in paths if not path.name.startswith("._")]


def reference_identity(reference_root: str | Path) -> dict[str, list[int]]:
    """Size and mtime of every reference file that can change outputs (GRCh38 dbSNP index and merge map, 1000 Genomes
    frequencies, build-detection sample, chains, FASTA)."""
    root = Path(reference_root)
    patterns = ["dbsnp/GRCh38/*.json", "dbsnp/GRCh38/merged_*.u32", "1000genomes/GRCh38/*", "build_detection/*",
                "genome/*.chain*", "genome/*.fa", "genome/*.fasta"]
    identity: dict[str, list[int]] = {}
    for pattern in patterns:
        for path in without_appledouble(sorted(root.glob(pattern))):
            if path.is_file():
                stat = path.stat()
                identity[str(path.relative_to(root))] = [int(stat.st_size), int(stat.st_mtime)]
    return identity


def settings_fingerprint(record: dict[str, Any], config: dict[str, Any], reference_root: str | Path) -> str:
    """Hash of everything that determines a row's outputs; a mismatch means the outputs are stale."""
    payload = {
        "schema": OUTPUT_SCHEMA_VERSION,
        "row": {col: str(record.get(col, "") or "").strip() for col in MANIFEST_COLUMNS + MANIFEST_OPTIONAL_COLUMNS},
        "config": {key: value for key, value in config.items() if key != "_config_path"},
        "references": reference_identity(reference_root),
    }
    return short_hash(json.dumps(payload, sort_keys=True, default=str), 16)


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


def download_url(url: str, out_dir: str | Path, filename: str | None = None, retries: int = 3, timeout: int = 600,
                 backoff: float = 2.0) -> Path:
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
            if attempt < retries and backoff:
                time.sleep(min(30.0, backoff**attempt))
    if last_error:
        raise last_error
    return out


def sanitize_filename(filename: str) -> str:
    filename = Path(filename).name
    filename = re.sub(r"[^A-Za-z0-9._+-]+", "_", filename)
    return filename or "gwas"


def pick_existing_column(columns: Iterable[str], aliases: list[str]) -> str | None:
    columns_list = list(columns)
    lower_to_real = {c.lower(): c for c in columns_list}
    for alias in aliases:
        if alias in columns_list:
            return alias
        if alias.lower() in lower_to_real:
            return lower_to_real[alias.lower()]
    return None


def _header_line_and_metadata_count(path: str | Path) -> tuple[str, int]:
    """Return the first non-empty, non-'##' line (line ending removed) and how many '##' lines precede it."""
    opener = gzip.open if str(path).endswith(".gz") else open
    metadata_lines = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.rstrip("\r\n")
            if stripped.startswith("##"):
                metadata_lines += 1
                continue
            if stripped.strip():
                return stripped, metadata_lines
            metadata_lines += 1
    return "", metadata_lines


def count_metadata_lines(path: str | Path) -> int:
    """Number of leading '##' metadata (or blank) lines before a table's header row."""
    return _header_line_and_metadata_count(path)[1]


def infer_delimiter(path: str | Path) -> str:
    """Tab if the header has a tab, else comma, else runs of whitespace (space-aligned files)."""
    header, _ = _header_line_and_metadata_count(path)
    if "\t" in header:
        return "\t"
    if "," in header:
        return ","
    if re.search(r"\s", header.strip()):
        return r"\s+"
    return "\t"


def read_table(path: str | Path, delimiter: str | None = None, comment: str | None = None, **kwargs: Any) -> pd.DataFrame:
    """Read a summary-statistics table; extra keyword arguments (usecols, nrows, ...) go to pandas.read_csv."""
    path = Path(path)
    sep = delimiter or infer_delimiter(path)
    _, skip = _header_line_and_metadata_count(path)
    compression = "gzip" if str(path).endswith(".gz") else None
    return pd.read_csv(
        path, sep=sep, compression=compression, comment=comment, low_memory=False,
        skiprows=skip, na_values=["."], keep_default_na=True, **kwargs,
    )


def resolve_ancestry(value: Any) -> str | None:
    """Map a manifest ancestry label to a 1000 Genomes super-population code, or None if not unambiguous."""
    text = str(value if value is not None else "").strip()
    if not text:
        return None
    if text.upper() in SUPER_POPULATIONS:
        return text.upper()
    return ANCESTRY_SYNONYMS.get(sanitize(text, default="").replace("_", "-"))


def ancestry_to_population(ancestry: str) -> str:
    return resolve_ancestry(ancestry) or "ALL"


def normalize_build_value(value: Any) -> str:
    """Return 'hg19', 'hg38', or '' for labels such as 'HG19/GRCh37', 'GRCh38.p13', 'b37' or '38'."""
    compact = re.sub(r"[\s_.\-]", "", str(value or "").strip().lower())
    if not compact:
        return ""
    is38 = bool(re.search(r"(hg38|grch38|b38|build38|ncbi38)", compact)) or compact == "38"
    is37 = bool(re.search(r"(hg19|grch37|b37|build37|ncbi37)", compact)) or compact in {"37", "19"}
    if is38 and not is37:
        return "hg38"
    if is37 and not is38:
        return "hg19"
    return ""


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
        "25": "X",  # PLINK XY: pseudo-autosomal region of X
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
            if integer == 25:
                return "X"
            if integer == 26:
                return "MT"
    except ValueError:
        pass
    if text in {"X", "Y", "MT"}:
        return text
    return text


def numeric_series(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


PRIMARY_CHROMS = tuple([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
_COMPLEMENT_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def _map_unique(series: pd.Series, function: Any, missing: Any) -> pd.Series:
    """Apply `function` once per distinct non-missing value (fast for low-cardinality columns)."""
    values = pd.Series(series).astype(object)
    present = values.notna().to_numpy()
    out = np.full(len(values), missing, dtype=object)
    if present.any():
        observed = values[present]
        mapping = {value: function(value) for value in pd.unique(observed)}
        out[present] = observed.map(mapping).to_numpy(dtype=object)
    return pd.Series(out, index=values.index, dtype=object)


def upper_alleles(series: pd.Series) -> pd.Series:
    """Upper-cased, stripped allele strings with missing values as empty strings (object dtype)."""
    return _map_unique(series, lambda v: str(v).strip().upper(), "")


def complement_alleles(series: pd.Series) -> pd.Series:
    """Base-by-base complement of allele strings (A<->T, C<->G); other characters unchanged."""
    return _map_unique(upper_alleles(series), lambda v: v.translate(_COMPLEMENT_TABLE), "")


def normalize_chromosome_series(series: pd.Series) -> pd.Series:
    """normalize_chromosome applied once per distinct value; missing values become pd.NA."""
    return _map_unique(series, normalize_chromosome, pd.NA)


def is_palindromic(a: pd.Series, b: pd.Series) -> np.ndarray:
    a, b = upper_alleles(a), upper_alleles(b)
    return ((a + b).isin(["AT", "TA", "CG", "GC"])).to_numpy()


def alleles_match(a: pd.Series, b: pd.Series, ref_a: pd.Series, ref_b: pd.Series, allow_complement: bool = True) -> np.ndarray:
    """True where {a, b} equals {ref_a, ref_b}, optionally after complementing a and b."""
    a, b, ref_a, ref_b = (upper_alleles(s).reset_index(drop=True) for s in (a, b, ref_a, ref_b))
    valid = (a != "") & (b != "") & (ref_a != "") & (ref_b != "")
    same = ((a == ref_a) & (b == ref_b)) | ((a == ref_b) & (b == ref_a))
    if allow_complement:
        ca, cb = complement_alleles(a), complement_alleles(b)
        same = same | ((ca == ref_a) & (cb == ref_b)) | ((ca == ref_b) & (cb == ref_a))
    return (valid & same).to_numpy()


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
