#!/usr/bin/env python3
"""GWAS-VCF writing: header contract, body formatting, numeric chromosome ordering, bgzip + tabix.

Every output declares GRCh38 (##reference plus ##contig lines with GRCh38 lengths), carries the
trait metadata lines downstream consumers parse, and is bgzip-compressed with a real tabix index.
"""

from __future__ import annotations

import os
import re
import shutil
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
    make_work_dir,
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
METADATA_KEYS = ("study_id", "phenotype", "ancestry", "sex", "PMID")
RESERVED_INFO = {"POS19", "SNPID", "EAF_SOURCE"}
_RANK = {chrom: rank for rank, chrom in enumerate(PRIMARY_CHROMS, start=1)}

FORMAT_HEADER = [
    '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">',
    '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">',
    '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">',
    '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency">',
    '##FORMAT=<ID=SS,Number=A,Type=Float,Description="Sample size used to estimate genetic effect">',
]
POS19_INFO = '##INFO=<ID=POS19,Number=1,Type=Integer,Description="GRCh37 position of the variant">'


def chrom_rank(chrom: Any) -> int:
    return _RANK.get(str(normalize_chromosome(chrom)), 99)


def sanitize_header_value(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value if value is not None else "")).strip()


def sample_name(record: dict[str, Any]) -> str:
    return re.sub(r"\s+", "_", sanitize_header_value(record.get("phenotype"))) or "GWAS"


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


def metadata_header_lines(record: dict[str, Any]) -> list[str]:
    return [f"##{key}={sanitize_header_value(record.get(key)) or 'NA'}" for key in METADATA_KEYS]


def contig_header_lines() -> list[str]:
    return [f"##contig=<ID={chrom},length={length},assembly=GRCh38>" for chrom, length in GRCH38_LENGTHS.items()]


def assembly_header_lines() -> list[str]:
    return ["##reference=GRCh38", *contig_header_lines()]


def table_header(record: dict[str, Any], extras: Iterable[str]) -> list[str]:
    extras_info = [f'##INFO=<ID={info_key(name)},Number=1,Type=String,Description="Source column {sanitize_header_value(name)}">'
                   for name in extras]
    return [
        "##fileformat=VCFv4.2",
        f"##fileDate={datetime.now(timezone.utc).strftime('%Y%m%d')}",
        "##source=GWASprocessing",
        f"##GWASprocessing_version={OUTPUT_SCHEMA_VERSION}",
        "##GWASVCF=1.0",
        *assembly_header_lines(),
        *metadata_header_lines(record),
        POS19_INFO,
        '##INFO=<ID=SNPID,Number=1,Type=String,Description="Source variant identifier">',
        '##INFO=<ID=EAF_SOURCE,Number=1,Type=String,Description="Origin of FORMAT/AF: input (reported by the study) '
        'or reference (1000 Genomes frequency of the effect allele)">',
        *extras_info,
        *FORMAT_HEADER,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name(record),
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


def table_body_lines(df: pd.DataFrame, extras: Iterable[str]) -> list[str]:
    snp = df["SNP"].astype("string").str.strip()
    ids = np.where(snp.str.match(r"(?i)^rs\d+$").fillna(False).to_numpy(dtype=bool), snp.str.lower().fillna(".").to_numpy(dtype=object), ".")
    pos19 = _fmt_int(df["POS19"])
    snpid = [info_value(v) for v in df["SNPID"]] if "SNPID" in df.columns else ["."] * len(df)
    source = df["EAF_SOURCE"].fillna("").astype(str).to_numpy() if "EAF_SOURCE" in df.columns else np.full(len(df), "")
    source = np.where(np.isin(source, ["input", "reference"]), source, ".")
    info = np.array([f"POS19={p};SNPID={s};EAF_SOURCE={e}" for p, s, e in zip(pos19, snpid, source)], dtype=object)
    for name in extras:
        values = [info_value(v) for v in df[name]] if name in df.columns else ["."] * len(df)
        info = info + np.array([f";{info_key(name)}={v}" for v in values], dtype=object)
    sample = (_fmt_float(df["BETA"]) + ":" + _fmt_float(df["SE"]) + ":" + _fmt_float(df["LP"]) + ":"
              + _fmt_float(df["EAF"]) + ":" + _fmt_int(df["N"]))
    chrom = normalize_chromosome_series(df["CHR"]).astype(str).to_numpy(dtype=object)
    pos = _fmt_int(df["POS38"])
    lines = (chrom + "\t" + pos + "\t" + ids + "\t" + _alleles(df["OA"]) + "\t" + _alleles(df["EA"])
             + "\t.\tPASS\t" + info + "\tES:SE:LP:AF:SS\t" + sample)
    return list(lines)


def sort_frame(df: pd.DataFrame) -> pd.DataFrame:
    rank = normalize_chromosome_series(df["CHR"]).map(lambda c: _RANK.get(str(c), 99)).to_numpy(dtype="int64")
    order = np.lexsort((numeric_series(df["POS38"]).to_numpy(), rank))
    return df.iloc[order]


def require_htslib() -> None:
    if not (command_exists("bgzip") and command_exists("tabix")):
        raise config_error(
            "bgzip and tabix (htslib) are required to write indexed GWAS-VCF files but were not found on PATH. "
            "Install htslib (e.g. `brew install htslib`, `conda install -c bioconda htslib`, or `module load htslib`) "
            "or run `bash src/source.sh --tools-only`."
        )


def bgzip_and_index(plain_path: Path, vcf_path: Path) -> None:
    """Compress `plain_path` to `vcf_path` (atomically) and create a tabix index; removes the plain file."""
    require_htslib()
    tmp = vcf_path.with_name(f".{vcf_path.name}.tmp-{os.getpid()}")
    try:
        with open(tmp, "wb") as handle:
            subprocess.run(["bgzip", "-c", str(plain_path)], stdout=handle, check=True)
        os.replace(tmp, vcf_path)
    finally:
        tmp.unlink(missing_ok=True)
    index = Path(str(vcf_path) + ".tbi")
    index.unlink(missing_ok=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", str(vcf_path)], check=True)
    if not index.exists():
        raise RuntimeError(f"tabix did not create {index}")
    plain_path.unlink(missing_ok=True)


def write_table_vcf(df: pd.DataFrame, vcf_path: Path, record: dict[str, Any], extras: Iterable[str]) -> None:
    extras = list(extras)
    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    work = make_work_dir()
    try:
        plain = work / "table.vcf"
        ordered = sort_frame(df)
        with open(plain, "w", encoding="utf-8") as handle:
            handle.write("\n".join(table_header(record, extras)) + "\n")
            chunk = 500_000
            for start in range(0, len(ordered), chunk):
                lines = table_body_lines(ordered.iloc[start:start + chunk], extras)
                if lines:
                    handle.write("\n".join(lines) + "\n")
        bgzip_and_index(plain, vcf_path)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def sort_vcf_body(plain_path: Path) -> None:
    """Sort a plain VCF body by (numeric chromosome rank, position) in place, keeping the header first."""
    plain_path = Path(plain_path)
    header_path = plain_path.with_name(plain_path.name + ".hdr")
    body_path = plain_path.with_name(plain_path.name + ".body")
    sorted_path = plain_path.with_name(plain_path.name + ".sorted")
    try:
        with open(plain_path, "r", encoding="utf-8") as src, open(header_path, "w", encoding="utf-8") as head, \
                open(body_path, "w", encoding="utf-8") as body:
            for line in src:
                if line.startswith("#"):
                    head.write(line)
                elif line.strip():
                    chrom = line.split("\t", 1)[0]
                    body.write(f"{chrom_rank(chrom)}\t{line}")
        env = {**os.environ, "LC_ALL": "C"}
        with open(sorted_path, "w", encoding="utf-8") as out:
            subprocess.run(["sort", "-t", "\t", "-k1,1n", "-k3,3n", "-T", str(plain_path.parent), str(body_path)],
                           check=True, stdout=out, env=env)
        with open(plain_path, "w", encoding="utf-8") as out, open(header_path, "r", encoding="utf-8") as head, \
                open(sorted_path, "r", encoding="utf-8") as body:
            out.writelines(head)
            for line in body:
                out.write(line.split("\t", 1)[1])
    finally:
        for path in (header_path, body_path, sorted_path):
            path.unlink(missing_ok=True)
