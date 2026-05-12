#!/usr/bin/env python3
"""Audit a GWAS-VCF and export the guide-style flat GWAS columns."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


DEFAULT_VCF = Path("GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783/standardised.gwas.vcf.gz")
DEFAULT_OUTPUT = Path("GWAS/locke_2015_25673413/standardise/GWAScatalog_GCST002783/standardised.flat.python.tsv")
OUTPUT_COLUMNS = ["CHR", "POS", "SNP", "REF", "ALT", "EAF", "BETA", "SE", "P", "N"]
REQUIRED_FORMAT = ["ES", "SE", "LP", "SS"]


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def parse_float(value: str) -> float | None:
    if value in {"", ".", "NA", "NaN", "nan"}:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def parse_integer(value: str) -> int | None:
    number = parse_float(value)
    if number is None:
        return None
    rounded = math.floor(number + 0.5)
    if abs(number - rounded) > 1e-8:
        raise ValueError(f"Sample size is not an integer after VCF standardisation: {value}")
    return int(rounded)


def format_value(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not math.isfinite(value):
            return "NA"
        return f"{value:.8g}"
    return str(value)


def read_flat_records(vcf_path: Path) -> tuple[list[str], list[dict[str, object]]]:
    samples: list[str] = []
    records: list[dict[str, object]] = []
    header_seen = False

    with open_text(vcf_path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                if len(fields) < 10:
                    raise ValueError("VCF header does not include a sample column")
                samples = fields[9:]
                header_seen = True
                continue
            if not header_seen:
                raise ValueError("VCF data encountered before #CHROM header")

            fields = line.split("\t")
            if len(fields) < 10:
                raise ValueError(f"VCF record has fewer than 10 columns: {line[:120]}")
            chrom, pos, snp, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            format_keys = fields[8].split(":")
            missing_format = [key for key in REQUIRED_FORMAT if key not in format_keys]
            if missing_format:
                raise ValueError(f"VCF record missing FORMAT keys: {', '.join(missing_format)}")
            values = fields[9].split(":")
            format_map = dict(zip(format_keys, values, strict=False))
            lp = parse_float(format_map.get("LP", "."))
            records.append(
                {
                    "CHR": chrom,
                    "POS": int(pos),
                    "SNP": snp,
                    "REF": ref,
                    "ALT": alt,
                    "EAF": parse_float(format_map.get("AF", ".")),
                    "BETA": parse_float(format_map.get("ES", ".")),
                    "SE": parse_float(format_map.get("SE", ".")),
                    "P": None if lp is None else 10 ** (-lp),
                    "N": parse_integer(format_map.get("SS", ".")),
                }
            )

    if not header_seen:
        raise ValueError("VCF #CHROM header not found")
    return samples, records


def validate(records: list[dict[str, object]]) -> None:
    if not records:
        raise ValueError("No VCF records were parsed")
    for column in OUTPUT_COLUMNS:
        if column not in records[0]:
            raise ValueError(f"Flat records missing column: {column}")
    required = ["CHR", "POS", "REF", "ALT", "BETA", "SE", "P"]
    for index, record in enumerate(records, start=1):
        missing = [column for column in required if record.get(column) in {None, "", "."}]
        if missing:
            raise ValueError(f"Record {index} has missing required values: {', '.join(missing)}")


def write_flat(output_path: Path, records: list[dict[str, object]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow({column: format_value(record[column]) for column in OUTPUT_COLUMNS})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", nargs="?", type=Path, default=DEFAULT_VCF)
    parser.add_argument("output", nargs="?", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if not args.vcf.exists():
        raise SystemExit(f"VCF not found: {args.vcf}")

    samples, records = read_flat_records(args.vcf)
    validate(records)
    write_flat(args.output, records)
    print(f"VCF rows: {len(records)}")
    print(f"Samples: {','.join(samples)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
