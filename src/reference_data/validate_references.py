#!/usr/bin/env python3
"""Check that the reference data GWASprocessing needs is present and complete.

Usage: python3 src/reference_data/validate_references.py [--reference-root references]

Required (exit 1 when missing): the GRCh38 and GRCh37 FASTAs with their .fai indexes, the hg19ToHg38 chain, the dbSNP
GRCh38 index, the 1000 Genomes GRCh38 frequencies for AFR, AMR, EAS, EUR, SAS and ALL, and the build-detection sample.
Recommended (reported as WARNING): the dbSNP merge map (without it, rsIDs that dbSNP has merged into another are not
placed) and the hg38ToHg19 chain (needed only to rebuild the build-detection sample).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

POPULATIONS = ["AFR", "AMR", "EAS", "EUR", "SAS", "ALL"]
FREQUENCY_FILES = ("key.i64", "ref.u1", "alt.u1", "af.f4")
DETECTION_FILES = [f"{build}.{kind}" for build in ("grch37", "grch38") for kind in ("key.i64", "ref.u1", "alt.u1")]
BUILD_SCRIPT = "src/reference_data/build_grch38_references.sh"


def nonempty(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def check(root: Path) -> tuple[list[str], list[str], list[str]]:
    """(errors, warnings, notes)."""
    errors: list[str] = []
    warnings: list[str] = []
    notes: list[str] = []
    genome = root / "genome"
    for name in ("GRCh38.fa", "GRCh37.fa"):
        if not nonempty(genome / name):
            errors.append(f"missing {genome / name} (run download_fasta_and_chains.sh)")
        elif not nonempty(genome / f"{name}.fai"):
            errors.append(f"missing {genome / name}.fai (samtools faidx, or rerun download_fasta_and_chains.sh)")
    if not nonempty(genome / "hg19ToHg38.over.chain.gz"):
        errors.append(f"missing {genome / 'hg19ToHg38.over.chain.gz'} (run download_fasta_and_chains.sh)")
    if not nonempty(genome / "hg38ToHg19.over.chain.gz"):
        warnings.append("hg38ToHg19 chain missing: needed only to rebuild the build-detection sample")
    dbsnp = root / "dbsnp" / "GRCh38"
    if not (dbsnp / "meta.json").exists() or not nonempty(dbsnp / "rs.u32"):
        errors.append(f"missing dbSNP GRCh38 index in {dbsnp} (run {BUILD_SCRIPT})")
    else:
        meta = json.loads((dbsnp / "meta.json").read_text(encoding="utf-8"))
        notes.append(f"dbSNP build {meta.get('dbsnp_build', '?')}: {meta.get('rsids', 0):,} rsIDs")
        if (dbsnp / "merged_meta.json").exists() and nonempty(dbsnp / "merged_from.u32"):
            merged = json.loads((dbsnp / "merged_meta.json").read_text(encoding="utf-8"))
            notes.append(f"dbSNP merge map: {merged.get('retired_rsids', 0):,} retired rsIDs")
        else:
            warnings.append(f"dbSNP merge map missing in {dbsnp}: rsIDs merged into another rsID will not be placed "
                            f"(run build_dbsnp_merges.py, or {BUILD_SCRIPT})")
    freqs = root / "1000genomes" / "GRCh38"
    for pop in POPULATIONS:
        missing = [kind for kind in FREQUENCY_FILES if not nonempty(freqs / f"{pop}.{kind}")]
        if missing:
            errors.append(f"missing 1000 Genomes GRCh38 {pop} frequencies ({', '.join(missing)}) in {freqs} (run {BUILD_SCRIPT})")
    detection = root / "build_detection"
    missing = [name for name in DETECTION_FILES if not nonempty(detection / name)]
    if missing or not (detection / "meta.json").exists():
        errors.append(f"missing build-detection sample in {detection} (run build_build_detection_sample.py)")
    return errors, warnings, notes


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--reference-root", default="references")
    args = parser.parse_args()
    errors, warnings, notes = check(Path(args.reference_root))
    for line in notes:
        print(f"OK: {line}")
    for line in warnings:
        print(f"WARNING: {line}")
    for line in errors:
        print(f"ERROR: {line}", file=sys.stderr)
    if errors:
        sys.exit(1)
    print("Reference data complete.")


if __name__ == "__main__":
    main()
