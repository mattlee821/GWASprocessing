#!/usr/bin/env python3
"""Coordinates: chr:pos variant-ID parsing, rsID numbers, position keys and the build-detection sample size.

Liftover is BCFtools/liftover (standardise_gwas.lift_with_bcftools); build detection matches the paired GRCh37/GRCh38
sample (gwas_reference.build_hits).
"""

from __future__ import annotations

import re
from typing import Any

import numpy as np
import pandas as pd

from gwas_pipeline_common import PRIMARY_CHROMS, normalize_chromosome, normalize_chromosome_series, numeric_series

CHROM_CODES = {chrom: code for code, chrom in enumerate(PRIMARY_CHROMS, start=1)}
KEY_FACTOR = 10_000_000_000  # position keys: chromosome code * 1e10 + position (fits int64)
DETECTION_SAMPLE = 20_000  # rows sampled for build detection, with a fixed seed
DETECTION_SEED = 821


def position_key(chrom: Any, pos: Any) -> np.ndarray:
    """int64 key per (chromosome, position); -1 where either is missing or not a primary chromosome."""
    codes = normalize_chromosome_series(pd.Series(chrom).reset_index(drop=True)).map(CHROM_CODES).astype("float64").to_numpy()
    positions = numeric_series(pd.Series(pos).reset_index(drop=True)).to_numpy(dtype="float64")
    valid = np.isfinite(codes) & np.isfinite(positions) & (positions > 0)
    key = np.full(len(positions), -1, dtype="int64")
    key[valid] = codes[valid].astype("int64") * KEY_FACTOR + positions[valid].astype("int64")
    return key


_VARIANT_ID = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|X|Y|MT|M|XY)[:_](?P<pos>[0-9]+)(?:[:_](?P<a1>[ACGTN]+)[:_](?P<a2>[ACGTN]+))?$",
    flags=re.I,
)


def parse_variant_ids(ids: pd.Series) -> pd.DataFrame:
    """Parse chr:pos[:a1:a2] / chr_pos_a1_a2 identifiers into CHR and POS (missing where unparseable)."""
    text = pd.Series(ids).astype("string").str.strip()
    parts = text.str.extract(_VARIANT_ID)
    chrom = parts["chrom"].map(lambda v: normalize_chromosome(v) if isinstance(v, str) else pd.NA)
    return pd.DataFrame({"CHR": chrom.astype(object), "POS": numeric_series(parts["pos"]).astype("float64")}, index=text.index)


def rs_number(ids: pd.Series) -> np.ndarray:
    """Numeric part of bare rsIDs (rs123 -> 123); -1 for anything else."""
    digits = pd.Series(ids).astype("string").str.strip().str.extract(r"^[rR][sS](\d+)$")[0]
    return numeric_series(digits).fillna(-1).astype("int64").to_numpy()
