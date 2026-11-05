#!/usr/bin/env python3
"""Association statistics helpers: -log10 P (LP) derivation and genomic inflation (lambda GC).

LP is carried natively so P-values beyond double precision (P < ~1e-308, common for cis-pQTLs)
are never lost: a LOG10P column is used as-is, P converts directly, and a P that underflowed to 0
is recovered from the Wald z-score |BETA/SE|.
"""

from __future__ import annotations

import math
from statistics import NormalDist
from typing import Any

import numpy as np

_SQRT2 = math.sqrt(2.0)
_LN10 = math.log(10.0)
_ERFC_LIMIT = 37.0  # two-sided P for |z| below this is representable as a double
_CHI2_MEDIAN_1DF = 0.454936423119572  # qchisq(0.5, df = 1)
_erfc = np.vectorize(math.erfc, otypes=[float])


def lp_from_z(z: Any) -> np.ndarray:
    """Two-sided -log10 P for a standard-normal statistic; exact below |z|=37, asymptotic above."""
    z = np.abs(np.asarray(z, dtype="float64"))
    out = np.full(z.shape, np.nan)
    finite = np.isfinite(z)
    small = finite & (z < _ERFC_LIMIT)
    if small.any():
        out[small] = -np.log10(_erfc(z[small] / _SQRT2))
    large = finite & ~small
    if large.any():
        x = z[large] / _SQRT2
        x2 = x * x
        # erfc(x) ~ exp(-x^2) / (x sqrt(pi)) * (1 - 1/(2x^2) + 3/(4x^4) - 15/(8x^6) + 105/(16x^8))
        series = 1.0 - 1.0 / (2.0 * x2) + 3.0 / (4.0 * x2**2) - 15.0 / (8.0 * x2**3) + 105.0 / (16.0 * x2**4)
        ln_p = -x2 - np.log(x * math.sqrt(math.pi)) + np.log(series)
        out[large] = -ln_p / _LN10
    return out


def _as_array(values: Any, n: int) -> np.ndarray:
    if values is None:
        return np.full(n, np.nan)
    return np.asarray(values, dtype="float64").reshape(-1)


def compute_lp(p: Any = None, log10p: Any = None, beta: Any = None, se: Any = None) -> tuple[np.ndarray, dict[str, Any]]:
    """Return (LP array, provenance counts).

    Order per row: a non-negative LOG10P value (a column that is mostly negative is read as log10 P and
    negated); else -log10 P for 0 < P <= 1; else, when P is 0 or absent, -log10 of the two-sided Wald P
    from BETA/SE (SE > 0). P outside [0, 1] stays missing.
    """
    lengths = [len(np.asarray(v).reshape(-1)) for v in (p, log10p, beta, se) if v is not None]
    n = lengths[0] if lengths else 0
    p_arr, l_arr, b_arr, s_arr = (_as_array(v, n) for v in (p, log10p, beta, se))
    lp = np.full(n, np.nan)
    info: dict[str, Any] = {"from_log10p": 0, "from_p": 0, "from_z": 0, "invalid": 0, "log10p_negated": False}

    finite_l = np.isfinite(l_arr)
    if finite_l.any() and float(np.median(l_arr[finite_l])) < 0:
        l_arr = -l_arr
        info["log10p_negated"] = True
    use_l = finite_l & (l_arr >= 0)
    lp[use_l] = l_arr[use_l]
    info["from_log10p"] = int(use_l.sum())

    need = np.isnan(lp)
    use_p = need & np.isfinite(p_arr) & (p_arr > 0) & (p_arr <= 1)
    lp[use_p] = -np.log10(p_arr[use_p])
    info["from_p"] = int(use_p.sum())

    zero_or_absent = np.isnan(lp) & ((p_arr == 0) | np.isnan(p_arr))
    valid_z = zero_or_absent & np.isfinite(b_arr) & np.isfinite(s_arr) & (s_arr > 0)
    if valid_z.any():
        lp[valid_z] = lp_from_z(b_arr[valid_z] / s_arr[valid_z])
    info["from_z"] = int(valid_z.sum())
    info["invalid"] = int(np.isnan(lp).sum())
    return lp, info


def lambda_gc(lp: Any) -> float | None:
    """Genomic inflation factor: median chi-square (1 df) implied by the LP values over its null median."""
    values = np.asarray(lp, dtype="float64").reshape(-1)
    values = values[np.isfinite(values) & (values >= 0)]
    if values.size == 0:
        return None
    p_median = 10.0 ** (-float(np.median(values)))
    if p_median <= 0.0:
        return math.inf
    if p_median >= 1.0:
        return 0.0
    z = -NormalDist().inv_cdf(p_median / 2.0)
    return z * z / _CHI2_MEDIAN_1DF
