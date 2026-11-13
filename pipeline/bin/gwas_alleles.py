#!/usr/bin/env python3
"""Allele operations, 1000 Genomes GRCh38 annotation and reference-genome alignment.

Comparisons with 1000 Genomes happen in GRCh38 orientation at (CHR, POS38). EAF filled from 1000 Genomes is the
reference frequency of the effect allele (not the minor allele) and is marked EAF_SOURCE=reference so it never
drives palindromic strand decisions.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gwas_coordinates import position_key
from gwas_pipeline_common import (
    alleles_match,
    complement_alleles,
    is_palindromic,
    normalize_chromosome_series,
    numeric_series,
    upper_alleles,
    without_appledouble,
)

# Palindromic (A/T, C/G) SNPs are resolved by allele frequency with TwoSampleMR's harmonise tolerance.
PALINDROMIC_TOLERANCE = 0.08
_ACGT = re.compile(r"^[ACGT]+$")


# ---------------------------------------------------------------------------------------------------
# 1000 Genomes annotation
# ---------------------------------------------------------------------------------------------------

def annotate_with_frequencies(frame: pd.DataFrame, freqs: Any) -> pd.DataFrame:
    """REF_A1 (1000 Genomes ALT), REF_A2 (REF) and REF_AF1 (ALT frequency in freqs.population) for rows matching a
    record at (CHR, POS38) with the same alleles on either strand."""
    out = frame.copy()
    n = len(out)
    ref_a1 = np.full(n, "", dtype=object)
    ref_a2 = np.full(n, "", dtype=object)
    ref_af1 = np.full(n, np.nan)
    if n:
        candidates = freqs.annotate(out["CHR"].to_numpy(dtype=object), out["POS38"].to_numpy(dtype="float64"))
        ea, oa = upper_alleles(out["EA"]).reset_index(drop=True), upper_alleles(out["OA"]).reset_index(drop=True)
        done = np.zeros(n, dtype=bool)
        for k in range(3):
            a1, a2 = candidates[f"A1_{k}"].reset_index(drop=True), candidates[f"A2_{k}"].reset_index(drop=True)
            match = ~done & alleles_match(ea, oa, a1, a2)
            ref_a1[match], ref_a2[match] = a1[match].to_numpy(), a2[match].to_numpy()
            ref_af1[match] = candidates[f"AF1_{k}"].to_numpy()[match]
            done |= match
    out["REF_A1"], out["REF_A2"], out["REF_AF1"] = ref_a1, ref_a2, ref_af1
    return out


def reference_af_of_effect_allele(frame: pd.DataFrame) -> np.ndarray:
    """Reference frequency of each row's effect allele (NaN when the alleles do not match the reference).

    Non-palindromic SNPs may match on the opposite strand; palindromic SNPs are read on the same strand.
    """
    if "REF_AF1" not in frame.columns:
        return np.full(len(frame), np.nan)
    ea37, oa37 = upper_alleles(frame["EA"]), upper_alleles(frame["OA"])
    ea37_arr = ea37.to_numpy()
    comp_ea = complement_alleles(ea37).to_numpy()
    a1 = upper_alleles(pd.Series(frame["REF_A1"])).to_numpy()
    a2 = upper_alleles(pd.Series(frame["REF_A2"])).to_numpy()
    af1 = numeric_series(pd.Series(frame["REF_AF1"])).to_numpy(dtype="float64")
    pal = is_palindromic(ea37, oa37)
    has = (a1 != "") & np.isfinite(af1)
    af = np.full(len(frame), np.nan)
    af = np.where(has & (ea37_arr == a1), af1, af)
    af = np.where(has & np.isnan(af) & (ea37_arr == a2), 1.0 - af1, af)
    af = np.where(has & np.isnan(af) & ~pal & (comp_ea == a1), af1, af)
    af = np.where(has & np.isnan(af) & ~pal & (comp_ea == a2), 1.0 - af1, af)
    return af


def eaf_source(frame: pd.DataFrame) -> np.ndarray:
    eaf = numeric_series(frame["EAF"]).to_numpy(dtype="float64")
    if "EAF_SOURCE" in frame.columns:
        source = frame["EAF_SOURCE"].fillna("").astype(str).to_numpy()
    else:
        source = np.full(len(frame), "", dtype=object)
    return np.where((source == "") & np.isfinite(eaf), "input", source).astype(object)


def fill_eaf_from_reference(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    out = frame.copy()
    out["EAF"] = numeric_series(out["EAF"]).astype("float64")
    out["EAF_SOURCE"] = eaf_source(out)
    ref_af = reference_af_of_effect_allele(out)
    fill = ~np.isfinite(out["EAF"].to_numpy()) & np.isfinite(ref_af)
    out.loc[fill, "EAF"] = ref_af[fill]
    out.loc[fill, "EAF_SOURCE"] = "reference"
    return out, int(fill.sum())


# ---------------------------------------------------------------------------------------------------
# Analysis QC
# ---------------------------------------------------------------------------------------------------

def resolve_palindromic_strand(frame: pd.DataFrame, qc: dict[str, Any], tolerance: float = PALINDROMIC_TOLERANCE) -> pd.DataFrame:
    """TwoSampleMR action=2 style strand inference for A/T and C/G SNPs, using study-reported EAF only.

    Adds PAL: '' for non-palindromic rows; 'forward' when study and reference frequencies of the effect allele
    fall on the same side of 0.5; 'complemented' (alleles complemented) when they fall on opposite sides;
    'ambiguous' when either frequency is within the tolerance of 0.5, there is no allele-matched reference,
    or the study reported no EAF. Nothing is removed: ambiguous SNPs are kept and marked (INFO/HM)."""
    out = frame.copy()
    out["PAL"] = ""
    pal = is_palindromic(out["EA"], out["OA"])
    detail = {"palindromic_total": int(pal.sum()), "kept_same_strand": 0, "flipped": 0, "ambiguous_frequency_near_half": 0,
              "ambiguous_no_reference": 0, "ambiguous_no_study_eaf": 0}
    if pal.any():
        low, high = 0.5 - tolerance, 0.5 + tolerance
        ref_af = reference_af_of_effect_allele(out)
        eaf = numeric_series(out["EAF"]).to_numpy(dtype="float64")
        study = (eaf_source(out) == "input") & np.isfinite(eaf)
        no_reference = pal & ~np.isfinite(ref_af)
        no_study = pal & np.isfinite(ref_af) & ~study
        usable = pal & np.isfinite(ref_af) & study
        near_half = usable & (((eaf > low) & (eaf < high)) | ((ref_af > low) & (ref_af < high)))
        resolvable = usable & ~near_half
        opposite = ((eaf < 0.5) & (ref_af > 0.5)) | ((eaf > 0.5) & (ref_af < 0.5))
        flip = resolvable & opposite
        if flip.any():
            out.loc[flip, "EA"] = complement_alleles(out.loc[flip, "EA"]).to_numpy()
            out.loc[flip, "OA"] = complement_alleles(out.loc[flip, "OA"]).to_numpy()
        out.loc[resolvable & ~flip, "PAL"] = "forward"
        out.loc[flip, "PAL"] = "complemented"
        out.loc[pal & ~resolvable, "PAL"] = "ambiguous"
        detail.update({"kept_same_strand": int((resolvable & ~flip).sum()), "flipped": int(flip.sum()),
                       "ambiguous_frequency_near_half": int(near_half.sum()), "ambiguous_no_reference": int(no_reference.sum()),
                       "ambiguous_no_study_eaf": int(no_study.sum())})
    qc["palindromic_resolution"] = detail
    return out


def deduplicate(frame: pd.DataFrame, action: np.ndarray | None = None) -> tuple[pd.DataFrame, int]:
    """Keep the first record per (CHR, POS38, allele pair from dedup_alleles); never keyed on rsID."""
    low, high = dedup_alleles(frame, action)
    key = pd.DataFrame({"KEY": position_key(frame["CHR"], frame["POS38"]), "LOW": low, "HIGH": high})
    duplicate = key.duplicated(keep="first").to_numpy()
    return frame[~duplicate], int(duplicate.sum())


ALIGNED_ACTIONS = ("keep", "swap", "complement", "complement_swap")


def dedup_alleles(frame: pd.DataFrame, action: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """The allele pair (LOW, HIGH) that deduplicate() keys on, with CHR and POS38.

    A row aligned to the GRCh38 reference (action from orient_to_reference) is keyed on (REF, ALT) in that order, so an
    insertion and a deletion at the same position stay apart (G>GCA and GCA>G). Any other row is keyed on the unordered
    pair, so the same variant reported either way round is still found."""
    ea, oa = upper_alleles(frame["EA"]).to_numpy(), upper_alleles(frame["OA"]).to_numpy()
    first = ea <= oa
    low, high = np.where(first, ea, oa), np.where(first, oa, ea)
    if action is not None:
        aligned = np.isin(action, ALIGNED_ACTIONS)
        low, high = np.where(aligned, oa, low), np.where(aligned, ea, high)
    return low, high


def af_moments(frame: pd.DataFrame) -> dict[str, float]:
    """Sufficient statistics for af_concordance over the rows of one frame (combine with combine_af_moments)."""
    ref_af = reference_af_of_effect_allele(frame)
    eaf = numeric_series(frame["EAF"]).to_numpy(dtype="float64")
    use = (eaf_source(frame) == "input") & np.isfinite(eaf) & np.isfinite(ref_af) & ~is_palindromic(frame["EA"], frame["OA"])
    x, y = eaf[use], ref_af[use]
    if not len(x):
        return {"n": 0, "mean_x": 0.0, "mean_y": 0.0, "m2_x": 0.0, "m2_y": 0.0, "c_xy": 0.0, "far": 0}
    mx, my = float(np.mean(x)), float(np.mean(y))
    return {"n": int(len(x)), "mean_x": mx, "mean_y": my, "m2_x": float(np.sum((x - mx) ** 2)),
            "m2_y": float(np.sum((y - my) ** 2)), "c_xy": float(np.sum((x - mx) * (y - my))),
            "far": int(np.sum(np.abs(x - y) > 0.2))}


def combine_af_moments(a: dict[str, float], b: dict[str, float]) -> dict[str, float]:
    """Pool two sets of moments (Chan et al. pairwise update), so pieces of a table give the whole-table result."""
    if not a["n"]:
        return dict(b)
    if not b["n"]:
        return dict(a)
    n = a["n"] + b["n"]
    dx, dy = b["mean_x"] - a["mean_x"], b["mean_y"] - a["mean_y"]
    w = a["n"] * b["n"] / n
    return {"n": n, "mean_x": a["mean_x"] + dx * b["n"] / n, "mean_y": a["mean_y"] + dy * b["n"] / n,
            "m2_x": a["m2_x"] + b["m2_x"] + dx * dx * w, "m2_y": a["m2_y"] + b["m2_y"] + dy * dy * w,
            "c_xy": a["c_xy"] + b["c_xy"] + dx * dy * w, "far": a["far"] + b["far"]}


def af_concordance_from_moments(m: dict[str, float]) -> dict[str, Any]:
    n = int(m["n"])
    result: dict[str, Any] = {"n": n, "pearson_r": None, "frac_abs_diff_gt_0.2": None}
    if n >= 3 and m["m2_x"] > 0 and m["m2_y"] > 0:
        result["pearson_r"] = float(max(-1.0, min(1.0, m["c_xy"] / np.sqrt(m["m2_x"] * m["m2_y"]))))
        result["frac_abs_diff_gt_0.2"] = float(m["far"] / n)
    return result


def af_concordance(frame: pd.DataFrame) -> dict[str, Any]:
    """Correlation of study EAF with the reference frequency of the same allele (non-palindromic SNPs)."""
    return af_concordance_from_moments(af_moments(frame))


# ---------------------------------------------------------------------------------------------------
# Reference genome (FASTA) access and alignment
# ---------------------------------------------------------------------------------------------------

def build_fai(fasta_path: str | Path) -> Path:
    """Write a samtools-compatible .fai index for an uncompressed FASTA."""
    fasta_path = Path(fasta_path)
    entries: list[str] = []
    name = None
    length = seq_offset = 0
    line_bases = line_width = 0
    offset = 0
    with open(fasta_path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    entries.append(f"{name}\t{length}\t{seq_offset}\t{line_bases}\t{line_width}")
                name = line[1:].split()[0].decode("ascii")
                seq_offset, length, line_bases, line_width = offset + len(line), 0, 0, 0
            else:
                stripped = line.rstrip(b"\r\n")
                if stripped and line_bases == 0:
                    line_bases, line_width = len(stripped), len(line)
                length += len(stripped)
            offset += len(line)
    if name is not None:
        entries.append(f"{name}\t{length}\t{seq_offset}\t{line_bases}\t{line_width}")
    fai = Path(str(fasta_path) + ".fai")
    tmp = fai.with_name(f".{fai.name}.tmp-{os.getpid()}")
    tmp.write_text("\n".join(entries) + "\n", encoding="utf-8")
    os.replace(tmp, fai)
    return fai


class Fasta:
    """Random access to an uncompressed GRCh38 FASTA via its .fai index (built on first use)."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        fai = Path(str(self.path) + ".fai")
        if not fai.exists() or fai.stat().st_mtime < self.path.stat().st_mtime:
            build_fai(self.path)
        self.index: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            fields = line.split("\t")
            if len(fields) >= 5:
                self.index[fields[0]] = (int(fields[1]), int(fields[2]), int(fields[3]), int(fields[4]))
        self._data: np.ndarray | None = None

    @classmethod
    def find(cls, reference_root: str | Path, build: str = "GRCh38") -> "Fasta | None":
        genome = Path(reference_root) / "genome"
        if not genome.is_dir():
            return None
        alias = {"GRCh38": "hg38", "GRCh37": "hg19"}[build]
        candidates = [genome / f"{build}.fa", genome / f"{build}.fasta"]
        for pattern in (f"*{build}*.fa", f"*{build}*.fasta", f"*{alias}*.fa", f"*{alias}*.fasta"):
            candidates.extend(without_appledouble(sorted(genome.glob(pattern))))
        for candidate in candidates:
            if candidate.is_file():
                return cls(candidate)
        return None

    def contig(self, chrom: str) -> str | None:
        names = [chrom, f"chr{chrom}"] + (["M", "chrM", "chrMT"] if chrom == "MT" else [])
        return next((name for name in names if name in self.index), None)

    def _bytes(self) -> np.ndarray:
        if self._data is None:
            self._data = np.memmap(self.path, dtype=np.uint8, mode="r")
        return self._data

    def bases(self, chrom: str, positions: Any, length: int = 1) -> np.ndarray:
        """Upper-case reference sequence of `length` bases starting at each 1-based position ('' if unavailable)."""
        positions = np.asarray(positions, dtype="int64")
        out = np.full(len(positions), "", dtype=object)
        name = self.contig(str(chrom))
        if name is None or len(positions) == 0:
            return out
        seq_len, offset, line_bases, line_width = self.index[name]
        valid = (positions >= 1) & (positions + length - 1 <= seq_len)
        if not valid.any():
            return out
        data = self._bytes()
        start = positions[valid] - 1
        columns = []
        for k in range(length):
            q = start + k
            columns.append(data[offset + (q // line_bases) * line_width + (q % line_bases)])
        matrix = np.stack(columns, axis=1).astype(np.uint8)
        strings = matrix.view(f"S{length}").reshape(-1).astype(str)
        out[valid] = np.char.upper(strings).astype(object)
        return out


def reference_alignment(chrom: pd.Series, pos: pd.Series, ea: pd.Series, oa: pd.Series, fasta: Fasta) -> pd.DataFrame:
    """Per-variant action needed to make REF (other allele) equal the GRCh38 reference sequence.

    Actions: keep, swap (effect allele is the reference base), complement (other allele's complement is the
    reference: reverse-strand report), complement_swap, drop (SNV whose alleles cannot match the reference),
    indel_unmatched (indel whose alleles do not match the reference sequence; removed, as REF must match the
    genome), unalignable (no reference sequence at the position, or an N in it; kept).
    """
    chroms = normalize_chromosome_series(pd.Series(chrom).reset_index(drop=True)).astype(str).to_numpy()
    positions = numeric_series(pd.Series(pos).reset_index(drop=True)).to_numpy(dtype="float64")
    ea_s = upper_alleles(pd.Series(ea).reset_index(drop=True))
    oa_s = upper_alleles(pd.Series(oa).reset_index(drop=True))
    ea_arr, oa_arr = ea_s.to_numpy(), oa_s.to_numpy()
    n = len(ea_arr)
    action = np.full(n, "unalignable", dtype=object)
    acgt_ea, acgt_oa = (ea_s.str.fullmatch(r"[ACGT]+").to_numpy(dtype=bool), oa_s.str.fullmatch(r"[ACGT]+").to_numpy(dtype=bool)) if n else (np.zeros(0, bool), np.zeros(0, bool))
    acgt = acgt_ea & acgt_oa
    acgt &= np.isfinite(positions)
    lengths_ea = ea_s.str.len().to_numpy(dtype=int) if n else np.zeros(0, int)
    lengths_oa = oa_s.str.len().to_numpy(dtype=int) if n else np.zeros(0, int)
    snv = acgt & (lengths_ea == 1) & (lengths_oa == 1)
    indel = acgt & ~snv
    comp_ea, comp_oa = complement_alleles(ea_s).to_numpy(), complement_alleles(oa_s).to_numpy()
    pal = is_palindromic(ea_s, oa_s)
    for c in np.unique(chroms[acgt]) if acgt.any() else []:
        rows = np.where((chroms == c) & snv)[0]
        if len(rows):
            ref = fasta.bases(c, positions[rows].astype("int64"))
            has = np.isin(ref, ["A", "C", "G", "T"])  # an N (or no sequence) leaves the variant unverified
            r_oa, r_ea = oa_arr[rows] == ref, ea_arr[rows] == ref
            r_coa = ~pal[rows] & (comp_oa[rows] == ref)
            r_cea = ~pal[rows] & (comp_ea[rows] == ref)
            decided = np.where(r_oa, "keep", np.where(r_ea, "swap", np.where(r_coa, "complement",
                               np.where(r_cea, "complement_swap", "drop"))))
            action[rows] = np.where(has, decided, "unalignable")
        rows = np.where((chroms == c) & indel)[0]
        if len(rows):
            seq_oa = np.full(len(rows), "", dtype=object)
            seq_ea = np.full(len(rows), "", dtype=object)
            for length_value in np.unique(lengths_oa[rows]):
                sel = lengths_oa[rows] == length_value
                seq_oa[sel] = fasta.bases(c, positions[rows][sel].astype("int64"), int(length_value))
            for length_value in np.unique(lengths_ea[rows]):
                sel = lengths_ea[rows] == length_value
                seq_ea[sel] = fasta.bases(c, positions[rows][sel].astype("int64"), int(length_value))
            matched_oa = seq_oa == oa_arr[rows]
            matched_ea = ~matched_oa & (seq_ea == ea_arr[rows])
            # a window that runs off the sequence or holds an N cannot be checked: kept as unverified
            unknown = ~(pd.Series(seq_oa, dtype=object).str.fullmatch(r"[ACGT]+").fillna(False).to_numpy(dtype=bool)
                        & pd.Series(seq_ea, dtype=object).str.fullmatch(r"[ACGT]+").fillna(False).to_numpy(dtype=bool))
            action[rows] = np.where(matched_oa, "keep", np.where(matched_ea, "swap",
                                                                 np.where(unknown, "unalignable", "indel_unmatched")))
    return pd.DataFrame({"action": action})


REMOVED_ACTIONS = ("drop", "indel_unmatched")  # alignment actions that remove the row


def normalise_indels(frame: pd.DataFrame, fasta: Fasta, action: np.ndarray) -> tuple[pd.DataFrame, int]:
    """Left-align and trim indels whose REF (OA) matches the GRCh38 sequence (Tan et al. 2015, as bcftools norm):
    while the alleles end in the same base, drop it (re-anchoring on the base to the left when an allele empties),
    then remove shared leading bases while both alleles are longer than 1. Returns the frame and how many moved."""
    out = frame.copy()
    ea, oa = upper_alleles(out["EA"]).to_numpy(), upper_alleles(out["OA"]).to_numpy()
    rows = np.where(np.isin(action, ["keep", "swap"]) & ((pd.Series(ea).str.len() > 1) | (pd.Series(oa).str.len() > 1)).to_numpy())[0]
    if not len(rows):
        return out, 0
    chroms = normalize_chromosome_series(out["CHR"]).astype(str).to_numpy()
    positions = numeric_series(out["POS38"]).to_numpy(dtype="float64")
    moved = 0
    new_pos, new_ea, new_oa = positions.copy(), ea.copy(), oa.copy()
    for i in rows:
        chrom, pos, ref, alt = chroms[i], int(positions[i]), oa[i], ea[i]
        start = (pos, ref, alt)
        for _ in range(1000):
            if ref and alt and ref[-1] == alt[-1]:
                ref, alt = ref[:-1], alt[:-1]
            elif not ref or not alt:
                base = fasta.bases(chrom, [pos - 1])[0] if pos > 1 else ""
                if not base:
                    break
                ref, alt, pos = base + ref, base + alt, pos - 1
            else:
                break
        while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
            ref, alt, pos = ref[1:], alt[1:], pos + 1
        if ref and alt and (pos, ref, alt) != start:
            new_pos[i], new_oa[i], new_ea[i] = pos, ref, alt
            moved += 1
    out["POS38"], out["EA"], out["OA"] = new_pos, new_ea, new_oa
    return out, moved


def alignment_summary(action: np.ndarray, fasta: Fasta) -> dict[str, Any]:
    counts = {str(k): int(v) for k, v in zip(*np.unique(action, return_counts=True))} if len(action) else {}
    return alignment_summary_from_counts(counts, fasta)


def alignment_summary_from_counts(counts: dict[str, int], fasta: Fasta) -> dict[str, Any]:
    return {
        "performed": True, "fasta": str(fasta.path), "kept": counts.get("keep", 0),
        "swapped": counts.get("swap", 0) + counts.get("complement_swap", 0),
        "strand_flipped": counts.get("complement", 0) + counts.get("complement_swap", 0),
        "dropped_mismatch": counts.get("drop", 0), "indels_unmatched_dropped": counts.get("indel_unmatched", 0),
        "unalignable_kept": counts.get("unalignable", 0),
    }


def align_table_to_reference(frame: pd.DataFrame, fasta: Fasta) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Orient each variant so OA (written as REF) is the GRCh38 base; flip BETA and EAF when alleles swap."""
    out, action = orient_to_reference(frame, fasta)
    return out[~np.isin(action, REMOVED_ACTIONS)], alignment_summary(action, fasta)


def orient_to_reference(frame: pd.DataFrame, fasta: Fasta) -> tuple[pd.DataFrame, np.ndarray]:
    """align_table_to_reference without dropping: every row, oriented, with its action ('drop' rows unchanged)."""
    action = reference_alignment(frame["CHR"], frame["POS38"], frame["EA"], frame["OA"], fasta)["action"].to_numpy()
    out = frame.copy()
    complement = np.isin(action, ["complement", "complement_swap"])
    if complement.any():
        out.loc[complement, "EA"] = complement_alleles(out.loc[complement, "EA"]).to_numpy()
        out.loc[complement, "OA"] = complement_alleles(out.loc[complement, "OA"]).to_numpy()
    swap = np.isin(action, ["swap", "complement_swap"])
    if swap.any():
        ea_before = out.loc[swap, "EA"].to_numpy(copy=True)  # a view under pandas 3 would see the next line
        out.loc[swap, "EA"] = out.loc[swap, "OA"].to_numpy()
        out.loc[swap, "OA"] = ea_before
        out.loc[swap, "BETA"] = -numeric_series(out.loc[swap, "BETA"])
        out.loc[swap, "EAF"] = 1.0 - numeric_series(out.loc[swap, "EAF"])
        if "EZ" in out.columns:
            out.loc[swap, "EZ"] = -numeric_series(out.loc[swap, "EZ"])
    return out, action
