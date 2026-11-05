#!/usr/bin/env python3
"""Reference lookup, allele operations and reference-genome alignment.

All comparisons with the 1000 Genomes lookup happen in GRCh37 orientation: alleles of a variant that
lifted onto the minus strand are complemented first. EAF filled from the lookup is the reference
frequency of the effect allele (not the minor allele) and is marked EAF_SOURCE=reference so it never
drives palindromic strand decisions.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from gwas_coordinates import position_key, rs_number
from gwas_pipeline_common import (
    alleles_match,
    complement_alleles,
    is_palindromic,
    normalize_chromosome_series,
    numeric_series,
    upper_alleles,
    without_appledouble,
)

LOOKUP_COLUMNS = ["CHR", "SNP", "POS19", "A1", "A2", "AF1"]  # columns of the lookup files on disk
PREPARED_COLUMNS = ["CHR", "RS", "POS19", "A1", "A2", "AF1", "KEY"]  # compact in-memory form
LOOKUP_CHUNK = 2_000_000
REFERENCE_ANNOTATION = ["REF_SNP", "REF_A1", "REF_A2", "REF_AF1"]

# Analysis-QC defaults. The stage is opt-in (default off); a study YAML `analysis_qc` block can toggle
# filters and override thresholds. `palindromic_tolerance` matches TwoSampleMR's harmonise tolerance.
ANALYSIS_QC_DEFAULTS = {
    "biallelic_snp": True,
    "remove_palindromic": True,
    "maf": True,
    "maf_threshold": 0.01,
    "deduplicate": True,
    "resolve_palindromic_strand": True,
    "palindromic_tolerance": 0.08,
}
_ACGT = re.compile(r"^[ACGT]+$")


# ---------------------------------------------------------------------------------------------------
# Reference lookup
# ---------------------------------------------------------------------------------------------------

def prepare_lookup(table: pd.DataFrame) -> pd.DataFrame:
    """Compact, normalised lookup: CHR (category), RS (rsID number, -1 if none), POS19, A1/A2 (category),
    AF1 and KEY (chromosome/position key). Accepts files' SNP column or an existing RS column; idempotent."""
    if list(table.columns) == PREPARED_COLUMNS:
        return table
    n = len(table)
    get = lambda column: table[column] if column in table.columns else pd.Series(np.nan, index=table.index)  # noqa: E731
    chrom = normalize_chromosome_series(get("CHR"))
    pos19 = numeric_series(get("POS19")).astype("float64")
    rs = table["RS"].to_numpy(dtype="int64") if "RS" in table.columns else rs_number(get("SNP").astype("string").str.strip())
    out = pd.DataFrame({
        "CHR": pd.Categorical(chrom.to_numpy(dtype=object)), "RS": np.asarray(rs, dtype="int64"), "POS19": pos19.to_numpy(),
        "A1": pd.Categorical(upper_alleles(get("A1")).to_numpy(dtype=object)),
        "A2": pd.Categorical(upper_alleles(get("A2")).to_numpy(dtype=object)),
        "AF1": numeric_series(get("AF1")).astype("float64").to_numpy(),
        "KEY": position_key(chrom, pos19) if n else np.array([], dtype="int64"),
    })
    return out


def read_reference_lookup(reference_root: str | Path, population: str) -> tuple[pd.DataFrame, str]:
    """Load <root>/1000genomes/phase3/lookup/<population>.tsv(.gz), falling back to ALL.

    Returns (prepared lookup, population actually used; "" when nothing was found). The file is read in
    chunks into the compact form of prepare_lookup; legacy columns (SNPID, EAF_<pop>, POS38) are ignored.
    """
    root = Path(reference_root) / "1000genomes" / "phase3" / "lookup"
    for pop in dict.fromkeys([str(population or "ALL").upper(), "ALL"]):
        for suffix in (".tsv.gz", ".tsv"):
            path = root / f"{pop}{suffix}"
            if not path.exists():
                continue
            header = pd.read_csv(path, sep="\t", nrows=0).columns
            usecols = [c for c in LOOKUP_COLUMNS if c in header]
            chunks = [prepare_lookup(chunk) for chunk in pd.read_csv(
                path, sep="\t", usecols=usecols, low_memory=False, chunksize=LOOKUP_CHUNK,
                dtype={c: str for c in ("CHR", "SNP", "A1", "A2") if c in usecols})]
            table = pd.concat(chunks, ignore_index=True) if chunks else prepare_lookup(pd.DataFrame(columns=usecols))
            for column in ("CHR", "A1", "A2"):
                table[column] = pd.Categorical(table[column].astype(object))
            return table, pop
    return prepare_lookup(pd.DataFrame(columns=LOOKUP_COLUMNS)), ""


def lookup_has_alleles(lookup: pd.DataFrame) -> bool:
    return (not lookup.empty) and bool((upper_alleles(lookup["A1"]) != "").any()) and bool(np.isfinite(lookup["AF1"]).any())


def grch37_alleles(frame: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Effect/other alleles in GRCh37 orientation (complemented where the liftover strand was '-')."""
    ea, oa = upper_alleles(frame["EA"]), upper_alleles(frame["OA"])
    if "LIFT_STRAND" not in frame.columns:
        return ea, oa
    minus = frame["LIFT_STRAND"].astype(str).eq("-").to_numpy()
    ea37 = pd.Series(np.where(minus, complement_alleles(ea).to_numpy(), ea.to_numpy()), index=frame.index, dtype=object)
    oa37 = pd.Series(np.where(minus, complement_alleles(oa).to_numpy(), oa.to_numpy()), index=frame.index, dtype=object)
    return ea37, oa37


def annotate_with_reference(frame: pd.DataFrame, lookup: pd.DataFrame) -> pd.DataFrame:
    """Add REF_SNP/REF_A1/REF_A2/REF_AF1 for rows matching the lookup on (CHR, POS19) and alleles."""
    out = frame.copy()
    n = len(out)
    ref_snp = np.full(n, None, dtype=object)
    ref_a1 = np.full(n, "", dtype=object)
    ref_a2 = np.full(n, "", dtype=object)
    ref_af1 = np.full(n, np.nan)
    if n and lookup_has_alleles(lookup):
        ref = prepare_lookup(lookup)
        ea37, oa37 = grch37_alleles(out)
        probe = pd.DataFrame({"_row": np.arange(n), "KEY": position_key(out["CHR"], out["POS19"]),
                              "A": ea37.to_numpy(), "B": oa37.to_numpy()})
        probe = probe[probe["KEY"] >= 0]
        candidates = ref.loc[ref["KEY"].isin(probe["KEY"]).to_numpy(), ["KEY", "RS", "A1", "A2", "AF1"]]
        merged = probe.merge(candidates, on="KEY", how="inner")
        if not merged.empty:
            merged = merged[alleles_match(merged["A"], merged["B"], merged["A1"].astype(object), merged["A2"].astype(object))]
            merged = merged.drop_duplicates("_row", keep="first")
            rows = merged["_row"].to_numpy()
            rs = merged["RS"].to_numpy(dtype="int64")
            ref_snp[rows] = np.where(rs >= 0, np.char.add("rs", rs.astype(str)).astype(object), None)
            ref_a1[rows] = merged["A1"].astype(object).to_numpy()
            ref_a2[rows] = merged["A2"].astype(object).to_numpy()
            ref_af1[rows] = merged["AF1"].to_numpy(dtype="float64")
    out["REF_SNP"], out["REF_A1"], out["REF_A2"], out["REF_AF1"] = ref_snp, ref_a1, ref_a2, ref_af1
    return out


def fill_missing_rsids(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    out = frame.copy()
    snp = out["SNP"].astype("string").str.strip()
    missing = (snp.isna() | ~snp.str.lower().str.startswith("rs", na=False)).to_numpy()
    available = pd.Series(out["REF_SNP"]).notna().to_numpy()
    fill = missing & available
    out["SNP"] = out["SNP"].astype(object)
    out.loc[fill, "SNP"] = out.loc[fill, "REF_SNP"]
    return out, int(fill.sum())


def reference_af_of_effect_allele(frame: pd.DataFrame) -> np.ndarray:
    """Reference frequency of each row's effect allele (NaN when the alleles do not match the reference).

    Non-palindromic SNPs may match on the opposite strand; palindromic SNPs are read on the same strand.
    """
    if "REF_AF1" not in frame.columns:
        return np.full(len(frame), np.nan)
    ea37, oa37 = grch37_alleles(frame)
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

def resolve_palindromic_strand(frame: pd.DataFrame, settings: dict[str, Any], qc: dict[str, Any]) -> tuple[pd.DataFrame, int]:
    """TwoSampleMR action=2 style strand inference for A/T and C/G SNPs, using study-reported EAF only.

    Kept when study and reference frequencies of the effect allele fall on the same side of 0.5,
    complemented when they fall on opposite sides, dropped when either is within the tolerance of 0.5,
    when there is no allele-matched reference, or when the EAF was not reported by the study.
    """
    out = frame.copy()
    pal = is_palindromic(out["EA"], out["OA"])
    detail = {"palindromic_total": int(pal.sum()), "kept_same_strand": 0, "flipped": 0, "dropped_ambiguous": 0,
              "dropped_no_reference": 0, "dropped_no_study_eaf": 0, "dropped_resolution_disabled": 0}
    if not pal.any():
        qc["palindromic_resolution"] = detail
        return out, 0
    if not settings.get("resolve_palindromic_strand", True):
        detail["dropped_resolution_disabled"] = int(pal.sum())
        qc["palindromic_resolution"] = detail
        return out[~pal], int(pal.sum())

    tolerance = float(settings.get("palindromic_tolerance", 0.08))
    low, high = 0.5 - tolerance, 0.5 + tolerance
    ref_af = reference_af_of_effect_allele(out)
    eaf = numeric_series(out["EAF"]).to_numpy(dtype="float64")
    study = (eaf_source(out) == "input") & np.isfinite(eaf)
    no_reference = pal & ~np.isfinite(ref_af)
    no_study = pal & np.isfinite(ref_af) & ~study
    usable = pal & np.isfinite(ref_af) & study
    ambiguous = usable & (((eaf > low) & (eaf < high)) | ((ref_af > low) & (ref_af < high)))
    resolvable = usable & ~ambiguous
    opposite = ((eaf < 0.5) & (ref_af > 0.5)) | ((eaf > 0.5) & (ref_af < 0.5))
    flip = resolvable & opposite
    if flip.any():
        out.loc[flip, "EA"] = complement_alleles(out.loc[flip, "EA"]).to_numpy()
        out.loc[flip, "OA"] = complement_alleles(out.loc[flip, "OA"]).to_numpy()
    drop = pal & ~resolvable
    detail.update({"kept_same_strand": int((resolvable & ~flip).sum()), "flipped": int(flip.sum()),
                   "dropped_ambiguous": int(ambiguous.sum()), "dropped_no_reference": int(no_reference.sum()),
                   "dropped_no_study_eaf": int(no_study.sum())})
    qc["palindromic_resolution"] = detail
    return out[~drop], int(drop.sum())


def deduplicate(frame: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Keep the first record per (CHR, POS38, unordered allele pair); never keyed on rsID."""
    ea, oa = upper_alleles(frame["EA"]).to_numpy(), upper_alleles(frame["OA"]).to_numpy()
    first = ea <= oa
    key = pd.DataFrame({"KEY": position_key(frame["CHR"], frame["POS38"]), "LOW": np.where(first, ea, oa),
                        "HIGH": np.where(first, oa, ea)})
    duplicate = key.duplicated(keep="first").to_numpy()
    return frame[~duplicate], int(duplicate.sum())


def apply_analysis_qc(frame: pd.DataFrame, config_settings: dict[str, Any] | None, qc: dict[str, Any]) -> pd.DataFrame:
    settings = {**ANALYSIS_QC_DEFAULTS, **(config_settings or {})}
    removed: dict[str, int] = {}
    out = frame
    if settings.get("biallelic_snp", True):
        ea, oa = upper_alleles(out["EA"]), upper_alleles(out["OA"])
        keep = (ea.isin(list("ACGT")) & oa.isin(list("ACGT")) & (ea != oa)).to_numpy()
        removed["non_biallelic"] = int((~keep).sum())
        out = out[keep]
    if settings.get("remove_palindromic", True):
        out, dropped = resolve_palindromic_strand(out, settings, qc)
        removed["palindromic"] = dropped
    if settings.get("maf", True):
        threshold = float(settings.get("maf_threshold", 0.01))
        eaf = numeric_series(out["EAF"]).to_numpy(dtype="float64")
        maf = np.minimum(eaf, 1.0 - eaf)
        keep = np.isfinite(eaf) & (eaf >= 0) & (eaf <= 1) & (maf >= threshold)
        removed[f"maf_below_{threshold}"] = int((~keep).sum())
        out = out[keep]
    if settings.get("deduplicate", True):
        out, removed["duplicate"] = deduplicate(out)
    qc["qc_removed"] = removed
    qc["analysis_qc_settings"] = settings
    return out


def af_concordance(frame: pd.DataFrame) -> dict[str, Any]:
    """Correlation of study EAF with the reference frequency of the same allele (non-palindromic SNPs)."""
    ref_af = reference_af_of_effect_allele(frame)
    eaf = numeric_series(frame["EAF"]).to_numpy(dtype="float64")
    use = (eaf_source(frame) == "input") & np.isfinite(eaf) & np.isfinite(ref_af) & ~is_palindromic(frame["EA"], frame["OA"])
    n = int(use.sum())
    result: dict[str, Any] = {"n": n, "pearson_r": None, "frac_abs_diff_gt_0.2": None}
    if n >= 3 and np.std(eaf[use]) > 0 and np.std(ref_af[use]) > 0:
        result["pearson_r"] = float(np.corrcoef(eaf[use], ref_af[use])[0, 1])
        result["frac_abs_diff_gt_0.2"] = float(np.mean(np.abs(eaf[use] - ref_af[use]) > 0.2))
    return result


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
    def find(cls, reference_root: str | Path) -> "Fasta | None":
        genome = Path(reference_root) / "genome"
        if not genome.is_dir():
            return None
        candidates = [genome / "GRCh38.fa", genome / "GRCh38.fasta"]
        for pattern in ("*GRCh38*.fa", "*GRCh38*.fasta", "*hg38*.fa", "*hg38*.fasta"):
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
    indel_unmatched (kept), unalignable (non-ACGT alleles or no reference sequence; kept).
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
            has = ref != ""
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
            action[rows] = np.where(matched_oa, "keep", np.where(matched_ea, "swap", "indel_unmatched"))
    return pd.DataFrame({"action": action})


def alignment_summary(action: np.ndarray, fasta: Fasta) -> dict[str, Any]:
    counts = {str(k): int(v) for k, v in zip(*np.unique(action, return_counts=True))} if len(action) else {}
    return {
        "performed": True, "fasta": str(fasta.path), "kept": counts.get("keep", 0),
        "swapped": counts.get("swap", 0) + counts.get("complement_swap", 0),
        "strand_flipped": counts.get("complement", 0) + counts.get("complement_swap", 0),
        "dropped_mismatch": counts.get("drop", 0), "indels_unmatched_kept": counts.get("indel_unmatched", 0),
        "unalignable_kept": counts.get("unalignable", 0),
    }


def align_table_to_reference(frame: pd.DataFrame, fasta: Fasta) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Orient each variant so OA (written as REF) is the GRCh38 base; flip BETA and EAF when alleles swap."""
    action = reference_alignment(frame["CHR"], frame["POS38"], frame["EA"], frame["OA"], fasta)["action"].to_numpy()
    out = frame.copy()
    complement = np.isin(action, ["complement", "complement_swap"])
    if complement.any():
        out.loc[complement, "EA"] = complement_alleles(out.loc[complement, "EA"]).to_numpy()
        out.loc[complement, "OA"] = complement_alleles(out.loc[complement, "OA"]).to_numpy()
    swap = np.isin(action, ["swap", "complement_swap"])
    if swap.any():
        ea_before = out.loc[swap, "EA"].to_numpy()
        out.loc[swap, "EA"] = out.loc[swap, "OA"].to_numpy()
        out.loc[swap, "OA"] = ea_before
        out.loc[swap, "BETA"] = -numeric_series(out.loc[swap, "BETA"])
        out.loc[swap, "EAF"] = 1.0 - numeric_series(out.loc[swap, "EAF"])
    return out[action != "drop"], alignment_summary(action, fasta)
