"""Step log: what each stage of standardisation did to the rows, how many were removed and why.

build_steps(qc) turns the pooled QC record into an ordered list of steps. Each step has rows_in, removed,
rows_out, the reasons rows were removed (disjoint, so they add up) and the changes made without removing
rows. The last rows_out must equal qc["output_rows"]; if it does not, the log says so and a warning is added.
format_steps() renders the list as plain text for the task log. Counts only: no variant-level content.
"""

from __future__ import annotations

from typing import Any

CAUSE_TEXT = {
    "rsid_not_in_coordinate_index": "no position in the file and the rsID is not in the coordinate index",
    "no_position_or_rsid": "no position in the file and no rsID to look one up",
    "liftover_to_grch38_failed": "GRCh37 position could not be lifted to GRCh38",
    "chromosome_not_primary": "chromosome is not 1-22, X, Y or MT",
    "invalid_position": "position is zero or negative",
    "missing_allele": "effect or other allele missing",
    "invalid_allele": "allele is not made of A, C, G and T (e.g. I/D codes or '-'), or both alleles are the same",
    "missing_beta": "effect size missing and not derivable",
    "missing_or_invalid_se": "standard error missing, zero or negative and not derivable",
    "missing_or_invalid_p": "P value missing or invalid and not derivable",
}
COORDINATE_CAUSES = ("rsid_not_in_coordinate_index", "no_position_or_rsid", "liftover_to_grch38_failed",
                     "chromosome_not_primary", "invalid_position")
STATISTIC_CAUSES = ("missing_allele", "invalid_allele", "missing_beta", "missing_or_invalid_se", "missing_or_invalid_p")
LIFTOVER_TEXT = {
    "unmapped": "position absent from the chain file",
    "other_chrom": "maps to a different chromosome",
    "multiple": "maps to more than 1 position",
    "minus_strand_non_snv_dropped": "maps to the minus strand and is not a SNV, so the alleles cannot be complemented",
}
RARE_THRESHOLD = 0.01


def _n(value: Any) -> int:
    try:
        return int(value or 0)
    except (TypeError, ValueError):
        return 0


def _step(steps: list[dict[str, Any]], rows: int | None, title: str, why: str, removed_by: dict[str, int] | None = None,
          changed: dict[str, int] | None = None, notes: list[str] | None = None) -> int | None:
    removed_by = {k: v for k, v in (removed_by or {}).items() if v}
    removed = sum(removed_by.values())
    rows_out = None if rows is None else rows - removed
    steps.append({"step": len(steps) + 1, "title": title, "why": why, "rows_in": rows, "removed": removed,
                  "rows_out": rows_out, "removed_by": removed_by, "changed": {k: v for k, v in (changed or {}).items() if v},
                  "notes": notes or []})
    return rows_out


def _rare_note(detail: dict[str, Any]) -> str:
    rare, no_eaf = _n(detail.get("study_maf_below_threshold")), _n(detail.get("no_study_eaf"))
    return f"of which study MAF < {RARE_THRESHOLD:g}: {rare:,}; no study EAF: {no_eaf:,}"


def build_steps(qc: dict[str, Any]) -> list[dict[str, Any]]:
    steps: list[dict[str, Any]] = []
    rows: int | None = _n(qc["input_rows"]) if "input_rows" in qc else None
    vcf = qc.get("input_format") == "gwas_vcf"
    ids = qc.get("variant_id_normalisation", {}) or {}
    p_source = qc.get("p_value_source", {}) or {}
    derivations = qc.get("derivations", {}) or {}
    rows = _step(steps, rows, "Read input",
                 "Records read from the source GWAS-VCF, 1 row per ALT allele." if vcf else "Rows read from the source table.",
                 None,
                 {"rsIDs extracted from compound identifiers": _n(ids.get("rsid_extracted")),
                  "LP read directly (-log10 P)": _n(p_source.get("from_log10p")),
                  "LP computed from BETA/SE because P was 0 or missing": _n(p_source.get("from_z")),
                  "effect sizes from odds ratios": _n(derivations.get("beta_from_or")),
                  "standard errors from confidence intervals": _n(derivations.get("se_from_ci")),
                  "standard errors from P and effect size": _n(derivations.get("se_from_p"))})
    exclusions = qc.get("exclusions") or []
    if exclusions:
        rows = _step(steps, rows, "Configured exclusions", "Variants listed in exclusion files named in the study YAML.",
                     {f"listed in {e.get('file')}": _n(e.get("removed")) for e in exclusions})
    joins = qc.get("joins") or []
    if any("rows_before" in j for j in joins):
        unmatched = sum(_n(j.get("rows_before")) - _n(j.get("rows_after")) for j in joins)
        rows = _step(steps, rows, "Configured joins", "Rows without a match in an inner join named in the study YAML.",
                     {"no match in a joined file": unmatched})

    causes = qc.get("incomplete_by_cause", {}) or {}
    fills = qc.get("reference_fills", {}) or {}
    lift = (qc.get("liftover", {}) or {}).get("hg19_to_hg38", {}) or {}
    lift_detail = ", ".join(f"{LIFTOVER_TEXT.get(k, k)}: {_n(v):,}" for k, v in lift.items()
                            if k not in ("ok", "minus_strand") and _n(v))
    notes = [f"{CAUSE_TEXT[c]}: {_rare_note(causes[c])}" for c in COORDINATE_CAUSES if c in causes]
    if lift_detail and "liftover_to_grch38_failed" in causes:
        notes.append(f"liftover failures: {lift_detail}")
    build = qc.get("build_detection", {}) or {}
    if build.get("build"):
        notes.append(f"input build: {build.get('build')} ({build.get('method', '')})")
    rows = _step(steps, rows, "GRCh38 coordinates",
                 "Every variant needs a GRCh38 chromosome and position: taken from the file, parsed from chr:pos identifiers, "
                 "or looked up by rsID, then lifted from GRCh37 where needed.",
                 {CAUSE_TEXT[c]: _n((causes.get(c) or {}).get("rows")) for c in COORDINATE_CAUSES},
                 {"positions parsed from chr:pos identifiers": _n(qc.get("positions_from_variant_ids")),
                  "position looked up by rsID": _n(fills.get("POS19")),
                  "lifted from GRCh37 to GRCh38": _n(lift.get("ok")),
                  "alleles complemented after a minus-strand liftover": _n(lift.get("minus_strand"))}, notes)
    rows = _step(steps, rows, "Alleles and statistics",
                 "Every variant needs 2 different alleles made of A, C, G and T, an effect size, a positive standard error "
                 "and a P value.",
                 {CAUSE_TEXT[c]: _n((causes.get(c) or {}).get("rows")) for c in STATISTIC_CAUSES}, None,
                 [f"{CAUSE_TEXT[c]}: {_rare_note(causes[c])}" for c in STATISTIC_CAUSES if c in causes])
    rows = _step(steps, rows, "1000 Genomes annotation",
                 "Allele frequencies are matched from the 1000 Genomes lookup; where the study reports none, FORMAT/AF takes "
                 "the 1000 Genomes frequency and FORMAT/AFS says so. Nothing is removed here.",
                 None, {"AF filled from 1000 Genomes": _n(fills.get("EAF")), "rsIDs filled from 1000 Genomes": _n(fills.get("rsid"))})
    pal = qc.get("palindromic_resolution", {}) or {}
    rows = _step(steps, rows, "Palindromic (A/T, C/G) SNPs",
                 "The strand of an A/T or C/G SNP cannot be read from its alleles, so it is inferred by comparing the study "
                 "and 1000 Genomes allele frequencies (tolerance 0.08). SNPs where that is not possible are kept as reported "
                 "and marked HM=palindrome_ambiguous. Nothing is removed here.",
                 None,
                 {"kept, same strand (palindrome_forward)": _n(pal.get("kept_same_strand")),
                  "complemented to the reference strand (palindrome_complemented)": _n(pal.get("flipped")),
                  "ambiguous: a frequency within 0.42-0.58": _n(pal.get("ambiguous_frequency_near_half")),
                  "ambiguous: no allele-matched 1000 Genomes frequency": _n(pal.get("ambiguous_no_reference")),
                  "ambiguous: no study-reported frequency": _n(pal.get("ambiguous_no_study_eaf"))},
                 [f"palindromic SNPs: {_n(pal.get('palindromic_total')):,}"])
    rows = _step(steps, rows, "Duplicates",
                 "Records with the same GRCh38 chromosome, position and allele pair; the first in input order is kept.",
                 {"duplicate of an earlier record": _n(qc.get("duplicates_removed"))})
    alignment = qc.get("reference_alignment", {}) or {}
    if alignment.get("performed"):
        rows = _step(steps, rows, "GRCh38 reference alignment",
                     "REF must match the GRCh38 genome: alleles are swapped or complemented to match it, and SNVs that match on "
                     "neither strand are removed.",
                     {"alleles match neither strand of the GRCh38 base": _n(alignment.get("dropped_mismatch"))},
                     {"already aligned": _n(alignment.get("kept")), "effect and other allele swapped": _n(alignment.get("swapped")),
                      "complemented to the forward strand": _n(alignment.get("strand_flipped")),
                      "indels kept without a REF check": _n(alignment.get("indels_unmatched_kept")),
                      "kept, reference base unavailable": _n(alignment.get("unalignable_kept"))})
    output = _n(qc.get("output_rows"))
    if rows is not None and rows != output:
        steps[-1]["notes"].append(f"UNRECONCILED: the steps leave {rows:,} rows but {output:,} were written")
        qc.setdefault("warnings", []).append(f"step log does not reconcile: {rows} rows expected, {output} written")
    return steps


def format_steps(qc: dict[str, Any], record: dict[str, Any]) -> str:
    """Plain-text summary for the task log (the HTML report renders the same steps in full)."""
    steps = qc.get("steps") or build_steps(qc)
    width = 12
    lines = [f"GWASprocessing step log: {record.get('study_id', '')} / {record.get('output_id', '')}",
             f"Source: {record.get('raw_path', '')}",
             f"Format: {qc.get('input_format', '')}; reference population {qc.get('reference_population_used', '')}",
             f"Started {qc.get('started_at', '')}, completed {qc.get('completed_at', '')}", ""]
    pieces = qc.get("table_pieces")
    if pieces:
        lines.insert(3, f"Pieces: {pieces.get('pieces')} of {pieces.get('rows_per_piece')} rows, {pieces.get('workers')} workers")
    lines.append(f"{'Step':<44}{'Rows in':>{width}}{'Removed':>{width}}{'Rows out':>{width}}")
    fmt = lambda v: "" if v is None else f"{v:,}"  # noqa: E731
    for s in steps:
        lines.append(f"{str(s['step']) + '. ' + s['title']:<44}{fmt(s['rows_in']):>{width}}{fmt(s['removed']):>{width}}"
                     f"{fmt(s['rows_out']):>{width}}")
    lines += ["", f"Written: {_n(qc.get('output_rows')):,} of {_n(qc.get('input_rows')):,} input rows", ""]
    for s in steps:
        lines.append(f"{s['step']}. {s['title']}")
        lines.append(f"   {s['why']}")
        for reason, count in s["removed_by"].items():
            lines.append(f"   removed {count:>12,}  {reason}")
        for change, count in s["changed"].items():
            lines.append(f"   changed {count:>12,}  {change}")
        for note in s["notes"]:
            lines.append(f"   note: {note}")
        lines.append("")
    timings = qc.get("timings_s")
    if timings:
        lines.append("Time (seconds)")
        for k, v in (timings.get("wall") or {}).items():
            lines.append(f"   {k:<40}{v:>10.1f}")
        lines.append("   per-piece steps, summed over pieces (run in parallel with more than 1 worker):")
        for k, v in (timings.get("piece_steps_summed_over_pieces") or {}).items():
            lines.append(f"      {k:<37}{v:>10.1f}")
        lines.append("")
    if qc.get("warnings"):
        lines.append("Warnings")
        lines += [f"   {w}" for w in qc["warnings"]]
    return "\n".join(lines).rstrip() + "\n"
