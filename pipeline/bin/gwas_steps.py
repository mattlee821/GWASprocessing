"""Step log: what each stage of standardisation did to the rows, how many were removed and why.

build_steps(qc) turns the pooled QC record into an ordered list of steps. Each step has rows_in, removed,
rows_out, the reasons rows were removed (disjoint, so they add up) and the changes made without removing
rows. The last rows_out must equal qc["output_rows"]; if it does not, the log says so and a warning is added.
format_steps() renders the list as the plain-text log written beside the output (<prefix>.steps.txt).
Counts only: no variant-level content is written.
"""

from __future__ import annotations

from typing import Any

CAUSE_TEXT = {
    "rsid_not_in_coordinate_index": "no position in the file and the rsID is not in the 1000 Genomes coordinate index",
    "no_position_or_rsid": "no position in the file and no rsID to look one up",
    "liftover_to_grch38_failed": "GRCh37 position could not be lifted to GRCh38",
    "chromosome_not_primary": "chromosome is not 1-22, X, Y or MT",
    "invalid_position": "position is zero or negative",
    "missing_allele": "effect or other allele missing",
    "missing_beta": "effect size missing and not derivable",
    "missing_or_invalid_se": "standard error missing, zero or negative and not derivable",
    "missing_or_invalid_p": "P value missing or invalid and not derivable",
}
LIFTOVER_TEXT = {
    "unmapped": "position absent from the chain file",
    "other_chrom": "maps to a different chromosome",
    "multiple": "maps to more than 1 position",
    "minus_strand_non_snv_dropped": "maps to the minus strand and is not a SNV, so the alleles cannot be complemented",
}
PALINDROMIC_TEXT = {
    "dropped_no_reference": "no allele-matched frequency in the 1000 Genomes lookup to infer the strand",
    "dropped_ambiguous": "study or reference frequency within the tolerance of 0.5, so the strand cannot be inferred",
    "dropped_no_study_eaf": "no study-reported allele frequency to infer the strand",
    "dropped_resolution_disabled": "strand resolution switched off in the settings",
}


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


def _rare_note(detail: dict[str, Any], threshold: float) -> str:
    rare, no_eaf = _n(detail.get("study_maf_below_threshold")), _n(detail.get("no_study_eaf"))
    return f"of which study MAF < {threshold:g}: {rare:,}; no study EAF: {no_eaf:,}"


def _analysis_steps(steps: list[dict[str, Any]], rows: int | None, qc: dict[str, Any]) -> int | None:
    settings = qc.get("analysis_qc_settings", {}) or {}
    removed = qc.get("qc_removed", {}) or {}
    threshold = float(settings.get("maf_threshold", 0.01))
    if settings.get("biallelic_snp", True):
        rows = _step(steps, rows, "Biallelic SNPs only",
                     "Analysis QC keeps single-base A/C/G/T alleles that differ; indels and other alleles are removed.",
                     {"not a biallelic SNP": _n(removed.get("non_biallelic"))})
    if settings.get("maf", True):
        maf_key = next((k for k in removed if k.startswith("maf_below_")), "")
        detail = qc.get("maf_detail") or {}
        by = {f"MAF < {threshold:g}": _n(detail.get("below_threshold")),
              "no valid allele frequency (study or 1000 Genomes)": _n(detail.get("no_valid_frequency"))}
        if not detail:
            by = {f"MAF < {threshold:g} or no frequency": _n(removed.get(maf_key))}
        rows = _step(steps, rows, f"Minor allele frequency >= {threshold:g}",
                     "Analysis QC removes rare variants, using the study EAF, or the 1000 Genomes frequency where the study gives none.",
                     by)
    if settings.get("remove_palindromic", True):
        pal = qc.get("palindromic_resolution", {}) or {}
        by = {PALINDROMIC_TEXT[k]: _n(pal.get(k)) for k in PALINDROMIC_TEXT}
        rows = _step(steps, rows, "Palindromic (A/T, C/G) SNPs",
                     "For A/T and C/G SNPs the strand cannot be read from the alleles, so it is inferred by comparing the study "
                     "and 1000 Genomes allele frequencies; SNPs where that is not possible are removed.",
                     by, {"palindromic SNPs kept, same strand": _n(pal.get("kept_same_strand")),
                          "palindromic SNPs kept, complemented to the reference strand": _n(pal.get("flipped"))},
                     [f"palindromic SNPs reaching this step: {_n(pal.get('palindromic_total')):,}"])
    if settings.get("deduplicate", True):
        rows = _step(steps, rows, "Duplicates",
                     "Records with the same GRCh38 chromosome, position and allele pair; the first in input order is kept.",
                     {"duplicate of an earlier record": _n(removed.get("duplicate"))})
    return rows


def build_steps(qc: dict[str, Any]) -> list[dict[str, Any]]:
    steps: list[dict[str, Any]] = []
    rows: int | None = _n(qc["input_rows"]) if "input_rows" in qc else None
    threshold = float((qc.get("analysis_qc_settings") or {}).get("maf_threshold", 0.01))
    vcf = qc.get("input_format") == "gwas_vcf"
    ids = qc.get("variant_id_normalisation", {}) or {}
    p_source = qc.get("p_value_source", {}) or {}
    derivations = qc.get("derivations", {}) or {}
    rows = _step(steps, rows, "Read input",
                 "Rows read from the source table." if not vcf else "Records read from the source GWAS-VCF.",
                 {"malformed VCF record (fewer than 8 fields)": _n((qc.get("liftover", {}).get("hg19_to_hg38") or
                                                                     qc.get("liftover", {}).get("hg38_to_hg19") or {}).get("malformed_record"))}
                 if vcf else None,
                 {"rsIDs extracted from compound identifiers": _n(ids.get("rsid_extracted")),
                  "P values from LOG10P": _n(p_source.get("from_log10p")), "P values from Z": _n(p_source.get("from_z")),
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

    lift = (qc.get("liftover", {}) or {}).get("hg19_to_hg38", {}) or {}
    lift_detail = ", ".join(f"{LIFTOVER_TEXT.get(k, k)}: {_n(v):,}" for k, v in lift.items()
                            if k not in ("ok", "minus_strand", "malformed_record", "non_primary_chromosome", "bad_position") and _n(v))
    if vcf:
        lifted = _n(qc.get("lifted_rows"))
        failed = sum(_n(v) for k, v in lift.items() if k not in ("ok", "minus_strand", "malformed_record", "non_primary_chromosome",
                                                                    "bad_position", "minus_strand_non_snv_dropped"))
        if rows is not None and "lifted_rows" in qc:
            failed = rows - lifted - _n(lift.get("non_primary_chromosome")) - _n(lift.get("bad_position")) \
                - _n(lift.get("minus_strand_non_snv_dropped"))
        rows = _step(steps, rows, "GRCh38 coordinates",
                     "Records must be on a primary chromosome with a valid position, and GRCh37 sources must lift to GRCh38.",
                     {CAUSE_TEXT["chromosome_not_primary"]: _n(lift.get("non_primary_chromosome")),
                      CAUSE_TEXT["invalid_position"]: _n(lift.get("bad_position")),
                      CAUSE_TEXT["liftover_to_grch38_failed"]: failed,
                      LIFTOVER_TEXT["minus_strand_non_snv_dropped"]: _n(lift.get("minus_strand_non_snv_dropped"))},
                     {"lifted from GRCh37 to GRCh38": _n(lift.get("ok")),
                      "alleles complemented after a minus-strand liftover": _n(lift.get("minus_strand"))},
                     [f"liftover failures: {lift_detail}"] if lift_detail else [])
    else:
        causes = qc.get("incomplete_by_cause", {}) or {}
        fills = qc.get("reference_fills", {}) or {}
        coordinate_causes = ("rsid_not_in_coordinate_index", "no_position_or_rsid", "liftover_to_grch38_failed",
                             "chromosome_not_primary", "invalid_position")
        notes = [f"{CAUSE_TEXT[c]}: {_rare_note(causes[c], threshold)}" for c in coordinate_causes if c in causes]
        if lift_detail and "liftover_to_grch38_failed" in causes:
            notes.append(f"liftover failures: {lift_detail}")
        rows = _step(steps, rows, "GRCh38 coordinates",
                     "Every variant needs a GRCh38 chromosome and position: taken from the file, parsed from chr:pos identifiers, "
                     "or looked up by rsID in the 1000 Genomes coordinate index, then lifted from GRCh37 where needed.",
                     {CAUSE_TEXT[c]: _n((causes.get(c) or {}).get("rows")) for c in coordinate_causes},
                     {"positions parsed from chr:pos identifiers": _n(qc.get("positions_from_variant_ids")),
                      "chromosome and GRCh37 position looked up by rsID": _n(fills.get("POS19")),
                      "lifted from GRCh37 to GRCh38": _n(lift.get("ok")),
                      "alleles complemented after a minus-strand liftover": _n(lift.get("minus_strand"))}, notes)
        statistic_causes = ("missing_allele", "missing_beta", "missing_or_invalid_se", "missing_or_invalid_p")
        rows = _step(steps, rows, "Alleles and statistics",
                     "Every variant needs both alleles, an effect size, a positive standard error and a P value.",
                     {CAUSE_TEXT[c]: _n((causes.get(c) or {}).get("rows")) for c in statistic_causes}, None,
                     [f"{CAUSE_TEXT[c]}: {_rare_note(causes[c], threshold)}" for c in statistic_causes if c in causes])
        rows = _step(steps, rows, "1000 Genomes annotation",
                     "Allele frequencies and rsIDs are matched from the 1000 Genomes lookup; nothing is removed here.",
                     None, {"EAF filled from 1000 Genomes": _n(fills.get("EAF")), "rsIDs filled from 1000 Genomes": _n(fills.get("rsid"))})

    if qc.get("analysis_qc"):
        rows = _analysis_steps(steps, rows, qc)
    else:
        steps.append({"step": len(steps) + 1, "title": "Analysis QC", "why": "Switched off (--analysis-qc false): no MAF, "
                      "palindromic, biallelic or duplicate filters.", "rows_in": rows, "removed": 0, "rows_out": rows,
                      "removed_by": {}, "changed": {}, "notes": []})
    alignment = qc.get("reference_alignment", {}) or {}
    if alignment.get("performed"):
        rows = _step(steps, rows, "GRCh38 reference alignment",
                     "REF must match the GRCh38 genome: alleles are swapped or complemented to match it, and SNVs that match on "
                     "neither strand are removed.",
                     {"alleles match neither strand of the GRCh38 base": _n(alignment.get("dropped_mismatch"))},
                     {"already aligned": _n(alignment.get("kept")), "effect and other allele swapped": _n(alignment.get("swapped")),
                      "complemented to the forward strand": _n(alignment.get("strand_flipped")),
                      "indels kept without a match": _n(alignment.get("indels_unmatched_kept")),
                      "kept, reference base unavailable": _n(alignment.get("unalignable_kept"))})
    output = _n(qc.get("output_rows"))
    if rows is not None and rows != output:
        steps[-1]["notes"].append(f"UNRECONCILED: the steps leave {rows:,} rows but {output:,} were written")
        qc.setdefault("warnings", []).append(f"step log does not reconcile: {rows} rows expected, {output} written")
    return steps


def format_steps(qc: dict[str, Any], record: dict[str, Any]) -> str:
    steps = qc.get("steps") or build_steps(qc)
    width = 12
    lines = [f"GWASprocessing step log: {record.get('study_id', '')} / {record.get('output_id', '')}",
             f"Source: {record.get('raw_path', '')}",
             f"Format: {qc.get('input_format', '')}; analysis QC {'on' if qc.get('analysis_qc') else 'off'}; "
             f"reference population {qc.get('reference_population_used', '')}",
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
