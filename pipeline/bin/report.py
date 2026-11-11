#!/usr/bin/env python3
"""Per-GWAS HTML report: step table, harmonisation counts, QC metrics and flags, warnings, timings and a raw vs
standardised QQ plot, in 1 self-contained file (no external files or links).

Usage, in the pipeline (after standardise_gwas.py):
    report.py --standardise-record standardise_record.json [--force false] [--dry-run false] --output row_record.json
After processing, for a spot check of 1 GWAS (its standardise folder or any file in it):
    report.py --dir GWAS/<STUDYID>/standardise/<SOURCE_ID> [--raw <raw file>]

QQ plots work in -log10 P (LP) space so P-values beyond double precision are plotted. Every point with LP >= 2
is drawn; the dense null region below is thinned to plotting resolution, which is visually lossless. The raw
panel reads the columns and delimiter the standardiser used. An existing report is reused only while it is newer
than the standardised VCF.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import html
import io
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import gwas_stats  # noqa: E402
from gwas_pipeline_common import (  # noqa: E402
    BETA_ALIASES,
    LOGP_ALIASES,
    P_ALIASES,
    SE_ALIASES,
    count_metadata_lines,
    infer_delimiter,
    is_gwas_vcf,
    numeric_series,
    parse_bool,
    pick_existing_column,
    read_json,
    write_json,
)

TABLE_CHUNK = 1_000_000
THIN_BELOW = 2.0
THIN_RESOLUTION = 0.005


def lp_from_vcf(path: Path) -> np.ndarray:
    opener = gzip.open if str(path).endswith(".gz") else open
    values: list[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            keys = parts[8].split(":")
            if "LP" in keys:
                sample = parts[9].split(":")
                index = keys.index("LP")
                values.append(sample[index] if index < len(sample) else ".")
    return numeric_series(pd.Series(values, dtype=object)).to_numpy(dtype="float64")


def lp_from_table(path: Path, source_columns: dict[str, str] | None = None, delimiter: str | None = None) -> np.ndarray:
    sep = delimiter or infer_delimiter(path)
    skip = count_metadata_lines(path)
    compression = "gzip" if str(path).endswith(".gz") else None
    header = list(pd.read_csv(path, sep=sep, compression=compression, nrows=0, skiprows=skip).columns)
    stripped = {str(c).lstrip("#"): c for c in header}
    source_columns = source_columns or {}

    def column(field: str, aliases: list[str]) -> str | None:
        chosen = source_columns.get(field)
        if chosen in stripped:
            return stripped[chosen]
        if source_columns:
            return None
        found = pick_existing_column(list(stripped), aliases)
        return stripped[found] if found else None

    cols = {field: column(field, aliases) for field, aliases in
            (("P", P_ALIASES), ("LOG10P", LOGP_ALIASES), ("BETA", BETA_ALIASES), ("SE", SE_ALIASES))}
    usecols = sorted({c for c in cols.values() if c})
    if not (cols["P"] or cols["LOG10P"]):
        return np.array([], dtype="float64")
    chunks = []
    reader = pd.read_csv(path, sep=sep, compression=compression, usecols=usecols, chunksize=TABLE_CHUNK, low_memory=False,
                         skiprows=skip, na_values=["."])
    for chunk in reader:
        get = lambda field: numeric_series(chunk[cols[field]]) if cols[field] else None  # noqa: E731
        lp, _ = gwas_stats.compute_lp(p=get("P"), log10p=get("LOG10P"), beta=get("BETA"), se=get("SE"))
        chunks.append(lp)
    return np.concatenate(chunks) if chunks else np.array([], dtype="float64")


def extract_lp(path: str | Path, source_columns: dict[str, str] | None = None, delimiter: str | None = None) -> np.ndarray:
    path = Path(path)
    if not path.exists():
        return np.array([], dtype="float64")
    if is_gwas_vcf(path):
        return lp_from_vcf(path)
    return lp_from_table(path, source_columns, delimiter)


def qq_points(lp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(expected, observed) LP pairs, largest first, with the LP < 2 bulk thinned to plotting resolution."""
    values = np.asarray(lp, dtype="float64")
    observed = np.sort(values[np.isfinite(values) & (values >= 0)])[::-1]
    n = observed.size
    if n == 0:
        return np.array([]), np.array([])
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    bulk = (observed < THIN_BELOW) & (expected < THIN_BELOW)
    cells = np.stack([np.round(expected[bulk] / THIN_RESOLUTION), np.round(observed[bulk] / THIN_RESOLUTION)], axis=1)
    _, first = np.unique(cells, axis=0, return_index=True)
    keep = np.concatenate([np.where(~bulk)[0], np.where(bulk)[0][np.sort(first)]])
    keep.sort()
    return expected[keep], observed[keep]


def draw_qq(ax: Any, lp: np.ndarray, title: str) -> float | None:
    expected, observed = qq_points(lp)
    if observed.size == 0:
        ax.text(0.5, 0.5, "No valid P values", ha="center", va="center")
        ax.set_title(title)
        return None
    lam = gwas_stats.lambda_gc(lp)
    n = int(np.isfinite(lp).sum())
    ax.scatter(expected, observed, s=5, alpha=0.6, linewidths=0, rasterized=True)
    limit = float(max(expected.max(), observed.max(), 1.0))
    ax.plot([0, limit], [0, limit], color="black", linewidth=0.8)
    ax.set_xlim(0, limit * 1.02)
    ax.set_ylim(0, limit * 1.02)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(f"{title} (n={n:,}, lambda={lam:.3f})" if lam is not None else f"{title} (n={n:,})")
    return lam


CSS = """
:root { --bg: #ffffff; --fg: #1b1f24; --muted: #5b6470; --line: #d9dee4; --flag: #a33a00; --head: #f3f5f7; }
@media (prefers-color-scheme: dark) { :root { --bg: #15181c; --fg: #e6e9ed; --muted: #9aa4b0; --line: #2e343b; --flag: #ff9b62; --head: #1d2126; } }
body { background: var(--bg); color: var(--fg); font: 14px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; margin: 0; padding: 24px 16px; }
main { max-width: 1040px; margin: 0 auto; }
h1 { font-size: 22px; margin: 0 0 4px; } h2 { font-size: 17px; margin: 28px 0 8px; } h3 { font-size: 14px; margin: 18px 0 4px; }
.muted { color: var(--muted); } .flag { color: var(--flag); font-weight: 600; }
table { border-collapse: collapse; width: 100%; margin: 8px 0; } th, td { border-bottom: 1px solid var(--line); padding: 4px 8px; text-align: left; vertical-align: top; }
th { background: var(--head); font-weight: 600; } td.n, th.n { text-align: right; font-variant-numeric: tabular-nums; }
ul { margin: 4px 0 4px 18px; padding: 0; } img { max-width: 100%; height: auto; }
"""


def _e(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def _n(value: Any) -> str:
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return _e(value)


def _fixed(value: Any, digits: int = 3) -> str:
    return f"{value:.{digits}f}" if isinstance(value, (int, float)) and value == value else "not available"


def _table(rows: list[tuple[Any, Any]]) -> str:
    return "<table>" + "".join(f"<tr><th>{_e(k)}</th><td>{v}</td></tr>" for k, v in rows) + "</table>"


def qq_png(raw_lp: np.ndarray, std_lp: np.ndarray, title: str) -> tuple[str, float | None, float | None]:
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)
    lam_raw = draw_qq(axes[0], raw_lp, "Raw QQ")
    lam_std = draw_qq(axes[1], std_lp, "Standardised QQ")
    fig.suptitle(title, fontsize=11)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=110)
    plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode("ascii"), lam_raw, lam_std


def render(record: dict[str, Any], qc: dict[str, Any], png: str) -> str:
    sample = qc.get("sample_header", {}) or {}
    build = qc.get("build_detection", {}) or {}
    flags = qc.get("qc_flags", {}) or {}
    metrics = qc.get("qc_metrics", {}) or {}
    parts = [f"<h1>{_e(record.get('study_id'))} / {_e(record.get('output_id'))}</h1>",
             f"<p class='muted'>{_e(record.get('phenotype'))} &middot; {_e(record.get('ancestry'))} &middot; "
             f"GWASprocessing {_e((qc.get('provenance') or {}).get('schema_version'))} &middot; completed {_e(qc.get('completed_at'))}</p>"]
    parts.append("<h2>Summary</h2>" + _table([
        ("Source file", _e(record.get("raw_path"))), ("Input format", _e(qc.get("input_format"))),
        ("Input build", _e(f"{build.get('build', '')} ({build.get('method', '')})" if build else "positions from rsIDs (dbSNP)")),
        ("Reference frequencies", _e(f"{qc.get('reference_population_used', '')} {qc.get('reference_frequencies', '')}")),
        ("Rows read", _n(qc.get("input_rows"))), ("Variants written", _n(qc.get("output_rows"))),
        ("Lambda GC", _e(_fixed(qc.get("lambda_gc")))),
        ("AF correlation with 1000 Genomes", _e(_fixed((qc.get("af_concordance") or {}).get("pearson_r")))),
        ("Settings hash", _e(qc.get("settings_hash")))]))
    parts.append("<h2>QC flags</h2>" + ("<ul>" + "".join(f"<li class='flag'>{_e(k)}: {_e(v)}</li>" for k, v in flags.items()) + "</ul>"
                                         if flags else "<p>No OpenGWAS QC flag raised.</p>")
                 + _table([(k, _e(f"{v:.4g}" if isinstance(v, float) else v)) for k, v in metrics.items()]))
    if sample:
        parts.append("<h2>##SAMPLE header</h2>" + _table([(k, _n(v) if isinstance(v, int) else _e(v)) for k, v in sample.items()]))
    steps = qc.get("steps") or []
    rows = "".join(f"<tr><td>{s['step']}. {_e(s['title'])}</td><td class='n'>{_n(s['rows_in'])}</td><td class='n'>{_n(s['removed'])}</td>"
                   f"<td class='n'>{_n(s['rows_out'])}</td></tr>" for s in steps)
    parts.append("<h2>Steps</h2><table><tr><th>Step</th><th class='n'>Rows in</th><th class='n'>Removed</th><th class='n'>Rows out</th></tr>"
                 + rows + "</table>")
    for s in steps:
        items = [f"<li>removed {_n(c)}: {_e(r)}</li>" for r, c in s["removed_by"].items()]
        items += [f"<li>changed {_n(c)}: {_e(r)}</li>" for r, c in s["changed"].items()]
        items += [f"<li class='muted'>{_e(n)}</li>" for n in s["notes"]]
        parts.append(f"<h3>{s['step']}. {_e(s['title'])}</h3><p class='muted'>{_e(s['why'])}</p>"
                     + (f"<ul>{''.join(items)}</ul>" if items else ""))
    if png:
        parts.append(f"<h2>QQ plots</h2><img alt='Raw and standardised QQ plots' src='data:image/png;base64,{png}'>")
    if qc.get("warnings"):
        parts.append("<h2>Warnings</h2><ul>" + "".join(f"<li>{_e(w)}</li>" for w in qc["warnings"]) + "</ul>")
    timings = qc.get("timings_s") or {}
    if timings:
        parts.append("<h2>Time (seconds)</h2>" + _table([(k, _e(v)) for k, v in {**(timings.get("wall") or {}),
                                                                                  **(timings.get("piece_steps_summed_over_pieces") or {})}.items()]))
    title = f"{record.get('output_id', '')} report"
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>"
            f"<title>{_e(title)}</title><style>{CSS}</style></head><body><main>{''.join(parts)}</main></body></html>\n")


def write_report(record: dict[str, Any], qc: dict[str, Any], report_path: Path, raw_path: str | Path | None) -> dict[str, Any]:
    raw_lp = extract_lp(raw_path, qc.get("source_columns"), qc.get("input_delimiter")) if raw_path else np.array([])
    std_lp = extract_lp(record["standard_vcf"])
    png, lam_raw, lam_std = qq_png(raw_lp, std_lp, f"{record.get('study_id', '')} | {record.get('phenotype', '')} | {record.get('ancestry', '')}")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = report_path.with_name(f".{report_path.name}.tmp")
    tmp.write_text(render(record, qc, png), encoding="utf-8")
    tmp.replace(report_path)
    return {"lambda_gc_raw": lam_raw, "lambda_gc_standardised": lam_std}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--standardise-record", default="")
    parser.add_argument("--dir", default="", help="a standardise folder (or a file in it), to report 1 GWAS after processing")
    parser.add_argument("--raw", default="", help="with --dir: the raw file for the raw QQ panel (default: the one recorded)")
    parser.add_argument("--force", default="false")
    parser.add_argument("--dry-run", default="false")
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    if args.dir:
        folder = Path(args.dir) if Path(args.dir).is_dir() else Path(args.dir).parent
        qc = read_json(folder / "standardised.qc.json")
        metadata = read_json(folder / "standardised.metadata.json")
        record = {**metadata, "standard_vcf": str(folder / "standardised.gwas.vcf.gz")}
        raw = args.raw or metadata.get("raw_path", "")
        out = folder / "standardised.report.html"
        write_report(record, qc, out, raw if raw and Path(raw).exists() else None)
        print(out)
        return

    record = read_json(args.standardise_record)
    report_path = Path(record["report_html"])
    vcf_path = Path(record["standard_vcf"])
    if parse_bool(args.dry_run):
        record["report_created"] = False
        write_json(args.output, record)
        return
    fresh = report_path.exists() and report_path.stat().st_size > 0 and vcf_path.exists() \
        and report_path.stat().st_mtime >= vcf_path.stat().st_mtime
    if fresh and not parse_bool(args.force):
        record["report_created"] = False
        write_json(args.output, record)
        return
    record.update(write_report(record, record.get("qc") or {}, report_path, record.get("raw_path")))
    record["report_created"] = True
    write_json(args.output, record)


if __name__ == "__main__":
    main()
