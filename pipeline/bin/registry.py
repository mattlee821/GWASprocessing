#!/usr/bin/env python3
"""Write <work root>/manifest.html: 1 self-contained page listing every GWAS this install has processed or tried to.

Usage:
    python3 pipeline/bin/registry.py --work-root GWAS [--log-root temp/logs]

src/standardise.sh runs it at the end of every run. Rows come from 3 places: the state file of completed rows
(<work root>/.standardise_state.tsv), the rows the last run tried (<log root>/run_manifest.tsv with
state_update_summary.json, which records which of them completed), and the registry kept beside the page
(<work root>/.manifest_registry.json), which remembers failed rows between runs. For each row the page shows the
manifest fields, its status, the variant counts and QC flags from its QC JSON, and links (relative to the work root)
to its output folder, VCF and HTML report.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
from pathlib import Path
from typing import Any

from gwas_pipeline_common import output_stem, read_json, read_state, utc_now, write_json

REGISTRY = ".manifest_registry.json"
PAGE = "manifest.html"
STATE = ".standardise_state.tsv"
MANIFEST_FIELDS = ("GWAS_location", "phenotype", "ancestry", "sex", "author", "year", "PMID", "source_type", "study_yaml")


def _read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _json(path: Any) -> dict[str, Any]:
    try:
        return read_json(path) if path and Path(path).exists() else {}
    except (OSError, ValueError):
        return {}


def update_registry(work_root: Path, log_root: Path | None) -> dict[str, dict[str, Any]]:
    """Merge the state file and the last run into the registry and return it, keyed by output folder (1 row per GWAS
    output, whichever manifest row last wrote it)."""
    path = work_root / REGISTRY
    registry: dict[str, dict[str, Any]] = _json(path).get("rows", {}) if path.exists() else {}
    state = sorted(read_state(work_root / STATE).values(), key=lambda r: r.get("updated_at", ""))
    for row in state:  # oldest first, so the latest completion of a folder wins
        folder = str(Path(row["standard_vcf"]).parent) if row.get("standard_vcf") else row["row_hash"]
        entry = registry.setdefault(folder, {})
        if entry.get("status") == "failed" and entry.get("failed_at", "") > row.get("updated_at", ""):
            continue  # a later attempt at this folder failed
        entry.update({"study_id": row.get("study_id", ""), "output_id": row.get("output_id", ""),
                      "GWAS_location": entry.get("GWAS_location") or row.get("gwas_location", ""), "standardise_dir": folder,
                      "status": "complete", "updated_at": max(entry.get("updated_at", ""), row.get("updated_at", ""))})
    if log_root is not None:
        summary = _json(log_root / "run_manifest_summary.json")
        done = set(_json(log_root / "state_update_summary.json").get("updated_row_hashes", []))
        if summary and not summary.get("dry_run"):
            for row in _read_tsv(log_root / "run_manifest.tsv"):
                entry = registry.setdefault(row.get("standardise_dir") or row["row_hash"], {})
                entry.update({key: row.get(key, "") for key in MANIFEST_FIELDS + ("study_id", "output_id", "standardise_dir")})
                if row["row_hash"] in done:
                    entry.update({"status": "complete", "updated_at": utc_now()})
                else:
                    entry.update({"status": "failed", "failed_at": utc_now()})
    write_json(path, {"updated_at": utc_now(), "rows": registry})
    return registry


def row_view(entry: dict[str, Any], work_root: Path) -> dict[str, Any]:
    """Everything the page shows for 1 row, read from its output folder."""
    folder = Path(entry.get("standardise_dir") or "")
    stem = output_stem(folder) if folder.parts else "standardised"
    vcf, qc_path, report = (folder / name for name in (f"{stem}.gwas.vcf.gz", f"{stem}.qc.json", f"{stem}.report.html"))
    metadata = _json(folder / f"{stem}.metadata.json")
    qc = _json(qc_path)
    status = entry.get("status", "")
    if status == "complete" and not vcf.exists():
        status = "outputs missing"
    sample = qc.get("sample_header", {}) or {}
    build = (qc.get("build_detection") or {}).get("build") or ("rsIDs" if qc.get("coordinate_source") == "dbsnp" else "")
    relative = (lambda p: os.path.relpath(p, work_root)) if folder.parts else (lambda p: "")
    return {
        "status": status, "study_id": entry.get("study_id", ""), "output_id": entry.get("output_id", ""),
        "phenotype": entry.get("phenotype") or metadata.get("phenotype", ""), "ancestry": entry.get("ancestry") or metadata.get("ancestry", ""),
        "sex": entry.get("sex") or metadata.get("sex", ""), "author": entry.get("author") or metadata.get("author", ""),
        "year": entry.get("year") or metadata.get("year", ""), "PMID": entry.get("PMID") or metadata.get("PMID", ""),
        "GWAS_location": entry.get("GWAS_location") or metadata.get("GWAS_location", ""), "study_yaml": entry.get("study_yaml", ""),
        "study_type": sample.get("StudyType", ""), "cases": sample.get("TotalCases", ""), "controls": sample.get("TotalControls", ""),
        "input_rows": qc.get("input_rows", ""), "output_rows": qc.get("output_rows", ""), "af_filled": sample.get("AFFilled", ""),
        "build": {"hg19": "GRCh37", "hg38": "GRCh38"}.get(build, build), "lambda_gc": qc.get("lambda_gc", ""),
        "flags": list((qc.get("qc_flags") or {}).values()), "warnings": qc.get("warnings") or [],
        "processed": (qc.get("completed_at") or entry.get("updated_at") or entry.get("failed_at") or "")[:16].replace("T", " "),
        "settings_hash": qc.get("settings_hash", ""),
        "folder": relative(folder) if folder.exists() else "", "vcf": relative(vcf) if vcf.exists() else "",
        "report": relative(report) if report.exists() else "",
    }


def _cell(value: Any) -> str:
    return html.escape("" if value is None else str(value))


def _number(value: Any) -> str:
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return _cell(value)


def render(rows: list[dict[str, Any]], work_root: Path) -> str:
    counts = {status: sum(r["status"] == status for r in rows) for status in ("complete", "failed", "outputs missing")}
    body = []
    for r in sorted(rows, key=lambda r: (r["study_id"], r["output_id"])):
        kept = ""
        if r["input_rows"] not in ("", None) and r["output_rows"] not in ("", None) and int(r["input_rows"]):
            kept = f'{_number(r["output_rows"])} of {_number(r["input_rows"])} ({100 * int(r["output_rows"]) / int(r["input_rows"]):.2f}%)'
        sample = ""
        if r["study_type"] == "CaseControl":
            sample = f'{_number(r["cases"])} cases, {_number(r["controls"])} controls'
        elif r["controls"] not in ("", None):
            sample = f'N {_number(r["controls"])}'
        location = r["GWAS_location"]
        source = f'<a href="{_cell(location)}">{_cell(location)}</a>' if location.startswith(("http://", "https://")) else _cell(location)
        pmid = f'<a href="https://pubmed.ncbi.nlm.nih.gov/{_cell(r["PMID"])}/">{_cell(r["PMID"])}</a>' if str(r["PMID"]).isdigit() else _cell(r["PMID"])
        links = " ".join(f'<a href="{_cell(r[key])}">{label}</a>' for key, label in (("folder", "folder"), ("report", "report"), ("vcf", "VCF")) if r[key])
        flags = "".join(f"<li>{_cell(f)}</li>" for f in r["flags"])
        warnings = f'<details><summary>{len(r["warnings"])}</summary><ul>{"".join(f"<li>{_cell(w)}</li>" for w in r["warnings"])}</ul></details>' if r["warnings"] else ""
        lam = f'{float(r["lambda_gc"]):.3f}' if r["lambda_gc"] not in ("", None) else ""
        status_class = r["status"].replace(" ", "-")
        body.append(
            f'<tr class="{status_class}"><td><span class="status {status_class}">{_cell(r["status"])}</span></td>'
            f'<td>{_cell(r["study_id"])}</td><td>{_cell(r["output_id"])}</td><td>{_cell(r["phenotype"])}</td>'
            f'<td>{_cell(r["ancestry"])}</td><td>{_cell(r["sex"])}</td><td>{_cell(r["author"])} {_cell(r["year"])}</td><td>{pmid}</td>'
            f'<td class="wrap">{source}</td><td>{_cell(r["build"])}</td><td>{sample}</td><td class="num">{kept}</td>'
            f'<td class="num">{_number(r["af_filled"])}</td><td class="num">{lam}</td><td><ul>{flags}</ul></td><td>{warnings}</td>'
            f'<td>{_cell(r["processed"])}</td><td>{links}</td></tr>')
    headers = ["Status", "Study", "GWAS", "Phenotype", "Ancestry", "Sex", "Author", "PMID", "Source", "Input build", "Sample",
               "Variants written", "AF filled", "Lambda GC", "QC flags", "Warnings", "Processed (UTC)", "Links"]
    head = "".join(f"<th>{h}</th>" for h in headers)
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>GWASprocessing manifest</title>
<style>
:root {{ --bg: #ffffff; --fg: #1f2328; --muted: #59636e; --line: #d1d9e0; --head: #f6f8fa; --ok: #1a7f37; --bad: #cf222e; --warn: #9a6700; --link: #0969da; }}
@media (prefers-color-scheme: dark) {{ :root {{ --bg: #0d1117; --fg: #e6edf3; --muted: #9198a1; --line: #3d444d; --head: #151b23; --ok: #3fb950; --bad: #f85149; --warn: #d29922; --link: #4493f8; }} }}
body {{ margin: 0; padding: 24px 16px; background: var(--bg); color: var(--fg); font: 14px/1.45 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }}
h1 {{ font-size: 20px; margin: 0 0 4px; }} p {{ margin: 4px 0 12px; color: var(--muted); }}
input {{ width: min(420px, 100%); padding: 6px 10px; border: 1px solid var(--line); border-radius: 6px; background: var(--bg); color: var(--fg); margin-bottom: 12px; }}
.scroll {{ overflow-x: auto; border: 1px solid var(--line); border-radius: 6px; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ padding: 6px 10px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; white-space: nowrap; }}
th {{ background: var(--head); position: sticky; top: 0; cursor: pointer; user-select: none; }}
td.wrap {{ white-space: normal; word-break: break-all; min-width: 220px; }} td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
ul {{ margin: 0; padding-left: 16px; }} a {{ color: var(--link); }}
.status {{ font-weight: 600; }} .status.complete {{ color: var(--ok); }} .status.failed {{ color: var(--bad); }} .status.outputs-missing {{ color: var(--warn); }}
</style></head><body>
<h1>GWASprocessing manifest</h1>
<p>{counts["complete"]} complete, {counts["failed"]} failed, {counts["outputs missing"]} with outputs missing. Work root: {_cell(work_root)}. Updated {_cell(utc_now()[:16].replace("T", " "))} UTC. Links are relative to this page.</p>
<input id="filter" type="search" placeholder="Filter rows" aria-label="Filter rows">
<div class="scroll"><table id="rows"><thead><tr>{head}</tr></thead><tbody>
{chr(10).join(body)}
</tbody></table></div>
<script>
const rows = Array.from(document.querySelectorAll('#rows tbody tr'));
document.getElementById('filter').addEventListener('input', e => {{
  const q = e.target.value.toLowerCase();
  rows.forEach(r => {{ r.style.display = r.textContent.toLowerCase().includes(q) ? '' : 'none'; }});
}});
document.querySelectorAll('#rows th').forEach((th, i) => th.addEventListener('click', () => {{
  const asc = th.dataset.asc !== 'true'; th.dataset.asc = asc;
  const key = r => r.children[i].textContent.trim();
  const num = s => parseFloat(s.replace(/,/g, ''));
  rows.sort((a, b) => {{ const x = key(a), y = key(b); const d = (!isNaN(num(x)) && !isNaN(num(y))) ? num(x) - num(y) : x.localeCompare(y); return asc ? d : -d; }});
  const tbody = document.querySelector('#rows tbody'); rows.forEach(r => tbody.appendChild(r));
}}));
</script>
</body></html>
"""


def write_page(work_root: Path, log_root: Path | None = None) -> Path:
    work_root.mkdir(parents=True, exist_ok=True)
    registry = update_registry(work_root, log_root)
    page = work_root / PAGE
    tmp = page.with_suffix(".html.tmp")
    tmp.write_text(render([row_view(entry, work_root) for entry in registry.values()], work_root), encoding="utf-8")
    os.replace(tmp, page)
    return page


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--work-root", required=True)
    parser.add_argument("--log-root", default="")
    args = parser.parse_args()
    page = write_page(Path(args.work_root), Path(args.log_root) if args.log_root else None)
    print(f"  Manifest page: {page}")


if __name__ == "__main__":
    main()
