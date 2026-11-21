#!/usr/bin/env python3
"""Split a study manifest into fixed-size batch manifests, check batch progress, and safely purge batch VCFs.

Subcommands:
    batches.py split     --manifest FULL.tsv --out-dir DIR [--size 100] [--force]
    batches.py status    --batches DIR [--work-root GWAS] [--export-dir DIR]
    batches.py clean-raw --batches DIR --batch ID [--work-root GWAS] [--clean-temp]
    batches.py purge     --batches DIR --batch ID [--work-root GWAS] [--export-dir DIR] [--allow-partial]
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any

from gwas_pipeline_common import (
    atomic_write_text,
    infer_source_type,
    normalize_source_type,
    output_id_from_row,
    read_manifest,
    study_id_from_row,
    utc_now,
    valid_output_name,
    write_json,
)
import registry as registry_mod


def _dir_bytes(path: Path) -> int:
    if not path.exists():
        return 0
    if path.is_file():
        return path.stat().st_size
    total = 0
    for item in path.rglob("*"):
        if item.is_file() and not item.is_symlink():
            try:
                total += item.stat().st_size
            except OSError:
                pass
    return total


def _row_paths(row: dict[str, str], work_root: Path) -> dict[str, Any]:
    study_id = study_id_from_row(row)
    source_type = normalize_source_type(row["source_type"]) if row.get("source_type") else infer_source_type(row["GWAS_location"])
    output_id = output_id_from_row(row, source_type)
    output_name = str(row.get("output_name", "") or "").strip()
    stem = output_name if output_name and valid_output_name(output_name) else "standardised"
    std_dir = work_root / study_id / "standardise" / output_id
    raw_dir = work_root / study_id / "raw" / output_id
    return {
        "study_id": study_id,
        "output_id": output_id,
        "stem": stem,
        "std_dir": std_dir,
        "raw_dir": raw_dir,
        "vcf": std_dir / f"{stem}.gwas.vcf.gz",
        "tbi": std_dir / f"{stem}.gwas.vcf.gz.tbi",
        "qc": std_dir / f"{stem}.qc.json",
        "meta": std_dir / f"{stem}.metadata.json",
        "report": std_dir / f"{stem}.report.html",
        "purged": std_dir / f"{stem}.purged.json",
    }


def _export_present(export_dir: Path | None, output_id: str) -> bool:
    if export_dir is None:
        return True
    if not export_dir.exists():
        return False
    direct = export_dir / output_id
    if direct.is_dir() and any(p.is_file() and p.stat().st_size > 0 for p in direct.rglob("*")):
        return True
    if direct.is_file() and direct.stat().st_size > 0:
        return True
    for match in export_dir.glob(f"*{output_id}*"):
        if match.is_file() and match.stat().st_size > 0:
            return True
        if match.is_dir() and any(p.is_file() and p.stat().st_size > 0 for p in match.rglob("*")):
            return True
    return False


def _batch_files(batches_dir: Path) -> list[Path]:
    if not batches_dir.is_dir():
        raise SystemExit(f"Batches folder not found: {batches_dir}")
    files = sorted(p for p in batches_dir.glob("batch_*.tsv") if p.is_file())
    if not files:
        raise SystemExit(f"No batch_*.tsv manifests found in {batches_dir}")
    return files


def _resolve_batch(batches_dir: Path, spec: str) -> Path:
    files = _batch_files(batches_dir)
    clean = spec.strip()
    for candidate in files:
        if candidate.name == clean or candidate.stem == clean:
            return candidate
    if re.fullmatch(r"\d+", clean):
        num = int(clean)
        for candidate in files:
            m = re.fullmatch(r"batch_0*(\d+)\.tsv", candidate.name)
            if m and int(m.group(1)) == num:
                return candidate
    raise SystemExit(f"Batch '{spec}' not found in {batches_dir} (available: {files[0].name} .. {files[-1].name})")


def cmd_split(args: argparse.Namespace) -> None:
    manifest = Path(args.manifest)
    out_dir = Path(args.out_dir)
    size = int(args.size)
    if size < 1:
        raise SystemExit("--size must be at least 1")

    parsed_rows = read_manifest(manifest)
    if not parsed_rows:
        raise SystemExit(f"Manifest has no rows: {manifest}")

    raw_lines = manifest.read_text(encoding="utf-8").splitlines()
    header = raw_lines[0]
    body_lines = [raw_lines[int(r["_line"]) - 1] for r in parsed_rows]

    seen_ids: set[tuple[str, str]] = set()
    for r in parsed_rows:
        p = _row_paths(r, Path("GWAS"))
        key = (p["study_id"], p["output_id"])
        if key in seen_ids:
            raise SystemExit(f"Duplicate study/output_id in manifest: {key[0]}/{key[1]}")
        seen_ids.add(key)

    if out_dir.exists() and any(out_dir.glob("batch_*.tsv")) and not args.force:
        raise SystemExit(f"{out_dir} already contains batch_*.tsv files; pass --force to overwrite them")
    out_dir.mkdir(parents=True, exist_ok=True)
    for old in out_dir.glob("batch_*.tsv"):
        old.unlink()

    n_batches = (len(body_lines) + size - 1) // size
    width = max(2, len(str(n_batches)))
    index_lines = ["batch\tmanifest\trows\tfirst_output\tlast_output\n"]

    for idx in range(n_batches):
        start, end = idx * size, min(len(body_lines), (idx + 1) * size)
        chunk_lines = body_lines[start:end]
        chunk_rows = parsed_rows[start:end]
        batch_name = f"batch_{idx + 1:0{width}d}.tsv"
        batch_path = out_dir / batch_name
        atomic_write_text(batch_path, "\n".join([header, *chunk_lines]) + "\n")
        first_out = _row_paths(chunk_rows[0], Path("GWAS"))["output_id"]
        last_out = _row_paths(chunk_rows[-1], Path("GWAS"))["output_id"]
        index_lines.append(f"{idx + 1:0{width}d}\t{batch_path}\t{len(chunk_lines)}\t{first_out}\t{last_out}\n")

    atomic_write_text(out_dir / "batches_index.tsv", "".join(index_lines))
    print(f"Split {len(body_lines):,} rows from {manifest} into {n_batches} batches (size {size}) in {out_dir}")


def cmd_status(args: argparse.Namespace) -> None:
    batches_dir = Path(args.batches)
    work_root = Path(args.work_root)
    export_dir = Path(args.export_dir) if args.export_dir else None
    files = _batch_files(batches_dir)

    headers = ["Batch", "Rows", "RawStaged", "VCFs", "QCDone"] + (["Exported"] if export_dir else []) + ["State"]
    table: list[list[str]] = []
    counts = {"DONE": 0, "STANDARDISED": 0, "PARTIAL": 0, "PENDING": 0}
    next_pending = ""

    for bf in files:
        rows = read_manifest(bf)
        n = len(rows)
        n_raw = n_vcf = n_qc = n_purged = n_exp = 0
        for r in rows:
            p = _row_paths(r, work_root)
            has_vcf = p["vcf"].is_file() and p["vcf"].stat().st_size > 0 and p["tbi"].is_file() and p["tbi"].stat().st_size > 0
            has_qc = p["qc"].is_file() and p["qc"].stat().st_size > 0
            has_purged = p["purged"].is_file() or (has_qc and not has_vcf)
            has_raw = p["raw_dir"].is_dir() and any(p["raw_dir"].iterdir())
            n_vcf += int(has_vcf)
            n_qc += int(has_qc)
            n_purged += int(has_qc and not has_vcf and has_purged)
            n_raw += int(has_raw)
            if export_dir:
                n_exp += int(_export_present(export_dir, p["output_id"]))

        if n > 0 and n_purged == n and n_vcf == 0:
            state = "DONE (purged)"
            counts["DONE"] += 1
        elif n > 0 and n_vcf == n and n_qc == n:
            if export_dir and n_exp == n:
                state = "EXPORTED (ready to purge)"
            elif n_raw > 0:
                state = "STANDARDISED (raw staged)"
            else:
                state = "STANDARDISED (ready)"
            counts["STANDARDISED"] += 1
        elif n_vcf == 0 and n_qc == 0 and n_raw == 0:
            state = "PENDING"
            counts["PENDING"] += 1
            if not next_pending:
                next_pending = str(bf)
        else:
            state = "PARTIAL / IN PROGRESS"
            counts["PARTIAL"] += 1

        row_cells = [bf.stem, str(n), f"{n_raw}/{n}", f"{n_vcf}/{n}", f"{n_qc}/{n}"]
        if export_dir:
            row_cells.append(f"{n_exp}/{n}")
        row_cells.append(state)
        table.append(row_cells)

    widths = [max(len(h), *(len(r[i]) for r in table)) for i, h in enumerate(headers)]
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    print("  ".join("-" * w for w in widths))
    for r in table:
        print(fmt.format(*r))
    print("")
    print(f"Summary ({len(files)} batches, work_root={work_root}):")
    print(f"  Done (purged): {counts['DONE']} | Standardised: {counts['STANDARDISED']} | "
          f"Partial: {counts['PARTIAL']} | Pending: {counts['PENDING']}")
    if next_pending:
        print(f"  Next pending batch: {next_pending}")
    if work_root.exists():
        usage = shutil.disk_usage(work_root)
        print(f"  Disk free at work_root: {usage.free / 1e9:.1f} GB (of {usage.total / 1e9:.1f} GB)")


def cmd_clean_raw(args: argparse.Namespace) -> None:
    bf = _resolve_batch(Path(args.batches), args.batch)
    work_root = Path(args.work_root)
    rows = read_manifest(bf)
    freed = 0
    cleaned = skipped = 0
    for r in rows:
        p = _row_paths(r, work_root)
        ready = p["vcf"].is_file() and p["vcf"].stat().st_size > 0 and p["tbi"].is_file() and p["qc"].is_file()
        if not ready:
            skipped += 1
            continue
        if p["raw_dir"].exists():
            freed += _dir_bytes(p["raw_dir"])
            shutil.rmtree(p["raw_dir"], ignore_errors=True)
            cleaned += 1

    if args.clean_temp:
        seen_archives: set[Path] = set()
        for r in rows:
            arch_dir = _row_paths(r, work_root)["raw_dir"].parent / "_archives"
            if arch_dir not in seen_archives and arch_dir.is_dir():
                seen_archives.add(arch_dir)
                freed += _dir_bytes(arch_dir)
                shutil.rmtree(arch_dir, ignore_errors=True)
        for env_key in ("GWAS_NEXTFLOW_WORK", "GWAS_TMPDIR"):
            val = os.environ.get(env_key, "").strip()
            if val and Path(val).is_dir():
                freed += _dir_bytes(Path(val))
                shutil.rmtree(val, ignore_errors=True)
                Path(val).mkdir(parents=True, exist_ok=True)

    print(f"{bf.name}: removed staged raw folders for {cleaned} standardised row(s) "
          f"({skipped} row(s) not ready/already clean); freed {freed / 1e9:.2f} GB")


def cmd_purge(args: argparse.Namespace) -> None:
    bf = _resolve_batch(Path(args.batches), args.batch)
    work_root = Path(args.work_root)
    export_dir = Path(args.export_dir) if args.export_dir else None
    rows = read_manifest(bf)

    missing_qc: list[str] = []
    missing_export: list[str] = []
    eligible_paths: list[dict[str, Any]] = []

    for r in rows:
        p = _row_paths(r, work_root)
        has_qc = p["qc"].is_file() and p["qc"].stat().st_size > 0 and p["meta"].is_file() and p["meta"].stat().st_size > 0
        if not has_qc:
            missing_qc.append(p["output_id"])
            continue
        if export_dir and not _export_present(export_dir, p["output_id"]):
            missing_export.append(p["output_id"])
            continue
        eligible_paths.append(p)

    if missing_export:
        preview = ", ".join(missing_export[:10])
        raise SystemExit(
            f"Refusing to purge {bf.name}: {len(missing_export)} row(s) have no exported output in {export_dir} "
            f"(e.g. {preview})"
        )
    if missing_qc and not args.allow_partial:
        preview = ", ".join(missing_qc[:10])
        raise SystemExit(
            f"Refusing to purge {bf.name}: {len(missing_qc)} row(s) are missing .qc.json/.metadata.json "
            f"(e.g. {preview}). Finish standardising the batch first, or pass --allow-partial."
        )

    freed = 0
    purged_count = 0
    stamp = utc_now()
    for p in eligible_paths:
        for target in (p["vcf"], p["tbi"], p["report"]):
            if target.is_file():
                freed += target.stat().st_size
                target.unlink()
        if p["raw_dir"].exists():
            freed += _dir_bytes(p["raw_dir"])
            shutil.rmtree(p["raw_dir"], ignore_errors=True)
        write_json(p["purged"], {"batch": bf.name, "output_id": p["output_id"], "purged_at": stamp})
        purged_count += 1

    if work_root.exists():
        try:
            with registry_mod.file_lock(work_root / registry_mod.PAGE):
                reg = registry_mod.update_registry(work_root, None)
                views = [registry_mod.row_view(entry, work_root) for entry in reg.values()]
                atomic_write_text(work_root / registry_mod.PAGE, registry_mod.render(views, work_root))
        except Exception as exc:  # noqa: BLE001
            print(f"Warning: could not refresh {work_root / 'manifest.html'}: {exc}", file=sys.stderr)

    print(f"{bf.name}: purged VCF/raw files for {purged_count}/{len(rows)} row(s) "
          f"(kept .qc.json and .metadata.json); freed {freed / 1e9:.2f} GB")


def main() -> None:
    parser = argparse.ArgumentParser(description="Manage fixed-size GWASprocessing manifest batches.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_split = sub.add_parser("split", help="Split a full manifest TSV into fixed-size batch_XX.tsv files")
    p_split.add_argument("--manifest", required=True, help="Full study manifest TSV")
    p_split.add_argument("--out-dir", required=True, help="Output folder for batch_XX.tsv files")
    p_split.add_argument("--size", type=int, default=100, help="Rows per batch (default: 100)")
    p_split.add_argument("--force", action="store_true", help="Overwrite existing batch_*.tsv files in --out-dir")
    p_split.set_defaults(func=cmd_split)

    p_status = sub.add_parser("status", help="Show status of every batch in a batches folder")
    p_status.add_argument("--batches", required=True, help="Folder containing batch_XX.tsv files")
    p_status.add_argument("--work-root", default=os.environ.get("GWAS_WORK_ROOT", "GWAS"), help="GWAS work root")
    p_status.add_argument("--export-dir", default="", help="Optional finemapping export folder to verify per protein")
    p_status.set_defaults(func=cmd_status)

    p_clean = sub.add_parser("clean-raw", help="Remove raw/<output_id> staging files for standardised rows in 1 batch")
    p_clean.add_argument("--batches", required=True, help="Folder containing batch_XX.tsv files")
    p_clean.add_argument("--batch", required=True, help="Batch number or name (e.g. 01, batch_01, batch_01.tsv)")
    p_clean.add_argument("--work-root", default=os.environ.get("GWAS_WORK_ROOT", "GWAS"), help="GWAS work root")
    p_clean.add_argument("--clean-temp", action="store_true", help="Also clear GWAS_NEXTFLOW_WORK and GWAS_TMPDIR")
    p_clean.set_defaults(func=cmd_clean_raw)

    p_purge = sub.add_parser("purge", help="Delete .gwas.vcf.gz* and raw staging for 1 batch while keeping .qc.json")
    p_purge.add_argument("--batches", required=True, help="Folder containing batch_XX.tsv files")
    p_purge.add_argument("--batch", required=True, help="Batch number or name (e.g. 01, batch_01, batch_01.tsv)")
    p_purge.add_argument("--work-root", default=os.environ.get("GWAS_WORK_ROOT", "GWAS"), help="GWAS work root")
    p_purge.add_argument("--export-dir", default="", help="Optional finemapping export folder that must contain each protein before purge")
    p_purge.add_argument("--allow-partial", action="store_true", help="Purge completed rows even if some rows in the batch lack .qc.json")
    p_purge.set_defaults(func=cmd_purge)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
