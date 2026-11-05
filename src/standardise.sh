#!/usr/bin/env bash
#SBATCH --job-name=GWASprocessing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=8G
# Run the GWASprocessing Nextflow workflow over a manifest.
#
# Usage: bash src/standardise.sh [options]      (from any directory)
#        sbatch src/standardise.sh [options]    (submit from the repository root; uses the slurm profile)
# Settings come from command-line options, then environment variables, then .env (see .env.example).
# Exit status: Nextflow's status if it failed, 1 if any manifest row failed, otherwise 0.

set -euo pipefail

find_root() {
  if [[ -n "${GWASPROCESSING_ROOT:-}" ]]; then
    (cd "${GWASPROCESSING_ROOT}" && pwd)
    return
  fi
  local script_root candidate
  script_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  if [[ -f "${script_root}/pipeline/main.nf" ]]; then
    printf '%s\n' "${script_root}"
    return
  fi
  # sbatch executes a spooled copy of this script, so fall back to the submission directory.
  for candidate in "${SLURM_SUBMIT_DIR:-}" "${PWD}"; do
    if [[ -n "${candidate}" && -f "${candidate}/pipeline/main.nf" ]]; then
      (cd "${candidate}" && pwd)
      return
    fi
  done
  echo "Cannot find the GWASprocessing repository (pipeline/main.nf). Run from the repository root or set GWASPROCESSING_ROOT." >&2
  exit 2
}

lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }

usage() {
  cat <<'USAGE'
Usage: bash src/standardise.sh [options]
       sbatch src/standardise.sh [options]

Options:
  --manifest PATH       Manifest TSV. Default: $GWAS_MANIFEST or manifests/main.tsv
  --study-yaml PATH     Optional study-specific YAML overrides applied to every row.
  --profile NAME        Nextflow profile. Default: slurm inside an sbatch job, else $GWAS_PROFILE, else local
  --only-study STUDYID  Run one derived study ID, e.g. locke_2015_25673413.
  --row-limit N         Eligible manifest rows to run. Default: 1 for local, 0 (all) for slurm.
  --force               Re-run rows even when complete outputs exist for the current settings.
  --dry-run             Build the run manifest without downloading or standardising.
  --qcplot true|false   Create the raw/standardised QQ plot TIFF. Default: true.
  --analysis-qc true|false
                        Apply analysis QC (biallelic SNPs, palindromic strand resolution, MAF>=0.01,
                        positional de-duplication) before writing. Default: false.
  -h, --help            Show this help.

Environment (or .env): GWAS_WORK_ROOT, GWAS_REFERENCE_ROOT, GWAS_CONFIG_ROOT, GWAS_TEMP_ROOT, GWAS_LOG_ROOT,
GWAS_NEXTFLOW_WORK, GWAS_NEXTFLOW_LOG, GWAS_PROFILE, GWAS_QCPLOT, GWAS_ANALYSIS_QC, GWAS_ROW_LIMIT,
GWAS_VENV (Python virtual environment to put first on PATH), NEXTFLOW (path to the nextflow launcher).
USAGE
}

ROOT_DIR="$(find_root)"
# shellcheck source=src/lib/env.sh
. "${ROOT_DIR}/src/lib/env.sh"
load_env_defaults "${ROOT_DIR}/.env"

resolve_path() {
  local value="$1"
  if [[ "${value}" =~ ^[a-zA-Z][a-zA-Z0-9+.-]*:// ]] || [[ "${value}" = /* ]]; then
    printf '%s\n' "${value}"
  else
    printf '%s\n' "${ROOT_DIR}/${value}"
  fi
}

MANIFEST="${GWAS_MANIFEST:-manifests/main.tsv}"
WORK_ROOT="${GWAS_WORK_ROOT:-GWAS}"
REFERENCE_ROOT="${GWAS_REFERENCE_ROOT:-references}"
CONFIG_ROOT="${GWAS_CONFIG_ROOT:-manifests}"
TEMP_ROOT="${GWAS_TEMP_ROOT:-temp}"
LOG_ROOT="${GWAS_LOG_ROOT:-temp/logs}"
QCPLOT="${GWAS_QCPLOT:-true}"
ANALYSIS_QC="${GWAS_ANALYSIS_QC:-false}"
ROW_LIMIT="${GWAS_ROW_LIMIT:-}"
NEXTFLOW_WORK="${GWAS_NEXTFLOW_WORK:-temp/nextflow_work}"
NEXTFLOW_LOG="${GWAS_NEXTFLOW_LOG:-temp/.nextflow.log}"
PROFILE=""
STUDY_YAML=""
ONLY_STUDY=""
FORCE="false"
DRY_RUN="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --study-yaml) STUDY_YAML="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --only-study) ONLY_STUDY="$2"; shift 2 ;;
    --row-limit) ROW_LIMIT="$2"; shift 2 ;;
    --force) FORCE="true"; shift ;;
    --dry-run) DRY_RUN="true"; shift ;;
    --qcplot) QCPLOT="$2"; shift 2 ;;
    --analysis-qc) ANALYSIS_QC="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z "${PROFILE}" ]]; then
  if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    PROFILE="slurm"
  else
    PROFILE="${GWAS_PROFILE:-local}"
  fi
fi
QCPLOT="$(lower "${QCPLOT}")"
ANALYSIS_QC="$(lower "${ANALYSIS_QC}")"
for pair in "--qcplot:${QCPLOT}" "--analysis-qc:${ANALYSIS_QC}"; do
  case "${pair#*:}" in
    true|false) ;;
    *) echo "${pair%%:*} must be true or false" >&2; exit 2 ;;
  esac
done
if [[ -z "${ROW_LIMIT}" ]]; then
  if [[ "${PROFILE}" == "local" ]]; then ROW_LIMIT="1"; else ROW_LIMIT="0"; fi
fi
if ! [[ "${ROW_LIMIT}" =~ ^[0-9]+$ ]]; then
  echo "--row-limit must be a non-negative integer" >&2
  exit 2
fi

MANIFEST="$(resolve_path "${MANIFEST}")"
WORK_ROOT="$(resolve_path "${WORK_ROOT}")"
REFERENCE_ROOT="$(resolve_path "${REFERENCE_ROOT}")"
CONFIG_ROOT="$(resolve_path "${CONFIG_ROOT}")"
TEMP_ROOT="$(resolve_path "${TEMP_ROOT}")"
LOG_ROOT="$(resolve_path "${LOG_ROOT}")"
NEXTFLOW_WORK="$(resolve_path "${NEXTFLOW_WORK}")"
NEXTFLOW_LOG="$(resolve_path "${NEXTFLOW_LOG}")"
if [[ -n "${STUDY_YAML}" ]]; then
  STUDY_YAML="$(resolve_path "${STUDY_YAML}")"
fi

for path in "${WORK_ROOT}" "${REFERENCE_ROOT}" "${TEMP_ROOT}" "${NEXTFLOW_WORK}"; do
  case "${path}" in
    */Library/CloudStorage/*|*OneDrive*|*Dropbox*|*"Google Drive"*|*iCloud*)
      echo "Warning: ${path} is in a cloud-synced folder. Large outputs and Nextflow work files will be synced and" \
           "cloud-only reference files re-downloaded on every run; set GWAS_WORK_ROOT, GWAS_REFERENCE_ROOT," \
           "GWAS_TEMP_ROOT and GWAS_NEXTFLOW_WORK to local or scratch storage." >&2
      break ;;
  esac
done

if [[ -n "${GWAS_VENV:-}" ]]; then
  VENV_DIR="$(resolve_path "${GWAS_VENV}")"
  if [[ -x "${VENV_DIR}/bin/python3" ]]; then
    export PATH="${VENV_DIR}/bin:${PATH}"
  else
    echo "Warning: GWAS_VENV=${GWAS_VENV} has no bin/python3; using python3 from PATH" >&2
  fi
fi

if [[ -n "${NEXTFLOW:-}" ]]; then
  NEXTFLOW_BIN="${NEXTFLOW}"
elif command -v nextflow >/dev/null 2>&1; then
  NEXTFLOW_BIN="$(command -v nextflow)"
elif [[ -x "${ROOT_DIR}/.tools/bin/nextflow" ]]; then
  NEXTFLOW_BIN="${ROOT_DIR}/.tools/bin/nextflow"
else
  echo "nextflow not found. Run 'bash src/source.sh --tools-only' or set NEXTFLOW=/path/to/nextflow in .env" >&2
  exit 2
fi

mkdir -p "${WORK_ROOT}" "${REFERENCE_ROOT}" "${TEMP_ROOT}" "${LOG_ROOT}" "${NEXTFLOW_WORK}" "$(dirname "${NEXTFLOW_LOG}")"
# Summaries from an earlier run must never be reported as this run's results.
rm -f "${LOG_ROOT}/run_manifest_summary.json" "${LOG_ROOT}/run_manifest.tsv" "${LOG_ROOT}/state_update_summary.json"

export GWAS_WORK_ROOT="${WORK_ROOT}" GWAS_REFERENCE_ROOT="${REFERENCE_ROOT}" GWAS_CONFIG_ROOT="${CONFIG_ROOT}" GWAS_LOG_ROOT="${LOG_ROOT}"

cmd=(
  "${NEXTFLOW_BIN}"
  -log "${NEXTFLOW_LOG}"
  run "${ROOT_DIR}/pipeline/main.nf"
  -profile "${PROFILE}"
  -work-dir "${NEXTFLOW_WORK}"
  --repo_root "${ROOT_DIR}"
  --pipeline_root "${ROOT_DIR}/pipeline"
  --manifest "${MANIFEST}"
  --work_root "${WORK_ROOT}"
  --reference_root "${REFERENCE_ROOT}"
  --config_root "${CONFIG_ROOT}"
  --log_root "${LOG_ROOT}"
  --force "${FORCE}"
  --dry_run "${DRY_RUN}"
  --qcplot "${QCPLOT}"
  --analysis_qc "${ANALYSIS_QC}"
  --row_limit "${ROW_LIMIT}"
)
if [[ -n "${STUDY_YAML}" ]]; then cmd+=(--study_yaml "${STUDY_YAML}"); fi
if [[ -n "${ONLY_STUDY}" ]]; then cmd+=(--only_study "${ONLY_STUDY}"); fi

cd "${TEMP_ROOT}"
status=0
"${cmd[@]}" || status=$?

report_status=0
python3 - "${LOG_ROOT}" "${DRY_RUN}" <<'PY' || report_status=$?
import csv
import json
import sys
from pathlib import Path

log_root, dry_run = Path(sys.argv[1]), sys.argv[2] == "true"
summary_path = log_root / "run_manifest_summary.json"
if not summary_path.exists():
    print("\nGWASprocessing: the run stopped before the manifest was prepared; see the Nextflow log.")
    raise SystemExit(0)
summary = json.loads(summary_path.read_text(encoding="utf-8"))
state_path = log_root / "state_update_summary.json"
state = json.loads(state_path.read_text(encoding="utf-8")) if state_path.exists() else {}
print("")
print("GWASprocessing summary")
print(f"  Manifest: {summary.get('manifest', '')}")
print(f"  Rows to run: {summary.get('rows_to_run', 0)}")
print(f"  Eligible rows: {summary.get('rows_eligible', 0)}")
print(f"  Deferred by row limit: {summary.get('rows_deferred_due_to_limit', 0)}")
print(f"  Rows skipped (complete for current settings): {summary.get('rows_skipped_complete', 0)}")
print(f"  Rows filtered: {summary.get('rows_filtered', 0)}")
print(f"  Rows completed: {state.get('updated', 0)}")
print(f"  State file: {summary.get('state_file', '')}")
if int(summary.get("rows_to_run", 0) or 0) == 0:
    if int(summary.get("rows_skipped_complete", 0) or 0) > 0:
        print("  Nothing launched: all selected rows already have complete outputs for the current settings.")
        print("  To rerun them anyway, use: bash src/standardise.sh --force")
    else:
        print("  Nothing launched: no selected rows were found in the manifest.")
    raise SystemExit(0)
if dry_run:
    raise SystemExit(0)
done = set(state.get("updated_row_hashes", []))
run_manifest = log_root / "run_manifest.tsv"
rows = list(csv.DictReader(open(run_manifest, encoding="utf-8"), delimiter="\t")) if run_manifest.exists() else []
failed = [r for r in rows if r.get("row_hash") not in done]
if failed:
    print(f"  {len(failed)} row(s) failed:")
    for r in failed:
        print(f"    {r.get('study_id')}/{r.get('output_id')}")
    print("  Details are in the Nextflow log; fix the inputs and re-run (completed rows are skipped).")
    raise SystemExit(1)
PY

if [[ "${status}" -ne 0 ]]; then
  exit "${status}"
fi
exit "${report_status}"
