#!/usr/bin/env bash
#SBATCH --job-name=GWASprocessing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=8G

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ -f "${ROOT_DIR}/.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  . "${ROOT_DIR}/.env"
  set +a
fi

usage() {
  cat <<'USAGE'
Usage: bash src/standardise.sh [options]
       sbatch src/standardise.sh [options]

Options:
  --manifest PATH       Manifest TSV. Default: $GWAS_MANIFEST or manifests/main.tsv
  --study-yaml PATH     Optional study-specific YAML overrides.
  --profile NAME        Nextflow profile. Default: $GWAS_PROFILE, slurm inside sbatch, otherwise local
  --only-study STUDYID  Run one derived study ID, e.g. locke_2015_25673413.
  --row-limit N         Eligible manifest rows to run. Default: 1 for local, 0/all for slurm.
  --force              Re-run rows even when persistent outputs already exist.
  --dry-run            Build the run manifest without downloading or standardising.
  --qcplot true|false  Create raw/standardised QQ and Manhattan TIFF. Default: true.
  -h, --help           Show this help.
USAGE
}

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
if [[ -n "${GWAS_PROFILE:-}" ]]; then
  PROFILE="${GWAS_PROFILE}"
elif [[ -n "${SLURM_JOB_ID:-}" ]]; then
  PROFILE="slurm"
else
  PROFILE="local"
fi
QCPLOT="${GWAS_QCPLOT:-true}"
ROW_LIMIT="${GWAS_ROW_LIMIT:-}"
NEXTFLOW_WORK="${GWAS_NEXTFLOW_WORK:-temp/nextflow_work}"
NEXTFLOW_LOG="${GWAS_NEXTFLOW_LOG:-temp/.nextflow.log}"
STUDY_YAML=""
ONLY_STUDY=""
FORCE="false"
DRY_RUN="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest)
      MANIFEST="$2"
      shift 2
      ;;
    --study-yaml)
      STUDY_YAML="$2"
      shift 2
      ;;
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    --only-study)
      ONLY_STUDY="$2"
      shift 2
      ;;
    --row-limit)
      ROW_LIMIT="$2"
      shift 2
      ;;
    --force)
      FORCE="true"
      shift
      ;;
    --dry-run)
      DRY_RUN="true"
      shift
      ;;
    --qcplot)
      QCPLOT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

case "${QCPLOT,,}" in
  true|false) ;;
  *) echo "--qcplot must be true or false" >&2; exit 2 ;;
esac
if [[ -z "${ROW_LIMIT}" ]]; then
  if [[ "${PROFILE}" == "local" ]]; then
    ROW_LIMIT="1"
  else
    ROW_LIMIT="0"
  fi
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

mkdir -p "${WORK_ROOT}" "${REFERENCE_ROOT}" "${TEMP_ROOT}" "${LOG_ROOT}" "${NEXTFLOW_WORK}" "$(dirname "${NEXTFLOW_LOG}")"

export GWAS_WORK_ROOT="${WORK_ROOT}"
export GWAS_REFERENCE_ROOT="${REFERENCE_ROOT}"
export GWAS_CONFIG_ROOT="${CONFIG_ROOT}"
export GWAS_LOG_ROOT="${LOG_ROOT}"

cmd=(
  nextflow
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
  --row_limit "${ROW_LIMIT}"
)

if [[ -n "${STUDY_YAML}" ]]; then
  cmd+=(--study_yaml "${STUDY_YAML}")
fi
if [[ -n "${ONLY_STUDY}" ]]; then
  cmd+=(--only_study "${ONLY_STUDY}")
fi

cd "${TEMP_ROOT}"
status=0
"${cmd[@]}" || status=$?

summary="${LOG_ROOT}/run_manifest_summary.json"
if [[ -f "${summary}" ]]; then
  python3 - "${summary}" <<'PY'
import json
import sys
from pathlib import Path

summary = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
print("")
print("GWASprocessing summary")
print(f"  Manifest: {summary.get('manifest', '')}")
print(f"  Rows to run: {summary.get('rows_to_run', 0)}")
print(f"  Eligible rows: {summary.get('rows_eligible', summary.get('rows_to_run', 0))}")
print(f"  Deferred by row limit: {summary.get('rows_deferred_due_to_limit', 0)}")
print(f"  Rows skipped complete: {summary.get('rows_skipped_complete', 0)}")
print(f"  Rows filtered: {summary.get('rows_filtered', 0)}")
print(f"  Row limit: {summary.get('row_limit', 0)}")
print(f"  State file: {summary.get('state_file', '')}")
if int(summary.get("rows_to_run", 0) or 0) == 0:
    if int(summary.get("rows_skipped_complete", 0) or 0) > 0:
        print("  Nothing launched: all selected manifest rows already have complete outputs.")
        print("  To rerun them, use: bash src/standardise.sh --force")
    else:
        print("  Nothing launched: no selected rows were found in the manifest.")
PY
fi

exit "${status}"
