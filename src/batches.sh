#!/usr/bin/env bash
# Split a large manifest into fixed batches, inspect batch progress, clean staged raw files, and purge batch VCFs.
#
# Usage:
#   bash src/batches.sh split     --manifest manifests/sun_combined.tsv --out-dir manifests/sun_batches [--size 100]
#   bash src/batches.sh status    --batches manifests/sun_batches [--export-dir /path/to/finemap_exports]
#   bash src/batches.sh clean-raw --batches manifests/sun_batches --batch 01 [--clean-temp]
#   bash src/batches.sh purge     --batches manifests/sun_batches --batch 01 [--export-dir /path/to/finemap_exports]

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=src/lib/env.sh
. "${ROOT_DIR}/src/lib/env.sh"
load_settings "${ROOT_DIR}"

for _var in GWAS_WORK_ROOT GWAS_NEXTFLOW_WORK GWAS_TMPDIR; do
  _val="${!_var:-}"
  if [[ -n "${_val}" && "${_val}" != /* ]]; then
    export "${_var}=${ROOT_DIR}/${_val}"
  fi
done
export GWAS_WORK_ROOT="${GWAS_WORK_ROOT:-${ROOT_DIR}/GWAS}"

exec python3 "${ROOT_DIR}/pipeline/bin/batches.py" "$@"
