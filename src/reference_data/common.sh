#!/usr/bin/env bash
# Shared setup for the reference-data scripts: repository root, settings, and GWAS_REFERENCE_ROOT.
set -euo pipefail

REFERENCE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${REFERENCE_SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=src/lib/env.sh
. "${REPO_ROOT}/src/lib/env.sh"
load_settings "${REPO_ROOT}"

resolve_path() {
  local value="$1"
  if [[ "${value}" = /* ]]; then
    printf '%s\n' "${value}"
  else
    printf '%s\n' "${REPO_ROOT}/${value}"
  fi
}

GWAS_REFERENCE_ROOT="$(resolve_path "${GWAS_REFERENCE_ROOT:-references}")"
require_real_path GWAS_REFERENCE_ROOT "${GWAS_REFERENCE_ROOT}"
export GWAS_REFERENCE_ROOT
mkdir -p "${GWAS_REFERENCE_ROOT}"
