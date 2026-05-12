#!/usr/bin/env bash
set -euo pipefail

REFERENCE_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${REFERENCE_SCRIPT_DIR}/../.." && pwd)"

if [[ -f "${REPO_ROOT}/.env" ]]; then
  set -a
  # shellcheck disable=SC1091
  . "${REPO_ROOT}/.env"
  set +a
fi

resolve_path() {
  local value="$1"
  if [[ "${value}" = /* ]]; then
    printf '%s\n' "${value}"
  else
    printf '%s\n' "${REPO_ROOT}/${value}"
  fi
}

GWAS_REFERENCE_ROOT="$(resolve_path "${GWAS_REFERENCE_ROOT:-references}")"
mkdir -p "${GWAS_REFERENCE_ROOT}"
