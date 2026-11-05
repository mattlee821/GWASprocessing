#!/usr/bin/env bash
# Install the tools GWASprocessing needs and build or download its reference data.
#
# Usage:
#   bash src/source.sh                  # tools + Python packages, then all reference data
#   bash src/source.sh --tools-only     # only install/verify tools and Python packages
#   bash src/source.sh --reference-only # only build/download reference data (tools must exist)
#   bash src/source.sh --force          # rebuild lookup tables and the rsID coordinate index
#
# Idempotent: every step skips work whose output already exists, so it is safe to re-run after an
# interruption. Tool locations (PLINK2, NEXTFLOW, GWAS_VENV) are written to server.env, which
# src/standardise.sh reads. Run it with bash; sourcing it only defines its functions.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFDATA="${ROOT_DIR}/src/reference_data"
TOOLS_DIR="${ROOT_DIR}/.tools/bin"
ENV_FILE="${ROOT_DIR}/server.env"
REQUIRED_MODULES="numpy, pandas, yaml, requests, matplotlib, pyliftover, certifi"
# shellcheck source=src/lib/env.sh
. "${ROOT_DIR}/src/lib/env.sh"

# set_env_var KEY VALUE: create or replace KEY=VALUE in server.env (other lines untouched) and export it.
set_env_var() {
  local key="$1" value="$2" tmp line found=0
  [[ -f "${ENV_FILE}" ]] || : > "${ENV_FILE}"
  tmp="${ENV_FILE}.tmp.$$"
  while IFS= read -r line || [[ -n "${line}" ]]; do
    if [[ "${line}" == "${key}="* ]]; then
      printf '%s=%s\n' "${key}" "${value}"
      found=1
    else
      printf '%s\n' "${line}"
    fi
  done < "${ENV_FILE}" > "${tmp}"
  if [[ "${found}" -eq 0 ]]; then
    printf '%s=%s\n' "${key}" "${value}" >> "${tmp}"
  fi
  mv "${tmp}" "${ENV_FILE}"
  export "${key}=${value}"
}

install_plink2() {
  if command -v "${PLINK2:-plink2}" >/dev/null 2>&1 && "${PLINK2:-plink2}" --version >/dev/null 2>&1; then
    echo "plink2: found ($(command -v "${PLINK2:-plink2}"))"
    return
  fi
  local url
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64) url="${PLINK2_URL:-https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260504.zip}" ;;
    Darwin-x86_64) url="${PLINK2_URL:-https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_avx2_20260504.zip}" ;;
    Linux-x86_64) url="${PLINK2_URL:-https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260504.zip}" ;;
    Linux-aarch64|Linux-arm64) url="${PLINK2_URL:-https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_arm64_20260504.zip}" ;;
    *) echo "No known plink2 build for $(uname -s)-$(uname -m). Set PLINK2=/path/to/plink2 in server.env." >&2; return 1 ;;
  esac
  echo "plink2: installing from ${url}"
  mkdir -p "${TOOLS_DIR}"
  curl -fsSL -o "${TOOLS_DIR}/plink2.zip" "${url}"
  (cd "${TOOLS_DIR}" && unzip -o plink2.zip plink2 >/dev/null)
  chmod +x "${TOOLS_DIR}/plink2"
  rm -f "${TOOLS_DIR}/plink2.zip"
  xattr -d com.apple.quarantine "${TOOLS_DIR}/plink2" 2>/dev/null || true
  set_env_var PLINK2 "${TOOLS_DIR}/plink2"
}

install_nextflow() {
  if [[ -n "${NEXTFLOW:-}" && -x "${NEXTFLOW}" ]]; then
    echo "nextflow: found (${NEXTFLOW})"
  elif command -v nextflow >/dev/null 2>&1; then
    echo "nextflow: found ($(command -v nextflow))"
  elif [[ -x "${TOOLS_DIR}/nextflow" ]]; then
    set_env_var NEXTFLOW "${TOOLS_DIR}/nextflow"
    echo "nextflow: found (${TOOLS_DIR}/nextflow)"
  else
    echo "nextflow: installing into ${TOOLS_DIR}"
    mkdir -p "${TOOLS_DIR}"
    (cd "${TOOLS_DIR}" && curl -fsSL https://get.nextflow.io | bash >/dev/null)
    chmod +x "${TOOLS_DIR}/nextflow"
    set_env_var NEXTFLOW "${TOOLS_DIR}/nextflow"
  fi
}

install_htslib() {
  if command -v bgzip >/dev/null 2>&1 && command -v tabix >/dev/null 2>&1; then
    echo "bgzip/tabix: found ($(command -v bgzip))"
    return
  fi
  echo "bgzip/tabix: not found, attempting install"
  if command -v brew >/dev/null 2>&1; then
    brew install htslib
  elif command -v mamba >/dev/null 2>&1; then
    mamba install -y -c bioconda -c conda-forge htslib
  elif command -v conda >/dev/null 2>&1; then
    conda install -y -c bioconda -c conda-forge htslib
  elif command -v apt-get >/dev/null 2>&1 && sudo -n true 2>/dev/null; then
    sudo apt-get update -y && sudo apt-get install -y tabix
  elif command -v yum >/dev/null 2>&1 && sudo -n true 2>/dev/null; then
    sudo yum install -y htslib
  fi
  if ! command -v bgzip >/dev/null 2>&1 || ! command -v tabix >/dev/null 2>&1; then
    echo "bgzip/tabix are required. Install htslib (e.g. 'module load htslib' on a cluster) and re-run." >&2
    return 1
  fi
}

install_python_packages() {
  local python="python3"
  if [[ -n "${GWAS_VENV:-}" && -x "${ROOT_DIR}/${GWAS_VENV#./}/bin/python3" ]]; then
    python="${ROOT_DIR}/${GWAS_VENV#./}/bin/python3"
  fi
  if "${python}" -c "import ${REQUIRED_MODULES}" 2>/dev/null; then
    echo "python packages: already available (${python})"
    return
  fi
  echo "python packages: installing from requirements.txt"
  if [[ -z "${GWAS_VENV:-}" ]] && python3 -m pip install --quiet -r "${ROOT_DIR}/requirements.txt"; then
    return
  fi
  # A system Python that refuses pip installs (PEP 668) or a configured GWAS_VENV: use a project virtualenv.
  local venv="${GWAS_VENV:-.venv}"
  venv="${venv#./}"
  if [[ ! -x "${ROOT_DIR}/${venv}/bin/python3" ]]; then
    python3 -m venv "${ROOT_DIR}/${venv}"
  fi
  "${ROOT_DIR}/${venv}/bin/python3" -m pip install --quiet -r "${ROOT_DIR}/requirements.txt"
  set_env_var GWAS_VENV "${venv}"
  echo "python packages: installed into ${ROOT_DIR}/${venv} (src/standardise.sh puts it first on PATH)"
}

build_reference_data() {
  local force_flag=()
  [[ "${FORCE}" -eq 1 ]] && force_flag=(--force)
  bash "${REFDATA}/download_1000genomes_phase3.sh"
  bash "${REFDATA}/build_1000genomes_phase3.sh"
  python3 "${REFDATA}/make_1000genomes_lookup_tables.py" --reference-root "${GWAS_REFERENCE_ROOT:-${ROOT_DIR}/references}" ${force_flag[@]+"${force_flag[@]}"}
  bash "${REFDATA}/download_fasta_and_chains.sh"
  python3 "${REFDATA}/validate_references.py" --reference-root "${GWAS_REFERENCE_ROOT:-${ROOT_DIR}/references}"
}

main() {
  set -euo pipefail
  TOOLS_ONLY=0
  REFERENCE_ONLY=0
  FORCE=0
  local arg
  for arg in "$@"; do
    case "${arg}" in
      --tools-only) TOOLS_ONLY=1 ;;
      --reference-only) REFERENCE_ONLY=1 ;;
      --force) FORCE=1 ;;
      -h|--help) sed -n '2,13p' "${BASH_SOURCE[0]}"; exit 0 ;;
      *) echo "Unknown option: ${arg}" >&2; exit 2 ;;
    esac
  done
  load_settings "${ROOT_DIR}"
  echo "== GWASprocessing setup =="
  if [[ "${REFERENCE_ONLY}" -eq 0 ]]; then
    echo "-- Tools --"
    install_plink2
    install_nextflow
    install_htslib
    install_python_packages
  fi
  if [[ "${TOOLS_ONLY}" -eq 0 ]]; then
    echo "-- Reference data --"
    build_reference_data
  fi
  echo "== Setup complete =="
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  main "$@"
fi
