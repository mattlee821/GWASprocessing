#!/usr/bin/env bash
# Install the tools GWASprocessing needs and build or download its reference data.
#
# Usage:
#   bash src/source.sh                  # tools + Python packages, then all reference data
#   bash src/source.sh --tools-only     # only install/verify tools and Python packages
#   bash src/source.sh --reference-only # only build/download reference data (tools must exist)
#   bash src/source.sh --force          # rebuild the built reference lookups (dbSNP index and merge map, 1000 Genomes, build sample)
#   bash src/source.sh --r-packages     # only install the R packages for the R audit test (gwasvcf, VariantAnnotation)
#   bash src/source.sh --env-only       # only create or complete .env
#
# Every run first creates .env from .env.example if it is missing and fills in the data paths; tool locations
# (PLINK2, NEXTFLOW, GWAS_VENV, GWAS_BCFTOOLS, BCFTOOLS_PLUGINS) are written to .env as they are installed. Values
# already in .env are kept. The only setting to add by hand is OPENGWAS_JWT.
# Idempotent: every step skips work whose output already exists, so it is safe to re-run after an
# interruption. Run it with bash; sourcing it only defines its functions.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFDATA="${ROOT_DIR}/src/reference_data"
TOOLS_DIR="${ROOT_DIR}/.tools/bin"
ENV_FILE="${ROOT_DIR}/.env"
REQUIRED_MODULES="numpy, pandas, yaml, requests, matplotlib, certifi"
# shellcheck source=src/lib/env.sh
. "${ROOT_DIR}/src/lib/env.sh"

# set_env_var KEY VALUE: create or replace KEY=VALUE in .env (other lines untouched) and export it.
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

# Paths filled into a new .env when they are still placeholders or empty, relative to the data home.
ENV_DEFAULT_PATHS=(GWAS_WORK_ROOT=GWAS GWAS_REFERENCE_ROOT=references GWAS_TEMP_ROOT=temp GWAS_LOG_ROOT=temp/logs
                   GWAS_NEXTFLOW_WORK=temp/nextflow_work GWAS_NEXTFLOW_LOG=temp/.nextflow.log)

# ensure_env_file: create .env from .env.example when it is missing (moving the values of an old server.env into
# it), then fill each path still a /WHERE/... placeholder or empty: from the environment when set there, otherwise
# under the data home (GWAS_DATA_HOME; default the repository, or ~/GWASprocessing-data when the repository is in
# a cloud-synced folder). Values already in .env are never changed.
ensure_env_file() {
  local template="${ROOT_DIR}/.env.example" old="${ROOT_DIR}/server.env" line key value home pair filled=() moved=0
  if [[ ! -f "${ENV_FILE}" ]]; then
    if [[ ! -f "${template}" ]]; then
      echo "${template} not found: cannot create .env" >&2
      return 1
    fi
    cp "${template}" "${ENV_FILE}"
    echo "Created .env from .env.example"
    if [[ -f "${old}" ]]; then
      while IFS= read -r line || [[ -n "${line}" ]]; do
        [[ "${line}" =~ ^[[:space:]]*(#|$) ]] && continue
        [[ "${line}" == *=* ]] || continue
        key="${line%%=*}"; key="${key#export }"; key="${key//[[:space:]]/}"
        [[ "${key}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
        value="${line#*=}"; value="${value%\"}"; value="${value#\"}"; value="${value%\'}"; value="${value#\'}"
        set_env_var "${key}" "${value}" >/dev/null
        moved=$((moved + 1))
      done < "${old}"
      echo "Moved ${moved} settings from server.env into .env; server.env is no longer read and can be deleted"
    fi
  fi
  home="${GWAS_DATA_HOME:-${ROOT_DIR}}"
  case "${home}" in
    */Library/CloudStorage/*|*OneDrive*|*Dropbox*|*"Google Drive"*|*iCloud*)
      [[ -n "${GWAS_DATA_HOME:-}" ]] || home="${HOME}/GWASprocessing-data" ;;
  esac
  for pair in "${ENV_DEFAULT_PATHS[@]}"; do
    key="${pair%%=*}"
    value="$(env_value "${ENV_FILE}" "${key}")"
    if [[ -z "${value}" || "${value}" == /WHERE/* ]]; then
      if [[ -n "${!key:-}" && "${!key}" != /WHERE/* ]]; then
        set_env_var "${key}" "${!key}"
      else
        set_env_var "${key}" "${home}/${pair#*=}"
      fi
      filled+=("${key}")
    fi
  done
  if [[ "${#filled[@]}" -gt 0 ]]; then
    echo "Filled in .env: ${filled[*]} (under ${home} unless set in the environment); edit .env to move them"
  fi
}

# remind_token: say where the OpenGWAS token goes while .env has none.
remind_token() {
  if [[ -z "$(env_value "${ENV_FILE}" OPENGWAS_JWT)" ]]; then
    echo "Next: add your OpenGWAS token to .env as OPENGWAS_JWT=<token> (needed only for OpenGWAS sources;" \
         "generate one at https://api.opengwas.io/profile/)."
  fi
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
    *) echo "No known plink2 build for $(uname -s)-$(uname -m). Set PLINK2=/path/to/plink2 in .env." >&2; return 1 ;;
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

# BCFtools/liftover (Genovese et al. 2024) is a bcftools plugin that is not packaged with bcftools: build bcftools with
# it from source into GWAS_TOOLS_ROOT (default .tools), then record GWAS_BCFTOOLS and BCFTOOLS_PLUGINS in .env.
# It is required: GRCh37 inputs are lifted with it, and the build-detection sample is built with it.
BCFTOOLS_VERSION="${BCFTOOLS_VERSION:-1.25}"
install_bcftools_liftover() {
  if [[ -n "${GWAS_BCFTOOLS:-}" && -x "${GWAS_BCFTOOLS}" ]] && BCFTOOLS_PLUGINS="${BCFTOOLS_PLUGINS:-}" "${GWAS_BCFTOOLS}" plugin -l 2>/dev/null | grep -qx liftover; then
    echo "bcftools +liftover: found (${GWAS_BCFTOOLS})"
    return
  fi
  if command -v bcftools >/dev/null 2>&1 && bcftools plugin -l 2>/dev/null | grep -qx liftover; then
    echo "bcftools +liftover: found on PATH ($(command -v bcftools))"
    return
  fi
  local root="${GWAS_TOOLS_ROOT:-${ROOT_DIR}/.tools}" dir
  dir="${root}/bcftools-${BCFTOOLS_VERSION}"
  echo "bcftools +liftover: building bcftools ${BCFTOOLS_VERSION} with the liftover plugin in ${dir}"
  if ! { mkdir -p "${root}" \
        && { [[ -d "${dir}" ]] || { curl -fsSL -o "${root}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
               "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
             && tar -xjf "${root}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" -C "${root}" \
             && rm -f "${root}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"; }; } \
        && curl -fsSL -o "${dir}/plugins/liftover.c" https://raw.githubusercontent.com/freeseek/score/master/liftover.c \
        && (cd "${dir}" && ./configure --disable-libcurl --disable-lzma --disable-bz2 --without-curses >/dev/null \
            && (cd "htslib-${BCFTOOLS_VERSION}" && ./configure --disable-libcurl --disable-lzma --disable-bz2 >/dev/null) \
            && make -j4 >/dev/null); } 2>/dev/null; then
    echo "Error: could not build bcftools with the liftover plugin (needs a C compiler, make and zlib). Install them" \
         "and re-run, or set GWAS_BCFTOOLS and BCFTOOLS_PLUGINS in .env to an existing build." >&2
    return 1
  fi
  set_env_var GWAS_BCFTOOLS "${dir}/bcftools"
  set_env_var BCFTOOLS_PLUGINS "${dir}/plugins"
  echo "bcftools +liftover: built (${dir}/bcftools)"
}

# R packages are needed only by pipeline/tests/test_audit.py, which checks that R's gwasvcf reads the output.
install_r_packages() {
  if ! command -v Rscript >/dev/null 2>&1; then
    echo "R packages: Rscript not found; install R (https://cran.r-project.org) first" >&2
    return 1
  fi
  if Rscript -e 'quit(status = !all(sapply(c("gwasvcf", "VariantAnnotation", "GenomeInfoDb"), requireNamespace, quietly = TRUE)))' >/dev/null 2>&1; then
    echo "R packages: gwasvcf and VariantAnnotation already available"
    return
  fi
  echo "R packages: installing VariantAnnotation (Bioconductor) and gwasvcf (github.com/MRCIEU/gwasvcf)"
  Rscript -e '
    repos <- "https://cloud.r-project.org"
    for (p in c("BiocManager", "remotes")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = repos)
    # gwasvcf imports GenomeInfoDb, which Bioconductor 3.23 no longer installs alongside VariantAnnotation
    BiocManager::install(c("VariantAnnotation", "GenomeInfoDb"), update = FALSE, ask = FALSE)
    remotes::install_github("MRCIEU/gwasvcf", upgrade = "never")
    stopifnot(requireNamespace("gwasvcf", quietly = TRUE))'
}

build_reference_data() {
  [[ "${FORCE}" -eq 1 ]] && export GWAS_REFERENCE_FORCE=1
  bash "${REFDATA}/download_fasta_and_chains.sh"
  bash "${REFDATA}/build_grch38_references.sh"
  python3 "${REFDATA}/validate_references.py" --reference-root "${GWAS_REFERENCE_ROOT:-${ROOT_DIR}/references}"
}

main() {
  set -euo pipefail
  TOOLS_ONLY=0
  REFERENCE_ONLY=0
  FORCE=0
  R_ONLY=0
  ENV_ONLY=0
  local arg
  for arg in "$@"; do
    case "${arg}" in
      --tools-only) TOOLS_ONLY=1 ;;
      --reference-only) REFERENCE_ONLY=1 ;;
      --force) FORCE=1 ;;
      --r-packages) R_ONLY=1 ;;
      --env-only) ENV_ONLY=1 ;;
      -h|--help) sed -n '2,18p' "${BASH_SOURCE[0]}"; exit 0 ;;
      *) echo "Unknown option: ${arg}" >&2; exit 2 ;;
    esac
  done
  ensure_env_file
  load_settings "${ROOT_DIR}"
  if [[ "${ENV_ONLY}" -eq 1 ]]; then
    remind_token
    return
  fi
  if [[ "${R_ONLY}" -eq 1 ]]; then
    install_r_packages
    return
  fi
  echo "== GWASprocessing setup =="
  if [[ "${REFERENCE_ONLY}" -eq 0 ]]; then
    echo "-- Tools --"
    install_plink2
    install_nextflow
    install_htslib
    install_python_packages
    install_bcftools_liftover
  fi
  if [[ "${TOOLS_ONLY}" -eq 0 ]]; then
    echo "-- Reference data --"
    build_reference_data
  fi
  echo "== Setup complete =="
  remind_token
}

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  main "$@"
fi
