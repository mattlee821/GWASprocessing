#!/usr/bin/env bash
# Install the tools GWASprocessing needs and build or download its reference data.
#
# Usage:
#   bash src/source.sh                       # tools, then all reference data
#   bash src/source.sh --tools-only          # only the tool environment (and bcftools with the liftover plugin)
#   bash src/source.sh --reference-only      # only download and build the reference data (tools must exist)
#   bash src/source.sh --container apptainer # run the pipeline steps in the published image (or: docker)
#   bash src/source.sh --force               # rebuild the built references (and re-pull the container image)
#   bash src/source.sh --r-packages          # only the R packages for the R reader test (gwasvcf, VariantAnnotation)
#   bash src/source.sh --env-only            # only create or complete .env
#
# Needs only bash, curl and tar. Every run first creates .env from .env.example if it is missing and fills in the
# data paths; the tool environment (GWAS_ENV), NEXTFLOW, PLINK2, GWAS_BCFTOOLS and BCFTOOLS_PLUGINS are written to
# .env as they are installed. Values already in .env are kept. The only setting to add by hand is OPENGWAS_JWT.
# Idempotent: every step skips work whose output already exists, so it is safe to re-run after an
# interruption. Run it with bash; sourcing it only defines its functions.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REFDATA="${ROOT_DIR}/src/reference_data"
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

# ensure_env_file: create .env from .env.example when it is missing, then fill each path still a /WHERE/...
# placeholder or empty: from the environment when set there, otherwise
# under the data home (GWAS_DATA_HOME; default the repository, or ~/GWASprocessing-data when the repository is in
# a cloud-synced folder). Values already in .env are never changed.
ensure_env_file() {
  local template="${ROOT_DIR}/.env.example" home pair filled=()
  if [[ ! -f "${ENV_FILE}" ]]; then
    if [[ ! -f "${template}" ]]; then
      echo "${template} not found: cannot create .env" >&2
      return 1
    fi
    cp "${template}" "${ENV_FILE}"
    echo "Created .env from .env.example"
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

# ---- Tools ------------------------------------------------------------------------------------------------------
# Everything the pipeline runs (Java, Nextflow, Python and its packages, bgzip/tabix, plink2, and bcftools with the
# BCFtools/liftover plugin) is installed into 1 private environment under GWAS_TOOLS_ROOT (default .tools), from
# environment.yml, with micromamba (1 static binary, no administrator rights). Nothing outside GWAS_TOOLS_ROOT is
# changed. The scripts put the environment first on PATH (src/lib/env.sh, activate_tools). A tool path the user has
# set in .env (NEXTFLOW, PLINK2, GWAS_BCFTOOLS, BCFTOOLS_PLUGINS) is kept.
BCFTOOLS_VERSION="${BCFTOOLS_VERSION:-1.25}"
LIFTOVER_URL="${LIFTOVER_URL:-https://raw.githubusercontent.com/freeseek/score/master/liftover.c}"
CONTAINER_REF="${GWAS_CONTAINER_REF:-ghcr.io/mattlee821/gwasprocessing:latest}"

tools_root() {
  local root="${GWAS_TOOLS_ROOT:-.tools}"
  [[ "${root}" == /* ]] || root="${ROOT_DIR}/${root}"
  printf '%s\n' "${root}"
}

mamba_platform() {
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64) echo osx-arm64 ;;
    Darwin-x86_64) echo osx-64 ;;
    Linux-x86_64) echo linux-64 ;;
    Linux-aarch64|Linux-arm64) echo linux-aarch64 ;;
    *) echo "Unsupported platform $(uname -s)-$(uname -m): install the tools in environment.yml yourself and set" \
            "GWAS_ENV, NEXTFLOW and GWAS_BCFTOOLS in .env." >&2; return 1 ;;
  esac
}

# set_tool KEY PATH: record a tool in .env unless the user already points KEY at a working program.
set_tool() {
  local current
  current="$(env_value "${ENV_FILE}" "$1")"
  if [[ -n "${current}" && "${current}" != "$2" && -x "${current}" ]]; then
    echo "$1: keeping ${current} from .env"
    return
  fi
  set_env_var "$1" "$2"
}

install_micromamba() {
  local root bin platform
  root="$(tools_root)"
  bin="${root}/bin/micromamba"
  if [[ ! -x "${bin}" ]]; then
    platform="$(mamba_platform)"
    echo "micromamba: downloading for ${platform} into ${root}/bin"
    mkdir -p "${root}"
    curl -fsSL "https://micro.mamba.pm/api/micromamba/${platform}/latest" | tar -xj -C "${root}" bin/micromamba
    chmod +x "${bin}"
  fi
  MAMBA_BIN="${bin}"
  export MAMBA_ROOT_PREFIX="${root}/mamba"  # package cache, kept inside the tools folder
}

file_hash() { (shasum -a 256 "$1" 2>/dev/null || sha256sum "$1") | awk '{print $1}'; }

# install_environment: create or update the tool environment from environment.yml; skipped while environment.yml is
# unchanged since the last install.
install_environment() {
  local root prefix spec stamp hash
  root="$(tools_root)"
  prefix="${root}/env"
  spec="${ROOT_DIR}/environment.yml"
  stamp="${prefix}/.environment.yml.sha256"
  hash="$(file_hash "${spec}")"
  install_micromamba
  if [[ -f "${stamp}" && "$(cat "${stamp}")" == "${hash}" ]]; then
    echo "tool environment: up to date (${prefix})"
  elif [[ -d "${prefix}/conda-meta" ]]; then
    echo "tool environment: updating ${prefix} from environment.yml"
    "${MAMBA_BIN}" install -y -q -p "${prefix}" -f "${spec}"
  else
    echo "tool environment: creating ${prefix} from environment.yml (several minutes the first time)"
    "${MAMBA_BIN}" create -y -q -p "${prefix}" -f "${spec}"
  fi
  printf '%s\n' "${hash}" > "${stamp}"
  ENV_PREFIX="${prefix}"
  set_env_var GWAS_ENV "${prefix}"
  activate_tools "${ROOT_DIR}"
  set_tool NEXTFLOW "${prefix}/bin/nextflow"
  if [[ -z "${GWAS_CONTAINER:-}" ]]; then
    ensure_plink2 "${prefix}"
    set_tool PLINK2 "${prefix}/bin/plink2"
  fi
}

# ensure_plink2 PREFIX: some conda plink2 builds link to libraries outside the environment (the osx-64 build looks for
# Homebrew's zstd); when the environment's plink2 does not run, put the official static build in its place.
ensure_plink2() {
  local prefix="$1" url
  ( "${prefix}/bin/plink2" --version ) >/dev/null 2>&1 && return 0  # subshell: a crashing binary prints nothing
  case "$(uname -s)-$(uname -m)" in
    Darwin-arm64) url="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260504.zip" ;;
    Darwin-x86_64) url="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_avx2_20260504.zip" ;;
    Linux-x86_64) url="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_avx2_20260504.zip" ;;
    Linux-aarch64|Linux-arm64) url="https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_arm64_20260504.zip" ;;
    *) echo "plink2 from the environment does not run; set PLINK2 in .env" >&2; return 1 ;;
  esac
  url="${PLINK2_URL:-${url}}"
  echo "plink2: the environment's build does not run here; installing the official build from ${url}"
  curl -fsSL -o "${prefix}/plink2.zip" "${url}"
  rm -f "${prefix}/bin/plink2"
  "${prefix}/bin/python3" -c 'import sys, zipfile; zipfile.ZipFile(sys.argv[1]).extract("plink2", sys.argv[2])' \
    "${prefix}/plink2.zip" "${prefix}/bin"
  rm -f "${prefix}/plink2.zip"
  chmod +x "${prefix}/bin/plink2"
  xattr -d com.apple.quarantine "${prefix}/bin/plink2" 2>/dev/null || true
  "${prefix}/bin/plink2" --version >/dev/null
}

# BCFtools/liftover (Genovese et al. 2024) is a bcftools plugin that no package provides: build bcftools from source
# into the tool environment, with the environment's compiler and its htslib, so `bcftools` there lists the plugin.
# It is required: GRCh37 inputs are lifted with it, and the build-detection sample is built with it.
install_bcftools_liftover() {
  local plugins="${BCFTOOLS_PLUGINS:-}" src
  if [[ -n "${GWAS_BCFTOOLS:-}" && -x "${GWAS_BCFTOOLS}" ]] \
      && BCFTOOLS_PLUGINS="${plugins}" "${GWAS_BCFTOOLS}" plugin -l 2>/dev/null | grep -qx liftover; then
    echo "bcftools +liftover: found (${GWAS_BCFTOOLS})"
    return
  fi
  if [[ -x "${ENV_PREFIX}/bin/bcftools" ]] \
      && BCFTOOLS_PLUGINS="${ENV_PREFIX}/libexec/bcftools" "${ENV_PREFIX}/bin/bcftools" plugin -l 2>/dev/null | grep -qx liftover; then
    echo "bcftools +liftover: found (${ENV_PREFIX}/bin/bcftools)"
  else
    src="$(tools_root)/src/bcftools-${BCFTOOLS_VERSION}"
    echo "bcftools +liftover: building bcftools ${BCFTOOLS_VERSION} with the liftover plugin into ${ENV_PREFIX}"
    mkdir -p "$(dirname "${src}")"
    if [[ ! -d "${src}" ]]; then
      curl -fsSL "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
        | tar -xj -C "$(dirname "${src}")"
    fi
    curl -fsSL -o "${src}/plugins/liftover.c" "${LIFTOVER_URL}"
    # The environment's compiler and its htslib, named explicitly (no activation needed); on macOS the compiler uses
    # the system SDK when the Command Line Tools provide one.
    local sdk=""
    [[ "$(uname -s)" == Darwin ]] && sdk="$(xcrun --show-sdk-path 2>/dev/null || true)"
    if ! (cd "${src}" && export PATH="${ENV_PREFIX}/bin:${PATH}" CC="${ENV_PREFIX}/bin/cc" \
            CPPFLAGS="-I${ENV_PREFIX}/include" LDFLAGS="-L${ENV_PREFIX}/lib -Wl,-rpath,${ENV_PREFIX}/lib" \
            ${sdk:+SDKROOT="${sdk}"} \
          && ./configure --prefix="${ENV_PREFIX}" --with-htslib=system --disable-libgsl --disable-perl-filters >build.log 2>&1 \
          && make -j4 >>build.log 2>&1 && make install >>build.log 2>&1); then
      echo "Error: could not build bcftools with the liftover plugin; see ${src}/build.log. Set GWAS_BCFTOOLS and" \
           "BCFTOOLS_PLUGINS in .env to use an existing build." >&2
      return 1
    fi
  fi
  set_env_var GWAS_BCFTOOLS "${ENV_PREFIX}/bin/bcftools"
  set_env_var BCFTOOLS_PLUGINS "${ENV_PREFIX}/libexec/bcftools"
  "${ENV_PREFIX}/bin/bcftools" plugin -l | grep -qx liftover
  echo "bcftools +liftover: ready (${ENV_PREFIX}/bin/bcftools)"
}

# ---- Container (optional) ---------------------------------------------------------------------------------------
# With --container apptainer|docker (or GWAS_CONTAINER in .env) the pipeline steps run inside the published image,
# which holds the same tools; the host keeps only the environment above for Nextflow, Java and the run wrapper, and
# bcftools is not built on the host. The reference build also runs inside the image.
install_container_image() {
  local root sif cmd
  root="$(tools_root)"
  case "${GWAS_CONTAINER}" in
    apptainer|singularity)
      cmd="$(command -v apptainer || command -v singularity || true)"
      if [[ -z "${cmd}" ]]; then
        echo "Error: GWAS_CONTAINER=${GWAS_CONTAINER} but neither apptainer nor singularity is on PATH (try 'module load apptainer')." >&2
        return 1
      fi
      sif="${root}/gwasprocessing.sif"
      if [[ ! -s "${sif}" || "${FORCE}" -eq 1 ]]; then
        echo "container: pulling ${CONTAINER_REF} into ${sif}"
        mkdir -p "${root}"
        "${cmd}" pull --force "${sif}" "docker://${CONTAINER_REF}"
      else
        echo "container: found ${sif}"
      fi
      set_env_var GWAS_CONTAINER_IMAGE "${sif}" ;;
    docker)
      if ! command -v docker >/dev/null 2>&1; then
        echo "Error: GWAS_CONTAINER=docker but docker is not on PATH." >&2
        return 1
      fi
      echo "container: pulling ${CONTAINER_REF}"
      docker pull "${CONTAINER_REF}"
      set_env_var GWAS_CONTAINER_IMAGE "${CONTAINER_REF}" ;;
    *) echo "Error: GWAS_CONTAINER must be apptainer or docker (got ${GWAS_CONTAINER})." >&2; return 1 ;;
  esac
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
  local run=()
  [[ "${FORCE}" -eq 1 ]] && export GWAS_REFERENCE_FORCE=1
  [[ -n "${GWAS_CONTAINER:-}" ]] && run=(in_container)
  "${run[@]+"${run[@]}"}" bash "${REFDATA}/download_fasta_and_chains.sh"
  "${run[@]+"${run[@]}"}" bash "${REFDATA}/build_grch38_references.sh"
  "${run[@]+"${run[@]}"}" python3 "${REFDATA}/validate_references.py" --reference-root "${GWAS_REFERENCE_ROOT:-${ROOT_DIR}/references}"
}

main() {
  set -euo pipefail
  TOOLS_ONLY=0
  REFERENCE_ONLY=0
  FORCE=0
  R_ONLY=0
  ENV_ONLY=0
  local container=""
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --tools-only) TOOLS_ONLY=1 ;;
      --reference-only) REFERENCE_ONLY=1 ;;
      --force) FORCE=1 ;;
      --r-packages) R_ONLY=1 ;;
      --env-only) ENV_ONLY=1 ;;
      --container) container="${2:-}"; shift ;;
      --container=*) container="${1#*=}" ;;
      -h|--help) sed -n '2,19p' "${BASH_SOURCE[0]}"; exit 0 ;;
      *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
    shift
  done
  ensure_env_file
  if [[ -n "${container}" ]]; then
    case "${container}" in
      apptainer|singularity|docker) set_env_var GWAS_CONTAINER "${container}" >/dev/null ;;
      *) echo "--container must be apptainer or docker" >&2; exit 2 ;;
    esac
  fi
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
    install_environment
    if [[ -n "${GWAS_CONTAINER:-}" ]]; then
      install_container_image
    else
      install_bcftools_liftover
    fi
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
