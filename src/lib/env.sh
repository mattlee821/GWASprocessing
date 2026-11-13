#!/usr/bin/env bash
# Shared helpers for GWASprocessing shell scripts. Source this file; it defines functions only.

# load_env_defaults FILE: export KEY=VALUE lines from FILE for keys that are not already set, so values
# from the environment or the command line always win over file defaults. Comments and blank lines are skipped.
load_env_defaults() {
  local file="$1" line key value
  [[ -f "${file}" ]] || return 0
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ "${line}" =~ ^[[:space:]]*(#|$) ]] && continue
    [[ "${line}" == *=* ]] || continue
    key="${line%%=*}"
    key="${key#export }"
    key="${key//[[:space:]]/}"
    [[ "${key}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] || continue
    value="${line#*=}"
    value="${value%\"}"; value="${value#\"}"; value="${value%\'}"; value="${value#\'}"
    if [[ -z "${!key+x}" ]]; then
      export "${key}=${value}"
    fi
  done < "${file}"
}

# load_settings ROOT: load the settings for a checkout from ROOT/.env (untracked; src/source.sh creates it from the
# tracked .env.example). The template itself is never read: its values are placeholders. Without .env, the scripts' built-in defaults and the
# environment apply, and a note says how to create .env.
load_settings() {
  if [[ -f "$1/.env" ]]; then
    load_env_defaults "$1/.env"
  else
    echo "Note: no .env in $1. Run 'bash src/source.sh --env-only' to create it from .env.example; using the" \
         "environment and built-in defaults." >&2
  fi
  activate_tools "$1"
}

# activate_tools ROOT: put the tool environment that src/source.sh created (GWAS_ENV, relative to ROOT or absolute)
# first on PATH, so python3, nextflow, bgzip, tabix, plink2 and bcftools all come from it, point bcftools at the
# plugins built into it, and point Java at the environment's (NXF_JAVA_HOME, if set, still wins for Nextflow). Inside the container image (GWAS_IN_CONTAINER=1) the image's own tools are used instead.
activate_tools() {
  local env="${GWAS_ENV:-}"
  # Python uses the tool environment's packages only, never ones installed in the user's ~/.local
  export PYTHONNOUSERSITE=1
  [[ -n "${GWAS_IN_CONTAINER:-}" || -z "${env}" ]] && return 0
  [[ "${env}" == /* ]] || env="$1/${env}"
  [[ -d "${env}/bin" ]] || return 0
  case ":${PATH}:" in
    *":${env}/bin:"*) ;;
    *) export PATH="${env}/bin:${PATH}" ;;
  esac
  if [[ -z "${BCFTOOLS_PLUGINS:-}" && -d "${env}/libexec/bcftools" ]]; then
    export BCFTOOLS_PLUGINS="${env}/libexec/bcftools"
  fi
  # The environment's Java (conda openjdk lives in lib/jvm and is normally found through activation), for Nextflow.
  if [[ -x "${env}/lib/jvm/bin/java" && -z "${NXF_JAVA_HOME:-}" ]]; then
    export JAVA_HOME="${env}/lib/jvm" JAVA_LD_LIBRARY_PATH="${env}/lib/jvm/lib/server"
  fi
}

# require_real_path NAME VALUE: stop while a path still holds a /WHERE/... placeholder.
require_real_path() {
  if [[ "$2" == /WHERE/* ]]; then
    echo "$1 is still the placeholder $2. Set it in .env (bash src/source.sh --env-only creates .env)." >&2
    exit 2
  fi
}

# env_value FILE KEY: the value of KEY in FILE (empty when absent).
env_value() {
  local line
  [[ -f "$1" ]] || return 0
  line="$(grep -E "^[[:space:]]*(export[[:space:]]+)?$2=" "$1" | tail -n 1)" || true
  line="${line#*=}"
  line="${line%\"}"; line="${line#\"}"; line="${line%\'}"; line="${line#\'}"
  printf '%s' "${line}"
}

# fetch URL OUT: download to OUT.part and rename on success, so an interrupted download is never
# mistaken for a complete file. Skips the download when OUT already exists and is non-empty.
fetch() {
  local url="$1" out="$2"
  if [[ -s "${out}" ]]; then
    echo "Skip: ${out} already exists"
    return 0
  fi
  if command -v wget >/dev/null 2>&1; then
    wget -c -O "${out}.part" "${url}"
  else
    curl -fL --retry 5 -C - -o "${out}.part" "${url}"
  fi
  mv "${out}.part" "${out}"
}

# in_container CMD...: run a command inside the GWASprocessing image (GWAS_CONTAINER, GWAS_CONTAINER_IMAGE) with the
# repository (ROOT_DIR) and reference folders mounted at the same paths, so .env and the scripts work unchanged.
in_container() {
  local folders=("${ROOT_DIR}") folder args=() resolved
  for folder in "${GWAS_REFERENCE_ROOT:-references}" "${GWAS_REFERENCE_RAW_ROOT:-}"; do
    [[ -n "${folder}" ]] || continue
    resolved="${folder}"
    [[ "${resolved}" == /* ]] || resolved="${ROOT_DIR}/${resolved}"
    mkdir -p "${resolved}"
    folders+=("${resolved}")
  done
  case "${GWAS_CONTAINER}" in
    apptainer|singularity)
      for folder in "${folders[@]}"; do args+=(--bind "${folder}"); done
      "$(command -v apptainer || command -v singularity)" exec "${args[@]}" --pwd "${ROOT_DIR}" "${GWAS_CONTAINER_IMAGE}" "$@" ;;
    docker)
      for folder in "${folders[@]}"; do args+=(-v "${folder}:${folder}"); done
      docker run --rm -u "$(id -u):$(id -g)" "${args[@]}" -w "${ROOT_DIR}" \
        -e GWAS_REFERENCE_ROOT -e GWAS_REFERENCE_RAW_ROOT -e GWAS_REFERENCE_FORCE "${GWAS_CONTAINER_IMAGE}" "$@" ;;
  esac
}
