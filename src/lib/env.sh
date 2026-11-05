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

# load_settings ROOT: load the settings for a checkout. The untracked server.env holds the real values and is
# read first, so it wins over the placeholders in the tracked .env.
load_settings() {
  load_env_defaults "$1/server.env"
  load_env_defaults "$1/.env"
}

# require_real_path NAME VALUE: stop while a path still holds a /WHERE/... placeholder from .env.
require_real_path() {
  if [[ "$2" == /WHERE/* ]]; then
    echo "$1 is still the placeholder $2. Set it in server.env." >&2
    exit 2
  fi
}

# fetch URL OUT: download to OUT.part and rename on success, so an interrupted download is never
# mistaken for a complete file. Skips the download when OUT already exists and is non-empty.
fetch() {
  local url="$1" out="$2"
  if [[ -s "${out}" ]]; then
    echo "Skip: ${out} already exists"
    return 0
  fi
  wget -c -O "${out}.part" "${url}"
  mv "${out}.part" "${out}"
}
