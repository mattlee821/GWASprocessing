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
# tracked .env.example). The template itself is never read: its values are placeholders. A checkout that still has
# only the old server.env is read from that, with a warning. Without either, the scripts' built-in defaults and the
# environment apply, and a note says how to create .env.
load_settings() {
  if [[ -f "$1/.env" ]]; then
    load_env_defaults "$1/.env"
  elif [[ -f "$1/server.env" ]]; then
    echo "Warning: reading settings from server.env, the old name. Run 'bash src/source.sh --env-only' to move them" \
         "into .env." >&2
    load_env_defaults "$1/server.env"
  else
    echo "Note: no .env in $1. Run 'bash src/source.sh --env-only' to create it from .env.example; using the" \
         "environment and built-in defaults." >&2
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
