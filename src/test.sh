#!/usr/bin/env bash
# Run the test suite with the tools that src/source.sh installed. Takes about 3 minutes and uses no network.
# Usage: bash src/test.sh [-v] [test_module[.Class[.test]] ...]
# On the container route (GWAS_CONTAINER in .env) the suite runs inside the GWASprocessing image, so it tests the
# image's tools; otherwise it runs with the tool environment. Tests whose tools are missing (R with gwasvcf, for
# example) are skipped and listed.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=src/lib/env.sh
. "${ROOT_DIR}/src/lib/env.sh"
load_settings "${ROOT_DIR}"

if [[ -n "${GWAS_CONTAINER:-}" && -z "${GWAS_IN_CONTAINER:-}" ]]; then
  if [[ -z "${GWAS_CONTAINER_IMAGE:-}" ]]; then
    echo "Error: GWAS_CONTAINER=${GWAS_CONTAINER} but GWAS_CONTAINER_IMAGE is not set; run 'bash src/source.sh --container ${GWAS_CONTAINER}'." >&2
    exit 1
  fi
  echo "Running the test suite inside ${GWAS_CONTAINER_IMAGE}"
  # Inside the image the steps run directly (GWAS_CONTAINER empty), with the image's Java and a writable Nextflow home.
  in_container env -u JAVA_HOME -u JAVA_LD_LIBRARY_PATH GWAS_CONTAINER= NXF_HOME=/tmp/gwasprocessing-nextflow-home \
    bash src/test.sh "$@"
  exit
fi

cd "${ROOT_DIR}/pipeline/tests"
echo "Running the GWASprocessing test suite with $(command -v python3)"
names=()
options=()
for arg in "$@"; do
  if [[ "${arg}" == -* ]]; then options+=("${arg}"); else names+=("${arg}"); fi
done
if [[ ${#names[@]} -gt 0 ]]; then
  python3 -m unittest ${options[@]+"${options[@]}"} "${names[@]}"
else
  python3 -m unittest discover -s . ${options[@]+"${options[@]}"}
fi
