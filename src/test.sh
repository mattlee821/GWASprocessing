#!/usr/bin/env bash
# Run the test suite with the tools that src/source.sh installed. Takes about 3 minutes and uses no network.
# Usage: bash src/test.sh [-v]
# Tests whose tools are missing (bcftools with the liftover plugin, Nextflow, R with gwasvcf) are skipped and listed.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=src/lib/env.sh
. "${ROOT_DIR}/src/lib/env.sh"
load_settings "${ROOT_DIR}"
activate_tools "${ROOT_DIR}"
if [[ -n "${GWAS_VENV:-}" && -x "${GWAS_VENV}/bin/python3" ]]; then
  export PATH="${GWAS_VENV}/bin:${PATH}"
fi
cd "${ROOT_DIR}"
echo "Running the GWASprocessing test suite with $(command -v python3)"
python3 -m unittest discover -s pipeline/tests "$@"
