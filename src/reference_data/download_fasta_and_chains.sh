#!/usr/bin/env bash
set -euo pipefail
# Download FASTA and chain files used for VCF conversion and liftover.

# shellcheck source=scripts/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

OUT="${GWAS_REFERENCE_ROOT}/genome"
mkdir -p "${OUT}"

echo "Download FASTA/chain references into: ${OUT}"
echo "Set FASTA_GRCH37_URL, FASTA_GRCH38_URL, HG19_TO_HG38_CHAIN_URL, and HG38_TO_HG19_CHAIN_URL in .env."

for var in FASTA_GRCH37_URL FASTA_GRCH38_URL HG19_TO_HG38_CHAIN_URL HG38_TO_HG19_CHAIN_URL; do
  url="${!var:-}"
  if [[ -n "${url}" ]]; then
    wget -c -P "${OUT}" "${url}"
  fi
done
