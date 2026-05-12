#!/usr/bin/env bash
set -euo pipefail
# Download dbSNP VCFs used for build detection and rsID lookup.

# shellcheck source=scripts/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

OUT="${GWAS_REFERENCE_ROOT}/dbsnp"
mkdir -p "${OUT}"

echo "Download dbSNP GRCh37/GRCh38 VCFs into: ${OUT}"
echo "Set DBSNP_GRCH37_URL and DBSNP_GRCH38_URL in .env to enable direct downloads."

if [[ -n "${DBSNP_GRCH37_URL:-}" ]]; then
  wget -c -P "${OUT}" "${DBSNP_GRCH37_URL}"
fi
if [[ -n "${DBSNP_GRCH38_URL:-}" ]]; then
  wget -c -P "${OUT}" "${DBSNP_GRCH38_URL}"
fi
