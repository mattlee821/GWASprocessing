#!/usr/bin/env bash
# Download the UCSC liftover chains and the GRCh38 reference FASTA into <reference root>/genome.
# Usage: bash src/reference_data/download_fasta_and_chains.sh
# The FASTA is decompressed to genome/GRCh38.fa and indexed (GRCh38.fa.fai) so standardisation can align
# REF/ALT to the GRCh38 reference. Override the sources with HG19_TO_HG38_CHAIN_URL, HG38_TO_HG19_CHAIN_URL
# and FASTA_GRCH38_URL in server.env. Idempotent: existing files are skipped.
set -euo pipefail
# shellcheck source=src/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

OUT="${GWAS_REFERENCE_ROOT}/genome"
mkdir -p "${OUT}"

HG19_TO_HG38_CHAIN_URL="${HG19_TO_HG38_CHAIN_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz}"
HG38_TO_HG19_CHAIN_URL="${HG38_TO_HG19_CHAIN_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz}"
FASTA_GRCH38_URL="${FASTA_GRCH38_URL:-https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"

fetch "${HG19_TO_HG38_CHAIN_URL}" "${OUT}/hg19ToHg38.over.chain.gz"
fetch "${HG38_TO_HG19_CHAIN_URL}" "${OUT}/hg38ToHg19.over.chain.gz"

if [[ -s "${OUT}/GRCh38.fa" ]]; then
  echo "Skip: ${OUT}/GRCh38.fa already exists"
else
  fetch "${FASTA_GRCH38_URL}" "${OUT}/GRCh38.fa.gz"
  gunzip -c "${OUT}/GRCh38.fa.gz" > "${OUT}/GRCh38.fa.part"
  mv "${OUT}/GRCh38.fa.part" "${OUT}/GRCh38.fa"
  rm -f "${OUT}/GRCh38.fa.gz"
fi
if [[ ! -s "${OUT}/GRCh38.fa.fai" ]]; then
  python3 - "${REPO_ROOT}/pipeline/bin" "${OUT}/GRCh38.fa" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from gwas_alleles import build_fai
print(f"Indexed {build_fai(sys.argv[2])}")
PY
fi
