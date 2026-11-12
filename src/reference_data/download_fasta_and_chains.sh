#!/usr/bin/env bash
# Download the UCSC liftover chains and the GRCh38 and GRCh37 reference FASTAs into <reference root>/genome.
# GRCh37 is used only for BCFtools/liftover: as the source genome when GRCh37 inputs are lifted (hg19ToHg38), and as the
# target when the build-detection sample is placed in GRCh37 (hg38ToHg19); outputs are GRCh38.
# Usage: bash src/reference_data/download_fasta_and_chains.sh
# The FASTA is decompressed to genome/GRCh38.fa and indexed (GRCh38.fa.fai) so standardisation can align
# REF/ALT to the GRCh38 reference. Override the sources with HG19_TO_HG38_CHAIN_URL, HG38_TO_HG19_CHAIN_URL
# and FASTA_GRCH38_URL in .env. Idempotent: existing files are skipped.
set -euo pipefail
# shellcheck source=src/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

OUT="${GWAS_REFERENCE_ROOT}/genome"
mkdir -p "${OUT}"

HG19_TO_HG38_CHAIN_URL="${HG19_TO_HG38_CHAIN_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz}"
HG38_TO_HG19_CHAIN_URL="${HG38_TO_HG19_CHAIN_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz}"
FASTA_GRCH38_URL="${FASTA_GRCH38_URL:-https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz}"
FASTA_GRCH37_URL="${FASTA_GRCH37_URL:-https://ftp.ensembl.org/pub/grch37/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz}"

fetch "${HG19_TO_HG38_CHAIN_URL}" "${OUT}/hg19ToHg38.over.chain.gz"
fetch "${HG38_TO_HG19_CHAIN_URL}" "${OUT}/hg38ToHg19.over.chain.gz"

for build in GRCh38 GRCh37; do
  url_var="FASTA_GRCH38_URL"; [[ "${build}" == GRCh37 ]] && url_var="FASTA_GRCH37_URL"
  if [[ -s "${OUT}/${build}.fa" ]]; then
    echo "Skip: ${OUT}/${build}.fa already exists"
  else
    fetch "${!url_var}" "${OUT}/${build}.fa.gz"
    gunzip -c "${OUT}/${build}.fa.gz" > "${OUT}/${build}.fa.part"
    mv "${OUT}/${build}.fa.part" "${OUT}/${build}.fa"
    rm -f "${OUT}/${build}.fa.gz"
  fi
  if [[ ! -s "${OUT}/${build}.fa.fai" ]]; then
    python3 - "${REPO_ROOT}/pipeline/bin" "${OUT}/${build}.fa" <<'PY'
import sys
sys.path.insert(0, sys.argv[1])
from gwas_alleles import build_fai
print(f"Indexed {build_fai(sys.argv[2])}")
PY
  fi
done
