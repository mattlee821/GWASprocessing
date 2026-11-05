#!/usr/bin/env bash
# Download the 1000 Genomes phase 3 PLINK 2 files and relatedness lists into <reference root>/1000genomes/phase3/raw.
# Usage: bash src/reference_data/download_1000genomes_phase3.sh
# Idempotent: complete files are skipped; interrupted downloads resume from their .part file.
set -euo pipefail
# shellcheck source=src/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

RAW="${GWAS_REFERENCE_ROOT}/1000genomes/phase3/raw"
mkdir -p "${RAW}"

fetch https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst "${RAW}/all_phase3.pgen.zst"
fetch https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst "${RAW}/all_phase3.pvar.zst"
fetch https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam "${RAW}/all_phase3.psam"
fetch "https://www.dropbox.com/s/0omyj2tyu7jmmw9/deg1_phase3.king.cutoff.out.id?dl=1" "${RAW}/deg1_phase3.king.cutoff.out.id"
fetch "https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1" "${RAW}/deg2_phase3.king.cutoff.out.id"
fetch https://genetics.ghpc.au.dk/doug/genetic_map_b37.zip "${RAW}/genetic_map_b37.zip"
