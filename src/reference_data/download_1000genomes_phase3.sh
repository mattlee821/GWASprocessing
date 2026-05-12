#!/usr/bin/env bash
set -euo pipefail
# Download raw 1000 Genomes phase 3 files and genetic maps.

# shellcheck source=scripts/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

RAW="${GWAS_REFERENCE_ROOT}/1000genomes/phase3/raw"
mkdir -p "${RAW}"

wget -c -P "${RAW}" https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst
wget -c -P "${RAW}" https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst
wget -c -O "${RAW}/all_phase3.psam" https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam
wget -c -O "${RAW}/deg1_phase3.king.cutoff.out.id" "https://www.dropbox.com/s/0omyj2tyu7jmmw9/deg1_phase3.king.cutoff.out.id?dl=1"
wget -c -O "${RAW}/deg2_phase3.king.cutoff.out.id" "https://www.dropbox.com/s/zj8d14vv9mp6x3c/deg2_phase3.king.cutoff.out.id?dl=1"
wget -c -P "${RAW}" https://genetics.ghpc.au.dk/doug/genetic_map_b37.zip
