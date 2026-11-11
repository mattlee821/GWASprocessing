#!/usr/bin/env bash
# Download and build the GRCh38 references: the dbSNP index, 1000 Genomes GRCh38 frequencies and the paired
# build-detection sample.
# Usage: bash src/reference_data/build_grch38_references.sh
#
# Raw downloads (about 34 GB: dbSNP 30 GB, 1000 Genomes 4 GB) go to GWAS_REFERENCE_RAW_ROOT (default
# <reference root>/raw); they are read once, start to finish, so they can live on a slow or external drive. The
# built lookups (about 25 GB) go under GWAS_REFERENCE_ROOT, which should be local disk: they are read at random
# positions on every run. Needs bcftools and plink2 on PATH (PLINK2 overrides), and python3 with numpy and pandas.
# Idempotent: existing downloads and builds are skipped (delete a build's meta.json to rebuild it).
set -euo pipefail
# shellcheck source=src/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

RAW="$(resolve_path "${GWAS_REFERENCE_RAW_ROOT:-${GWAS_REFERENCE_ROOT}/raw}")"
DBSNP_URL="${DBSNP_URL:-https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz}"
ASSEMBLY_REPORT_URL="${ASSEMBLY_REPORT_URL:-https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt}"
HG38_PSAM_URL="${HG38_PSAM_URL:-https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x&dl=1}"
PHASE3_PSAM_URL="${PHASE3_PSAM_URL:-https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1}"
PYTHON="${PYTHON:-python3}"

mkdir -p "${RAW}/dbsnp" "${RAW}/1000genomes_hg38/by_chromosome"
fetch "${ASSEMBLY_REPORT_URL}" "${RAW}/dbsnp/GCF_000001405.40_assembly_report.txt"
fetch "${DBSNP_URL}.md5" "${RAW}/dbsnp/GCF_000001405.40.gz.md5"
fetch "${DBSNP_URL}" "${RAW}/dbsnp/GCF_000001405.40.gz"
expected="$(awk '{print $1}' "${RAW}/dbsnp/GCF_000001405.40.gz.md5")"
observed="$( (md5sum "${RAW}/dbsnp/GCF_000001405.40.gz" 2>/dev/null || md5 -r "${RAW}/dbsnp/GCF_000001405.40.gz") | awk '{print $1}')"
if [[ "${expected}" != "${observed}" ]]; then
  echo "dbSNP checksum mismatch (expected ${expected}, got ${observed}); delete the file and re-run" >&2
  exit 1
fi
fetch "${HG38_PSAM_URL}" "${RAW}/1000genomes_hg38/all_hg38.psam"
fetch "${PHASE3_PSAM_URL}" "${RAW}/1000genomes_hg38/phase3_corrected.psam"
while IFS=$'\t' read -r chrom pgen pvar; do
  [[ "${chrom}" == chrom ]] && continue
  fetch "${pgen}" "${RAW}/1000genomes_hg38/by_chromosome/chr${chrom}_hg38.pgen.zst"
  fetch "${pvar}" "${RAW}/1000genomes_hg38/by_chromosome/chr${chrom}_hg38.pvar.zst"
done < "${REFERENCE_SCRIPT_DIR}/1000genomes_grch38_urls.tsv"

"${PYTHON}" "${REFERENCE_SCRIPT_DIR}/build_dbsnp_index.py" --dbsnp "${RAW}/dbsnp/GCF_000001405.40.gz" \
  --assembly-report "${RAW}/dbsnp/GCF_000001405.40_assembly_report.txt" --out "${GWAS_REFERENCE_ROOT}/dbsnp/GRCh38"
"${PYTHON}" "${REFERENCE_SCRIPT_DIR}/build_1000genomes_grch38.py" --raw "${RAW}/1000genomes_hg38" \
  --out "${GWAS_REFERENCE_ROOT}/1000genomes/GRCh38"
"${PYTHON}" "${REFERENCE_SCRIPT_DIR}/build_build_detection_sample.py" --reference-root "${GWAS_REFERENCE_ROOT}"
