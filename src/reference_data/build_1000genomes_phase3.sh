#!/usr/bin/env bash
# Build per-population 1000 Genomes phase 3 PLINK references (unrelated samples; common biallelic SNPs on
# chromosomes 1-22 and X) and their allele frequencies, under <reference root>/1000genomes/phase3/processed.
# Usage: bash src/reference_data/build_1000genomes_phase3.sh   (needs plink2; see src/source.sh)
# Idempotent: each output is skipped if it already exists.
set -euo pipefail
# shellcheck source=src/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

RAW="${GWAS_REFERENCE_ROOT}/1000genomes/phase3/raw"
PROCESSED="${GWAS_REFERENCE_ROOT}/1000genomes/phase3/processed"
mkdir -p "${PROCESSED}"
PLINK2="${PLINK2:-plink2}"
if ! command -v "${PLINK2}" >/dev/null 2>&1; then
  echo "plink2 not found. Run 'bash src/source.sh --tools-only' first, or set PLINK2=/path/to/plink2 in server.env" >&2
  exit 1
fi

decompress() {
  local src="$1" dest="$2"
  if [[ -s "${dest}" ]]; then
    echo "Skip: ${dest} already exists"
    return
  fi
  "${PLINK2}" --zst-decompress "${src}" > "${dest}.part"
  mv "${dest}.part" "${dest}"
}
decompress "${RAW}/all_phase3.pgen.zst" "${PROCESSED}/all_phase3.pgen"
decompress "${RAW}/all_phase3.pvar.zst" "${PROCESSED}/all_phase3.pvar"
[[ -s "${PROCESSED}/all_phase3.psam" ]] || cp "${RAW}/all_phase3.psam" "${PROCESSED}/all_phase3.psam"

cd "${PROCESSED}"
if [[ ! -s "all_phase3_unrelated.pgen" ]]; then
  "${PLINK2}" --pfile all_phase3 --remove "${RAW}/deg2_phase3.king.cutoff.out.id" --make-pgen --out all_phase3_unrelated
else
  echo "Skip: all_phase3_unrelated.pgen already exists"
fi

FILTERS=(--maf 0.01 --chr 1-22,X --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all)
for POP in AFR AMR EAS EUR SAS ALL; do
  mkdir -p "${POP}"
  if [[ -s "${POP}/${POP}.bim" && -s "${POP}/stats.afreq" ]]; then
    echo "Skip: ${POP} already built"
    continue
  fi
  keep=()
  if [[ "${POP}" != "ALL" ]]; then
    awk -v pop="${POP}" '($5 == pop) {print 0, $1}' all_phase3_unrelated.psam > "${POP}.keep"
    keep=(--keep "${POP}.keep")
  fi
  "${PLINK2}" --make-bed --out "${POP}/${POP}" --pfile all_phase3_unrelated "${FILTERS[@]}" ${keep[@]+"${keep[@]}"}
  # stats.afreq (#CHROM ID REF ALT ALT_FREQS OBS_CT): ALT_FREQS is the ALT-allele frequency used as AF1.
  "${PLINK2}" --bfile "${POP}/${POP}" --freq --out "${POP}/stats"
done

echo "Population PLINK references written under ${PROCESSED}"
echo "Next: python3 src/reference_data/make_1000genomes_lookup_tables.py --reference-root ${GWAS_REFERENCE_ROOT}"
