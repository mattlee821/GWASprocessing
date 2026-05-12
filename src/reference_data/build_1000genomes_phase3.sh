#!/usr/bin/env bash
set -euo pipefail
# Build population-specific 1000 Genomes phase 3 PLINK references.

# shellcheck source=scripts/reference_data/common.sh
. "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

RAW="${GWAS_REFERENCE_ROOT}/1000genomes/phase3/raw"
PROCESSED="${GWAS_REFERENCE_ROOT}/1000genomes/phase3/processed"
mkdir -p "${PROCESSED}"

PLINK2="${PLINK2:-plink2}"
PLINK19="${PLINK19:-plink1.9}"
LDAK="${LDAK:-ldak6.1.linux}"

if ! command -v "${PLINK2}" >/dev/null 2>&1; then
  echo "plink2 not found. Set PLINK2=/path/to/plink2 in .env" >&2
  exit 1
fi

"${PLINK2}" --zst-decompress "${RAW}/all_phase3.pgen.zst" > "${PROCESSED}/all_phase3.pgen"
"${PLINK2}" --zst-decompress "${RAW}/all_phase3.pvar.zst" > "${PROCESSED}/all_phase3.pvar"
cp "${RAW}/all_phase3.psam" "${PROCESSED}/all_phase3.psam"

cd "${PROCESSED}"
"${PLINK2}" --pfile all_phase3 --remove "${RAW}/deg2_phase3.king.cutoff.out.id" --make-pgen --out all_phase3_unrelated
mv all_phase3_unrelated.pgen all_phase3.pgen
mv all_phase3_unrelated.pvar all_phase3.pvar
mv all_phase3_unrelated.psam all_phase3.psam

awk '($5=="AFR"){print 0, $1}' all_phase3.psam > AFR.keep
awk '($5=="AMR"){print 0, $1}' all_phase3.psam > AMR.keep
awk '($5=="EAS"){print 0, $1}' all_phase3.psam > EAS.keep
awk '($5=="EUR"){print 0, $1}' all_phase3.psam > EUR.keep
awk '($5=="SAS"){print 0, $1}' all_phase3.psam > SAS.keep

for POP in AFR AMR EAS EUR SAS; do
  mkdir -p "${POP}"
  "${PLINK2}" --make-bed --out "${POP}/${POP}" --pgen all_phase3.pgen --pvar all_phase3.pvar --psam all_phase3.psam \
    --maf 0.01 --autosome --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all --keep "${POP}.keep"
done

mkdir -p ALL
"${PLINK2}" --make-bed --out ALL/ALL --pgen all_phase3.pgen --pvar all_phase3.pvar --psam all_phase3.psam \
  --maf 0.01 --autosome --snps-only just-acgt --max-alleles 2 --rm-dup exclude-all

echo "Population PLINK references written under ${PROCESSED}"
echo "Run scripts/reference_data/make_1000genomes_lookup_tables.py to create rsID/SNPID/EAF lookup tables."
