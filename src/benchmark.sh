#!/usr/bin/env bash
# Check an install: standardise the benchmark GWAS in manifests/test.tsv and compare the results with the values
# stored in src/benchmark/expected/.
# Usage: bash src/benchmark.sh [--force]
#
# Rows already complete for the current settings are not rerun, so a second run only repeats the checks; --force reruns
# them. The 2 GWAS Catalog rows are compared variant by variant with the stored values (10,008 Locke 2015 BMI and
# 10,016 O'Mara 2018 endometrial cancer variants: GRCh38 position, alleles, ES, SE, LP, AF and AFS). The 2 OpenGWAS rows
# must complete; without OPENGWAS_JWT in .env they are skipped. Prints PASS or FAIL for each check and exits 1 on any FAIL.
set -euo pipefail
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=src/lib/env.sh
. "${ROOT_DIR}/src/lib/env.sh"
load_settings "${ROOT_DIR}"
activate_tools "${ROOT_DIR}"
if [[ -n "${GWAS_VENV:-}" && -x "${GWAS_VENV}/bin/python3" ]]; then
  export PATH="${GWAS_VENV}/bin:${PATH}"
fi

FORCE=()
case "${1:-}" in
  --force) FORCE=(--force) ;;
  "") ;;
  -h|--help) sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
  *) echo "Unknown option: $1 (use --force or nothing)" >&2; exit 2 ;;
esac

WORK_ROOT="${GWAS_WORK_ROOT:-GWAS}"
[[ "${WORK_ROOT}" == /* ]] || WORK_ROOT="${ROOT_DIR}/${WORK_ROOT}"
MANIFEST="${ROOT_DIR}/manifests/test.tsv"
EXPECTED="${ROOT_DIR}/src/benchmark/expected"

run_manifest="${MANIFEST}"
skipped_opengwas=0
if [[ -z "${OPENGWAS_JWT:-}" ]]; then
  run_manifest="$(mktemp "${TMPDIR:-/tmp}/gwasprocessing_benchmark.XXXXXX")"
  trap 'rm -f "${run_manifest}"' EXIT
  grep -v "opengwas.io" "${MANIFEST}" > "${run_manifest}"
  skipped_opengwas=1
  echo "OPENGWAS_JWT is not set in .env: the 2 OpenGWAS rows are skipped (add a token to include them)."
fi

echo "Standardising the benchmark GWAS (the first run downloads about 630 MB and takes about 50 minutes with 3 workers)"
run_status=0
bash "${ROOT_DIR}/src/standardise.sh" --manifest "${run_manifest}" --row-limit 0 --report false ${FORCE[@]+"${FORCE[@]}"} \
  || run_status=$?

results=()
failed=0
check() {  # check NAME EXPECTED_FILE OUTPUT_FOLDER
  local vcf="${WORK_ROOT}/$3/standardised.gwas.vcf.gz" line
  if [[ ! -s "${vcf}" ]]; then
    results+=("FAIL  $1: no output (${vcf})"); failed=1; return
  fi
  if line="$(python3 "${ROOT_DIR}/src/benchmark/check_install.py" --expected "$2" --output "${vcf}" 2>&1)"; then
    results+=("PASS  $1: $(printf '%s\n' "${line}" | grep 'expected variants' | head -1)")
  else
    results+=("FAIL  $1: $(printf '%s\n' "${line}" | grep 'expected variants' | head -1)"); failed=1
    printf '%s\n' "${line}" | sed -n '2,11p' | sed 's/^/      /' >&2
  fi
}
exists() {  # exists NAME OUTPUT_FOLDER
  if [[ -s "${WORK_ROOT}/$2/standardised.gwas.vcf.gz" ]]; then
    results+=("PASS  $1: completed")
  else
    results+=("FAIL  $1: no output"); failed=1
  fi
}

check "Locke 2015 BMI, GWAS Catalog" "${EXPECTED}/locke_2015.tsv" "locke_2015_25673413/standardise/GWAScatalog_GCST002783"
check "O'Mara 2018 endometrial cancer, GWAS Catalog" "${EXPECTED}/omara_2018.tsv" "omara_2018_30093612/standardise/GWAScatalog_GCST006464"
if [[ "${skipped_opengwas}" -eq 0 ]]; then
  exists "Locke 2015 BMI, OpenGWAS" "locke_2015_25673413/standardise/openGWAS_ebi-a-GCST002783"
  exists "O'Mara 2018 endometrial cancer, OpenGWAS" "omara_2018_30093612/standardise/openGWAS_ebi-a-GCST006464"
fi

echo ""
echo "GWASprocessing benchmark"
printf '  %s\n' "${results[@]}"
if [[ "${skipped_opengwas}" -eq 1 ]]; then
  echo "  SKIP  OpenGWAS rows: no OPENGWAS_JWT"
fi
if [[ "${failed}" -ne 0 || "${run_status}" -ne 0 ]]; then
  echo "  Result: FAIL. The step log of each GWAS is in ${WORK_ROOT}/manifest.html."
  exit 1
fi
echo "  Result: PASS. This install reproduces the stored benchmark."
