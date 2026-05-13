#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash cojo.sh PROJECT_NAME PLINK_FILE_PATTERN P_CUTOFF WINDOW [BASE_DIR]
  bash cojo.sh -h|--help

Arguments:
  PROJECT_NAME         Project name, for example: TWB1_LAA_test1
  PLINK_FILE_PATTERN   Suffix after chrN, for example: _DR2_0.7_QCFiltered
  P_CUTOFF             COJO p-value cutoff, for example: 1e-6
  WINDOW               COJO window size, for example: 10000
  BASE_DIR             Base directory. Default: /mnt/data/ai_agent/gene_analysis

Example:
  bash cojo.sh TWB1_LAA_test1 _DR2_0.7_QCFiltered 1e-6 10000
  bash cojo.sh TWB1_LAA_test1 _DR2_0.7_QCFiltered 1e-6 10000 /mnt/data/ai_agent/gene_analysis
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 4 ]]; then
  usage
  exit 1
fi

if ! command -v gcta64 >/dev/null 2>&1; then
  echo "ERROR: gcta64 command not found in PATH."
  exit 1
fi

PROJECT_NAME=$1
PLINK_FILE_PATTERN=$2
P_CUTOFF=$3
WINDOW=$4
BASE_DIR=${5:-/mnt/data/ai_agent/gene_analysis}

PLINK_INPUT_DIR="${BASE_DIR}/PLINK/imputed/plink_binary_files/${PROJECT_NAME}"
COJO_INPUT_DIR="${BASE_DIR}/COJO/${PROJECT_NAME}"
COJO_RESULT_DIR="${BASE_DIR}/COJO/${PROJECT_NAME}/result"
COJO_LOG_DIR="${COJO_RESULT_DIR}/log"

mkdir -p "${COJO_RESULT_DIR}"
mkdir -p "${COJO_LOG_DIR}"

LOG_FILE="${COJO_LOG_DIR}/log_$(date +"%Y%m%d_%H%M%S").txt"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========== COJO Start =========="
echo "Time: $(date)"
echo "Project name: ${PROJECT_NAME}"
echo "Base dir: ${BASE_DIR}"
echo "PLINK pattern: ${PLINK_FILE_PATTERN}"
echo "P cutoff: ${P_CUTOFF}"
echo "Window: ${WINDOW}"
echo "Log file: ${LOG_FILE}"
echo

# Convert scientific notation such as 1e-6 into a filename-safe token like 1e6.
P_CUTOFF_TAG=$(echo "${P_CUTOFF}" | tr -d '+-' )

CHR_COMPLETED=()
CHR_MISSING_PLINK=()
CHR_MISSING_COJO=()

for CHR in {1..22}; do
  echo "Running COJO for chr${CHR}"

  BFILE_PREFIX="${PLINK_INPUT_DIR}/chr${CHR}${PLINK_FILE_PATTERN}"
  COJO_FILE="${COJO_INPUT_DIR}/chr${CHR}.cojo.ma"
  OUT_PREFIX="${COJO_RESULT_DIR}/chr${CHR}_p${P_CUTOFF_TAG}"

  if [[ ! -f "${BFILE_PREFIX}.bed" || ! -f "${BFILE_PREFIX}.bim" || ! -f "${BFILE_PREFIX}.fam" ]]; then
    echo "WARNING: PLINK bfile missing for chr${CHR}. Skip."
    CHR_MISSING_PLINK+=("${CHR}")
    continue
  fi

  if [[ ! -f "${COJO_FILE}" ]]; then
    echo "WARNING: COJO input file missing for chr${CHR}. Skip."
    CHR_MISSING_COJO+=("${CHR}")
    continue
  fi

  gcta64 \
    --bfile "${BFILE_PREFIX}" \
    --maf 0.01 \
    --cojo-file "${COJO_FILE}" \
    --cojo-slct \
    --cojo-p "${P_CUTOFF}" \
    --cojo-wind "${WINDOW}" \
    --diff-freq 0.2 \
    --out "${OUT_PREFIX}"

  CHR_COMPLETED+=("${CHR}")
done

echo
echo "========== Summary =========="
echo "Chromosomes completed: ${CHR_COMPLETED[*]:-}"
echo "Chromosomes missing PLINK bfile: ${CHR_MISSING_PLINK[*]:-}"
echo "Chromosomes missing COJO input: ${CHR_MISSING_COJO[*]:-}"
echo "Result dir: ${COJO_RESULT_DIR}"
echo "Log file: ${LOG_FILE}"
echo "COJO finished at: $(date)"
echo "============================="
