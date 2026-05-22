#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash spredixcan_v7_twas.sh SUBJECT_NAME [BASE_DIR]
  bash spredixcan_v7_twas.sh -h|--help

Arguments:
  SUBJECT_NAME  Subject/project name, for example: TWB1_LAA
  BASE_DIR      Base directory. Default: /mnt/data/ai_agent/gene_analysis

Example:
  bash spredixcan_v7_twas.sh TWB1_LAA
  bash spredixcan_v7_twas.sh TWB1_LAA /mnt/data/ai_agent/gene_analysis
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 1 ]]; then
  usage
  exit 1
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda command not found in PATH."
  exit 1
fi

SUBJECT_NAME=$1
BASE_DIR=${2:-/mnt/data/ai_agent/gene_analysis}

if [[ ! -d "${BASE_DIR}" ]]; then
  echo "ERROR: base_dir not found: ${BASE_DIR}"
  exit 1
fi

SPREDIXCAN_V7_GWAS_FILE="${BASE_DIR}/S_PrediXcan/GWAS/${SUBJECT_NAME}/chr1tochr22.sumstats"
SPREDIXCAN_V7_TISSUE_FILE="${BASE_DIR}/S_PrediXcan/tissues.txt"
SPREDIXCAN_V7_MODEL_DIR="${BASE_DIR}/S_PrediXcan/model"
SPREDIXCAN_V7_OUT_DIR="${BASE_DIR}/S_PrediXcan/result_v7/${SUBJECT_NAME}"
SPREDIXCAN_V7_LOG_DIR="${SPREDIXCAN_V7_OUT_DIR}/log"
SPREDIXCAN_V7_SCRIPT="${BASE_DIR}/tools/MetaXcan/software/SPrediXcan.py"

mkdir -p "${SPREDIXCAN_V7_OUT_DIR}"
mkdir -p "${SPREDIXCAN_V7_LOG_DIR}"

if [[ ! -f "${SPREDIXCAN_V7_GWAS_FILE}" ]]; then
  echo "ERROR: GWAS file not found: ${SPREDIXCAN_V7_GWAS_FILE}"
  exit 1
fi

if [[ ! -f "${SPREDIXCAN_V7_TISSUE_FILE}" ]]; then
  echo "ERROR: Tissue file not found: ${SPREDIXCAN_V7_TISSUE_FILE}"
  exit 1
fi

if [[ ! -d "${SPREDIXCAN_V7_MODEL_DIR}" ]]; then
  echo "ERROR: Model directory not found: ${SPREDIXCAN_V7_MODEL_DIR}"
  exit 1
fi

if [[ ! -f "${SPREDIXCAN_V7_SCRIPT}" ]]; then
  echo "ERROR: SPrediXcan.py not found: ${SPREDIXCAN_V7_SCRIPT}"
  exit 1
fi

LOG_FILE=${SPREDIXCAN_V7_LOG_DIR}/log_$(date +"%Y%m%d_%H%M%S").txt
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========== S-PrediXcan v7 TWAS Start =========="
echo "Time: $(date)"
echo "Subject name: ${SUBJECT_NAME}"
echo "GWAS file: ${SPREDIXCAN_V7_GWAS_FILE}"
echo "Tissue file: ${SPREDIXCAN_V7_TISSUE_FILE}"
echo "Model dir: ${SPREDIXCAN_V7_MODEL_DIR}"
echo "Output dir: ${SPREDIXCAN_V7_OUT_DIR}"
echo "Log file: ${LOG_FILE}"
echo

TISSUES_RUN=()
TISSUES_MODEL_MISSING=()
TISSUES_COV_MISSING=()

while read -r tissue; do
  [[ -z "${tissue}" ]] && continue

  model_db="${SPREDIXCAN_V7_MODEL_DIR}/gtex_v7_${tissue}_imputed_europeans_tw_0.5_signif.db"
  covariance_file="${SPREDIXCAN_V7_MODEL_DIR}/gtex_v7_${tissue}_imputed_eur_covariances.txt.gz"
  output_file="${SPREDIXCAN_V7_OUT_DIR}/${tissue}_SPrediXcan_v7_en.csv"

  if [[ ! -f "${model_db}" ]]; then
    echo "WARNING: model db not found, skip ${tissue}: ${model_db}"
    TISSUES_MODEL_MISSING+=("${tissue}")
    continue
  fi

  if [[ ! -f "${covariance_file}" ]]; then
    echo "WARNING: covariance file not found, skip ${tissue}: ${covariance_file}"
    TISSUES_COV_MISSING+=("${tissue}")
    continue
  fi

  echo "Running SPrediXcan v7 for ${tissue}"

  conda run -n metaxcan python "${SPREDIXCAN_V7_SCRIPT}" \
    --model_db_path "${model_db}" \
    --covariance "${covariance_file}" \
    --gwas_file "${SPREDIXCAN_V7_GWAS_FILE}" \
    --snp_column SNP \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --zscore_column Z \
    --output_file "${output_file}"

  TISSUES_RUN+=("${tissue}")
done < "${SPREDIXCAN_V7_TISSUE_FILE}"

echo
echo "========== Summary =========="
echo "Tissues completed: ${TISSUES_RUN[*]:-}"
echo "Tissues missing model db: ${TISSUES_MODEL_MISSING[*]:-}"
echo "Tissues missing covariance file: ${TISSUES_COV_MISSING[*]:-}"
echo "Log file: ${LOG_FILE}"
echo "S-PrediXcan v7 TWAS finished."
echo "============================="
