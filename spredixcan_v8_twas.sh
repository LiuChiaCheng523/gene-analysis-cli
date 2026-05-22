#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash spredixcan_v8_twas.sh SUBJECT_NAME [BASE_DIR]
  bash spredixcan_v8_twas.sh -h|--help

Arguments:
  SUBJECT_NAME  Subject/project name, for example: TWB1_LAA
  BASE_DIR      Base directory. Default: /mnt/data/ai_agent/gene_analysis

Example:
  bash spredixcan_v8_twas.sh TWB1_LAA
  bash spredixcan_v8_twas.sh TWB1_LAA /mnt/data/ai_agent/gene_analysis
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

SPREDIXCAN_V8_GWAS_FILE="${BASE_DIR}/S_PrediXcan/GWAS/${SUBJECT_NAME}/chr1tochr22.sumstats"
SPREDIXCAN_V8_TISSUE_FILE="${BASE_DIR}/S_PrediXcan/tissues_v8_elastic_net.txt"
SPREDIXCAN_V8_MODEL_DIR="${BASE_DIR}/S_PrediXcan/model/GTEx_v8/elastic_net_models"
SPREDIXCAN_V8_OUT_DIR="${BASE_DIR}/S_PrediXcan/result_v8/${SUBJECT_NAME}"
SPREDIXCAN_V8_SCRIPT="${BASE_DIR}/tools/MetaXcan/software/SPrediXcan.py"

mkdir -p "${SPREDIXCAN_V8_OUT_DIR}"

if [[ ! -f "${SPREDIXCAN_V8_GWAS_FILE}" ]]; then
  echo "ERROR: GWAS file not found: ${SPREDIXCAN_V8_GWAS_FILE}"
  exit 1
fi

if [[ ! -f "${SPREDIXCAN_V8_TISSUE_FILE}" ]]; then
  echo "ERROR: Tissue file not found: ${SPREDIXCAN_V8_TISSUE_FILE}"
  exit 1
fi

if [[ ! -d "${SPREDIXCAN_V8_MODEL_DIR}" ]]; then
  echo "ERROR: Model directory not found: ${SPREDIXCAN_V8_MODEL_DIR}"
  exit 1
fi

if [[ ! -f "${SPREDIXCAN_V8_SCRIPT}" ]]; then
  echo "ERROR: SPrediXcan.py not found: ${SPREDIXCAN_V8_SCRIPT}"
  exit 1
fi

while read -r tissue; do
  [[ -z "${tissue}" ]] && continue

  model_db="${SPREDIXCAN_V8_MODEL_DIR}/en_${tissue}.db"
  covariance_file="${SPREDIXCAN_V8_MODEL_DIR}/en_${tissue}.txt.gz"
  output_file="${SPREDIXCAN_V8_OUT_DIR}/${tissue}_SPrediXcan_v8_en.csv"

  if [[ ! -f "${model_db}" ]]; then
    echo "WARNING: model db not found, skip ${tissue}: ${model_db}"
    continue
  fi

  if [[ ! -f "${covariance_file}" ]]; then
    echo "WARNING: covariance file not found, skip ${tissue}: ${covariance_file}"
    continue
  fi

  echo "Running SPrediXcan v8 (elastic net) for ${tissue}"

  conda run -n metaxcan python "${SPREDIXCAN_V8_SCRIPT}" \
    --model_db_path "${model_db}" \
    --covariance "${covariance_file}" \
    --gwas_file "${SPREDIXCAN_V8_GWAS_FILE}" \
    --snp_column SNP \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --zscore_column Z \
    --output_file "${output_file}"
done < "${SPREDIXCAN_V8_TISSUE_FILE}"

echo "completed successfully"
