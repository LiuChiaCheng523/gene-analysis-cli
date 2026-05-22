#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash fusion_v7_twas.sh SUBJECT_NAME [BASE_DIR]
  bash fusion_v7_twas.sh -h|--help

Arguments:
  SUBJECT_NAME  Subject/project name, for example: TWB1_LAA
  BASE_DIR      Base directory. Default: /mnt/data/ai_agent/gene_analysis

Example:
  bash fusion_v7_twas.sh TWB1_LAA
  bash fusion_v7_twas.sh TWB1_LAA /mnt/data/ai_agent/gene_analysis
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

SUBJECT_NAME=$1

BASE_DIR=${2:-/mnt/data/ai_agent/gene_analysis}
if [[ ! -d "${BASE_DIR}" ]]; then
  echo "ERROR: base_dir not found: ${BASE_DIR}"
  exit 1
fi

FUSION_V7_GWAS_DIR=${BASE_DIR}/FUSION/GWAS/${SUBJECT_NAME}
FUSION_V7_WEIGHT_DIR=${BASE_DIR}/FUSION/WEIGHTS_v7/GTEx.ALL
FUSION_V7_LDREF=${BASE_DIR}/FUSION/LDREF/LDREF/1000G.EUR.
FUSION_V7_OUT_DIR=${BASE_DIR}/FUSION/result_v7/${SUBJECT_NAME}
FUSION_V7_LOG_DIR=${FUSION_V7_OUT_DIR}/log
FUSION_V7_TISSUE_LIST=${FUSION_V7_WEIGHT_DIR}/tissue_list_v7.txt
FUSION_V7_SCRIPT=${BASE_DIR}/FUSION/fusion_twas/FUSION.assoc_test.R

mkdir -p "${FUSION_V7_OUT_DIR}"
mkdir -p "${FUSION_V7_LOG_DIR}"

if [[ ! -d "${FUSION_V7_GWAS_DIR}" ]]; then
  echo "ERROR: FUSION GWAS directory not found: ${FUSION_V7_GWAS_DIR}"
  exit 1
fi

if [[ ! -d "${FUSION_V7_WEIGHT_DIR}" ]]; then
  echo "ERROR: FUSION weight directory not found: ${FUSION_V7_WEIGHT_DIR}"
  exit 1
fi

if [[ ! -f "${FUSION_V7_TISSUE_LIST}" ]]; then
  echo "ERROR: FUSION tissue list not found: ${FUSION_V7_TISSUE_LIST}"
  exit 1
fi

if [[ ! -f "${FUSION_V7_SCRIPT}" ]]; then
  echo "ERROR: FUSION.assoc_test.R not found: ${FUSION_V7_SCRIPT}"
  exit 1
fi

if [[ ! -d "${BASE_DIR}/FUSION/LDREF/LDREF" ]]; then
  echo "ERROR: FUSION LDREF directory not found: ${BASE_DIR}/FUSION/LDREF/LDREF"
  exit 1
fi

LOG_FILE=${FUSION_V7_LOG_DIR}/log_$(date +"%Y%m%d_%H%M%S").txt
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========== FUSION v7 TWAS Start =========="
echo "Time: $(date)"
echo "Subject name: ${SUBJECT_NAME}"
echo "GWAS dir: ${FUSION_V7_GWAS_DIR}"
echo "Weight dir: ${FUSION_V7_WEIGHT_DIR}"
echo "Output dir: ${FUSION_V7_OUT_DIR}"
echo "Log file: ${LOG_FILE}"
echo

for CHR in {1..22}; do
  mkdir -p "${FUSION_V7_OUT_DIR}/chr${CHR}"
done

CHR_WITH_DATA=()
CHR_WITHOUT_DATA=()
TISSUE_WITHOUT_POS=()

for CHR in {1..22}; do
  SUMSTATS=${FUSION_V7_GWAS_DIR}/chr${CHR}.sumstats
  CHR_OUTDIR=${FUSION_V7_OUT_DIR}/chr${CHR}

  if [[ ! -s "${SUMSTATS}" ]]; then
    echo "WARNING: ${SUMSTATS} does not exist or is empty. Skip chr${CHR}."
    CHR_WITHOUT_DATA+=("${CHR}")
    continue
  fi

  CHR_WITH_DATA+=("${CHR}")

  echo "Running TWAS (GTEx v7) for chr${CHR}"

  while read -r TISSUE; do
    [[ -z "${TISSUE}" ]] && continue

    POS=${FUSION_V7_WEIGHT_DIR}/GTEx.${TISSUE}.pos

    if [[ ! -f "${POS}" ]]; then
      echo "WARNING: POS file not found for tissue ${TISSUE}: ${POS}. Skip."
      TISSUE_WITHOUT_POS+=("${TISSUE}")
      continue
    fi

    echo "  Tissue: ${TISSUE}"

    Rscript "${FUSION_V7_SCRIPT}" \
      --sumstats "${SUMSTATS}" \
      --weights "${POS}" \
      --weights_dir "${FUSION_V7_WEIGHT_DIR}" \
      --ref_ld_chr "${FUSION_V7_LDREF}" \
      --chr "${CHR}" \
      --out "${CHR_OUTDIR}/${TISSUE}_FUSION_v7.dat"

  done < "${FUSION_V7_TISSUE_LIST}"
done

echo
echo "========== Summary =========="
echo "Chromosomes with data: ${CHR_WITH_DATA[*]:-}"
echo "Chromosomes without data: ${CHR_WITHOUT_DATA[*]:-}"
echo "Tissues without POS file: ${TISSUE_WITHOUT_POS[*]:-}"
echo "Log file: ${LOG_FILE}"
echo "FUSION v7 TWAS finished."
echo "============================="
