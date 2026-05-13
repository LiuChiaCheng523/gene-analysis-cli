#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash fusion_v8_twas.sh SUBJECT_NAME [BASE_DIR]
  bash fusion_v8_twas.sh -h|--help

Arguments:
  SUBJECT_NAME  Subject/project name, for example: TWB1_LAA
  BASE_DIR      Base directory. Default: /mnt/data/ai_agent/gene_analysis

Example:
  bash fusion_v8_twas.sh TWB1_LAA
  bash fusion_v8_twas.sh TWB1_LAA /mnt/data/ai_agent/gene_analysis
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
FUSION_V8_GWAS_DIR=${BASE_DIR}/FUSION/GWAS/${SUBJECT_NAME}
FUSION_V8_WEIGHT_DIR=${BASE_DIR}/FUSION/WEIGHTS
FUSION_V8_LDREF=${BASE_DIR}/FUSION/LDREF/LDREF/1000G.EUR.
FUSION_V8_OUT_DIR=${BASE_DIR}/FUSION/result_v8/${SUBJECT_NAME}
FUSION_V8_LOG_DIR=${FUSION_V8_OUT_DIR}/log
FUSION_V8_TISSUE_LIST=${FUSION_V8_WEIGHT_DIR}/tissue_list.txt
FUSION_V8_SCRIPT=${BASE_DIR}/FUSION/fusion_twas/FUSION.assoc_test.R

mkdir -p "${FUSION_V8_OUT_DIR}"
mkdir -p "${FUSION_V8_LOG_DIR}"

LOG_FILE=${FUSION_V8_LOG_DIR}/log_$(date +"%Y%m%d_%H%M%S").txt
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "========== FUSION v8 TWAS Start =========="
echo "Time: $(date)"
echo "Subject name: ${SUBJECT_NAME}"
echo "GWAS dir: ${FUSION_V8_GWAS_DIR}"
echo "Output dir: ${FUSION_V8_OUT_DIR}"
echo "Log file: ${LOG_FILE}"
echo

for CHR in {1..22}; do
  mkdir -p "${FUSION_V8_OUT_DIR}/chr${CHR}"
done

CHR_WITH_DATA=()
CHR_WITHOUT_DATA=()

for CHR in {1..22}; do
  SUMSTATS=${FUSION_V8_GWAS_DIR}/chr${CHR}.sumstats
  CHR_OUTDIR=${FUSION_V8_OUT_DIR}/chr${CHR}

  if [[ ! -s "${SUMSTATS}" ]]; then
    echo "WARNING: ${SUMSTATS} does not exist or is empty. Skip chr${CHR}."
    CHR_WITHOUT_DATA+=("${CHR}")
    continue
  fi

  CHR_WITH_DATA+=("${CHR}")

  echo "Running TWAS (GTEx v8) for chr${CHR}"

  while read -r TISSUE; do
    [[ -z "${TISSUE}" ]] && continue

    POS=${FUSION_V8_WEIGHT_DIR}/GTExv8.EUR.${TISSUE}.pos

    echo "  Tissue: ${TISSUE}"

    Rscript "${FUSION_V8_SCRIPT}" \
      --sumstats "${SUMSTATS}" \
      --weights "${POS}" \
      --weights_dir "${FUSION_V8_WEIGHT_DIR}" \
      --ref_ld_chr "${FUSION_V8_LDREF}" \
      --chr "${CHR}" \
      --out "${CHR_OUTDIR}/${TISSUE}_FUSION_v8.dat"

  done < "${FUSION_V8_TISSUE_LIST}"
done

echo
echo "========== Summary =========="
echo "Chromosomes with data: ${CHR_WITH_DATA[*]:-}"
echo "Chromosomes without data: ${CHR_WITHOUT_DATA[*]:-}"
echo "Log file: ${LOG_FILE}"
echo "FUSION v8 TWAS finished."
echo "============================="
