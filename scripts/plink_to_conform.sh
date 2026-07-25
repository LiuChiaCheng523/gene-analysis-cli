#!/bin/bash
set -euo pipefail

SUBJECT_NAME=${1:?Usage: bash plink_to_conform.sh SUBJECT_NAME PLINK_INPUT_FILE}
PLINK_INPUT_FILE=${2:?Usage: bash plink_to_conform.sh SUBJECT_NAME PLINK_INPUT_FILE}

BASE_DIR=/mnt/data/ai_agent/gene_analysis
PLINK_DIR=${BASE_DIR}/PLINK/original
VCF_DIR=${BASE_DIR}/VCF/original
BEAGLE_DIR=${BASE_DIR}/tools/beagle
REF_DIR=${BASE_DIR}/tools/EAS_clean

SUBJECT_PLINK_DIR=${PLINK_DIR}/${SUBJECT_NAME}
SUBJECT_VCF_DIR=${VCF_DIR}/${SUBJECT_NAME}

mkdir -p "${SUBJECT_VCF_DIR}"

echo "[1/4] Convert PLINK to VCF"
cd "${SUBJECT_PLINK_DIR}"

plink \
  --bfile "${PLINK_INPUT_FILE}" \
  --recode vcf bgz \
  --out "${SUBJECT_VCF_DIR}/final"


echo "[2/4] Normalize and sort VCF"
cd "${SUBJECT_VCF_DIR}"

bcftools norm -m -both -d exact final.vcf.gz \
  | bcftools sort -Oz -o input.vcf.gz

bcftools index -f input.vcf.gz


echo "[3/4] Split VCF by chromosome"
for CHR in {1..22}; do
  mkdir -p "chr${CHR}"

  bcftools view \
    -r "${CHR}" \
    input.vcf.gz \
    -Oz \
    -o "chr${CHR}/chr${CHR}_input.vcf.gz"

  bcftools index -f "chr${CHR}/chr${CHR}_input.vcf.gz"
done


echo "[4/4] Conform genotype to reference"
for CHR in {1..22}; do
  java -jar "${BEAGLE_DIR}/conform-gt.24May16.cee.jar" \
    ref="${REF_DIR}/chr${CHR}_EAS_clean.vcf.gz" \
    gt="chr${CHR}/chr${CHR}_input.vcf.gz" \
    chrom="${CHR}" \
    match=POS \
    out="chr${CHR}/chr${CHR}_conformed"
done

echo "Pipeline completed successfully: ${SUBJECT_NAME}"
