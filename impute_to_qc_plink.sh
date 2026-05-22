#!/bin/bash
set -euo pipefail

SUBJECT_NAME=${1:?Usage: bash impute_to_qc_plink.sh SUBJECT_NAME}

BASE_DIR=/mnt/data/ai_agent/gene_analysis
VCF_DIR=${BASE_DIR}/VCF/original
BEAGLE_DIR=${BASE_DIR}/tools/beagle
REF_DIR=${BASE_DIR}/tools/EAS_clean
IMPUTED_VCF_DIR=${BASE_DIR}/VCF/imputed

SUBJECT_VCF_DIR=${VCF_DIR}/${SUBJECT_NAME}
SUBJECT_IMPUTED_VCF_DIR=${IMPUTED_VCF_DIR}/${SUBJECT_NAME}

mkdir -p "${SUBJECT_IMPUTED_VCF_DIR}"

for CHR in {1..22}; do
  mkdir -p "${SUBJECT_IMPUTED_VCF_DIR}/chr${CHR}"
done


echo "[1/3] Run Beagle imputation"

cd "${SUBJECT_VCF_DIR}"

for CHR in {1..22}; do
  java -Xmx110g -jar "${BEAGLE_DIR}/beagle.27Feb25.75f.jar" \
    gt="chr${CHR}/chr${CHR}_conformed.vcf.gz" \
    ref="${REF_DIR}/chr${CHR}_EAS_clean.vcf.gz" \
    map="${BEAGLE_DIR}/plink.GRCh37.map/plink.chr${CHR}.GRCh37.map" \
    out="${SUBJECT_IMPUTED_VCF_DIR}/chr${CHR}/chr${CHR}_imputed" \
    nthreads=6 \
    window=30 \
    overlap=2
done


echo "[2/3] Filter imputed VCFs by DR2 >= 0.7"

for CHR in {1..22}; do
  cd "${SUBJECT_IMPUTED_VCF_DIR}/chr${CHR}"

  tabix -p vcf "chr${CHR}_imputed.vcf.gz"

  mkdir -p DR2_0.7

  bcftools view \
    -i 'INFO/DR2>=0.7' \
    -Oz \
    -o "DR2_0.7/chr${CHR}_imputed_DR2_0.7.vcf.gz" \
    "chr${CHR}_imputed.vcf.gz"

  bcftools index -f "DR2_0.7/chr${CHR}_imputed_DR2_0.7.vcf.gz"
done


echo "[3/3] Convert filtered VCFs to PLINK format"

for CHR in {1..22}; do
  cd "${SUBJECT_IMPUTED_VCF_DIR}/chr${CHR}"

  plink2 \
    --vcf "DR2_0.7/chr${CHR}_imputed_DR2_0.7.vcf.gz" \
    --make-bed \
    --out "DR2_0.7/stroke_chr${CHR}_DR2_0.7" \
    --max-alleles 2
done

echo "Imputation to QC PLINK completed successfully: ${SUBJECT_NAME}"
