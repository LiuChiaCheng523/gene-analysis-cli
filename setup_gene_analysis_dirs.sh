#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  bash setup_gene_analysis_dirs.sh BASE_DIR [PROJECT_NAME ...]
  bash setup_gene_analysis_dirs.sh -h|--help

Description:
  Create the expected folder skeleton for this repository under BASE_DIR.
  If one or more PROJECT_NAME values are provided, project-specific folders
  are also created under COJO, FUSION, and S_PrediXcan.

Examples:
  bash setup_gene_analysis_dirs.sh /mnt/data/ai_agent/gene_analysis
  bash setup_gene_analysis_dirs.sh /mnt/data/ai_agent/gene_analysis TWB1_LAA TWB2_VERTIGO
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

BASE_DIR=$1
shift || true

PROJECT_NAMES=("$@")

mkdir -p "${BASE_DIR}"

COMMON_DIRS=(
  "COJO"
  "FUSION"
  "FUSION/fusion_twas"
  "FUSION/GWAS"
  "FUSION/LDREF"
  "FUSION/LDREF/LDREF"
  "FUSION/result_v7"
  "FUSION/result_v8"
  "FUSION/WEIGHTS"
  "FUSION/WEIGHTS_v7"
  "FUSION/WEIGHTS_v7/GTEx.ALL"
  "overlap"
  "pathway"
  "pathway/GOBP"
  "pathway/KEGG"
  "pathway/log"
  "PLINK"
  "PLINK/imputed"
  "PLINK/imputed/glm_logistic"
  "PLINK/imputed/plink_binary_files"
  "PLINK/original"
  "scripts"
  "S_PrediXcan"
  "S_PrediXcan/GWAS"
  "S_PrediXcan/model"
  "S_PrediXcan/model/GTEx_v8"
  "S_PrediXcan/model/GTEx_v8/elastic_net_models"
  "S_PrediXcan/result_v7"
  "S_PrediXcan/result_v8"
  "tools"
  "tools/beagle"
  "tools/EAS_clean"
  "tools/ensembl"
  "tools/gcta"
  "tools/GENCODE"
  "tools/hg19toHg38"
  "tools/MetaXcan"
  "tools/ncbi_dbsnp"
  "VCF"
  "VCF/imputed"
  "VCF/original"
)

for rel_dir in "${COMMON_DIRS[@]}"; do
  mkdir -p "${BASE_DIR}/${rel_dir}"
done

for project_name in "${PROJECT_NAMES[@]}"; do
  mkdir -p "${BASE_DIR}/COJO/${project_name}"
  mkdir -p "${BASE_DIR}/COJO/${project_name}/result"
  mkdir -p "${BASE_DIR}/COJO/${project_name}/table"
  mkdir -p "${BASE_DIR}/COJO/${project_name}/manhattan"
  mkdir -p "${BASE_DIR}/COJO/${project_name}/log"

  mkdir -p "${BASE_DIR}/FUSION/GWAS/${project_name}"
  mkdir -p "${BASE_DIR}/FUSION/result_v7/${project_name}"
  mkdir -p "${BASE_DIR}/FUSION/result_v8/${project_name}"

  mkdir -p "${BASE_DIR}/S_PrediXcan/GWAS/${project_name}"
  mkdir -p "${BASE_DIR}/S_PrediXcan/result_v7/${project_name}"
  mkdir -p "${BASE_DIR}/S_PrediXcan/result_v8/${project_name}"

  mkdir -p "${BASE_DIR}/overlap/${project_name}"
  mkdir -p "${BASE_DIR}/pathway/GOBP/${project_name}"
  mkdir -p "${BASE_DIR}/pathway/KEGG/${project_name}"
  mkdir -p "${BASE_DIR}/pathway/log/${project_name}"
done

echo "========== Gene Analysis Directory Setup =========="
echo "Base dir: ${BASE_DIR}"
if [[ ${#PROJECT_NAMES[@]} -gt 0 ]]; then
  echo "Project folders created for: ${PROJECT_NAMES[*]}"
else
  echo "Project folders created for: none"
fi
echo
echo "Created common skeleton:"
printf '  - %s\n' "${COMMON_DIRS[@]}"
echo "==================================================="
