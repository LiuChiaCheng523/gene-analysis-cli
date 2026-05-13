# Gene Analysis Usage Guide

This document summarizes the scripts currently in:

`C:\Users\bruce\codex_code\gene_analysis`

Excluded from this guide:

- `cojo_annotation_manhattan_debug.R`

## Recommended Order

The current workflow is best understood as 4 stages:

1. Prepare phenotype files when needed
2. Convert GWAS results into downstream formats
3. Run association tools
4. Generate summary tables and figures

Recommended order:

1. `make_twb2_pheno_cli.R`
2. External PLINK GWAS step
3. `process_gwas_sumstats_cli.R`
4. `cojo.sh`
5. `cojo_annotation_manhattan_cli.R`
6. `fusion_v8_twas.sh`
7. `fusion_manhattan_heatmap_cli.R`
8. `spredixcan_v8_twas.sh`
9. `spredixcan_manhattan_heatmap_cli.R`

## Script Overview

### 1. `make_twb2_pheno_cli.R`

Purpose:
Create TWB2 phenotype and keep files from survey/lab data for downstream GWAS or phenotype-based analysis.

Usage:

```bash
Rscript make_twb2_pheno_cli.R VERTIGO_SELF categorical 1 0
Rscript make_twb2_pheno_cli.R BMI continuous
```

Inputs:

- `release_list_survey.csv`
- `lab_info.csv`
- `TWB2.hg38.impu.v4.fam`

Outputs:

- `TWB2_<var_name>_pheno.txt`
- `TWB2_<var_name>_keep.txt`

Output folder:

```text
/mnt/SP-siliconpower/TWB20250806download/Imputed.120161.TWB2/imputed_120161
```

Notes:

- For `categorical`, output coding is `1 = control`, `2 = case`
- For `continuous`, the original numeric value is kept

### 2. External GWAS step

Purpose:
Run your GWAS outside this folder, typically with PLINK logistic/linear regression.

Expected output for the scripts below:

```text
/mnt/data/ai_agent/gene_analysis/PLINK/imputed/glm_logistic/<project_name>/chrN.PHENO1.glm.logistic.hybrid
```

This step is not implemented as a script in the current folder.

### 3. `process_gwas_sumstats_cli.R`

Purpose:
Convert GWAS chromosome-level results into:

- COJO input files
- FUSION input files
- one combined S-PrediXcan GWAS file

Usage:

```bash
Rscript process_gwas_sumstats_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis
```

Outputs:

- `COJO/<project_name>/chrN.cojo.ma`
- `FUSION/GWAS/<project_name>/chrN.sumstats`
- `S_PrediXcan/GWAS/<project_name>/chr1tochr22.sumstats`

Notes:

- Missing chromosome files are skipped with warnings
- A process log is written under the FUSION GWAS log folder

### 4. `cojo.sh`

Purpose:
Run GCTA-COJO chromosome by chromosome using PLINK binary files and `.cojo.ma` files.

Usage:

```bash
bash cojo.sh TWB1_LAA_test1 _DR2_0.7_QCFiltered 1e-6 10000
```

Arguments:

- `PROJECT_NAME`
- `PLINK_FILE_PATTERN`
- `P_CUTOFF`
- `WINDOW`
- optional `BASE_DIR`

Expected PLINK input:

```text
/mnt/data/ai_agent/gene_analysis/PLINK/imputed/plink_binary_files/<project_name>/chrN<PLINK_FILE_PATTERN>.bed/.bim/.fam
```

Expected COJO input:

```text
/mnt/data/ai_agent/gene_analysis/COJO/<project_name>/chrN.cojo.ma
```

Outputs:

```text
/mnt/data/ai_agent/gene_analysis/COJO/<project_name>/result/chrN_p<P_CUTOFF_TAG>.*
```

Logs:

```text
/mnt/data/ai_agent/gene_analysis/COJO/<project_name>/result/log/
```

Notes:

- Missing PLINK or COJO input files are skipped
- Terminal summary is shown at the end

### 5. `cojo_annotation_manhattan_cli.R`

Purpose:
Annotate COJO hits using local GTF and generate:

- annotation tables
- nearest protein-coding gene table
- GWAS Manhattan plot

Usage:

```bash
Rscript cojo_annotation_manhattan_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis
```

Inputs:

- `COJO/<project_name>/result/chrN_p1e5.jma.cojo`
- `PLINK/imputed/glm_logistic/<project_name>/chrN.PHENO1.glm.logistic.hybrid`
- `tools/ensembl/Homo_sapiens.GRCh37.87.gtf`

Outputs:

```text
/mnt/data/ai_agent/gene_analysis/COJO/<project_name>/table
/mnt/data/ai_agent/gene_analysis/COJO/<project_name>/manhattan
/mnt/data/ai_agent/gene_analysis/COJO/<project_name>/log
```

Notes:

- This script uses local GTF, not live Ensembl queries
- Missing chromosome COJO files are skipped with warnings

### 6. `fusion_v8_twas.sh`

Purpose:
Run FUSION TWAS using GTEx v8 models.

Usage:

```bash
bash fusion_v8_twas.sh TWB1_LAA_test
```

Optional:

```bash
bash fusion_v8_twas.sh TWB1_LAA_test /mnt/data/ai_agent/gene_analysis
```

Expected input:

```text
/mnt/data/ai_agent/gene_analysis/FUSION/GWAS/<project_name>/chrN.sumstats
```

Outputs:

```text
/mnt/data/ai_agent/gene_analysis/FUSION/result_v8/<project_name>/chrN/
```

Logs:

```text
/mnt/data/ai_agent/gene_analysis/FUSION/result_v8/<project_name>/log/
```

Notes:

- Missing chromosome sumstats are skipped
- Log file is written and also streamed to terminal

### 7. `fusion_manhattan_heatmap_cli.R`

Purpose:
Summarize FUSION results into:

- merged TWAS tables
- FDR-filtered tables
- Manhattan plot
- heatmap

Usage:

```bash
Rscript fusion_manhattan_heatmap_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis v8 0.15
```

Also supports:

```bash
Rscript fusion_manhattan_heatmap_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis v7 0.15
```

Outputs:

```text
/mnt/data/ai_agent/gene_analysis/FUSION/result_v8/<project_name>/table
/mnt/data/ai_agent/gene_analysis/FUSION/result_v8/<project_name>/manhanttan
/mnt/data/ai_agent/gene_analysis/FUSION/result_v8/<project_name>/heatmap
```

or for v7:

```text
/mnt/data/ai_agent/gene_analysis/FUSION/result_v7/<project_name>/
```

Notes:

- Manhattan y-axis is dynamic, not fixed
- Extremely small P values are handled safely

### 8. `spredixcan_v8_twas.sh`

Purpose:
Run S-PrediXcan GTEx v8 elastic net models.

Usage:

```bash
bash spredixcan_v8_twas.sh TWB1_LAA_hg38
```

Optional:

```bash
bash spredixcan_v8_twas.sh TWB1_LAA_hg38 /mnt/data/ai_agent/gene_analysis
```

Expected input:

```text
/mnt/data/ai_agent/gene_analysis/S_PrediXcan/GWAS/<project_name>/chr1tochr22.sumstats
```

Outputs:

```text
/mnt/data/ai_agent/gene_analysis/S_PrediXcan/result_v8/<project_name>/
```

Notes:

- Uses `conda run -n metaxcan`
- Missing tissue model files are skipped with warnings

### 9. `spredixcan_manhattan_heatmap_cli.R`

Purpose:
Summarize S-PrediXcan results into:

- merged TWAS tables
- FDR-filtered tables
- Manhattan plot
- heatmap

Usage:

```bash
Rscript spredixcan_manhattan_heatmap_cli.R TWB2_VERTIGO /mnt/data/ai_agent/gene_analysis v8 0.15
```

Also supports:

```bash
Rscript spredixcan_manhattan_heatmap_cli.R TWB2_VERTIGO /mnt/data/ai_agent/gene_analysis v7 0.15
```

Outputs:

```text
/mnt/data/ai_agent/gene_analysis/S_PrediXcan/result_v8/<project_name>/table
/mnt/data/ai_agent/gene_analysis/S_PrediXcan/result_v8/<project_name>/manhanttan
/mnt/data/ai_agent/gene_analysis/S_PrediXcan/result_v8/<project_name>/heatmap
```

or for v7:

```text
/mnt/data/ai_agent/gene_analysis/S_PrediXcan/result_v7/<project_name>/
```

Notes:

- Manhattan y-axis is dynamic, not fixed
- Extremely small P values are handled safely

## Minimal Practical Pipelines

### Pipeline A: COJO only

1. Run external GWAS
2. Run `process_gwas_sumstats_cli.R`
3. Run `cojo.sh`
4. Run `cojo_annotation_manhattan_cli.R`

### Pipeline B: FUSION TWAS

1. Run external GWAS
2. Run `process_gwas_sumstats_cli.R`
3. Run `fusion_v8_twas.sh`
4. Run `fusion_manhattan_heatmap_cli.R`

### Pipeline C: S-PrediXcan TWAS

1. Run external GWAS
2. Run `process_gwas_sumstats_cli.R`
3. Run `spredixcan_v8_twas.sh`
4. Run `spredixcan_manhattan_heatmap_cli.R`

## Current Script List

- `cojo.sh`
- `cojo_annotation_manhattan_cli.R`
- `fusion_manhattan_heatmap_cli.R`
- `fusion_v8_twas.sh`
- `make_twb2_pheno_cli.R`
- `process_gwas_sumstats_cli.R`
- `spredixcan_manhattan_heatmap_cli.R`
- `spredixcan_v8_twas.sh`
