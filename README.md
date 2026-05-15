# Gene Analysis README

This repository contains the operational scripts currently used in:

`/mnt/data/ai_agent/gene_analysis/scripts`

This README focuses on:

- what each script does
- the recommended execution order
- all CLI / shell parameters
- the overall folder structure expected by the scripts

## Overall Workflow

There are 3 major downstream analysis branches after GWAS:

1. `COJO`
2. `FUSION`
3. `S-PrediXcan`

Recommended order:

1. Prepare phenotype files if needed
2. Run GWAS outside this repository
3. Convert GWAS output into downstream formats
4. Run COJO / FUSION / S-PrediXcan
5. Generate tables and plots

## Recommended Execution Order

### A. Phenotype preparation

1. `make_twb2_pheno_cli.R`

### B. GWAS preprocessing

2. External PLINK GWAS step
3. `process_gwas_sumstats_cli.R`

### C. COJO branch

4. `cojo.sh`
5. `cojo_annotation_manhattan_cli.R`

### D. FUSION branch

6. `fusion_v8_twas.sh`
7. `fusion_manhattan_heatmap_cli.R`

### E. S-PrediXcan branch

8. `spredixcan_v8_twas.sh`
9. `spredixcan_manhattan_heatmap_cli.R`

## Folder Structure

Below is the main directory structure assumed by the current code.

```text
/mnt/data/ai_agent/gene_analysis
├─ PLINK/
│  └─ imputed/
│     ├─ glm_logistic/
│     │  └─ <project_name>/
│     │     ├─ chr1.PHENO1.glm.logistic.hybrid
│     │     ├─ chr2.PHENO1.glm.logistic.hybrid
│     │     └─ ...
│     └─ plink_binary_files/
│        └─ <project_name>/
│           ├─ chr1<PLINK_FILE_PATTERN>.bed
│           ├─ chr1<PLINK_FILE_PATTERN>.bim
│           ├─ chr1<PLINK_FILE_PATTERN>.fam
│           └─ ...
├─ COJO/
│  └─ <project_name>/
│     ├─ chr1.cojo.ma
│     ├─ chr2.cojo.ma
│     ├─ ...
│     ├─ result/
│     │  ├─ chr1_p1e6.*
│     │  ├─ chr2_p1e6.*
│     │  └─ log/
│     ├─ table/
│     ├─ manhattan/
│     └─ log/
├─ FUSION/
│  ├─ GWAS/
│  │  └─ <project_name>/
│  │     ├─ chr1.sumstats
│  │     ├─ chr2.sumstats
│  │     ├─ ...
│  │     └─ log/
│  ├─ LDREF/
│  ├─ WEIGHTS/
│  ├─ WEIGHTS_v7/
│  ├─ fusion_twas/
│  ├─ result_v8/
│  │  └─ <project_name>/
│  │     ├─ chr1/
│  │     ├─ chr2/
│  │     ├─ ...
│  │     ├─ table/
│  │     ├─ manhanttan/
│  │     ├─ heatmap/
│  │     └─ log/
│  └─ result_v7/
│     └─ <project_name>/
│        ├─ chr1/
│        ├─ ...
│        ├─ table/
│        ├─ manhanttan/
│        └─ heatmap/
├─ S_PrediXcan/
│  ├─ GWAS/
│  │  └─ <project_name>/
│  │     └─ chr1tochr22.sumstats
│  ├─ model/
│  │  └─ GTEx_v8/
│  │     └─ elastic_net_models/
│  ├─ tissues.txt
│  ├─ tissues_v8_elastic_net.txt
│  ├─ result_v8/
│  │  └─ <project_name>/
│  │     ├─ <tissue>_SPrediXcan_v8_en.csv
│  │     ├─ table/
│  │     ├─ manhanttan/
│  │     └─ heatmap/
│  └─ result_v7/
│     └─ <project_name>/
│        ├─ <tissue>_SPrediXcan_v7_en.csv
│        ├─ table/
│        ├─ manhanttan/
│        └─ heatmap/
├─ scripts/
|  ├─ cojo.sh
|  ├─ cojo_annotation_manhattan_cli.R
|  ├─ fusion_manhattan_heatmap_cli.R
|  ├─ fusion_v8_twas.sh
|  ├─ make_twb2_pheno_cli.R
|  ├─ process_gwas_sumstats_cli.R
|  ├─ spredixcan_manhattan_heatmap_cli.R
|  ├─spredixcan_v8_twas.sh
└─ tools/
   └─ ensembl/
      ├─ Homo_sapiens.GRCh37.87.gtf
      └─ Homo_sapiens.GRCh38.115.gtf
```

## Environment and Required Software

The current scripts assume a Linux environment, especially for:

- file paths under `/mnt/...`
- `bash` scripts
- `wget`
- `conda run`
- `gcta64`

Recommended environment:

- Ubuntu or compatible Linux distribution
- `bash`
- `Rscript`
- `python`
- `conda`

External command-line tools required:

- `PLINK` or `PLINK2`
  Used outside this repository for GWAS
- `gcta64`
  Required by `cojo.sh`
- `python` inside the `metaxcan` conda environment
  Required by `spredixcan_v8_twas.sh`
- `wget`
  Useful for downloading GTF and reference files

Required local reference/data resources:

- `tools/ensembl/Homo_sapiens.GRCh37.87.gtf`
- `tools/ensembl/Homo_sapiens.GRCh38.115.gtf`
- `FUSION/WEIGHTS/`
- `FUSION/WEIGHTS_v7/`
- `FUSION/LDREF/`
- `S_PrediXcan/model/GTEx_v8/elastic_net_models/`

Recommended R packages used across the repository:

- `data.table`
- `dplyr`
- `ggplot2`
- `ggrepel`
- `rtracklayer`
- `stringr`
- `tidyr`
- `tibble`
- `circlize`
- `ComplexHeatmap`

Some original code snippets also referenced these packages:

- `RSQLite`
- `R.utils`
- `reshape2`
- `qqman`
- `biomaRt`
- `ggvenn`
- `tableone`
- `broom`
- `RColorBrewer`

Not all of them are required by the current finalized CLI scripts, but they may still be useful in older notebooks or exploratory code.

Minimal R install example:

```r
install.packages(c(
  "data.table",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "stringr",
  "tidyr",
  "tibble",
  "circlize"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "rtracklayer",
  "ComplexHeatmap"
))
```

Conda environment note for S-PrediXcan:

- `spredixcan_v8_twas.sh` currently runs:
  `conda run -n metaxcan python ...`
- You therefore need a conda environment named `metaxcan`
- That environment must be able to run `SPrediXcan.py`

## Citation

Core software and methods:

- PLINK:
  Purcell S, et al. PLINK: a tool set for whole-genome association and population-based linkage analyses. *American Journal of Human Genetics* (2007).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/17701901/)

- PLINK 2:
  Chang CC, et al. Second-generation PLINK: rising to the challenge of larger and richer datasets. *GigaScience* (2015).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/25722852/)

- GCTA / COJO:
  Yang J, et al. Conditional and joint multiple-SNP analysis of GWAS summary statistics identifies additional variants influencing complex traits. *Nature Genetics* (2012).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/22426310/)

- FUSION TWAS:
  Gusev A, et al. Integrative approaches for large-scale transcriptome-wide association studies. *Nature Genetics* (2016).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/26854917/)

- PrediXcan:
  Gamazon ER, et al. A gene-based association method for mapping traits using reference transcriptome data. *Nature Genetics* (2015).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/26258848/)

- MetaXcan / S-PrediXcan:
  Barbeira AN, et al. Integrating predicted transcriptome from multiple tissues improves association detection. *PLoS Genetics* (2019).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/30694772/)

Reference data and annotations:

- GTEx Consortium:
  The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* (2020).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/32913098/)

- Ensembl:
  Yates AD, et al. Ensembl 2020. *Nucleic Acids Research* (2020).
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/31691826/)

Project websites and resources:

- GCTA:
  [Official site](https://yanglab.westlake.edu.cn/software/gcta/)

- FUSION:
  [Project page](http://gusevlab.org/projects/fusion/)

- MetaXcan / S-PrediXcan:
  [GitHub](https://github.com/hakyimlab/MetaXcan)

- PredictDB GTEx models:
  [PredictDB](https://predictdb.org/)

- Ensembl FTP:
  [FTP download page](https://www.ensembl.org/info/data/ftp/index.html)

## Script-by-Script Guide

### 1. `make_twb2_pheno_cli.R`

Purpose:
Create TWB2 phenotype and keep files for downstream analysis.

Usage:

```bash
Rscript make_twb2_pheno_cli.R <var_name> <var_type> [case_label] [control_label]
```

Named usage:

```bash
Rscript make_twb2_pheno_cli.R --var_name <var_name> --var_type <categorical|continuous> [--case_label <value>] [--control_label <value>]
```

Parameters:

- `var_name`
  Survey variable name, for example `VERTIGO_SELF` or `BMI`
- `var_type`
  Must be `categorical` or `continuous`
- `case_label`
  Required only for `categorical`
- `control_label`
  Required only for `categorical`

Examples:

```bash
Rscript make_twb2_pheno_cli.R VERTIGO_SELF categorical 1 0
Rscript make_twb2_pheno_cli.R BMI continuous
```

Input files:

- `/mnt/SP-siliconpower/TWB20250806download/survey/release_list_colnames.txt`
- `/mnt/SP-siliconpower/TWB20250806download/survey/release_list_survey.csv`
- `/mnt/SP-siliconpower/TWB20250806download/lab_info/lab_info.csv`
- `/mnt/SP-siliconpower/TWB20250806download/Imputed.120161.TWB2/imputed_120161/TWB2.hg38.impu.v4.fam`

Output files:

- `TWB2_<var_name>_pheno.txt`
- `TWB2_<var_name>_keep.txt`

Output folder:

```text
/mnt/SP-siliconpower/TWB20250806download/Imputed.120161.TWB2/imputed_120161
```

Notes:

- For `categorical`, output coding is `1 = control`, `2 = case`
- For `continuous`, the original numeric value is written

### 2. External GWAS step

Purpose:
Run GWAS outside this repository, typically with PLINK logistic or linear regression.

Expected output folder:

```text
/mnt/data/ai_agent/gene_analysis/PLINK/imputed/glm_logistic/<project_name>/
```

Expected files:

```text
chr1.PHENO1.glm.logistic.hybrid
chr2.PHENO1.glm.logistic.hybrid
...
chr22.PHENO1.glm.logistic.hybrid
```

Expected PLINK binary folder for COJO:

```text
/mnt/data/ai_agent/gene_analysis/PLINK/imputed/plink_binary_files/<project_name>/
```

### 3. `process_gwas_sumstats_cli.R`

Purpose:
Convert GWAS results into downstream input formats for COJO, FUSION, and S-PrediXcan.

Usage:

```bash
Rscript process_gwas_sumstats_cli.R <project_name> [base_dir]
```

Named usage:

```bash
Rscript process_gwas_sumstats_cli.R --project_name <project_name> [--base_dir <base_dir>]
```

Parameters:

- `project_name`
  Project name used in folder naming
- `base_dir`
  Base analysis directory
  Default: `/mnt/data/ai_agent/gene_analysis`

Example:

```bash
Rscript process_gwas_sumstats_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis
```

Input:

- `PLINK/imputed/glm_logistic/<project_name>/chrN.PHENO1.glm.logistic.hybrid`

Outputs:

- `COJO/<project_name>/chrN.cojo.ma`
- `FUSION/GWAS/<project_name>/chrN.sumstats`
- `S_PrediXcan/GWAS/<project_name>/chr1tochr22.sumstats`

Logs:

- `FUSION/GWAS/<project_name>/log/`

### 4. `cojo.sh`

Purpose:
Run GCTA-COJO chromosome by chromosome.

Usage:

```bash
bash cojo.sh PROJECT_NAME PLINK_FILE_PATTERN P_CUTOFF WINDOW [BASE_DIR]
```

Parameters:

- `PROJECT_NAME`
  Project name
- `PLINK_FILE_PATTERN`
  Suffix after `chrN`, for example `_DR2_0.7_QCFiltered`
- `P_CUTOFF`
  COJO p-value cutoff, for example `1e-6`
- `WINDOW`
  COJO window size, for example `10000`
- `BASE_DIR`
  Default: `/mnt/data/ai_agent/gene_analysis`

Example:

```bash
bash cojo.sh TWB1_LAA_test1 _DR2_0.7_QCFiltered 1e-6 10000
```

Expected inputs:

- `PLINK/imputed/plink_binary_files/<project_name>/chrN<PLINK_FILE_PATTERN>.bed`
- `PLINK/imputed/plink_binary_files/<project_name>/chrN<PLINK_FILE_PATTERN>.bim`
- `PLINK/imputed/plink_binary_files/<project_name>/chrN<PLINK_FILE_PATTERN>.fam`
- `COJO/<project_name>/chrN.cojo.ma`

Outputs:

- `COJO/<project_name>/result/chrN_p<P_CUTOFF_TAG>.*`

Logs:

- `COJO/<project_name>/result/log/`

Terminal summary:

- completed chromosomes
- chromosomes missing PLINK bfiles
- chromosomes missing COJO input files

### 5. `cojo_annotation_manhattan_cli.R`

Purpose:
Annotate COJO results with local GTF and generate tables plus GWAS Manhattan plot.

Usage:

```bash
Rscript cojo_annotation_manhattan_cli.R <project_name> [base_dir] [genome_build]
```

Named usage:

```bash
Rscript cojo_annotation_manhattan_cli.R --project_name <project_name> [--base_dir <base_dir>] [--genome_build <grch37|grch38>] [--gtf_file <path>]
```

Parameters:

- `project_name`
- `base_dir`
  Default: `/mnt/data/ai_agent/gene_analysis`
- `genome_build`
  `grch37` or `grch38`
  Default: `grch37`
- `gtf_file`
  Optional explicit GTF path
  If omitted, the script selects a default GTF by `genome_build`

Examples:

```bash
Rscript cojo_annotation_manhattan_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis grch37
Rscript cojo_annotation_manhattan_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis grch38
Rscript cojo_annotation_manhattan_cli.R --project_name TWB1_LAA_test --base_dir /mnt/data/ai_agent/gene_analysis --genome_build grch38
```

Inputs:

- `COJO/<project_name>/result/chrN_p1e5.jma.cojo`
- `PLINK/imputed/glm_logistic/<project_name>/chrN.PHENO1.glm.logistic.hybrid`
- `tools/ensembl/Homo_sapiens.GRCh37.87.gtf` for `grch37`
- `tools/ensembl/Homo_sapiens.GRCh38.115.gtf` for `grch38`

Outputs:

- `COJO/<project_name>/table/`
- `COJO/<project_name>/manhattan/`
- `COJO/<project_name>/log/`

Notes:

- Use `grch37` GTF for GRCh37/hg19-like coordinates
- Use `grch38` GTF for GRCh38/hg38 coordinates

### 6. `fusion_v8_twas.sh`

Purpose:
Run FUSION GTEx v8 TWAS.

Usage:

```bash
bash fusion_v8_twas.sh SUBJECT_NAME [BASE_DIR]
```

Parameters:

- `SUBJECT_NAME`
  Should match the project name used in FUSION GWAS folder
- `BASE_DIR`
  Default: `/mnt/data/ai_agent/gene_analysis`

Example:

```bash
bash fusion_v8_twas.sh TWB1_LAA_test
```

Expected input:

- `FUSION/GWAS/<project_name>/chrN.sumstats`

Outputs:

- `FUSION/result_v8/<project_name>/chrN/`

Logs:

- `FUSION/result_v8/<project_name>/log/`

### 7. `fusion_manhattan_heatmap_cli.R`

Purpose:
Summarize FUSION results into combined tables, Manhattan plot, and heatmap.

Usage:

```bash
Rscript fusion_manhattan_heatmap_cli.R <project_name> [base_dir] [fusion_ver] [fdr_cutoff]
```

Named usage:

```bash
Rscript fusion_manhattan_heatmap_cli.R --project_name <project_name> [--base_dir <base_dir>] [--fusion_ver <v7|v8>] [--fdr_cutoff <value>] [--gencode_v26_gtf <path>] [--gencode_v19_gtf <path>] [--weight_dir_v7 <path>]
```

Parameters:

- `project_name`
- `base_dir`
  Default: `/mnt/data/ai_agent/gene_analysis`
- `fusion_ver`
  `v7` or `v8`
  Default: `v8`
- `fdr_cutoff`
  Default: `0.15`
- `gencode_v26_gtf`
  Default: `/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/GENCODE/gencode.v26.annotation.gtf`
- `gencode_v19_gtf`
  Default: `/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/GENCODE/gencode.v19.annotation.gtf`
- `weight_dir_v7`
  Default: `/home/sysadmin/Desktop/dyc_lab/FUSION/WEIGHTS_v7/GTEx.ALL`

Examples:

```bash
Rscript fusion_manhattan_heatmap_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis v8 0.15
Rscript fusion_manhattan_heatmap_cli.R TWB1_LAA_test /mnt/data/ai_agent/gene_analysis v7 0.15
```

Outputs:

- `FUSION/result_v8/<project_name>/table/`
- `FUSION/result_v8/<project_name>/manhanttan/`
- `FUSION/result_v8/<project_name>/heatmap/`

or:

- `FUSION/result_v7/<project_name>/table/`
- `FUSION/result_v7/<project_name>/manhanttan/`
- `FUSION/result_v7/<project_name>/heatmap/`

Notes:

- Manhattan y-axis is dynamic
- Extremely small P values are handled safely

### 8. `spredixcan_v8_twas.sh`

Purpose:
Run S-PrediXcan GTEx v8 elastic net models.

Usage:

```bash
bash spredixcan_v8_twas.sh SUBJECT_NAME [BASE_DIR]
```

Parameters:

- `SUBJECT_NAME`
  Project name used in `S_PrediXcan/GWAS/<project_name>`
- `BASE_DIR`
  Default: `/mnt/data/ai_agent/gene_analysis`

Example:

```bash
bash spredixcan_v8_twas.sh TWB1_LAA_hg38
```

Expected input:

- `S_PrediXcan/GWAS/<project_name>/chr1tochr22.sumstats`

Outputs:

- `S_PrediXcan/result_v8/<project_name>/`

Notes:

- Uses `conda run -n metaxcan`
- Missing model files are skipped with warnings

### 9. `spredixcan_manhattan_heatmap_cli.R`

Purpose:
Summarize S-PrediXcan results into combined tables, Manhattan plot, and heatmap.

Usage:

```bash
Rscript spredixcan_manhattan_heatmap_cli.R <project_name> [base_dir] [spredixcan_ver] [fdr_cutoff]
```

Named usage:

```bash
Rscript spredixcan_manhattan_heatmap_cli.R --project_name <project_name> [--base_dir <base_dir>] [--spredixcan_ver <v7|v8>] [--fdr_cutoff <value>] [--gencode_v26_gtf <path>] [--gencode_v19_gtf <path>]
```

Parameters:

- `project_name`
- `base_dir`
  Default: `/mnt/data/ai_agent/gene_analysis`
- `spredixcan_ver`
  `v7` or `v8`
  Default: `v8`
- `fdr_cutoff`
  Default: `0.15`
- `gencode_v26_gtf`
  Default: `/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/GENCODE/gencode.v26.annotation.gtf`
- `gencode_v19_gtf`
  Default: `/home/sysadmin/Desktop/dyc_lab/TWB1_stroke_raw_data/GENCODE/gencode.v19.annotation.gtf`

Examples:

```bash
Rscript spredixcan_manhattan_heatmap_cli.R TWB2_VERTIGO /mnt/data/ai_agent/gene_analysis v8 0.15
Rscript spredixcan_manhattan_heatmap_cli.R TWB2_VERTIGO /mnt/data/ai_agent/gene_analysis v7 0.15
```

Outputs:

- `S_PrediXcan/result_v8/<project_name>/table/`
- `S_PrediXcan/result_v8/<project_name>/manhanttan/`
- `S_PrediXcan/result_v8/<project_name>/heatmap/`

or:

- `S_PrediXcan/result_v7/<project_name>/table/`
- `S_PrediXcan/result_v7/<project_name>/manhanttan/`
- `S_PrediXcan/result_v7/<project_name>/heatmap/`

Notes:

- Manhattan y-axis is dynamic
- Extremely small P values are handled safely

## Minimal Pipelines

### COJO pipeline

1. Prepare phenotype if needed using `make_twb2_pheno_cli.R`
2. Run external GWAS
3. Run `process_gwas_sumstats_cli.R`
4. Run `cojo.sh`
5. Run `cojo_annotation_manhattan_cli.R`

### FUSION pipeline

1. Prepare phenotype if needed using `make_twb2_pheno_cli.R`
2. Run external GWAS
3. Run `process_gwas_sumstats_cli.R`
4. Run `fusion_v8_twas.sh`
5. Run `fusion_manhattan_heatmap_cli.R`

### S-PrediXcan pipeline

1. Prepare phenotype if needed using `make_twb2_pheno_cli.R`
2. Run external GWAS
3. Run `process_gwas_sumstats_cli.R`
4. Run `spredixcan_v8_twas.sh`
5. Run `spredixcan_manhattan_heatmap_cli.R`

## Current Script List

- `cojo.sh`
- `cojo_annotation_manhattan_cli.R`
- `fusion_manhattan_heatmap_cli.R`
- `fusion_v8_twas.sh`
- `make_twb2_pheno_cli.R`
- `process_gwas_sumstats_cli.R`
- `spredixcan_manhattan_heatmap_cli.R`
- `spredixcan_v8_twas.sh`
