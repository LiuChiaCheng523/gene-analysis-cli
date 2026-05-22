# Gene Analysis README

This repository contains the operational scripts currently used in:

`C:\Users\bruce\codex_code\gene_analysis`

This README focuses on:

- what each script does
- the recommended execution order
- all CLI / shell parameters
- the overall folder structure expected by the scripts

Excluded from normal workflow:

- `cojo_annotation_manhattan_debug.R`

## Quick Setup

This repository assumes that users can choose any `base_dir`, but the folder
structure under that `base_dir` must follow this repository's expected layout.

### One-Shot Checklist

For a new user, the clearest setup order is:

1. Prepare a Linux environment
   - recommended: `Ubuntu`
   - Windows users: `WSL2 + Ubuntu`
2. Install core system tools
   - `bash`
   - `R` and `Rscript`
   - `python`
   - `conda`
   - `wget`
   - `tar`
   - `dos2unix` (recommended if scripts were edited/transferred from Windows)
3. Install system libraries required to compile R packages
   - see `System Libraries For R Packages` below
   - critical: `libglpk-dev` and `libgmp-dev` (otherwise `clusterProfiler` fails via `igraph`)
4. Install external workflow tools
   - `gcta64` for `cojo.sh`
   - `MetaXcan` under `BASE_DIR/tools/MetaXcan/`
5. Install required R packages
   - recommended: run `install_r_packages.sh` (handles CRAN + Bioconductor + plink2R)
   - or follow `Minimal R install example`
6. Choose a base directory
   - example: `/mnt/data/ai_agent/gene_analysis`
7. Create the folder skeleton
   - run `setup_gene_analysis_dirs.sh`
8. Place reference data into the expected locations
   - `tools/ensembl/`
   - `FUSION/LDREF/`
   - `FUSION/WEIGHTS/`
   - `FUSION/WEIGHTS_v7/`
   - `S_PrediXcan/model/`
9. Place project-specific GWAS inputs into the expected folders
10. Run the workflow scripts in the order listed below

Minimal command sequence:

```bash
BASE_DIR=/your/base_dir

# 1) Create folder skeleton (now also creates tools/gcta and tools/MetaXcan)
bash setup_gene_analysis_dirs.sh "${BASE_DIR}"

# 2) Activate GCTA (make executable + symlink into PATH)
sudo chmod +x "${BASE_DIR}/tools/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64"
sudo ln -sf "${BASE_DIR}/tools/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64" /usr/local/bin/gcta64

# 3) Create the MetaXcan conda environment (named "metaxcan" to match the scripts)
#    NOTE: newer conda may ask you to accept channel Terms of Service first:
#    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
#    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
cd "${BASE_DIR}/tools/MetaXcan/software"
conda env create -n metaxcan -f conda_env.yaml

# 4) IMPORTANT: pin numpy back to 1.26.4
#    The pip stage (bgen-reader/cbgen) upgrades numpy to 2.x, which breaks
#    cyvcf2/scipy with "numpy.dtype size changed" errors. Downgrade fixes it.
#    (bgen-reader stays incompatible, but it is NOT needed for S-PrediXcan.)
conda run -n metaxcan pip install "numpy==1.26.4"

# 5) Install R packages (recommended: use the provided helper script)
bash install_r_packages.sh
```

Alternatively, install R packages manually (see `Minimal R install example` below).
The helper script `install_r_packages.sh` is preferred because it also installs the
required system libraries, installs `plink2R` from GitHub, and puts `KEGGREST` in the
correct (Bioconductor) channel.

To create the required directory skeleton, run:

```bash
bash setup_gene_analysis_dirs.sh /your/base_dir
```

If you already know one or more project names, you can create the common
folders and the project-specific folders at the same time:

```bash
bash setup_gene_analysis_dirs.sh /your/base_dir PROJECT_A PROJECT_B
```

Example:

```bash
bash setup_gene_analysis_dirs.sh /mnt/data/ai_agent/gene_analysis TWB1_LAA TWB2_VERTIGO
```

Important notes:

- `base_dir` is configurable.
- The folder structure inside `base_dir` is not arbitrary and should match this repository.
- Most scripts will stop with an error if required input folders or files are missing.
- `make_twb2_pheno_cli.R` uses its own separate `base_dir` because it reads from the TWB2 release bundle.

## What To Put On GitHub vs Cloud Storage

For shared use, this project is easiest to distribute in 2 layers:

1. GitHub repository for code and setup files
2. Cloud storage or shared drive for large reference data

Recommended to keep in GitHub:

- all `.R` CLI scripts
- all `.sh` workflow scripts
- `README.md`
- `setup_gene_analysis_dirs.sh`
- small example text files

Recommended to distribute separately through cloud storage, shared drive, or institutional storage:

- `FUSION/LDREF/`
- `FUSION/WEIGHTS/`
- `FUSION/WEIGHTS_v7/`
- `S_PrediXcan/model/`
- `tools/ensembl/`
- `tools/GENCODE/`
- any large GWAS intermediate files
- any project-specific result folders

Recommended sharing approach:

1. Users clone the GitHub repository
2. Users run `setup_gene_analysis_dirs.sh` to create the expected folder skeleton
3. Users download large reference data from your cloud storage
4. Users place those reference folders into the correct locations under their chosen `base_dir`

For folders such as `FUSION/fusion_twas`, `FUSION/LDREF`, `FUSION/WEIGHTS`, and `FUSION/WEIGHTS_v7`:

- `fusion_twas` can stay in the repository if it is code you want to ship with the workflow
- `LDREF`, `WEIGHTS`, and `WEIGHTS_v7` are better provided separately as large reference data

If you package large resources for sharing, it is usually cleaner to upload only the required reference subfolders rather than the entire `FUSION/` folder, because `FUSION/GWAS/`, `FUSION/result_v7/`, and `FUSION/result_v8/` are usually user- or project-specific outputs.

## Transferring Reference Data From Windows Into WSL2

If your reference data and model weights currently live on Windows (for example on
the Desktop), you can copy them into your WSL2 `BASE_DIR`. WSL automatically mounts
the Windows `C:` drive at `/mnt/c/`.

### Recommended: keep data inside the WSL filesystem

Put `BASE_DIR` under your Linux home (e.g. `~/gene_analysis_workflow`), NOT under
`/mnt/c/...`. Running PLINK / GCTA / S-PrediXcan against files on `/mnt/c/` is 5-10x
slower because of the Windows-Linux filesystem bridge.

### Set source and destination variables

To save typing, define a Windows source (`SRC`) and a WSL destination (`DST`).
These variables only exist in the current terminal session; if you open a new
terminal, set them again.

```bash
# Windows source (note: C:\Users\you\... becomes /mnt/c/Users/you/...)
SRC=/mnt/c/Users/<your_windows_user>/Desktop/set_up_data
# WSL destination (your chosen BASE_DIR)
DST=~/gene_analysis_workflow
```

### Copy folders with rsync

`rsync` is preferred over `cp` for large transfers: it shows overall progress,
can resume if interrupted, and skips files already copied.

```bash
# General form. KEEP the trailing slash on BOTH paths.
rsync -a --info=progress2 "$SRC/<subpath>/" "$DST/<subpath>/"
```

Trailing-slash rule (important, avoids accidental nested folders):

- `rsync -a SRC/folder/ DST/folder/`  -> copies the *contents* of `folder` into the
  destination `folder`. This is what you almost always want.
- `rsync -a SRC/folder  DST/folder/`  -> copies `folder` *itself* inside, producing
  `DST/folder/folder/` (a nested duplicate). Avoid this.

Examples used to populate this repository's expected layout:

```bash
rsync -a --info=progress2 "$SRC/FUSION/fusion_twas/"          "$DST/FUSION/fusion_twas/"
rsync -a --info=progress2 "$SRC/FUSION/LDREF/LDREF/"          "$DST/FUSION/LDREF/LDREF/"
rsync -a --info=progress2 "$SRC/FUSION/WEIGHTS/"              "$DST/FUSION/WEIGHTS/"
rsync -a --info=progress2 "$SRC/FUSION/WEIGHTS_v7/GTEx.ALL/"  "$DST/FUSION/WEIGHTS_v7/GTEx.ALL/"
rsync -a --info=progress2 "$SRC/S_PrediXcan/model/"           "$DST/S_PrediXcan/model/"
rsync -a --info=progress2 "$SRC/tools/ensembl/"              "$DST/tools/ensembl/"
rsync -a --info=progress2 --exclude='__MACOSX' "$SRC/tools/gcta/"     "$DST/tools/gcta/"
rsync -a --info=progress2 "$SRC/tools/MetaXcan/"             "$DST/tools/MetaXcan/"
```

### Copy single files with cp

```bash
# The small tissue list files belong directly under S_PrediXcan/
cp "$SRC/S_PrediXcan/tissues.txt" "$SRC/S_PrediXcan/tissues_v8_elastic_net.txt" "$DST/S_PrediXcan/"
```

### Verify a transfer

```bash
du -sh "$SRC/<subpath>/" "$DST/<subpath>/"   # sizes should be close
ls -lh "$DST/<subpath>/"                     # spot-check files
```

### When to convert line endings (CRLF -> LF)

Files created or edited on Windows often have CRLF (`\r\n`) line endings. Linux tools
expect LF (`\n`). A stray `\r` causes errors such as
`#!/bin/bash^M: bad interpreter: No such file or directory`.

Convert ONLY plain-text files. NEVER run `dos2unix` on binary files, because it
corrupts them.

| File type | Convert with dos2unix? |
| --- | --- |
| `.sh`, `.R`, `.py` scripts | Yes |
| `.txt` config / lists (e.g. `tissues.txt`) | Yes |
| `.csv`, `.tsv`, `.ma`, `.sumstats` (text) | Usually fine to convert if edited on Windows |
| `.bed` (PLINK binary) | NO - binary |
| `.db` (S-PrediXcan models) | NO - binary |
| `.RDat` / `.wgt.RDat` (FUSION weights) | NO - binary |
| `.png` and other images | NO - binary |
| `.gtf` (Ensembl, text but huge) | Not needed; downloaded files are already LF |

```bash
# Install dos2unix once
sudo apt install -y dos2unix

# Convert specific text files
dos2unix "$DST/S_PrediXcan/tissues.txt" "$DST/S_PrediXcan/tissues_v8_elastic_net.txt"

# Convert all scripts in a folder (text only - safe because the folder has no binaries)
find "$DST/scripts" -type f \( -name '*.sh' -o -name '*.R' -o -name '*.txt' \) -exec dos2unix {} \;

# After converting, give shell scripts execute permission
chmod +x "$DST"/scripts/*.sh
```

Tip: to check whether a file has CRLF, run `file <path>`. Output containing
`CRLF line terminators` means it should be converted; plain `ASCII text` is already LF.

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
6. Run overlap analysis if you want cross-method gene intersections
7. Run pathway enrichment if you want GO BP / KEGG interpretation

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

### F. Overlap branch

10. `overlap_venn_cli.R`

### G. Pathway branch

11. `pathway_enrichment_cli.R`

## Folder Structure

Below is the main directory structure assumed by the current code.

```text
/mnt/data/ai_agent/gene_analysis
|- PLINK/
|  |- imputed/
|  |  |- glm_logistic/
|  |  |  `- <project_name>/
|  |  |     |- chr1.PHENO1.glm.logistic.hybrid
|  |  |     |- chr2.PHENO1.glm.logistic.hybrid
|  |  |     `- ...
|  |  `- plink_binary_files/
|  |     `- <project_name>/
|  |        |- chr1<PLINK_FILE_PATTERN>.bed
|  |        |- chr1<PLINK_FILE_PATTERN>.bim
|  |        |- chr1<PLINK_FILE_PATTERN>.fam
|  |        `- ...
|  `- original/
|- COJO/
|  `- <project_name>/
|     |- chr1.cojo.ma
|     |- chr2.cojo.ma
|     |- ...
|     |- result/
|     |- table/
|     |- manhattan/
|     `- log/
|- FUSION/
|  |- fusion_twas/
|  |- GWAS/
|  |  `- <project_name>/
|  |     |- chr1.sumstats
|  |     |- chr2.sumstats
|  |     `- ...
|  |- LDREF/
|  |  `- LDREF/
|  |- WEIGHTS/
|  |- WEIGHTS_v7/
|  |  `- GTEx.ALL/
|  |- result_v8/
|  |  `- <project_name>/
|  |     |- chr1/
|  |     |- chr2/
|  |     |- ...
|  |     |- table/
|  |     |- manhanttan/
|  |     |- heatmap/
|  |     `- log/
|  `- result_v7/
|     `- <project_name>/
|        |- chr1/
|        |- ...
|        |- table/
|        |- manhanttan/
|        |- heatmap/
|        `- log/
|- S_PrediXcan/
|  |- GWAS/
|  |  `- <project_name>/
|  |     `- chr1tochr22.sumstats
|  |- model/
|  |  |- GTEx_v8/
|  |  |  `- elastic_net_models/
|  |  |- gtex_v7_<tissue>_imputed_europeans_tw_0.5_signif.db
|  |  `- gtex_v7_<tissue>_imputed_eur_covariances.txt.gz
|  |- tissues.txt
|  |- tissues_v8_elastic_net.txt
|  |- result_v8/
|  |  `- <project_name>/
|  |     |- <tissue>_SPrediXcan_v8_en.csv
|  |     |- table/
|  |     |- manhanttan/
|  |     |- heatmap/
|  |     `- log/
|  `- result_v7/
|     `- <project_name>/
|        |- <tissue>_SPrediXcan_v7_en.csv
|        |- table/
|        |- manhanttan/
|        |- heatmap/
|        `- log/
|- overlap/
|  `- <project_name_overlap>/
|     |- overlap_summary.csv
|     |- overlap_gene_name.csv
|     |- overlap_gene_name_pairwise_union.csv
|     |- twas_detected_gene_name.csv
|     `- venn_<methods>_FDRxx.png
|- pathway/
|  |- GOBP/
|  |  `- <project_name_overlap>/
|  |     |- <project_name_overlap>_GOBP_ORA_FDRxx.png
|  |     |- <project_name_overlap>_GOBP_ORA_FDRxx_summary.csv
|  |     `- <project_name_overlap>_GOBP_ORA_FDRxx_full.csv
|  |- KEGG/
|  |  `- <project_name_overlap>/
|  |     |- <project_name_overlap>_KEGG_ORA_FDRxx.png
|  |     |- <project_name_overlap>_KEGG_ORA_FDRxx_summary.csv
|  |     `- <project_name_overlap>_KEGG_ORA_FDRxx_full.csv
|  `- log/
|     `- <project_name_overlap>/
|        `- process_log_*.txt
`- tools/
   |- ensembl/
   |  |- Homo_sapiens.GRCh37.87.gtf
   |  `- Homo_sapiens.GRCh38.115.gtf
   |- GENCODE/
   |- gcta/
   `- MetaXcan/
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
- Windows 10/11 with WSL2 + Ubuntu
- `bash`
- `Rscript`
- `python`
- `conda`

Recommended runtime choice for Windows users:

- Do not run this workflow directly in native Windows PowerShell or Command Prompt unless you are prepared to adapt paths and shell scripts manually.
- The easiest and most reproducible option is:
  - install `WSL2`
  - install an `Ubuntu` distribution inside WSL
  - clone this repository inside the WSL environment
  - run all `.sh` and `.R` scripts from Ubuntu inside WSL

Typical setup for Windows users:

1. Install WSL2
2. Install Ubuntu from Microsoft Store
3. Install `R`, `python`, and `conda` inside Ubuntu
4. Install required external tools such as `gcta64`
5. Clone this repository inside WSL
6. Run:

```bash
bash setup_gene_analysis_dirs.sh /your/base_dir
```

This project is primarily designed around Linux-style paths and shell tools such as:

- `/mnt/...`
- `bash`
- `wget`
- `find`
- `tar`

Using WSL2 helps keep the runtime environment close to the original development environment.

External command-line tools required:

- `PLINK` or `PLINK2`
  Used outside this repository for GWAS
- `gcta64`
  Required by `cojo.sh`
- `python` inside the `metaxcan` conda environment
  Required by `spredixcan_v8_twas.sh`
- `wget`
  Useful for downloading GTF and reference files

Example installation commands for tools used by `cojo.sh`, `spredixcan_v7_twas.sh`, and `spredixcan_v8_twas.sh`:

```bash
BASE_DIR=/your/base_dir

sudo chmod +x "${BASE_DIR}/tools/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64"
sudo ln -sf "${BASE_DIR}/tools/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64" /usr/local/bin/gcta64

# Newer conda (>=26) may require accepting channel Terms of Service first:
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

cd "${BASE_DIR}/tools/MetaXcan/software"
conda env create -n metaxcan -f conda_env.yaml

# REQUIRED post-install fix: pin numpy back to 1.26.4
# The pip stage (bgen-reader/cbgen) pulls numpy 2.x, which breaks cyvcf2/scipy
# binary compatibility ("numpy.dtype size changed"). This restores a working env.
conda run -n metaxcan pip install "numpy==1.26.4"

which gcta64
gcta64 --help
conda run -n metaxcan python "${BASE_DIR}/tools/MetaXcan/software/SPrediXcan.py" --help
conda run -n metaxcan python -c "import pandas, scipy, numpy, h5py, cyvcf2; print('metaxcan core OK', numpy.__version__)"
```

Notes:

- These commands assume users have already copied or downloaded the official `gcta` and `MetaXcan` packages into `BASE_DIR/tools/`.
- `cojo.sh` expects `gcta64` to be available in `PATH`.
- `spredixcan_v7_twas.sh` and `spredixcan_v8_twas.sh` expect a conda environment named `metaxcan`.
- The bundled `conda_env.yaml` declares `name: imlabtools`; passing `-n metaxcan` overrides it so the environment name matches what the scripts call (`conda run -n metaxcan`).
- The numpy downgrade leaves `bgen-reader`/`cbgen` in a broken state (they want numpy >= 2.0), but these are only needed for BGEN genotype input and are NOT used by the S-PrediXcan summary-statistics workflow in this repository.
- `gcta64` does not support `--version`; running it prints the version banner before reporting that the option is invalid, which still confirms it works.

Required local reference/data resources:

- `tools/ensembl/Homo_sapiens.GRCh37.87.gtf`
- `tools/ensembl/Homo_sapiens.GRCh38.115.gtf`
- `FUSION/WEIGHTS/`
- `FUSION/WEIGHTS_v7/`
- `FUSION/LDREF/`
- `FUSION/fusion_twas/`
- `tools/MetaXcan/software/`
- `S_PrediXcan/model/GTEx_v8/elastic_net_models/`

FUSION-specific requirements for `fusion_v7_twas.sh` and `fusion_v8_twas.sh`:

- `Rscript` must be available
- official FUSION script must exist at:
  - `BASE_DIR/FUSION/fusion_twas/FUSION.assoc_test.R`
- required R packages:
  - `plink2R`
  - `optparse`
- optional R packages used only for advanced flags inside the official script:
  - `coloc`
  - `jlimR`

Minimal install example for the standard FUSION workflow:

```bash
Rscript -e "install.packages(c('plink2R','optparse'), repos='https://cloud.r-project.org')"
Rscript -e "library(plink2R); library(optparse); cat('FUSION R packages OK\n')"
```

Placement expectation for FUSION:

- `fusion_v7_twas.sh` and `fusion_v8_twas.sh` do not require a globally installed `fusion` command
- they directly call:
  - `BASE_DIR/FUSION/fusion_twas/FUSION.assoc_test.R`
- users therefore need to place the official FUSION repository files under `BASE_DIR/FUSION/fusion_twas/`

Recommended R packages used across the repository:

- `data.table`
- `dplyr`
- `ggplot2`
- `ggrepel`
- `rtracklayer`
- `stringr`
- `tidyr`
- `tibble`
- `ggvenn`
- `circlize`
- `ComplexHeatmap`
- `clusterProfiler`
- `optparse`
- `plink2R`
- `GO.db`
- `org.Hs.eg.db`
- `AnnotationDbi`
- `KEGGREST`

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
- `GO.db`
- `org.Hs.eg.db`
- `AnnotationDbi`

Not all of them are required by the current finalized CLI scripts, but they may still be useful in older notebooks or exploratory code. The package list above covers the current formal workflow scripts in this repository.

### System Libraries For R Packages

Several R packages compile from source and need system libraries. On Ubuntu/WSL2,
install these BEFORE installing the R packages:

```bash
sudo apt update
sudo apt install -y \
  libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libcurl4-openssl-dev libssl-dev libgit2-dev gfortran \
  libglpk-dev libgmp-dev cmake make g++
```

- `libglpk-dev` and `libgmp-dev` are required by `igraph`, a dependency of
  `clusterProfiler`. Without them, `clusterProfiler` installation fails.
- `libgit2-dev` is needed to install `plink2R` from GitHub via `remotes`.
- The `lib*-dev` graphics libraries are needed by `ggplot2`/`ComplexHeatmap` backends.

### Recommended: one-shot R package installer

The repository includes `install_r_packages.sh`, which installs the system
libraries above, all CRAN packages, `plink2R` (from GitHub), and all Bioconductor
packages. It is idempotent (already-installed packages are skipped):

```bash
bash install_r_packages.sh
```

### Minimal R install example (manual)

Note two corrections compared to a naive package list:

- `plink2R` is NOT on CRAN; install it from GitHub.
- `KEGGREST` is a Bioconductor package, NOT a CRAN package.

```r
# CRAN packages
install.packages(c(
  "remotes",
  "data.table",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "stringr",
  "tidyr",
  "tibble",
  "circlize",
  "ggvenn",
  "optparse"
), repos = "https://cloud.r-project.org")

# plink2R (FUSION dependency) — from GitHub, not CRAN
remotes::install_github("gabraham/plink2R/plink2R", upgrade = "never")

# Bioconductor packages (KEGGREST belongs here, not in CRAN)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(c(
  "rtracklayer",
  "ComplexHeatmap",
  "clusterProfiler",
  "GO.db",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "KEGGREST"
), update = FALSE, ask = FALSE)
```

Conda environment note for S-PrediXcan:

- `spredixcan_v8_twas.sh` currently runs:
  `conda run -n metaxcan python ...`
- You therefore need a conda environment named `metaxcan`
- That environment must be able to run `SPrediXcan.py`

## Dependency Matrix

This section summarizes what must already be installed before running each formal workflow script.

### Shell scripts

- `cojo.sh`
  - `bash`
  - `gcta64` available in `PATH`
  - PLINK binary files already prepared under `PLINK/imputed/plink_binary_files/<project_name>/`

- `fusion_v7_twas.sh`
  - `bash`
  - `Rscript`
  - official FUSION script at `BASE_DIR/FUSION/fusion_twas/FUSION.assoc_test.R`
  - R packages:
    - `plink2R`
    - `optparse`
  - reference folders:
    - `BASE_DIR/FUSION/WEIGHTS_v7/GTEx.ALL/`
    - `BASE_DIR/FUSION/LDREF/LDREF/`

- `fusion_v8_twas.sh`
  - `bash`
  - `Rscript`
  - official FUSION script at `BASE_DIR/FUSION/fusion_twas/FUSION.assoc_test.R`
  - R packages:
    - `plink2R`
    - `optparse`
  - reference folders:
    - `BASE_DIR/FUSION/WEIGHTS/`
    - `BASE_DIR/FUSION/LDREF/LDREF/`

- `spredixcan_v7_twas.sh`
  - `bash`
  - `conda`
  - conda environment named `metaxcan`
  - official MetaXcan script at `BASE_DIR/tools/MetaXcan/software/SPrediXcan.py`
  - model files under `BASE_DIR/S_PrediXcan/model/`

- `spredixcan_v8_twas.sh`
  - `bash`
  - `conda`
  - conda environment named `metaxcan`
  - official MetaXcan script at `BASE_DIR/tools/MetaXcan/software/SPrediXcan.py`
  - model files under `BASE_DIR/S_PrediXcan/model/GTEx_v8/elastic_net_models/`

- `setup_gene_analysis_dirs.sh`
  - `bash`
  - no extra package installation required

### R CLI scripts

- `process_gwas_sumstats_cli.R`
  - `Rscript`
  - R packages:
    - `data.table`

- `cojo_annotation_manhattan_cli.R`
  - `Rscript`
  - R packages:
    - `data.table`
    - `dplyr`
    - `ggplot2`
    - `ggrepel`
    - `rtracklayer`
  - local GTF under `BASE_DIR/tools/ensembl/`

- `fusion_manhattan_heatmap_cli.R`
  - `Rscript`
  - R packages:
    - `data.table`
    - `dplyr`
    - `stringr`
    - `ggplot2`
    - `tidyr`
    - `ggrepel`
    - `rtracklayer`
    - `tibble`
    - `circlize`
    - `ComplexHeatmap`

- `spredixcan_manhattan_heatmap_cli.R`
  - `Rscript`
  - R packages:
    - `data.table`
    - `dplyr`
    - `ggplot2`
    - `tidyr`
    - `ggrepel`
    - `rtracklayer`
    - `tibble`
    - `circlize`
    - `ComplexHeatmap`

- `overlap_venn_cli.R`
  - `Rscript`
  - R packages:
    - `data.table`
    - `dplyr`
    - `stringr`
    - `ggplot2`
    - `ggvenn`

- `pathway_enrichment_cli.R`
  - `Rscript`
  - R packages:
    - `data.table`
    - `dplyr`
    - `ggplot2`
    - `clusterProfiler`
    - `org.Hs.eg.db`
    - `AnnotationDbi`
    - `GO.db`
    - `KEGGREST`

- `make_twb2_pheno_cli.R`
  - `Rscript`
  - R packages:
    - `data.table`
    - `dplyr`

Excluded debug script:

- `cojo_annotation_manhattan_debug.R`
  - not part of the formal workflow

## Citation

If you use this workflow in analysis, manuscript writing, presentation, or reporting, you should cite the underlying methods and resources rather than only this repository.

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
  Default: `<base_dir>/tools/GENCODE/gencode.v26.annotation.gtf`
- `gencode_v19_gtf`
  Default: `<base_dir>/tools/GENCODE/gencode.v19.annotation.gtf`
- `weight_dir_v7`
  Default: `<base_dir>/FUSION/WEIGHTS_v7/GTEx.ALL`

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
  Default: `<base_dir>/tools/GENCODE/gencode.v26.annotation.gtf`
- `gencode_v19_gtf`
  Default: `<base_dir>/tools/GENCODE/gencode.v19.annotation.gtf`

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

### 10. `overlap_venn_cli.R`

Purpose:
Calculate gene overlap across available methods and generate overlap summary tables plus Venn plot.

Usage:

```bash
Rscript overlap_venn_cli.R --project_name_overlap <name> --base_dir <dir> [options]
```

Named parameters:

- `--project_name_v7 <name|None>`
  Project name for v7 TWAS tables
- `--project_name_v8 <name|None>`
  Project name for v8 TWAS tables
- `--project_name_cojo <name|None>`
  Project name for COJO annotation table
- `--project_name_overlap <name>`
  Output folder name under `overlap/`
- `--base_dir <dir>`
  Default: `/mnt/data/ai_agent/gene_analysis`
- `--ver <v7|v8|all>`
  Version selection for TWAS tables
  Default: `all`
- `--twas_method <fusion|spredixcan|all>`
  Method selection for TWAS tables
  Default: `all`
- `--fdr_cutoff <value>`
  Default: `0.15`
- `--cojo_p_cutoff <value>`
  Default: `1e-6`

Examples:

```bash
Rscript overlap_venn_cli.R \
  --project_name_v7 None \
  --project_name_v8 TWB1_LAA_hg38 \
  --project_name_cojo TWB1_LAA \
  --project_name_overlap TWB1_LAA \
  --base_dir /mnt/data/ai_agent/gene_analysis \
  --ver v8 \
  --twas_method all \
  --fdr_cutoff 0.2 \
  --cojo_p_cutoff 1e-5

Rscript overlap_venn_cli.R \
  --project_name_v7 None \
  --project_name_v8 TWB1_LAA_hg38 \
  --project_name_cojo TWB1_LAA \
  --project_name_overlap TWB1_LAA \
  --base_dir /mnt/data/ai_agent/gene_analysis \
  --ver v8 \
  --twas_method spredixcan \
  --fdr_cutoff 0.2 \
  --cojo_p_cutoff 1e-5
```

Inputs:

- `FUSION/result_v7/<project_name>/table/FUSION_v7_TWAS.csv`
- `FUSION/result_v8/<project_name>/table/FUSION_v8_TWAS.csv`
- `S_PrediXcan/result_v7/<project_name>/table/SPrediXcan_v7_TWAS.csv`
- `S_PrediXcan/result_v8/<project_name>/table/SPrediXcan_v8_TWAS.csv`
- `COJO/<project_name>/table/p1e5_annotation_genes.csv` or `p1e05_annotation_genes.csv`

Outputs:

- `overlap/<project_name_overlap>/overlap_summary.csv`
- `overlap/<project_name_overlap>/overlap_gene_name.csv`
- `overlap/<project_name_overlap>/overlap_gene_name_<METHOD1>_<METHOD2>.csv`
- `overlap/<project_name_overlap>/overlap_gene_name_pairwise_union.csv`
- `overlap/<project_name_overlap>/twas_detected_gene_name.csv`
- `overlap/<project_name_overlap>/union_gene_name_<METHOD1>_<METHOD2>.csv`
- `overlap/<project_name_overlap>/venn_<methods>_FDRxx.png`

Notes:

- The script uses whatever methods are actually available
- At least two valid methods are required
- If `twas_method=all` and all inputs exist, the Venn plot includes `FUSION`, `SPrediXcan`, and `GWAS`
- If only one TWAS method plus COJO is available, it generates a two-set Venn plot
- Pairwise intersection and pairwise union gene tables are also exported for every available method pair
- `overlap_gene_name_pairwise_union.csv` means the union of all pairwise overlaps, i.e. genes appearing in at least two methods
- `twas_detected_gene_name.csv` contains all unique genes detected by the selected TWAS method(s) and version(s)

### 11. `pathway_enrichment_cli.R`

Purpose:
Run GO Biological Process and KEGG over-representation analysis from overlap genes and TWAS universe genes.

Usage:

```bash
Rscript pathway_enrichment_cli.R <project_name_overlap> [base_dir] [fdr_cutoff]
```

Named usage:

```bash
Rscript pathway_enrichment_cli.R --project_name_overlap <name> [--base_dir <dir>] [--fdr_cutoff <value>] [--gene_file <path>] [--universe_file <path>]
```

Parameters:

- `project_name_overlap`
  Overlap project name under `overlap/`
- `base_dir`
  Default: `/mnt/data/ai_agent/gene_analysis`
- `fdr_cutoff`
  Used in output file naming
  Default: `0.15`
- `gene_file`
  Default: `overlap/<project_name_overlap>/overlap_gene_name_pairwise_union.csv`
- `universe_file`
  Default: `overlap/<project_name_overlap>/twas_detected_gene_name.csv`

Examples:

```bash
Rscript pathway_enrichment_cli.R TWB1_LAA /mnt/data/ai_agent/gene_analysis 0.2

Rscript pathway_enrichment_cli.R \
  --project_name_overlap TWB1_LAA \
  --base_dir /mnt/data/ai_agent/gene_analysis \
  --fdr_cutoff 0.2
```

Inputs:

- `overlap/<project_name_overlap>/overlap_gene_name_pairwise_union.csv`
- `overlap/<project_name_overlap>/twas_detected_gene_name.csv`

Outputs:

- `pathway/GOBP/<project_name_overlap>/<project_name_overlap>_GOBP_ORA_FDRxx.png`
- `pathway/GOBP/<project_name_overlap>/<project_name_overlap>_GOBP_ORA_FDRxx_summary.csv`
- `pathway/GOBP/<project_name_overlap>/<project_name_overlap>_GOBP_ORA_FDRxx_full.csv`
- `pathway/KEGG/<project_name_overlap>/<project_name_overlap>_KEGG_ORA_FDRxx.png`
- `pathway/KEGG/<project_name_overlap>/<project_name_overlap>_KEGG_ORA_FDRxx_summary.csv`
- `pathway/KEGG/<project_name_overlap>/<project_name_overlap>_KEGG_ORA_FDRxx_full.csv`
- `pathway/log/<project_name_overlap>/process_log_*.txt`

Output folder structure:

```text
pathway/
|- GOBP/
|  `- <project_name_overlap>/
|     |- <project_name_overlap>_GOBP_ORA_FDRxx.png
|     |- <project_name_overlap>_GOBP_ORA_FDRxx_summary.csv
|     `- <project_name_overlap>_GOBP_ORA_FDRxx_full.csv
|- KEGG/
|  `- <project_name_overlap>/
|     |- <project_name_overlap>_KEGG_ORA_FDRxx.png
|     |- <project_name_overlap>_KEGG_ORA_FDRxx_summary.csv
|     `- <project_name_overlap>_KEGG_ORA_FDRxx_full.csv
`- log/
   `- <project_name_overlap>/
      `- process_log_*.txt
```

Notes:

- GO analysis uses `enrichGO(..., ont = "BP")`
- KEGG analysis uses `enrichKEGG()`
- The universe is the selected TWAS-detected genes
- The candidate set defaults to the union of all pairwise overlap genes
- The top 20 pathways by `RichFactor` are plotted when results are available

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

### Overlap pipeline

1. Finish at least two of the following:
   `cojo_annotation_manhattan_cli.R`, `fusion_manhattan_heatmap_cli.R`, `spredixcan_manhattan_heatmap_cli.R`
2. Run `overlap_venn_cli.R`

### Pathway pipeline

1. Finish overlap analysis with `overlap_venn_cli.R`
2. Run `pathway_enrichment_cli.R`

## Current Script List

- `cojo.sh`
- `cojo_annotation_manhattan_cli.R`
- `fusion_manhattan_heatmap_cli.R`
- `fusion_v8_twas.sh`
- `make_twb2_pheno_cli.R`
- `overlap_venn_cli.R`
- `pathway_enrichment_cli.R`
- `process_gwas_sumstats_cli.R`
- `spredixcan_manhattan_heatmap_cli.R`
- `spredixcan_v8_twas.sh`
