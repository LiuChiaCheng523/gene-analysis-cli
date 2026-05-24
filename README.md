# Gene Analysis CLI

End-to-end TWAS pipeline for human GWAS summary statistics. Three downstream
analysis branches — **COJO** (GCTA), **FUSION TWAS**, and **S-PrediXcan** — followed
by cross-method overlap and pathway enrichment.

Designed for **Ubuntu 22.04** or **WSL2 + Ubuntu 22.04** on Windows.
Large reference data (FUSION weights, GTEx models, GENCODE GTFs) is distributed
separately from this repository; see
[Distribution: Code vs Reference Data](#distribution-code-vs-reference-data).

## Contents

- [Quick Bootstrap (New Machine)](#quick-bootstrap-new-machine)
- [Transferring Files from Windows](#transferring-files-from-windows)
- [Distribution: Code vs Reference Data](#distribution-code-vs-reference-data)
- [Folder Structure](#folder-structure)
- [Workflow Overview](#workflow-overview)
- [CLI Script Reference](#cli-script-reference)
- [Pipeline Recipes](#pipeline-recipes)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## Quick Bootstrap (New Machine)

Setting up a fresh machine takes seven steps. Run them in order. Every script in
this repository is **idempotent** — re-running it skips work that is already done.

### Step 0: Install WSL2 + Ubuntu (Windows users only)

In **PowerShell run as Administrator**:

```powershell
wsl --install -d Ubuntu-22.04
```

Reboot when prompted. Launch Ubuntu from the Start menu and create your Linux
user account + password on first launch. Linux users can skip this step.

### Step 1: Clone this repository

```bash
sudo apt update && sudo apt install -y git
cd ~
git clone https://github.com/LiuChiaCheng523/gene-analysis-cli.git
cd ~/gene-analysis-cli
chmod +x *.sh
```

### Step 2: Install base tools + R + Miniconda

```bash
bash setup_genomics_package.sh
# Reload shell so 'conda' becomes available (Miniconda init wrote to ~/.bashrc):
source ~/.bashrc
```

Installs: PLINK 1.9/2.0, HTSlib, Samtools, Bcftools, VCFtools, R 4.5.1 (via rig),
Miniconda (to `~/miniconda3`). Takes ~10 min on a fresh machine.

### Step 3: Install R packages

```bash
bash install_r_packages.sh
```

Installs system libraries required for compilation, all CRAN packages, `plink2R`
from GitHub, and all Bioconductor packages (`rtracklayer`, `ComplexHeatmap`,
`clusterProfiler`, `GO.db`, `org.Hs.eg.db`, `AnnotationDbi`, `KEGGREST`). Takes
20–40 min on first run because Bioconductor packages compile from source.

If the verification table shows `[MISSING]` for a Bioconductor package (usually
`clusterProfiler` due to parallel-build ordering or a transient download timeout),
simply re-run the script — its dependencies are now in place and the second pass
finishes cleanly.

### Step 4: Create the analysis directory skeleton

```bash
BASE_DIR=~/gene_analysis_workflow
bash setup_gene_analysis_dirs.sh "$BASE_DIR"
```

Optionally pass one or more project names to also create their per-project
sub-folders under COJO / FUSION / S_PrediXcan / overlap / pathway:

```bash
bash setup_gene_analysis_dirs.sh "$BASE_DIR" TWB1_LAA TWB2_VERTIGO
```

See [Folder Structure](#folder-structure) for what gets created.

### Step 5: Transfer reference data and external tools into BASE_DIR

The ~38 GB reference bundle (FUSION weights, GTEx models, GENCODE GTFs, GCTA
binary, MetaXcan code) is distributed separately. Set `SRC` to wherever your
bundle lives — typically a Windows folder mounted at `/mnt/c/...` or an external
drive at `/mnt/<letter>/...`.

```bash
SRC=/mnt/c/Users/<your_user>/Desktop/set_up_data   # adjust to your data source
DST="$BASE_DIR"

rsync -a --info=progress2 "$SRC/FUSION/fusion_twas/"          "$DST/FUSION/fusion_twas/"
rsync -a --info=progress2 "$SRC/FUSION/LDREF/LDREF/"          "$DST/FUSION/LDREF/LDREF/"
rsync -a --info=progress2 "$SRC/FUSION/WEIGHTS/"              "$DST/FUSION/WEIGHTS/"
rsync -a --info=progress2 "$SRC/FUSION/WEIGHTS_v7/GTEx.ALL/"  "$DST/FUSION/WEIGHTS_v7/GTEx.ALL/"
rsync -a --info=progress2 "$SRC/S_PrediXcan/model/"           "$DST/S_PrediXcan/model/"
cp "$SRC/S_PrediXcan/tissues.txt" "$SRC/S_PrediXcan/tissues_v8_elastic_net.txt" "$DST/S_PrediXcan/"
rsync -a --info=progress2 "$SRC/tools/ensembl/"               "$DST/tools/ensembl/"
rsync -a --info=progress2 "$SRC/tools/GENCODE/"               "$DST/tools/GENCODE/"
rsync -a --info=progress2 --exclude='__MACOSX' "$SRC/tools/gcta/" "$DST/tools/gcta/"
rsync -a --info=progress2 "$SRC/tools/MetaXcan/"              "$DST/tools/MetaXcan/"
```

`FUSION/WEIGHTS/` alone is ~22 GB across ~2.2 million small files; expect it to
take a long time when read from `/mnt/c/...`. See
[Transferring Files from Windows](#transferring-files-from-windows) for the
trailing-slash rule and the line-ending policy.

### Step 6: Activate GCTA

```bash
sudo chmod +x "$BASE_DIR/tools/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64"
sudo ln -sf "$BASE_DIR/tools/gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64" /usr/local/bin/gcta64
```

`gcta64` does not support `--version`; running it prints the version banner and
then reports `invalid option`, which still confirms it works.

### Step 7: Create the MetaXcan conda environment

```bash
# Newer conda (>= 26) may require accepting channel Terms of Service first:
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r

cd "$BASE_DIR/tools/MetaXcan/software"
conda env create -n metaxcan -f conda_env.yaml

# REQUIRED: pin numpy back to 1.26.4.  The pip stage pulls numpy 2.x via
# bgen-reader/cbgen, which breaks cyvcf2/scipy with
# "numpy.dtype size changed, may indicate binary incompatibility".
conda run -n metaxcan pip install "numpy==1.26.4"
```

Verify:

```bash
conda run -n metaxcan python -c "import pandas, scipy, numpy, h5py, cyvcf2; print('metaxcan core OK', numpy.__version__)"
conda run -n metaxcan python "$BASE_DIR/tools/MetaXcan/software/SPrediXcan.py" --help | head -3
```

The `metaxcan` environment leaves `bgen-reader` / `cbgen` in a broken state
(they require `numpy >= 2.0`), but those are only needed for BGEN genotype input
and are **not** used by the S-PrediXcan summary-statistics workflow in this
repository.

The bundled `conda_env.yaml` declares `name: imlabtools`; passing `-n metaxcan`
overrides it so the environment name matches what `spredixcan_v*_twas.sh` calls
with `conda run -n metaxcan`.

**Bootstrap complete.** See [Pipeline Recipes](#pipeline-recipes) to run analyses.

---

## Transferring Files from Windows

WSL2 mounts the Windows `C:` drive at `/mnt/c/` (other drives at `/mnt/d/`,
`/mnt/f/`, etc.). Keep `BASE_DIR` inside the WSL filesystem (for example
`~/gene_analysis_workflow`); analyses against `/mnt/c/...` run 5-10x slower
because of the Windows-Linux filesystem bridge.

### rsync for folders, cp for one-off files

```bash
SRC=/mnt/c/Users/<your_user>/Desktop/set_up_data
DST=~/gene_analysis_workflow

# General form. KEEP the trailing slash on BOTH paths.
rsync -a --info=progress2 "$SRC/<subpath>/" "$DST/<subpath>/"

# cp for a small text file
cp "$SRC/some_file.txt" "$DST/destination/"
```

Trailing-slash rule (avoids accidental nested folders):

| Form | Result |
| --- | --- |
| `rsync -a SRC/folder/ DST/folder/` | Copies *contents* of `folder` into the destination. Usually what you want. |
| `rsync -a SRC/folder DST/folder/` | Creates `DST/folder/folder/` (a nested duplicate). Avoid. |

### Verify a transfer

```bash
du -sh "$SRC/<subpath>/" "$DST/<subpath>/"   # sizes should be close
ls -lh "$DST/<subpath>/"                     # spot-check the files
```

### Line endings (CRLF → LF)

Windows-edited text files often have CRLF line endings; Linux tools expect LF.
A stray `\r` causes errors like `#!/bin/bash^M: bad interpreter`.
**Never run `dos2unix` on binary files — it corrupts them.**

| File type | Convert with dos2unix? |
| --- | --- |
| `.sh`, `.R`, `.py` scripts | Yes |
| `.txt` config / lists (e.g. `tissues.txt`) | Yes (only if edited on Windows) |
| `.csv`, `.tsv`, `.ma`, `.sumstats` text data | Yes only if edited on Windows |
| `.bed` (PLINK binary) | NO – binary |
| `.db` (S-PrediXcan models) | NO – binary |
| `.RDat` / `.wgt.RDat` (FUSION weights) | NO – binary |
| `.png` and other images | NO – binary |
| `.gtf` (Ensembl / GENCODE) | Not needed; downloaded files are already LF |

```bash
sudo apt install -y dos2unix
dos2unix "$DST/S_PrediXcan/tissues.txt" "$DST/S_PrediXcan/tissues_v8_elastic_net.txt"

# Convert all script-type text in a folder (no binaries in scripts/)
find "$DST/scripts" -type f \( -name '*.sh' -o -name '*.R' -o -name '*.txt' \) \
    -exec dos2unix {} \;
chmod +x "$DST"/scripts/*.sh

# Inspect line endings of an individual file:
file <path>      # output containing "CRLF line terminators" → convert
```

Files cloned with `git` on Linux are already LF, so `dos2unix` is not needed for
this repo's own scripts.

---

## Distribution: Code vs Reference Data

| Layer | Where it lives | Examples |
| --- | --- | --- |
| Code, setup scripts, documentation | This GitHub repo | `*.sh`, `*.R`, `README.md` |
| Small configs / tissue lists (~1 KB each) | Bundled with the reference data | `tissues.txt`, `tissues_v8_elastic_net.txt` |
| Large reference data | Cloud storage / shared drive / external drive | `FUSION/WEIGHTS/`, `FUSION/WEIGHTS_v7/`, `FUSION/LDREF/`, `S_PrediXcan/model/`, `tools/ensembl/`, `tools/GENCODE/`, `tools/gcta/`, `tools/MetaXcan/` |
| Project-specific inputs and outputs | Per-user, never committed | `GWAS/<project>/`, `result_v8/<project>/`, etc. |

Recommended sharing approach for collaborators:

1. Clone this repository.
2. Run the bootstrap (steps 1-4) to install tools and create the folder skeleton.
3. Obtain the large reference bundle separately and copy it into place
   (Step 5 of the bootstrap).
4. Place project-specific GWAS inputs into the expected sub-folders.

---

## Folder Structure

```text
BASE_DIR (e.g. ~/gene_analysis_workflow)
├── PLINK/
│   ├── imputed/
│   │   ├── glm_logistic/<project>/<chrN>.PHENO1.glm.logistic.hybrid
│   │   └── plink_binary_files/<project>/<chrN><PATTERN>.{bed,bim,fam}
│   └── original/
├── COJO/<project>/
│   ├── <chrN>.cojo.ma                     # input for cojo.sh
│   ├── result/                            # raw COJO output (.jma.cojo etc.)
│   ├── table/                             # annotated tables (p1e5_*.csv)
│   ├── manhattan/                         # PNG plot
│   └── log/
├── FUSION/
│   ├── fusion_twas/                       # official FUSION code + FUSION.assoc_test.R
│   ├── GWAS/<project>/<chrN>.sumstats
│   ├── LDREF/LDREF/                       # 1000G EUR reference (~180 MB)
│   ├── WEIGHTS/                           # GTEx v8 weights (~22 GB, many small .wgt.RDat)
│   ├── WEIGHTS_v7/GTEx.ALL/               # GTEx v7 weights
│   ├── result_v8/<project>/<chrN>/  +  table/ manhanttan/ heatmap/ log/
│   └── result_v7/<project>/<chrN>/  +  table/ manhanttan/ heatmap/ log/
├── S_PrediXcan/
│   ├── GWAS/<project>/chr1tochr22.sumstats
│   ├── model/
│   │   ├── GTEx_v8/elastic_net_models/<en_*.db,en_*.txt.gz>
│   │   ├── gtex_v7_<tissue>_imputed_europeans_tw_0.5_signif.db
│   │   └── gtex_v7_<tissue>_imputed_eur_covariances.txt.gz
│   ├── tissues.txt                        # v7 tissue list
│   ├── tissues_v8_elastic_net.txt         # v8 tissue list
│   ├── result_v8/<project>/<tissue>_SPrediXcan_v8_en.csv  +  table/ manhanttan/ heatmap/ log/
│   └── result_v7/<project>/<tissue>_SPrediXcan_v7_en.csv  +  table/ manhanttan/ heatmap/ log/
├── overlap/<project_overlap>/
│   ├── overlap_summary.csv
│   ├── overlap_gene_name.csv
│   ├── overlap_gene_name_pairwise_union.csv
│   ├── twas_detected_gene_name.csv
│   └── venn_<methods>_FDRxx.png
├── pathway/
│   ├── GOBP/<project_overlap>/<*_GOBP_ORA_FDRxx>.{png,csv}
│   ├── KEGG/<project_overlap>/<*_KEGG_ORA_FDRxx>.{png,csv}
│   └── log/<project_overlap>/process_log_*.txt
├── scripts/                               # empty; copy working scripts here if desired
├── tools/
│   ├── gcta/gcta-1.95.0-linux-kernel-3-x86_64/gcta64
│   ├── MetaXcan/software/SPrediXcan.py + conda_env.yaml
│   ├── ensembl/Homo_sapiens.GRCh{37,38}.*.gtf
│   └── GENCODE/gencode.v{19,26}.annotation.gtf
└── VCF/{imputed,original}/
```

`make_twb2_pheno_cli.R` is the only script that reads from a different base
directory: the TWB2 release bundle, default `/mnt/SP-siliconpower/TWB20250806download/`.

---

## Workflow Overview

```
phenotype prep  →  external PLINK GWAS  →  process_gwas_sumstats_cli.R
                                                  │
                ┌─────────────────────────────────┼─────────────────────────────────┐
                ▼                                 ▼                                 ▼
         cojo.sh                          fusion_v8_twas.sh                spredixcan_v8_twas.sh
                ▼                                 ▼                                 ▼
   cojo_annotation_manhattan_cli.R   fusion_manhattan_heatmap_cli.R   spredixcan_manhattan_heatmap_cli.R
                └─────────────────────────────────┼─────────────────────────────────┘
                                                  ▼
                                       overlap_venn_cli.R
                                                  ▼
                                     pathway_enrichment_cli.R
```

Genome build alignment — the GWAS sumstats coordinate system must match the
model's coordinate system:

| Analysis | Coordinate system |
| --- | --- |
| COJO | matches your GWAS |
| FUSION v8 / S-PrediXcan v8 | GRCh38 / hg38 |
| FUSION v7 / S-PrediXcan v7 | GRCh37 / hg19 |

---

## CLI Script Reference

Default `base_dir` for every script is `/mnt/data/ai_agent/gene_analysis`.

Replace `BASE_DIR` with your actual project root directory (e.g., `~/gene_analysis_workflow`).

### 1. `make_twb2_pheno_cli.R` — phenotype preparation

Create TWB2 phenotype + keep files for downstream PLINK GWAS.

```bash
Rscript make_twb2_pheno_cli.R <var_name> <var_type> [case_label] [control_label]
```
or
```bash
Rscript make_twb2_pheno_cli.R \
  --var_name <name> \
  --var_type <categorical|continuous> \
  --case_label <v>] \
  --control_label <v>]
```

| Argument | Notes |
| --- | --- |
| `var_name` | Survey variable name, e.g. `VERTIGO_SELF`, `BMI` |
| `var_type` | `categorical` or `continuous` |
| `case_label` / `control_label` | Required only for `categorical` |

Rscript will read Inputs: 

-`survey/release_list_*.{txt,csv}`
-`lab_info/lab_info.csv`
-`Imputed.120161.TWB2/imputed_120161/TWB2.hg38.impu.v4.fam`

Outputs files: 

-`TWB2_<var_name>_pheno.txt`
-`TWB2_<var_name>_keep.txt`

For step2 External PLINK GWAS (written next to the input `.fam`). For `categorical`, coding is `1 = control`,`2 = case`.

### 2. External PLINK GWAS (out of scope)

Run GWAS outside this repository, typically with PLINK logistic / linear
regression, and place outputs as:

```
PLINK/imputed/glm_logistic/<project>/chrN.PHENO1.glm.logistic.hybrid
PLINK/imputed/plink_binary_files/<project>/chrN<PATTERN>.{bed,bim,fam}
```

### 3. `process_gwas_sumstats_cli.R` — convert GWAS outputs

```bash
Rscript process_gwas_sumstats_cli.R <project> [base_dir]
# or named: --project_name --base_dir
```

Inputs: `PLINK/imputed/glm_logistic/<project>/chrN.PHENO1.glm.logistic.hybrid`

Outputs:
- `COJO/<project>/chrN.cojo.ma` (COJO input, per chromosome)
- `FUSION/GWAS/<project>/chrN.sumstats`
- `S_PrediXcan/GWAS/<project>/chr1tochr22.sumstats` (single concatenated file)

### 4. `cojo.sh` — GCTA-COJO

```bash
bash cojo.sh PROJECT_NAME PLINK_FILE_PATTERN P_CUTOFF WINDOW [BASE_DIR]
# example: bash cojo.sh TWB1_LAA _DR2_0.7_QCFiltered 1e-6 10000
```

| Argument | Notes |
| --- | --- |
| `PLINK_FILE_PATTERN` | Suffix after `chrN`, e.g. `_DR2_0.7_QCFiltered` |
| `P_CUTOFF` | COJO p-value cutoff, e.g. `1e-6` |
| `WINDOW` | COJO window size (in bp), e.g. `10000` |

Inputs: PLINK binary files and `COJO/<project>/chrN.cojo.ma`.
Outputs: `COJO/<project>/result/chrN_p<TAG>.*` (`.jma.cojo`, `.cma.cojo`, etc.).
Terminal summary lists completed / missing chromosomes.

### 5. `cojo_annotation_manhattan_cli.R` — annotate COJO + Manhattan plot

```bash
Rscript cojo_annotation_manhattan_cli.R <project> [base_dir] [grch37|grch38]
# or named: --project_name --base_dir --genome_build --gtf_file
```

| Argument | Default |
| --- | --- |
| `genome_build` | `grch37` |
| `gtf_file` | Auto-selected from `BASE_DIR/tools/ensembl/` by genome_build |

Outputs in `COJO/<project>/`:
- `table/p1e5_annotation_genes.csv`
- `table/p1e5_protein_coding_nearest_gene.csv`
- `manhattan/cojo_protein_coding_manhanttan.png`

### 6. `fusion_v8_twas.sh` — FUSION GTEx v8 TWAS

```bash
bash fusion_v8_twas.sh SUBJECT_NAME [BASE_DIR]
```

Requires `BASE_DIR/FUSION/fusion_twas/FUSION.assoc_test.R`, `WEIGHTS/`, and
`LDREF/LDREF/`. R packages: `plink2R`, `optparse`.
Input: `FUSION/GWAS/<project>/chrN.sumstats`. Output: `FUSION/result_v8/<project>/chrN/`.
Sister script `fusion_v7_twas.sh` runs the v7 weights (`WEIGHTS_v7/GTEx.ALL/`).

### 7. `fusion_manhattan_heatmap_cli.R` — FUSION tables + plots

```bash
Rscript fusion_manhattan_heatmap_cli.R <project> [base_dir] [v7|v8] [fdr_cutoff]

Rscript fusion_manhattan_heatmap_cli.R  \
  --project_name <project> \
  --base_dir [base_dir] \
  --fusion_ver [v7|v8] \
  --fdr_cutoff [fdr_cutoff]

```

| Argument | Default |
| --- | --- |
| `fusion_ver` | `v8` |
| `fdr_cutoff` | `0.15` |
| `gencode_v26_gtf` | `<base_dir>/tools/GENCODE/gencode.v26.annotation.gtf` |
| `gencode_v19_gtf` | `<base_dir>/tools/GENCODE/gencode.v19.annotation.gtf` |
| `weight_dir_v7` | `<base_dir>/FUSION/WEIGHTS_v7/GTEx.ALL` |

Outputs in `FUSION/result_v{7,8}/<project>/{table,manhanttan,heatmap}/`.

### 8. `spredixcan_v8_twas.sh` — S-PrediXcan GTEx v8

```bash
bash spredixcan_v8_twas.sh <SUBJECT_NAME> [BASE_DIR]
```

Uses `conda run -n metaxcan`. Reads
`S_PrediXcan/GWAS/<project>/chr1tochr22.sumstats` and the v8 elastic-net models;
writes per-tissue CSVs to `S_PrediXcan/result_v8/<project>/`. Sister script
`spredixcan_v7_twas.sh` runs the v7 models.

### 9. `spredixcan_manhattan_heatmap_cli.R` — S-PrediXcan tables + plots

```bash
Rscript spredixcan_manhattan_heatmap_cli.R <project> [base_dir] [v7|v8] [fdr_cutoff]
# named: --project_name --base_dir --spredixcan_ver --fdr_cutoff
#        --gencode_v26_gtf --gencode_v19_gtf
```

GENCODE defaults are the same as for `fusion_manhattan_heatmap_cli.R`.
Outputs in `S_PrediXcan/result_v{7,8}/<project>/{table,manhanttan,heatmap}/`.

### 10. `overlap_venn_cli.R` — cross-method gene overlap + Venn

```bash
Rscript overlap_venn_cli.R \
    --project_name_v7 <name|None> \
    --project_name_v8 <name|None> \
    --project_name_cojo <name|None> \
    --project_name_overlap <name> \
    --base_dir <dir> \
    --ver [v7|v8|all] \
    --twas_method [fusion|spredixcan|all] \
    --fdr_cutoff \
    --cojo_p_cutoff 
```
| Argument | Default |
| --- | --- |
| `fdr_cutoff` | `0.15` |
| `cojo_p_cutoff` | `1e-6` |

Reads `FUSION/result_v{7,8}/<project>/table/FUSION_v{7,8}_TWAS.csv`,
`S_PrediXcan/result_v{7,8}/<project>/table/SPrediXcan_v{7,8}_TWAS.csv`, and
`COJO/<project>/table/p1e{5,05}_annotation_genes.csv`.

Outputs in `overlap/<project_overlap>/`:
`overlap_summary.csv`, `overlap_gene_name.csv`,
`overlap_gene_name_pairwise_union.csv`, `twas_detected_gene_name.csv`, and
`venn_<methods>_FDRxx.png`. Needs at least two valid methods.

`--cojo_p_cutoff` must match the COJO table's filename tag. The COJO annotation
table is currently named `p1e5_annotation_genes.csv` regardless of the COJO
p-value cutoff actually used, so pass `--cojo_p_cutoff 1e-5` to pick it up.


### 11. `pathway_enrichment_cli.R` — GO:BP + KEGG

```bash
Rscript pathway_enrichment_cli.R <project_overlap> [base_dir] [fdr_cutoff]
```
or
```bash
Rscript pathway_enrichment_cli.R \
  --project_name_overlap <project_overlap> \
  --base_dir [base_dir] \
  --fdr_cutoff [fdr_cutoff] 
```

Rscript will read:

- `<base_dir>/overlap/<project_overlap>/overlap_gene_name_pairwise_union.csv` as candidate genes
- `<base_dir>/overlap/<project_overlap>/twas_detected_gene_name.csv` as universe genes

Outputs files:

- `<base_dir>/pathway/GOBP/<project_overlap>/*_GOBP_ORA_FDRxx.{png,csv}`
- `<base_dir>/pathway/KEGG/<project_overlap>/*_KEGG_ORA_FDRxx.{png,csv}`
- `<base_dir>/pathway/log/<project_overlap>/process_log_*.txt`

GO uses `enrichGO(ont="BP")`; KEGG uses `enrichKEGG(organism="hsa")` (needsnetwork access). 

Empty enrichment results are handled safely (no plot written, empty summary CSV).

### Helper scripts (not part of the analysis chain)

| Script | Purpose |
| --- | --- |
| `setup_genomics_package.sh` | Install base tools + R + Miniconda (Step 2). |
| `install_r_packages.sh` | Install all required R packages (Step 3). |
| `setup_gene_analysis_dirs.sh` | Create the BASE_DIR folder skeleton (Step 4). |
| `impute_to_qc_plink.sh` | Standalone helper for imputation-to-QC steps. |
| `plink_to_conform.sh` | Standalone helper for PLINK file conformance. |
| `external_plink_gwas_step_example.txt` | Example command lines for the external GWAS step. |

---

## Pipeline Recipes

Replace `<P>` with your project name. Pass the project name that matches the
input genome build (v8 ⇒ hg38, v7 ⇒ hg19).

### Full pipeline (single project, all three branches)

```bash
BASE_DIR=~/gene_analysis_workflow
cd ~/gene-analysis-cli   # or your scripts folder

# 1. Convert PLINK GWAS to downstream formats
Rscript process_gwas_sumstats_cli.R <P> "$BASE_DIR"

# 2. COJO branch
bash cojo.sh <P> _QC 1e-6 10000 "$BASE_DIR"
Rscript cojo_annotation_manhattan_cli.R <P> "$BASE_DIR" grch37

# 3. FUSION v8 branch
bash fusion_v8_twas.sh <P_hg38> "$BASE_DIR"
Rscript fusion_manhattan_heatmap_cli.R <P_hg38> "$BASE_DIR" v8 0.15

# 4. S-PrediXcan v8 branch
bash spredixcan_v8_twas.sh <P_hg38> "$BASE_DIR"
Rscript spredixcan_manhattan_heatmap_cli.R <P_hg38> "$BASE_DIR" v8 0.15

# 5. Overlap (need at least two of: COJO, FUSION TWAS, S-PrediXcan)
Rscript overlap_venn_cli.R \
    --project_name_cojo <P> \
    --project_name_v8 <P_hg38> \
    --project_name_overlap <P> \
    --base_dir "$BASE_DIR" \
    --ver v8 --twas_method all \
    --fdr_cutoff 0.2 --cojo_p_cutoff 1e-5

# 6. Pathway enrichment from overlap
Rscript pathway_enrichment_cli.R <P> "$BASE_DIR" 0.2
```

### Single branch only

| Branch | After `process_gwas_sumstats_cli.R`, run |
| --- | --- |
| COJO only | `cojo.sh` → `cojo_annotation_manhattan_cli.R` |
| FUSION only | `fusion_v8_twas.sh` → `fusion_manhattan_heatmap_cli.R` |
| S-PrediXcan only | `spredixcan_v8_twas.sh` → `spredixcan_manhattan_heatmap_cli.R` |

## Citation

If you use this workflow, cite the underlying methods and resources, not only
this repository.

### Core software and methods

- **PLINK** — Purcell S, *et al.* PLINK: a tool set for whole-genome association and
  population-based linkage analyses. *Am J Hum Genet* 2007.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/17701901/)
- **PLINK 2** — Chang CC, *et al.* Second-generation PLINK: rising to the challenge
  of larger and richer datasets. *GigaScience* 2015.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/25722852/)
- **GCTA / COJO** — Yang J, *et al.* Conditional and joint multiple-SNP analysis of
  GWAS summary statistics. *Nat Genet* 2012.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/22426310/)
- **FUSION TWAS** — Gusev A, *et al.* Integrative approaches for large-scale
  transcriptome-wide association studies. *Nat Genet* 2016.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/26854917/)
- **PrediXcan** — Gamazon ER, *et al.* A gene-based association method for mapping
  traits using reference transcriptome data. *Nat Genet* 2015.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/26258848/)
- **MetaXcan / S-PrediXcan** — Barbeira AN, *et al.* Integrating predicted
  transcriptome from multiple tissues improves association detection.
  *PLoS Genet* 2019. [PubMed](https://pubmed.ncbi.nlm.nih.gov/30694772/)

### Reference data

- **GTEx Consortium** — *Science* 2020.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/32913098/)
- **Ensembl** — Yates AD, *et al.* Ensembl 2020. *NAR* 2020.
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/31691826/)

### Project pages

- GCTA: <https://yanglab.westlake.edu.cn/software/gcta/>
- FUSION: <http://gusevlab.org/projects/fusion/>
- MetaXcan / S-PrediXcan: <https://github.com/hakyimlab/MetaXcan>
- PredictDB GTEx models: <https://predictdb.org/>
- Ensembl FTP: <https://www.ensembl.org/info/data/ftp/index.html>
