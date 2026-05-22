#!/usr/bin/env bash
# =====================================================================
# install_r_packages.sh
# Install all R packages required by the gene-analysis-cli pipeline.
#
# Covers:
#   - System libraries needed to compile R packages (apt)
#   - CRAN packages
#   - plink2R (from GitHub, used by FUSION)
#   - Bioconductor packages
#
# Usage:
#   chmod +x install_r_packages.sh
#   ./install_r_packages.sh
#
# Notes:
#   - Idempotent: already-installed packages are skipped.
#   - Installs into the user R library (no sudo needed for R itself).
#   - Bioconductor compilation can take 20-40 min on first run.
# =====================================================================

set -euo pipefail

GREEN='\033[0;32m'; BLUE='\033[0;34m'; YELLOW='\033[1;33m'; RED='\033[0;31m'; NC='\033[0m'
info()  { echo -e "${BLUE}[INFO]${NC} $*"; }
ok()    { echo -e "${GREEN}[ OK ]${NC} $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
err()   { echo -e "${RED}[FAIL]${NC} $*"; }

# ---------------------------------------------------------------------
# Step 1: System libraries (needed to compile rtracklayer, ComplexHeatmap,
#         clusterProfiler, ggplot2 graphics backends, git-based installs)
# ---------------------------------------------------------------------
info "Installing system build dependencies (sudo required)..."
sudo apt update
sudo apt install -y \
    libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libcurl4-openssl-dev libssl-dev libgit2-dev gfortran \
    libglpk-dev libgmp-dev \
    cmake make g++
ok "System dependencies installed"
# Note: libglpk-dev and libgmp-dev are REQUIRED by igraph (a clusterProfiler
# dependency). Without them clusterProfiler fails to compile.

# ---------------------------------------------------------------------
# Step 2: Run the R installer
# ---------------------------------------------------------------------
info "Starting R package installation (this may take a while)..."

Rscript - <<'RSCRIPT'
# ---- Setup ----------------------------------------------------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
ncpu <- max(1L, parallel::detectCores())
options(Ncpus = ncpu)
cat(sprintf(">>> Using %d cores for compilation\n", ncpu))

installed_ok <- character(0)
failed       <- character(0)

need <- function(p) !requireNamespace(p, quietly = TRUE)

# ---- CRAN packages --------------------------------------------------
cran_pkgs <- c(
  "remotes", "data.table", "dplyr", "ggplot2", "ggrepel",
  "stringr", "tidyr", "tibble", "circlize", "ggvenn", "optparse"
)

cat("\n========== CRAN packages ==========\n")
for (p in cran_pkgs) {
  if (need(p)) {
    cat(sprintf(">>> Installing CRAN: %s\n", p))
    tryCatch({
      install.packages(p)
      if (!need(p)) installed_ok <- c(installed_ok, p) else failed <- c(failed, p)
    }, error = function(e) { failed <<- c(failed, p); cat(sprintf("    ERROR: %s\n", conditionMessage(e))) })
  } else {
    cat(sprintf("[skip] %s already installed\n", p))
  }
}

# ---- plink2R (GitHub, used by FUSION) -------------------------------
cat("\n========== plink2R (GitHub) ==========\n")
if (need("plink2R")) {
  cat(">>> Installing plink2R from gabraham/plink2R\n")
  tryCatch({
    remotes::install_github("gabraham/plink2R/plink2R", upgrade = "never")
    if (!need("plink2R")) installed_ok <- c(installed_ok, "plink2R") else failed <- c(failed, "plink2R")
  }, error = function(e) { failed <<- c(failed, "plink2R"); cat(sprintf("    ERROR: %s\n", conditionMessage(e))) })
} else {
  cat("[skip] plink2R already installed\n")
}

# ---- Bioconductor packages ------------------------------------------
cat("\n========== Bioconductor packages ==========\n")
if (need("BiocManager")) install.packages("BiocManager")

bioc_pkgs <- c(
  "rtracklayer", "ComplexHeatmap", "clusterProfiler",
  "GO.db", "org.Hs.eg.db", "AnnotationDbi", "KEGGREST"
)
for (p in bioc_pkgs) {
  if (need(p)) {
    cat(sprintf(">>> Installing Bioconductor: %s\n", p))
    tryCatch({
      BiocManager::install(p, update = FALSE, ask = FALSE)
      if (!need(p)) installed_ok <- c(installed_ok, p) else failed <- c(failed, p)
    }, error = function(e) { failed <<- c(failed, p); cat(sprintf("    ERROR: %s\n", conditionMessage(e))) })
  } else {
    cat(sprintf("[skip] %s already installed\n", p))
  }
}

# ---- Summary --------------------------------------------------------
cat("\n========== Verification ==========\n")
all_pkgs <- c(cran_pkgs, "plink2R", bioc_pkgs)
for (p in all_pkgs) {
  status <- if (!need(p)) "OK   " else "MISSING"
  cat(sprintf("  [%s] %s\n", status, p))
}

if (length(failed) > 0) {
  cat(sprintf("\n>>> %d package(s) FAILED: %s\n", length(failed), paste(unique(failed), collapse = ", ")))
  quit(status = 1)
} else {
  cat("\n>>> All R packages installed successfully!\n")
}
RSCRIPT

rc=$?
echo ""
if [[ $rc -eq 0 ]]; then
  ok "==========================================="
  ok "  R package installation complete!"
  ok "==========================================="
else
  err "Some packages failed to install. Check the [MISSING] list above."
  err "Re-run this script (it skips already-installed packages) or share the error log."
fi
exit $rc
