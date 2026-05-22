#!/usr/bin/env bash
# =====================================================================
# setup-genomics.sh
# Bioinformatics tools setup for Ubuntu 22.04 (WSL2 or native)
#
# Installs:
#   - System updates & build dependencies
#   - PLINK 1.9 + PLINK 2.0
#   - HTSlib 1.21 (includes tabix, bgzip)
#   - Samtools 1.21
#   - Bcftools 1.21
#   - VCFtools (apt)
#   - R 4.5.1 (via rig, R Installation Manager)
#
# Usage:
#   chmod +x setup-genomics.sh
#   ./setup-genomics.sh
#
# Idempotent: safe to re-run; already-installed tools are skipped.
# =====================================================================

set -euo pipefail

# ---------- Tool versions (edit if you want different versions) ----------
HTSLIB_VERSION="1.21"
SAMTOOLS_VERSION="1.21"
BCFTOOLS_VERSION="1.21"
PLINK1_BUILD="20231211"
PLINK2_BUILD="20241222"
R_VERSION="4.5.1"

# ---------- Colors ----------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info()    { echo -e "${BLUE}[INFO]${NC} $*"; }
log_success() { echo -e "${GREEN}[ OK ]${NC} $*"; }
log_warn()    { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error()   { echo -e "${RED}[FAIL]${NC} $*"; }

section() {
    echo ""
    echo -e "${BLUE}=====================================================${NC}"
    echo -e "${BLUE}  $*${NC}"
    echo -e "${BLUE}=====================================================${NC}"
}

check_cmd() { command -v "$1" &> /dev/null; }

# ---------- Pre-flight checks ----------
if [[ $EUID -eq 0 ]]; then
    log_error "Do NOT run this script as root. Run as your normal user; sudo is used internally."
    exit 1
fi

log_info "Checking sudo access (you may be prompted for your password)..."
sudo -v
log_success "sudo OK"

START_TIME=$(date +%s)

# ============================================================
# Step 1: System update & build dependencies
# ============================================================
section "Step 1/7  System update & build dependencies"

sudo apt update
sudo apt upgrade -y

sudo apt install -y \
    build-essential autoconf automake pkg-config \
    libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    zlib1g-dev libncurses5-dev libdeflate-dev \
    wget curl git unzip ca-certificates \
    software-properties-common dirmngr gdebi-core \
    vcftools

log_success "Build dependencies installed"

# ============================================================
# Step 2: PLINK 1.9
# ============================================================
section "Step 2/7  PLINK 1.9"

if check_cmd plink; then
    log_warn "plink already installed: $(plink --version 2>&1 | head -1)"
else
    cd /tmp
    wget -q "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK1_BUILD}.zip"
    unzip -o "plink_linux_x86_64_${PLINK1_BUILD}.zip" -d plink1_pkg > /dev/null
    sudo mv plink1_pkg/plink /usr/local/bin/
    sudo chmod +x /usr/local/bin/plink
    log_success "PLINK 1.9 installed: $(plink --version 2>&1 | head -1)"
fi

# ============================================================
# Step 3: PLINK 2.0
# ============================================================
section "Step 3/7  PLINK 2.0"

if check_cmd plink2; then
    log_warn "plink2 already installed: $(plink2 --version 2>&1 | head -1)"
else
    cd /tmp
    wget -q "https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_${PLINK2_BUILD}.zip"
    unzip -o "plink2_linux_x86_64_${PLINK2_BUILD}.zip" -d plink2_pkg > /dev/null
    sudo mv plink2_pkg/plink2 /usr/local/bin/
    sudo chmod +x /usr/local/bin/plink2
    log_success "PLINK 2.0 installed: $(plink2 --version 2>&1 | head -1)"
fi

# ============================================================
# Step 4: HTSlib
# ============================================================
section "Step 4/7  HTSlib ${HTSLIB_VERSION}"

if check_cmd tabix && tabix --version 2>&1 | grep -q "${HTSLIB_VERSION}"; then
    log_warn "HTSlib ${HTSLIB_VERSION} already installed"
else
    cd /tmp
    wget -q "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"
    tar -xjf "htslib-${HTSLIB_VERSION}.tar.bz2"
    cd "htslib-${HTSLIB_VERSION}"
    ./configure --quiet
    make -j"$(nproc)" -s
    sudo make install -s
    log_success "HTSlib ${HTSLIB_VERSION} installed (provides tabix, bgzip)"
fi

# ============================================================
# Step 5: Samtools
# ============================================================
section "Step 5/7  Samtools ${SAMTOOLS_VERSION}"

if check_cmd samtools && samtools --version 2>&1 | head -1 | grep -q "${SAMTOOLS_VERSION}"; then
    log_warn "Samtools ${SAMTOOLS_VERSION} already installed"
else
    cd /tmp
    wget -q "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
    tar -xjf "samtools-${SAMTOOLS_VERSION}.tar.bz2"
    cd "samtools-${SAMTOOLS_VERSION}"
    ./configure --quiet
    make -j"$(nproc)" -s
    sudo make install -s
    log_success "Samtools ${SAMTOOLS_VERSION} installed"
fi

# ============================================================
# Step 6: Bcftools
# ============================================================
section "Step 6/7  Bcftools ${BCFTOOLS_VERSION}"

if check_cmd bcftools && bcftools --version 2>&1 | head -1 | grep -q "${BCFTOOLS_VERSION}"; then
    log_warn "Bcftools ${BCFTOOLS_VERSION} already installed"
else
    cd /tmp
    wget -q "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"
    tar -xjf "bcftools-${BCFTOOLS_VERSION}.tar.bz2"
    cd "bcftools-${BCFTOOLS_VERSION}"
    ./configure --quiet
    make -j"$(nproc)" -s
    sudo make install -s
    log_success "Bcftools ${BCFTOOLS_VERSION} installed"
fi

# Refresh shared library cache
sudo ldconfig

# ============================================================
# Step 7: R via rig
# ============================================================
section "Step 7/7  R ${R_VERSION} (via rig)"

if ! check_cmd rig; then
    log_info "Installing rig (R Installation Manager)..."
    curl -Ls https://github.com/r-lib/rig/releases/download/latest/rig-linux-latest.tar.gz \
        | sudo tar xzf - -C /usr/local
    log_success "rig installed: $(rig --version)"
else
    log_warn "rig already installed: $(rig --version)"
fi

if R --version 2>/dev/null | head -1 | grep -q "${R_VERSION}"; then
    log_warn "R ${R_VERSION} already installed"
else
    log_info "Installing R ${R_VERSION}... (this may take a few minutes)"
    sudo rig add "${R_VERSION}"
    sudo rig default "${R_VERSION}"
    log_success "R ${R_VERSION} installed and set as default"
fi

log_info "Configuring user R library..."
rig system setup-user-lib

# ============================================================
# Final verification
# ============================================================
section "Verification"

echo ""
printf "  %-10s : " "plink"    ; plink   --version 2>&1 | head -1
printf "  %-10s : " "plink2"   ; plink2  --version 2>&1 | head -1
printf "  %-10s : " "samtools" ; samtools --version 2>&1 | head -1
printf "  %-10s : " "bcftools" ; bcftools --version 2>&1 | head -1
printf "  %-10s : " "vcftools" ; vcftools --version 2>&1
printf "  %-10s : " "tabix"    ; tabix   --version 2>&1 | head -1
printf "  %-10s : " "bgzip"    ; bgzip   --version 2>&1 | head -1
printf "  %-10s : " "R"        ; R       --version 2>&1 | head -1
printf "  %-10s : " "Rscript"  ; Rscript --version 2>&1

# ---------- Cleanup ----------
log_info "Cleaning up /tmp build files..."
cd ~
rm -rf /tmp/plink1_pkg /tmp/plink2_pkg /tmp/htslib-* /tmp/samtools-* /tmp/bcftools-*
rm -f  /tmp/plink_linux_*.zip /tmp/plink2_linux_*.zip /tmp/htslib-*.tar.bz2 /tmp/samtools-*.tar.bz2 /tmp/bcftools-*.tar.bz2

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
log_success "==========================================="
log_success "  All tools installed successfully!"
log_success "  Elapsed: $((ELAPSED / 60))m $((ELAPSED % 60))s"
log_success "==========================================="
echo ""
echo "Next steps (optional):"
echo ""
echo "  # Install common R packages for genomics:"
echo "  Rscript -e 'install.packages(c(\"data.table\",\"qqman\",\"ggplot2\",\"dplyr\"), repos=\"https://cloud.r-project.org\")'"
echo ""
echo "  # Create standard working directories:"
echo "  mkdir -p ~/projects ~/data ~/scripts"
echo ""
