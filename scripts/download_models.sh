#!/bin/bash

# RFdiffusion, AlphaFold2, and Boltz-2 Model Downloader
# This script downloads all required models for protein prediction pipelines
# Total download size: ~15.0 GB

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check dependencies
print_status "Checking dependencies..."
if ! command_exists wget; then
    print_error "wget is required but not installed. Please install wget first."
    exit 1
fi

if ! command_exists tar; then
    print_error "tar is required but not installed. Please install tar first."
    exit 1
fi

# Function to download with retry
download_with_retry() {
    local url="$1"
    local filename="$2"
    local max_attempts=3
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        print_status "Downloading $filename (attempt $attempt/$max_attempts)..."
        if wget --progress=bar:force:noscroll --timeout=30 --tries=3 "$url" -O "$filename"; then
            print_success "Downloaded $filename"
            return 0
        else
            print_warning "Failed to download $filename (attempt $attempt/$max_attempts)"
            if [ $attempt -eq $max_attempts ]; then
                print_error "Failed to download $filename after $max_attempts attempts"
                return 1
            fi
            attempt=$((attempt + 1))
            sleep 5
        fi
    done
}

# Function to verify file size (optional check)
verify_file() {
    local file="$1"
    if [ -f "$file" ]; then
        local size=$(du -h "$file" | cut -f1)
        print_success "Verified: $file ($size)"
        return 0
    else
        print_error "Missing: $file"
        return 1
    fi
}

print_status "Starting model downloads..."
print_status "Total expected download size: ~15.0 GB. Final size: ~6.8 GB"
print_status "This may take a while depending on your internet connection."

# ============================================================================
# RFdiffusion Models (~3.7 GB)
# ============================================================================

print_status "Downloading RFdiffusion models (~3.7 GB)..."
mkdir -p models/rfd && cd models/rfd

# Array of RFdiffusion model URLs and filenames
declare -a rfd_models=(
    "http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt"
    "http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt"
    "http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt"
    "http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt"
    "http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt"
    "http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt"
    "http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt"
    "http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt"
)

for url in "${rfd_models[@]}"; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        print_success "$filename already exists, skipping..."
    else
        download_with_retry "$url" "$filename" || exit 1
    fi
done

cd ../..
print_success "RFdiffusion models download completed!"

# ============================================================================
# AlphaFold2 Models (~5.2 GB)
# ============================================================================

print_status "Downloading AlphaFold2 models (~5.2 GB)..."
mkdir -p models/af2 && cd models/af2

AF2_URL="https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar"
AF2_FILE="alphafold_params_2022-12-06.tar"

if [ -f "$AF2_FILE" ]; then
    print_success "$AF2_FILE already exists, skipping download..."
else
    download_with_retry "$AF2_URL" "$AF2_FILE" || exit 1
fi

# Extract if not already extracted
if [ -d "params" ]; then
    print_success "AlphaFold2 params already extracted, skipping..."
else
    print_status "Extracting AlphaFold2 parameters..."
    tar --extract --verbose --file="$AF2_FILE" params_model_1.npz params_model_1_multimer_v3.npz params_model_1_ptm.npz LICENSE|| exit 1
    rm -f alphafold_params_2022-12-06.tar

    print_success "AlphaFold2 parameters extracted!"
fi

cd ../..
print_success "AlphaFold2 models download completed!"

# ============================================================================
# Boltz-2 Models (~6.0 GB)
# ============================================================================

print_status "Downloading Boltz-2 models (~6.0 GB)..."
mkdir -p models/boltz && cd models/boltz

# Array of Boltz-2 model URLs and filenames
declare -a boltz_models=(
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt"
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar"
)

for url in "${boltz_models[@]}"; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        print_success "$filename already exists, skipping..."
    else
        download_with_retry "$url" "$filename" || exit 1
    fi
done

# Extract mols.tar if not already extracted
if [ -d "mols" ]; then
    print_success "Boltz-2 mols directory already extracted, skipping..."
else
    print_status "Extracting Boltz-2 molecular data... (only amino acids are needed for protein design)"
    tar -tf mols.tar | grep -E 'mols/(ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|UNK)\.pkl$' | xargs -d '\n' -I {} tar -xvf mols.tar "{}" || exit 1
    rm -f mols.tar
    print_success "Boltz-2 molecular data extracted!"
fi

print_status "Creating dummy affinity checkpoint file (not used for protein design)"
touch boltz2_aff.ckpt

cd ../..
print_success "Boltz-2 models download completed!"

# ============================================================================
# Verification
# ============================================================================

print_status "Verifying downloads..."

# Check RFdiffusion models
print_status "Checking RFdiffusion models..."
rfd_count=$(ls -1 models/rfd/*.pt 2>/dev/null | wc -l)
if [ "$rfd_count" -eq 8 ]; then
    print_success "RFdiffusion: Found all 8 required .pt files"
else
    print_error "RFdiffusion: Expected 8 .pt files, found $rfd_count"
fi

# Check AlphaFold2 models
print_status "Checking AlphaFold2 models..."
af2_count=$(ls -1 models/af2/params_model_1*.npz 2>/dev/null | wc -l)
if [ "$af2_count" -eq 3 ]; then
    print_success "AlphaFold2: Found all 3 params files"
else
    print_error "AlphaFold2: Expected 3 .npz files, found $af2_count"
fi

# Check Boltz-2 models
print_status "Checking Boltz-2 models..."
boltz_ckpt_count=$(ls -1 models/boltz/*.ckpt 2>/dev/null | wc -l)
if [ "$boltz_ckpt_count" -ge 2 ] && [ -d "models/boltz/mols" ]; then
    print_success "Boltz-2: Found .ckpt files and mols directory"
else
    print_error "Boltz-2: Missing required files (.ckpt files or mols directory)"
fi

# Final summary
print_status "Download Summary:"
echo "  - RFdiffusion models: $(du -sh models/rfd 2>/dev/null | cut -f1 || echo 'N/A')"
echo "  - AlphaFold2 models:  $(du -sh models/af2 2>/dev/null | cut -f1 || echo 'N/A')"
echo "  - Boltz-2 models:     $(du -sh models/boltz 2>/dev/null | cut -f1 || echo 'N/A')"
echo "  - Total size:         $(du -sh models/ 2>/dev/null | tail -1 | cut -f1 || echo 'N/A')"

print_success "Model download script completed!"
print_status "Remember to update your nextflow.config file with the model paths:"
echo "  - models/rfd = './models/rfd'"
echo "  - models/af2 = './models/af2'"
echo "  - models/boltz = './models/boltz'"
