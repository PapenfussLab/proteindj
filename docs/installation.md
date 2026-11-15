[🏠 ProteinDJ](../README.md) > Installation Guide

# ProteinDJ Installation Guide

This guide will walk you through installing ProteinDJ and all its dependencies and is intended for system admins. The installation process involves several steps that need to be completed once per cluster or system.

> **Upgrading from v2?** v3 is not backwards compatible — containers and models have been updated (e.g. FreeBindCraft replaces BindCraft, and new `mpnn_models`/`bg_models` paths are required), so a fresh installation following this guide is required rather than reusing an existing v2 setup.

## Prerequisites

Before starting, ensure you have:

- **Linux/Unix system** with internet access
- **Sufficient storage space**: ~18 GB for downloaded models and ~50 GB for containers
- **Administrative access** or ability to install software
- **SLURM cluster** (if using HPC environment)

## System Requirements

| Component   | Minimum                    | Recommended     |
| ----------- | -------------------------- | --------------- |
| **Storage** | 70 GB                      | 100 GB+         |
| **RAM**     | 32 GB                      | 48 GB+          |
| **GPU**     | NVIDIA GPU with 16GB+ VRAM | NVIDIA A30/A100 |
| **CPU**     | 8 cores                    | 24+ cores       |

---

## Step 1: Clone Repository

First, clone the ProteinDJ repository:

```bash
git clone https://github.com/PapenfussLab/proteindj
cd proteindj
```

---

## Step 2: Install Dependencies

ProteinDJ requires two key dependencies to be installed and accessible in your PATH. These are common software packages so they may already be implemented in your HPC environment (e.g. module load apptainer nextflow):

### Apptainer (Required)

Install [Apptainer](https://apptainer.org/docs/admin/main/installation.html) for containerization.

**Ubuntu/Debian:**

```bash
# Add repository and install
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

**CentOS/RHEL:**

```bash
# Install from EPEL
sudo yum install -y epel-release
sudo yum install -y apptainer
```

### Nextflow (Required)

Install [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ v24.04):

```bash
# Install Nextflow
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

**Verify installations:**

```bash
apptainer --version
nextflow -version
```

---

## Step 3: Download Required Models

> **💡 Tip:** This step downloads ~18 GB of model files. Consider doing this in a shared location so that other users can access the files.

### RFD4-Proteina (optional, for `rfd4_denovo`)

The `RFD4` process label uses the `${params.container_registry}/rfd4:${params.container_version}` image, built from `apptainer/rfd4.def`. `scripts/download_models.sh` downloads the pinned `nvidia/RFD4-Proteina` model bundle into `models/rfd4`, the default `rfd4_models` directory. Set `rfd4_models` to another prepopulated directory when models are stored centrally. cuEquivariance tuning data and compiled Triton kernels are shared between tasks through `rfd4_triton_cache`, which defaults to `$HOME/.cache/proteindj/rfd4` and is created automatically. See the [upstream installation guide](https://rfd4-proteina.org/latest/installation.html) for driver and model requirements.

To populate a shared cache directly, install the Hugging Face CLI and download the same pinned snapshot:

```bash
RFD4_MODELS=/path/to/shared/rfd4/models
mkdir -p "$RFD4_MODELS"
hf download nvidia/RFD4-Proteina \
    config.json model.safetensors \
    --revision a06e2509c7033b010fc43614e9fd840303588502 \
    --local-dir "$RFD4_MODELS"
```

Set `rfd4_models` to the value of `RFD4_MODELS`. The pipeline bind-mounts that directory and passes it to RFD4 as a local model bundle, so inference does not contact Hugging Face.

### RFdiffusion Models (~3.7 GB)

RFdiffusion requires several diffusion model checkpoints (~3.7 GB). If you have not already downloaded the models, use the commands below, and update the `rfd_models` variable in `nextflow.config` to the location of the model directory (e.g. './models/rfd'):

```bash
mkdir -p models/rfd && cd models/rfd

# Download all required checkpoints
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt

cd ../..
```

### AlphaFold2 Models (~5.3 GB)

To perform AlphaFold2 predictions, you will need to download the AF2 models from DeepMind's repository. If you have not already downloaded the models, use the commands below, and update the `af2_models` variable in `nextflow.config` to the location of the model directory (e.g. './models/af2'):

```bash
mkdir -p models/af2 && cd models/af2

# Download and extract AF2 parameters (only need the first model for AlphaFold2 Initial-Guess)
wget https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar
tar -xf alphafold_params_2022-12-06.tar 
rm -f alphafold_params_2022-12-06.tar

cd ../..
```

### Boltz-2 Models (~2.5 GB)

To perform Boltz-2 predictions, you will need to download the models (download ~3.8 GB, final size ~2.5 GB). If you have not already downloaded the models, use the commands below, and update the `boltz_models` variable in `nextflow.config` to the location of the model directory (e.g. './boltz_models'):

```bash
mkdir -p models/boltz && cd models/boltz

# Download Boltz-2 checkpoints and molecular data
wget https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt
wget https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar
# We only need the amino acid CCD files for protein design
tar -tf mols.tar | grep -E 'mols/(ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|UNK)\.pkl$' | xargs -d '\n' -I {} tar -xvf mols.tar "{}"
rm -f mols.tar
# Create dummy affinity checkpoint file (required by Boltz2 but not used for protein design)
touch boltz2_aff.ckpt

cd ../..
```

### BoltzGen Models (~4 GB)

ProteinDJ only runs BoltzGen's design step (sequences are re-designed downstream via MPNN/FAMPNN), so only the diffusion checkpoints and the molecule/CCD dataset need to be downloaded — not the inverse-folding, folding, or affinity checkpoints that `boltzgen run` would otherwise also fetch. Update the `bg_models` variable in `nextflow.config` to the location of the model directory (e.g. `'./models/boltzgen'`):

```bash
mkdir -p models/boltzgen

# These files are served from Hugging Face Hub's Xet storage backend. Using wget/curl is
# very slow so we recommend using the `hf` CLI instead (fetches chunks in parallel)
pip install -U "huggingface_hub[hf_xet]"
hf download boltzgen/boltzgen-1 boltzgen1_diverse.ckpt boltzgen1_adherence.ckpt --local-dir models/boltzgen
hf download boltzgen/inference-data mols.zip --repo-type dataset --local-dir models/boltzgen
```

### ProteinMPNN Model Weights (~0.1 GB)

ProteinMPNN weights are required for binder sequence design. Three sets of weights are provided: vanilla (standard), soluble (optimised for solubility), and HyperMPNN (improved accuracy). Update the `mpnn_models` variable in `nextflow.config` to the location of the model directory (e.g. `'./models/mpnn'`):

```bash
mkdir -p models/mpnn && cd models/mpnn

# Vanilla ProteinMPNN weights
mkdir vanilla_model_weights && cd vanilla_model_weights
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/vanilla_model_weights/v_48_002.pt
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/vanilla_model_weights/v_48_010.pt
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/vanilla_model_weights/v_48_020.pt
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/vanilla_model_weights/v_48_030.pt
cd ..

# Soluble ProteinMPNN weights
mkdir soluble_model_weights && cd soluble_model_weights
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/soluble_model_weights/v_48_002.pt
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/soluble_model_weights/v_48_010.pt
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/soluble_model_weights/v_48_020.pt
wget https://github.com/dauparas/ProteinMPNN/raw/refs/heads/main/soluble_model_weights/v_48_030.pt
cd ..

# HyperMPNN weights
mkdir hyper_model_weights && cd hyper_model_weights
wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_002_epoch240_hyper.pt -O v_48_002.pt
wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_010_epoch300_hyper.pt -O v_48_010.pt
wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_020_epoch300_hyper.pt -O v_48_020.pt
wget https://github.com/meilerlab/HyperMPNN/raw/refs/heads/main/retrained_models/v48_030_epoch300_hyper.pt -O v_48_030.pt
cd ../../../
```

**Verify downloads:**

```bash
# Check that all required files exist
ls -la models/rfd/                         # Should contain 8 .pt files
ls -la models/af2/                         # Should contain params directory
ls -la models/boltz/                       # Should contain .ckpt files and mols directory
ls -la models/boltzgen/                    # Should contain boltzgen1_diverse.ckpt, boltzgen1_adherence.ckpt and mols.zip
ls -la models/mpnn/vanilla_model_weights/  # Should contain 4 .pt files
ls -la models/mpnn/soluble_model_weights/  # Should contain 4 .pt files
ls -la models/mpnn/hyper_model_weights/    # Should contain 4 .pt files
```

---

## Step 4: Build Containers

ProteinDJ requires several containers for the different dependencies. By default, these will be fetched during execution of Nextflow and cached by Apptainer to the location specified by the environment variable `NXF_APPTAINER_CACHEDIR`. If you are on a shared environment, we recommend creating a profile in nextflow.config with `apptainer.cacheDir` set to a shared location with read/write permissions for all users (see the Milton/WEHI profile as an example). 

We recommend using our containers from the cloud, but if you would like to build/modify containers locally, we have provide def files for containers in `proteindj/apptainer`. When updating ProteinDJ make sure to also update your local containers to maintain compatibility with the pipeline. You may already have similar containers for some of these programs, but we have made changes to the source code and environment so we do not recommend using other containers with ProteinDJ. Once the containers have been built, you need to change the path to each built container in `nextflow.config` e.g. for BindCraft: 
```
withLabel: 'BC' {
   container         = "/path/to/containers/bindcraft.sif"
   containerOptions  = """--nv \
   --bind ${params.af2_models}:/af2params \
   """
}     
```

We have provided a script for building the containers in a series of sbatch jobs (`apptainer/build_containers.sh`). You may need to tweak the SLURM parameters and enviroment settings for your cluster.

> **⚠️ Important:** Container building requires significant resources and may take several hours. We recommend using our pre-compiled containers unless you have reason to build your own.

### Option A: Automated Parallel Build using SLURM script

1. **Navigate to the apptainer directory:**

   ```bash
   cd apptainer
   ```

2. **Edit build script:**

   ```bash
   # Edit the BUILD_DIRECTORY in build_containers.sh
   nano build_containers.sh

   # Update this line to your desired location:
   BUILD_DIRECTORY="/path/to/your/containers"
   ```

3. **Submit build jobs:**

   ```bash
   ./build_containers.sh
   ```

4. **Monitor progress:**
   ```bash
   # Check SLURM queue
   squeue -u $USER
   ```

### Option B: Manual Build

If you don't have SLURM or prefer manual building:

```bash
cd apptainer

# Build each container individually
apptainer build --fakeroot bindsweeper.sif bindsweeper.def
apptainer build --fakeroot boltz2.sif boltz2.def
apptainer build --fakeroot dl_binder_design.sif dl_binder_design.def
apptainer build --fakeroot fampnn.sif fampnn.def
apptainer build --fakeroot python_tools.sif python_tools.def
apptainer build --fakeroot rfdiffusion.sif rfdiffusion.def
cd ..
```

**Verify container builds:**

```bash
ls -la containers/  # Should contain 7 .sif files
```

---

## Step 5: Configure Cluster Parameters

Edit the `nextflow.config` file to match your system configuration:

```bash
nano nextflow.config
```

**Key parameters to update:**

| Parameter       | Description                                 | Examples                       |
| --------------- | ------------------------------------------- | ------------------------------ |
| `rfd_models`    | Path to RFdiffusion models                  | `"${projectDir}/models/rfd"`   |
| `mpnn_models`   | Path to ProteinMPNN models                  | `"${projectDir}/models/mpnn"`  |
| `af2_models`    | Path to AlphaFold2 models                   | `"${projectDir}/models/af2"`   |
| `boltz_models`  | Path to Boltz-2 models                      | `"${projectDir}/models/boltz"` |
| `bg_models`     | Path to BoltzGen models                     | `"${projectDir}/models/boltzgen"` |
| `gpu_model`     | Your GPU type                               | `'A30'`, `'V100'`, `'A100'`    |
| `gpus`          | Number of GPUs to request                   | `1`, `2`, `4`, `8`             |
| `cpus_per_gpu`  | Number of CPUs to request per GPU           | `8`, `12`                      |
| `memory_gpu`    | Memory to request for GPU jobs              | `'24GB'`, `'48GB'`             |
| `cpus`          | Number of CPUs to request for CPU-only jobs | `12`, `24`                     |
| `memory_cpu`    | Memory for request for CPU-only jobs        | `'24GB'`, `'48GB'`             |

---

## Step 6: Test Installation

Before running production workloads, verify your installation works correctly.

### Test 1: Monomer Design (5-10 minutes)

```bash
# From the proteindj root directory
nextflow run main.nf -profile test,rfd_denovo_monomer
```

**Expected output:**

- Uses RFdiffusion, Full-atom MPNN, and Boltz-2
- Generates small number of de novo monomers
- Creates output directory with results

### Test 2: Binder Design (5-10 minutes)

```bash
nextflow run main.nf -profile test,rfd_denovo_binder
```

**Expected output:**

- Uses RFdiffusion, ProteinMPNN, and AlphaFold2 Initial-Guess
- Generates small number of binders
- Creates output directory with results

### Comprehensive Testing (Optional)

For thorough validation, we have made an end-to-end testing script that performs multiple runs with different modes: `scripts/end2end_test.sh`. See our [testing documentation](docs/testing.md) for more details.:

```bash
# Run full end-to-end tests
./scripts/end2end_test.sh apptainer /home/$USER/test_outputs
```

---

## Troubleshooting

### Common Issues

**Container build failures:**

```bash
# Check available space
df -h

# Check Apptainer version
apptainer --version

# Try building with more verbose output
apptainer build --fakeroot --verbose <container>.sif <container>.def
```

**Configuration issues:**

```bash
# Test Nextflow configuration
nextflow config

# Verify file paths exist
ls -la /path/to/models
ls -la /path/to/containers
```

**Test failures:**

```bash
# Check Nextflow work directory
ls -la work/

# Review error logs
cat .nextflow.log

# Clean and retry
nextflow clean -f
```

### Performance Optimization

**For better performance:**

- Store models on fast SSD storage
- Use shared filesystem for multi-user setups
- Adjust CPU/memory requests based on available resources
- Consider using local scratch space for intermediate files

---

## Next Steps

After successful installation:

1. **Read the [Getting Started Guide](GETTING-STARTED.md)** for your first protein design
2. **Review [ProteinDJ Modes](MODES.md)** to understand available options
3. **Configure [Parameters](docs/PARAMETERS.md)** for your specific needs

---

## Getting Help

If you encounter issues:

1. **Check the [Troubleshooting Guide](TROUBLESHOOTING.md)**
2. **Search existing [GitHub Issues](https://github.com/PapenfussLab/proteindj/issues)**
3. **Create a new issue** with detailed error information

---

**🎉 Congratulations!** You're now ready to DJ some proteins with ProteinDJ!

[⬅️ Back to Main README](../README.md)
