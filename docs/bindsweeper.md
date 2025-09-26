[🏠 ProteinDJ](../README.md) > Bindsweeper User Guide

# BindSweeper User Guide

BindSweeper is a Python CLI tool that automates parameter sweeps for ProteinDJ/RFdiffusion workflows. It allows you to systematically explore different parameter combinations for protein binding analysis and design generation.

## Overview

BindSweeper works by:
1. Reading sweep configuration from YAML files
2. Generating parameter combinations to test
3. Creating Nextflow profiles for each combination
4. Executing the ProteinDJ pipeline with different parameters
5. Processing and organizing results

<img src="../img/bindsweeper_workflow.png" height="500">

## Quick Start

### Installation

> Please note: BindSweeper requires that ProteinDJ is installed and configured. If you have not installed ProteinDJ yet, follow the instructions [here](../docs/installation.md) 

1. **Install uv** (if not already installed):
   ```bash
   # On Linux/macOS
   curl -LsSf https://astral.sh/uv/install.sh | sh
   
   # Or using pip
   pip install uv
   ```

2. **Install BindSweeper**:
   ```bash
   # Navigate to the bindsweeper subdirectory inside the ProteinDJ installation directory
   cd <path-to-proteindj>/bindsweeper
   
   # Install BindSweeper globally
   uv tool install .
   cd ..
   
   # Verify installation
   bindsweeper --help
   ```

3. **Updating BindSweeper**:
   ```bash
   # To update bindsweeper to the latest version use the following commands
   git pull
   uv tool update bindsweeper --force-reinstall
   ```

### Basic Usage

1. Create a sweep configuration file (e.g., `sweep.yaml`)
2. Run the sweep:
```bash
bindsweeper
```

### Running on WEHI Systems

When running BindSweeper on WEHI systems, use `screen` to prevent disconnections during long-running jobs:

1. **Start a screen session**:
   ```bash
   screen -S bindsweeper_run
   ```

2. **Run BindSweeper within the screen session**:
   ```bash
   cd /path/to/your/work/directory
   bindsweeper --config sweep.yaml --output-dir results/
   ```

3. **Detach from screen** (keeps the job running):
   - Press `Ctrl+A`, then `D`

4. **Reattach to the screen session**:
   ```bash
   screen -r bindsweeper_run
   ```

5. **List all screen sessions**:
   ```bash
   screen -ls
   ```

6. **Terminate a screen session** (when job is complete):
   ```bash
   # From within the screen session
   exit
   
   # Or kill from outside
   screen -X -S bindsweeper_run quit
   ```


## Configuration Files

BindSweeper uses YAML configuration files to define parameter sweeps. Here are the key components:

### Basic Structure

```yaml
mode: binder_denovo  # or binder_foldcond
profile: milton      # Nextflow profile to use

# Parameters that remain constant across all runs
fixed_params:
  rfd_input_pdb: "/path/to/protein.pdb"
  rfd_noise_scale: 0.0

# Parameters to sweep across different values
sweep_params:
  parameter_name:
    type: range        # or values
    min: 0.0          # for range type
    max: 1.0
    step: 0.1
  # or
  other_parameter:
    values:           # for values type
      - "value1"
      - "value2"

# Results processing configuration
results_config:
  rank_dirname: results
  results_dirname: best_designs
  csv_filename: best_designs.csv
  output_csv: sweep_results.csv
  pdb_output_dir: sweep_designs
  zip_results: true
```

### Supported Modes

1. **`binder_denovo`**: De novo binder design
2. **`binder_foldcond`**: Fold-conditioned binder design

### Parameter Types

#### Range Parameters
```yaml
sweep_params:
  rfd_noise_scale:
    type: range
    min: 0.0
    max: 0.1
    step: 0.05
```

#### Value List Parameters
```yaml
sweep_params:
  rfd_hotspots:
    values:
      - "[A56,A115,A123]"
      - "[A56,A115]"
      - "[A56]"
```

## Example Configurations

### 1. Hotspots Sweep
Test different hotspot combinations for binder design:

```yaml
mode: binder_denovo
profile: milton

fixed_params:
  rfd_contigs: "[A17-131/0 60-100]"
  rfd_input_pdb: "./benchmarkdata/5o45_pd-l1.pdb"
  rfd_noise_scale: 1.0
  rfd_ckpt_override: "complex_beta"

  af2_max_pae_interaction: 5
  af2_max_rmsd_binder_bndaln: 1
  af2_min_plddt_overall: 80

sweep_params:
  rfd_hotspots:
    values:
      - null
      - "[A56]"
      - "[A56,A115,A123]"
```

### 2. Noise Scale Sweep
Explore different noise levels:

```yaml
mode: binder_denovo
profile: milton

fixed_params:
  rfd_contigs: "[A17-131/0 60-100]"
  rfd_input_pdb: "./benchmarkdata/5o45_pd-l1.pdb"
  rfd_hotspots: "[A56,A115,A123]"
  rfd_ckpt_override: "complex_beta"

sweep_params:
  rfd_noise_scale:
    type: range
    min: 0.0
    max: 0.1
    step: 0.05
```

### 3. Scaffold Directory Sweep
Test different scaffold sets:

```yaml
mode: binder_foldcond
profile: milton

fixed_params:
  rfd_input_pdb: "./benchmarkdata/5o45_pd-l1.pdb"
  rfd_hotspots: "[A56,A115,A123]"
  rfd_noise_scale: 0.0
  rfd_target_ss: "./benchmarkdata/5o45_pd-l1_ss.pt"
  rfd_target_adj:  "./benchmarkdata/5o45_pd-l1_adj.pt"

sweep_params:
  rfd_scaffold_dir:
    values:
      - "./binderscaffolds/scaffolds_100_EHEEHE"
      - "./binderscaffolds/scaffolds_100_HEEHE"
      - "./binderscaffolds/scaffolds_100_HHH"
      - "./binderscaffolds/scaffolds_100_HHHH"
```

## Command Line Options

### Basic Options
- `--help`: Show help message
- `--version`: Show version information
- `--debug`: Enable debug logging
- `--dry-run`: Print commands without executing them

### File and Directory Options
- `--config PATH`: Path to sweep configuration YAML file
- `--output-dir PATH`: Output directory for results
- `--pipeline-path PATH`: Path to the ProteinDJ main.nf file
- `--nextflow-config PATH`: Path to nextflow.config file

### Execution Options
- `--skip-sweep`: Skip parameter sweep and only process results
- `--continue-on-error`: Continue if individual parameter sweeps fail
- `--quick-test`: Run quick test with reduced parameters first
- `--auto-update`: Automatically sync/update dependencies

### Automation Options
- `-y, --yes-to-all`: Skip confirming files automatically

## Usage Examples

### Basic Sweep
```bash
bindsweeper --config sweep.yaml --output-dir ./results
```

### Quick Test First
```bash
bindsweeper --quick-test --config sweep.yaml
```

#### Quick Test Design Counts

When using `--quick-test`, BindSweeper reduces the number of designs to validate configurations efficiently:

- **Designs per combination**: 2 designs
- **Sequences per design**: 2 sequences
- **Total sequences per combination**: 4 sequences

For the standard quick test configurations:
- **Noise sweep**: 3 combinations → 12 total sequences (3 × 4)
- **Hotspots sweep**: 3 combinations → 12 total sequences (3 × 4) 
- **Scaffold sweep**: 3 combinations → 12 total sequences (3 × 4)
- **Multi-dimensional sweep**: 4 combinations → 16 total sequences (4 × 4)

This allows rapid validation of your parameter sweep configuration before committing to a full run with the default number of designs.

### Dry Run (Preview Commands)
```bash
bindsweeper --dry-run --config sweep.yaml
```

### Debug Mode
```bash
bindsweeper --debug --config sweep.yaml
```

### Continue on Errors
```bash
bindsweeper --continue-on-error --config sweep.yaml
```

## File Structure

BindSweeper generates organised output directories:

```
results/
├── combination1_param1_val1_param2_val2/
│   ├── results/
│   ├── best_designs/
│   └── logs/
├── combination2_param1_val3_param2_val4/
│   └── ...
├── sweep_results.csv        # Combined results
├── sweep_designs/          # Best designs from all runs
└── bindsweeper.log        # Execution log
```

## Results Processing

BindSweeper automatically:
1. Collects results from all parameter combinations
2. Ranks designs based on metrics
3. Copies best designs to a central directory
4. Generates summary CSV files
5. Optionally creates ZIP archives


## Tips and Best Practices

1. **Start with Quick Tests**: Use `--quick-test` to validate configuration
2. **Use Dry Runs**: Preview commands with `--dry-run` before execution
3. **Monitor Resources**: Large parameter sweeps can be resource-intensive
4. **Organize Results**: Use descriptive output directory names
5. **Check Dependencies**: Ensure ProteinDJ and required tools are installed

## Troubleshooting

### Common Issues

1. **Missing nextflow.config**: Ensure the file exists in current or parent directory
2. **Invalid parameters**: Check YAML syntax and parameter names
3. **Resource constraints**: Monitor system resources during execution
4. **Path issues**: Use absolute paths for input files

### Debug Information

Enable debug logging to see detailed execution information:
```bash
bindsweeper --debug
```

### Getting Help

For issues or questions:
- Check the log files in the output directory
- Use `--debug` flag for detailed information
- Review configuration file syntax
- Ensure all dependencies are installed

## Advanced Usage

### Custom Pipeline Paths
```bash
bindsweeper --pipeline-path /custom/path/main.nf
```

### Processing Existing Results
```bash
bindsweeper --skip-sweep --output-dir existing_results/
```

[⬅️ Back to Main README](../README.md)