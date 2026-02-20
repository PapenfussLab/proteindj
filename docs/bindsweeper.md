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
   cd bindsweeper && uv tool install --reinstall-package bindsweeper . && cd ..
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
  input_pdb: "/path/to/protein.pdb"
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
  hotspot_residues:
    values:
      - "A56,A115,A123"
      - "A56,A115"
      - "A56"
```

#### Paired Parameters
Paired parameters allow you to sweep multiple parameters in lock-step (zipped), rather than as a Cartesian product. This is useful when parameters are inherently linked — for example, each target PDB has a corresponding MSA file.

```yaml
sweep_params:
  uncropped_target_pdb:
    values:
      - "input/protein1.pdb"
      - "input/protein2.pdb"
      - "input/protein3.pdb"
    paired_with:
      boltz_msa_path:
        - "input/msas/protein1.a3m"
        - "input/msas/protein2.a3m"
        - "input/msas/protein3.a3m"
```

**Key behaviours:**
- All lists in `paired_with` must have the same length as the primary `values` list
- Paired values are **zipped** (not crossed): the first PDB always runs with the first MSA, etc.
- You can pair multiple secondary parameters at once — just add more keys under `paired_with`
- Paired parameters are combined via **Cartesian product** with any other (non-paired) sweep parameters
- A paired parameter cannot also appear as a separate sweep parameter or a fixed parameter

**Example with paired + unpaired:**

With 3 paired targets and 2 noise scale values, BindSweeper generates 3 × 2 = 6 combinations:

```yaml
sweep_params:
  uncropped_target_pdb:
    values: ["protein1.pdb", "protein2.pdb", "protein3.pdb"]
    paired_with:
      boltz_msa_path: ["protein1.a3m", "protein2.a3m", "protein3.a3m"]
  rfd_noise_scale:
    values: [0.0, 0.1]
```

This produces:
| Combination | `uncropped_target_pdb` | `boltz_msa_path` | `rfd_noise_scale` |
|:-----------:|:----------------------:|:----------------:|:-----------------:|
| 1           | protein1.pdb              | protein1.a3m        | 0.0               |
| 2           | protein1.pdb              | protein1.a3m        | 0.1               |
| 3           | protein2.pdb              | protein2.a3m        | 0.0               |
| 4           | protein2.pdb              | protein2.a3m        | 0.1               |
| 5           | protein3.pdb              | protein3.a3m        | 0.0               |
| 6           | protein3.pdb              | protein3.a3m        | 0.1               |

## Example Configurations

### 1. Hotspots Sweep
Test different hotspot combinations for binder design:

```yaml
mode: binder_denovo
profile: milton

fixed_params:
  design_length: "60-100"
  input_pdb: "./benchmarkdata/5o45_pd-l1.pdb"
  rfd_noise_scale: 1.0
  rfd_ckpt_override: "complex_beta"

  af2_max_pae_interaction: 5
  af2_max_rmsd_binder_bndaln: 1
  af2_min_plddt_overall: 80

sweep_params:
  hotspot_residues:
    values:
      - null
      - "A56"
      - "A56,A115,A123"
```

### 2. Noise Scale Sweep
Explore different noise levels:

```yaml
mode: binder_denovo
profile: milton

fixed_params:
  design_length: "60-100"
  input_pdb: "./benchmarkdata/5o45_pd-l1.pdb"
  hotspot_residues: "A56,A115,A123"
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
  input_pdb: "./benchmarkdata/5o45_pd-l1.pdb"
  hotspot_residues: "A56,A115,A123"
  rfd_noise_scale: 0.0

sweep_params:
  rfd_scaffold_dir:
    values:
      - "./binderscaffolds/scaffolds_100_EHEEHE"
      - "./binderscaffolds/scaffolds_100_HEEHE"
      - "./binderscaffolds/scaffolds_100_HHH"
      - "./binderscaffolds/scaffolds_100_HHHH"
```

### 4. Multi-Target Paired Sweep
Sweep across multiple targets, each with a corresponding MSA file:

```yaml
mode: bindcraft_denovo
profile: milton

fixed_params:
  skip_fold_seq: true
  pred_method: "boltz"

sweep_params:
  uncropped_target_pdb:
    values:
      - "input/protein1.pdb"
      - "input/protein2.pdb"
      - "input/protein3.pdb"
    paired_with:
      boltz_msa_path:
        - "input/msas/protein1.a3m"
        - "input/msas/protein2.a3m"
        - "input/msas/protein3.a3m"
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
- `--resume`: Add -resume flag to Nextflow commands to use cached tasks where inputs haven't changed
- `--parallel`: Execute parameter combinations in parallel (each with isolated Nextflow cache)
- `--max-parallel N`: Maximum number of parallel Nextflow runs (default: 4)
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
- **Paired target sweep**: N paired targets → N × 4 total sequences (e.g., 3 targets → 12 sequences)
- **Paired + unpaired sweep**: N paired × M unpaired → N × M × 4 total sequences

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

### Resume Interrupted Sweeps
```bash
bindsweeper --resume --config sweep.yaml
```

When using `--resume`, BindSweeper adds the `-resume` flag to all Nextflow commands. This enables Nextflow's caching mechanism, which:
- Skips tasks that have already completed successfully
- Re-runs only tasks where inputs, parameters, or scripts have changed
- Automatically detects parameter changes and re-executes affected tasks
- Preserves computational resources by avoiding redundant work

**Use cases for `--resume`:**
- **Interrupted runs**: Cluster timeouts, manual cancellation, or system failures
- **Iterative development**: Testing bug fixes in later pipeline stages while reusing early stage results
- **Parameter refinement**: Re-running with modified filtering thresholds while keeping expensive fold/sequence generation cached

**Important notes:**
- Nextflow determines what to cache based on task hashes (inputs, scripts, parameters, containers)
- If you modify any parameters (fixed or swept), Nextflow will automatically detect this and re-run affected tasks
- The `.nextflow/cache/` and `work/` directories must be preserved for resume to work
- Resume works at the task level within each parameter combination, not at the combination level

### Parallel Execution
```bash
# Execute up to 4 combinations in parallel (default)
bindsweeper --parallel --config sweep.yaml

# Control the maximum number of parallel runs
bindsweeper --parallel --max-parallel 8 --config sweep.yaml

# Combine with resume for robust parallel execution
bindsweeper --parallel --resume --max-parallel 6 --config sweep.yaml
```

When using `--parallel`, BindSweeper executes multiple parameter combinations concurrently, which:
- Runs multiple independent Nextflow pipelines simultaneously
- Provides isolated cache directories for each combination (prevents cache conflicts)
- Improves overall throughput on systems with available GPU and CPU resources
- Each Nextflow run still internally parallelizes tasks as normal
- Maintains proper resource allocation through the cluster scheduler

**Benefits of parallel execution:**
- **Faster completion**: Leverage multiple GPUs and CPUs concurrently across different parameter combinations
- **Better resource utilization**: Keep GPUs busy while other combinations process CPU tasks
- **Fault tolerance**: Failed combinations don't block others from completing
- **Natural batching**: Combinations complete and release resources as they finish

**Use cases for `--parallel`:**
- **Large parameter sweeps**: When testing 8+ parameter combinations
- **GPU-rich clusters**: Systems with multiple GPUs available for concurrent use
- **Mixed GPU/CPU workloads**: Combinations naturally interleave GPU and CPU-intensive stages
- **Time-sensitive projects**: Need results faster than sequential execution allows

**Resource considerations:**
- `--max-parallel` should be ≤ number of available GPUs to prevent GPU contention
- Each combination spawns its own Nextflow session with internal task parallelization
- Monitor cluster queue status to ensure combinations get scheduled efficiently
- Disk I/O can become a bottleneck with too many parallel runs
- Consider available memory: multiple AF2/Boltz runs require significant RAM per GPU

**Important notes:**
- Each combination uses isolated cache in `<output_dir>/.nextflow_cache/`
- Combinations are truly independent - no shared state or locks
- Works seamlessly with `--resume` flag for robust parallel execution
- Log output is captured per combination in respective output directories
- Works with `--quick-test` flag for rapid validation of large parameter sweeps

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
3. **Use Resume for Long Runs**: Always use `--resume` for multi-hour sweeps to recover from interruptions
4. **Leverage Parallel Execution**: Use `--parallel` for large sweeps on GPU-rich clusters to improve throughput
5. **Monitor Resources**: Large parameter sweeps can be resource-intensive, especially when running in parallel
6. **Optimize `--max-parallel`**: Set to match available GPUs (e.g., `--max-parallel 8` on systems with 8+ GPUs)
7. **Combine Resume and Parallel**: Use both `--resume --parallel` for robust and efficient large-scale sweeps
8. **Organize Results**: Use descriptive output directory names
9. **Check Dependencies**: Ensure ProteinDJ and required tools are installed
10. **Preserve Cache Directories**: Keep `.nextflow/` and `work/` directories to enable resume functionality
11. **Monitor Parallel Runs**: Check cluster queue to ensure combinations are being scheduled appropriately

## Troubleshooting

### Common Issues

1. **Missing nextflow.config**: Ensure the file exists in current or parent directory
2. **Invalid parameters**: Check YAML syntax and parameter names
3. **Resource constraints**: Monitor system resources during execution
4. **Path issues**: Use absolute paths for input files
5. **Resume not working**: Ensure `.nextflow/cache/` and `work/` directories exist and haven't been cleaned
6. **Unexpected re-execution with resume**: Nextflow detects input/parameter/script changes and correctly re-runs affected tasks

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