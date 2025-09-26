[🏠 ProteinDJ](../README.md) > End-To-End Testing

# ProteinDJ End-To-End Testing

## Overview

The `end2end_test.sh` script is a comprehensive testing framework for the ProteinDJ pipeline that automatically tests all supported pipeline modes and generates detailed reports. It helps ensure pipeline reliability and provides insights into the performance of different protein design modes. A full test takes ~40 minutes on 4 A30 GPUs.

## Features

- ✅ **Automated Testing**: Tests all 8 pipeline modes automatically
- 📊 **Comprehensive Reporting**: 

      - Generates detailed text reports
      - Tracks execution time and success rates
      - Analyzes generated files and directory structures
- 🔍 **Error Capture**: Captures and reports detailed error information


## Basic Usage

```bash
# Load nextflow module/environment if needed
module load nextflow/24.10.5

# Run tests with compute profile and output directory
./scripts/end2end_test.sh <compute_profile> <output_directory>

# Examples with different compute profiles
./scripts/end2end_test.sh milton /vast/scratch/users/$USER     # WEHI HPC cluster environment
./scripts/end2end_test.sh apptainer /home/$USER/test_outputs   # Apptainer containers
./scripts/end2end_test.sh singularity /tmp/proteindj_tests     # Singularity containers
```

### Prerequisites

1. **Nextflow**: Version 24.10.5 or later
2. **Environment**: Access to the compute environment specified
3. **Storage**: Sufficient space in the specified output directory for test outputs 

## Pipeline Modes Tested

The script tests the following ProteinDJ pipeline modes:

### Monomer Modes
- `monomer_denovo` - De novo protein design
- `monomer_foldcond` - Conditional folding design
- `monomer_motifscaff` - Motif-based scaffolding
- `monomer_partialdiff` - Partial diffusion design

### Binder-Specific Modes
- `binder_denovo` - De novo binder design
- `binder_foldcond` - Conditional binder folding
- `binder_motifscaff` - Binder motif scaffolding
- `binder_partialdiff` - Partial diffusion binder design

## Command Line Arguments

| Argument | Required | Description | Example |
|----------|----------|-------------|---------|
| `compute_profile` | Yes | The compute profile to use for testing | e.g. `milton`, `apptainer`, `singularity` |
| `output_directory` | Yes | Base directory where test outputs will be stored (< 50 MB) | e.g. `/vast/scratch/user/$USER`, `/home/$USER/tests` |

**Available Compute Profiles:**
- `milton` - WEHI HPC cluster environment
- `apptainer` - Uses Apptainer containers for execution
- `singularity` - Uses Singularity containers for execution

## Output Structure

The script creates a timestamped directory with the following structure:

```
test_results_YYYYMMDD_HHMMSS/
├── logs/
│   ├── test_execution.log          # Main execution log
│   ├── denovo.log                  # Individual mode logs
│   ├── binder_denovo.log
│   └── ...
└── analysis.txt                    # Comprehensive analysis and summary
```

## Report Types

### Text Analysis (`analysis.txt`)
- **Comprehensive report** combining all mode results
- **Summary statistics** with pass/fail counts and success rates
- **Output analysis** for each successful run
- **Performance metrics** and timing data

### Performance Metrics
- **Duration**: Time taken for each mode (in seconds)
- **Success Rate**: Percentage of modes that passed
- **Output Analysis**: Count of generated files (PDB, CSV) and number of designs in CSV files

### Output Directory

Test outputs are saved to:
```
<output_directory>/pdj_test_${mode}_${TIMESTAMP}
```

### Error Recovery

If tests fail:
1. Check the text analysis report (`analysis.txt`) for detailed information
2. Review individual log files in the `logs/` directory
3. Verify input data and configuration files
4. Ensure compute resources are available

---

**Note**: This testing framework is designed to validate the ProteinDJ pipeline functionality across all supported modes. Regular testing helps ensure pipeline reliability and catches regressions early in the development process.

[⬅️ Back to Main README](../README.md)