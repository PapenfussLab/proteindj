"""Pytest configuration and shared fixtures."""

import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir


@pytest.fixture
def sample_sweep_yaml():
    """Sample sweep.yaml content."""
    return """
mode: binder_denovo

fixed_params:
  rfd_contigs: "[A17-145/0 50-100]"
  input_pdb: "/path/to/target.pdb"
  mpnn_temperature: 0.0001
  num_designs: 4
  seqs_per_design: 2

sweep_params:
  rfd_noise_scale:
    type: range
    min: 0.0
    max: 1.0
    step: 0.5

  rfd_ckpt_override:
    values:
      - complex_base
      - complex_beta

  hotspot_residues:
    values:
      - "B208,B232,B239"
      - "B210,B234"
      - null

results_config:
  rank_dirname: rank
  extract_dirname: extract
  results_dirname: results
  csv_filename: best.csv
  output_csv: merged_best.csv
  pdb_output_dir: merged_best_designs
  zip_results: true
"""


@pytest.fixture
def sample_nextflow_config():
    """Sample nextflow.config content."""
    return """
params {
    // Essential parameters
    design_mode = null
    num_designs = 8
    seqs_per_design = 8
    out_dir = "/vast/scratch/users/$USER/pdj"

    // Mode-specific parameters
    rfd_contigs = null
    input_pdb = null
    hotspot_residues = null
    rfd_noise_scale = null

    // Advanced parameters
    mpnn_temperature = 0.1
    fampnn_temperature = 0.1
    rfd_ckpt_override = null
}

profiles {
    test {
        params {
            num_designs = 4
            seqs_per_design = 2
        }
    }
}
"""


@pytest.fixture
def sample_nextflow_schema():
    """Sample nextflow schema JSON."""
    return {
        "definitions": {
            "essential_parameters": {
                "properties": {
                    "design_mode": {
                        "type": "string",
                        "enum": ["denovo", "binder_denovo", "foldconditioning"],
                    },
                    "num_designs": {
                        "type": "integer",
                        "minimum": 1,
                        "maximum": 1000,
                    },
                    "rfd_noise_scale": {
                        "type": "number",
                        "minimum": 0.0,
                        "maximum": 2.0,
                    },
                }
            }
        }
    }


@pytest.fixture
def config_files(temp_dir, sample_sweep_yaml, sample_nextflow_config):
    """Create temporary config files."""
    sweep_path = Path(temp_dir) / "sweep.yaml"
    nextflow_path = Path(temp_dir) / "nextflow.config"

    sweep_path.write_text(sample_sweep_yaml)
    nextflow_path.write_text(sample_nextflow_config)

    return {
        "sweep_yaml": str(sweep_path),
        "nextflow_config": str(nextflow_path),
        "temp_dir": temp_dir,
    }
