"""Tests for sweep configuration loading and validation."""

import json
from pathlib import Path

import pytest

from bindsweeper.sweep_config import (
    ResultsConfig,
    SweepConfig,
    parse_nextflow_config,
    parse_nextflow_value,
    validate_param_value,
    validate_params_against_schema,
)
from bindsweeper.sweep_types import ListSweep, RangeSweep


class TestResultsConfig:
    """Test ResultsConfig dataclass."""

    def test_default_values(self):
        """Test default configuration values."""
        config = ResultsConfig()
        assert config.rank_dirname == "rank"
        assert config.extract_dirname == "extract"
        assert config.results_dirname == "results"
        assert config.csv_filename == "best.csv"
        assert config.output_csv == "merged_best.csv"
        assert config.pdb_output_dir == "merged_best_designs"
        assert config.zip_results is True

    def test_custom_values(self):
        """Test custom configuration values."""
        config = ResultsConfig(rank_dirname="custom_rank", zip_results=False)
        assert config.rank_dirname == "custom_rank"
        assert config.zip_results is False
        assert config.csv_filename == "best.csv"  # Default preserved


class TestSweepConfig:
    """Test SweepConfig loading and validation."""

    def test_from_yaml_basic(self, config_files):
        """Test loading basic YAML configuration."""
        config = SweepConfig.from_yaml(config_files["sweep_yaml"])

        assert config.mode == "binder_denovo"
        assert "rfd_contigs" in config.fixed_params
        assert "rfd_noise_scale" in config.sweep_params
        assert "rfd_ckpt_override" in config.sweep_params
        assert "hotspot_residues" in config.sweep_params

    def test_from_yaml_missing_mode(self, temp_dir):
        """Test loading YAML without required mode."""
        yaml_content = """
fixed_params:
  num_designs: 4
        """
        yaml_path = Path(temp_dir) / "invalid.yaml"
        yaml_path.write_text(yaml_content)

        with pytest.raises(ValueError, match="'mode' is required"):
            SweepConfig.from_yaml(str(yaml_path))

    def test_from_yaml_empty_sweep_params(self, temp_dir):
        """Test loading YAML with no sweep parameters."""
        yaml_content = """
mode: denovo
fixed_params:
  num_designs: 4
        """
        yaml_path = Path(temp_dir) / "no_sweep.yaml"
        yaml_path.write_text(yaml_content)

        config = SweepConfig.from_yaml(str(yaml_path))
        assert config.mode == "denovo"
        assert len(config.sweep_params) == 0
        assert "num_designs" in config.fixed_params

    def test_from_yaml_custom_results_config(self, temp_dir):
        """Test loading YAML with custom results configuration."""
        yaml_content = """
mode: denovo
results_config:
  rank_dirname: custom_rank
  zip_results: false
        """
        yaml_path = Path(temp_dir) / "custom_results.yaml"
        yaml_path.write_text(yaml_content)

        config = SweepConfig.from_yaml(str(yaml_path))
        assert config.results_config.rank_dirname == "custom_rank"
        assert config.results_config.zip_results is False

    def test_sweep_param_types(self, config_files):
        """Test that sweep parameters are created with correct types."""
        config = SweepConfig.from_yaml(config_files["sweep_yaml"])

        # Check range sweep
        assert isinstance(config.sweep_params["rfd_noise_scale"], RangeSweep)

        # Check list sweeps
        assert isinstance(config.sweep_params["rfd_ckpt_override"], ListSweep)
        assert isinstance(config.sweep_params["hotspot_residues"], ListSweep)


class TestValidateParamsAgainstSchema:
    """Test parameter validation against schema."""

    def test_validate_with_schema(self, temp_dir, sample_nextflow_schema):
        """Test validation against schema file."""
        schema_path = Path(temp_dir) / "schema.json"
        schema_path.write_text(json.dumps(sample_nextflow_schema))

        config_data = {
            "fixed_params": {"design_mode": "denovo", "num_designs": 8},
            "sweep_params": {"rfd_noise_scale": {"min": 0.0, "max": 1.0, "step": 0.5}},
        }

        # Should not raise exception
        validate_params_against_schema(config_data, str(schema_path))

    def test_validate_invalid_enum(self, temp_dir, sample_nextflow_schema):
        """Test validation with invalid enum value."""
        schema_path = Path(temp_dir) / "schema.json"
        schema_path.write_text(json.dumps(sample_nextflow_schema))

        config_data = {
            "fixed_params": {
                "design_mode": "invalid_mode"  # Not in enum
            }
        }

        with pytest.raises(ValueError, match="must be one of"):
            validate_params_against_schema(config_data, str(schema_path))

    def test_validate_out_of_range(self, temp_dir, sample_nextflow_schema):
        """Test validation with out-of-range values."""
        schema_path = Path(temp_dir) / "schema.json"
        schema_path.write_text(json.dumps(sample_nextflow_schema))

        config_data = {
            "sweep_params": {
                "num_designs": {
                    "type": "range",
                    "min": 100,
                    "max": 2000,
                    "step": 500,
                }  # Max exceeds 1000
            }
        }

        with pytest.raises(ValueError, match="must be <= 1000"):
            validate_params_against_schema(config_data, str(schema_path))


class TestValidateParamValue:
    """Test individual parameter value validation."""

    def test_validate_string_type(self):
        """Test string type validation."""
        param_def = {"type": "string"}

        # Valid
        validate_param_value("test_param", "valid_string", param_def)

        # Invalid
        with pytest.raises(ValueError, match="must be a string"):
            validate_param_value("test_param", 123, param_def)

    def test_validate_integer_type(self):
        """Test integer type validation."""
        param_def = {"type": "integer"}

        # Valid
        validate_param_value("test_param", 42, param_def)

        # Invalid
        with pytest.raises(ValueError, match="must be an integer"):
            validate_param_value("test_param", 3.14, param_def)

        with pytest.raises(ValueError, match="must be an integer"):
            validate_param_value("test_param", "not_int", param_def)

    def test_validate_number_type(self):
        """Test number type validation."""
        param_def = {"type": "number"}

        # Valid - int and float
        validate_param_value("test_param", 42, param_def)
        validate_param_value("test_param", 3.14, param_def)

        # Invalid
        with pytest.raises(ValueError, match="must be a number"):
            validate_param_value("test_param", "not_number", param_def)

    def test_validate_enum(self):
        """Test enum validation."""
        param_def = {"enum": ["option1", "option2", "option3"]}

        # Valid
        validate_param_value("test_param", "option1", param_def)

        # Invalid
        with pytest.raises(ValueError, match="must be one of"):
            validate_param_value("test_param", "invalid_option", param_def)

    def test_validate_range(self):
        """Test range validation."""
        param_def = {"type": "number", "minimum": 0, "maximum": 10}

        # Valid
        validate_param_value("test_param", 5, param_def)
        validate_param_value("test_param", 0, param_def)  # Boundary
        validate_param_value("test_param", 10, param_def)  # Boundary

        # Invalid - below minimum
        with pytest.raises(ValueError, match="must be >= 0"):
            validate_param_value("test_param", -1, param_def)

        # Invalid - above maximum
        with pytest.raises(ValueError, match="must be <= 10"):
            validate_param_value("test_param", 11, param_def)

    def test_validate_no_definition(self):
        """Test validation with no parameter definition."""
        # Should not raise exception
        validate_param_value("test_param", "any_value", None)


class TestParseNextflowConfig:
    """Test parsing Nextflow configuration files."""

    def test_parse_basic_params(self, config_files):
        """Test parsing basic parameters."""
        params = parse_nextflow_config(config_files["nextflow_config"])

        # Check that some expected params were parsed
        # The exact params may vary based on the regex parsing
        assert isinstance(params, dict)
        # Check for any numeric params that should be parsed
        numeric_params = [v for v in params.values() if isinstance(v, (int, float))]
        assert len(numeric_params) > 0  # Should parse at least some numbers

    def test_parse_nonexistent_file(self):
        """Test parsing non-existent file."""
        params = parse_nextflow_config("/nonexistent/path")
        assert params == {}

    def test_parse_empty_file(self, temp_dir):
        """Test parsing empty file."""
        empty_config = Path(temp_dir) / "empty.config"
        empty_config.write_text("")

        params = parse_nextflow_config(str(empty_config))
        assert params == {}


class TestParseNextflowValue:
    """Test parsing individual Nextflow values."""

    def test_parse_null(self):
        """Test parsing null values."""
        assert parse_nextflow_value("null") is None

    def test_parse_boolean(self):
        """Test parsing boolean values."""
        assert parse_nextflow_value("true") is True
        assert parse_nextflow_value("false") is False

    def test_parse_integer(self):
        """Test parsing integer values."""
        assert parse_nextflow_value("42") == 42
        assert parse_nextflow_value("-10") == -10

    def test_parse_float(self):
        """Test parsing float values."""
        assert parse_nextflow_value("3.14") == 3.14
        assert parse_nextflow_value("-2.5") == -2.5

    def test_parse_string_quoted(self):
        """Test parsing quoted strings."""
        assert parse_nextflow_value("'hello'") == "hello"
        assert parse_nextflow_value('"world"') == "world"

    def test_parse_string_unquoted(self):
        """Test parsing unquoted strings."""
        assert parse_nextflow_value("unquoted") == "unquoted"

    def test_parse_list_empty(self):
        """Test parsing empty list."""
        assert parse_nextflow_value("[]") == []

    def test_parse_list_simple(self):
        """Test parsing simple list."""
        result = parse_nextflow_value("['a', 'b', 'c']")
        assert result == ["a", "b", "c"]

    def test_parse_with_comments(self):
        """Test parsing values with trailing comments."""
        assert parse_nextflow_value("42 // comment") == 42
        assert parse_nextflow_value("'test' // another comment") == "test"
