#!/usr/bin/env python3
"""Configuration parsing and validation for sweep tool."""

import json
import os
import re
from dataclasses import dataclass, field
from typing import Any, Optional
from pathlib import Path

import yaml

from .binder_validator import BinderValidator
from .sweep_types import SweepType, create_sweep


@dataclass
class ResultsConfig:
    """Configuration for results processing."""

    rank_dirname: str = "results"
    results_dirname: str = "best_designs"
    csv_filename: str = "best_designs.csv"
    output_csv: str = "sweep_results.csv"
    pdb_output_dir: str = "sweep_designs"
    zip_results: bool = True


@dataclass
class SweepConfig:
    """Main configuration for parameter sweep."""

    mode: str
    fixed_params: dict[str, Any]
    sweep_params: dict[str, SweepType]
    results_config: ResultsConfig = field(default_factory=ResultsConfig)
    pipeline_path: str = "./main.nf"  # Default for proteinDJ root directory
    profile: Optional[str] = None  # Additional nextflow profile to include

    @classmethod
    def from_yaml(
        cls, yaml_path: str, schema_path: Optional[str] = None
    ) -> "SweepConfig":
        """Load configuration from YAML file."""
        with open(yaml_path) as f:
            data = yaml.safe_load(f)

        # Extract mode first
        mode = data.get("mode")
        if not mode:
            raise ValueError("'mode' is required in configuration")

        # Validate against schema if provided
        if schema_path and os.path.exists(schema_path):
            validate_params_against_schema(data, schema_path)
        else:
            # Use binder validator for binder modes as fallback
            if mode.startswith("binder_"):
                try:
                    validator = BinderValidator()
                    validator.validate_config(
                        mode, data.get("sweep_params", {}), data.get("fixed_params", {})
                    )
                except ValueError as e:
                    raise ValueError(str(e)) from e

        # Extract sweep parameters
        sweep_params = {}
        for param_name, param_def in data.get("sweep_params", {}).items():
            sweep_params[param_name] = create_sweep(param_def)

        # Extract fixed parameters
        fixed_params = data.get("fixed_params", {})

        # Validate that out_dir is not in fixed_params
        if "out_dir" in fixed_params:
            raise ValueError(
                "out_dir should not be specified in fixed_params. "
                "Use --output-dir command line option or modify the out_dir parameter in nextflow.config instead."
            )

        # Load results config from default location or user override
        results_config = cls._load_results_config(data.get("results_config"))

        # Extract pipeline path
        pipeline_path = data.get("pipeline_path", "./main.nf")

        # Extract profile
        profile = data.get("profile")

        return cls(
            mode=mode,
            fixed_params=fixed_params,
            sweep_params=sweep_params,
            results_config=results_config,
            pipeline_path=pipeline_path,
            profile=profile,
        )

    @classmethod
    def _load_results_config(cls, user_config: Optional[dict] = None) -> ResultsConfig:
        """Load results config from default file, with optional user overrides."""
        # Get path to default config file (in same package)
        default_config_path = Path(__file__).parent / "default_results_config.yaml"
        
        # Load default config
        default_config = {}
        if default_config_path.exists():
            with open(default_config_path) as f:
                default_data = yaml.safe_load(f)
                default_config = default_data.get("results_config", {})
        
        # Merge with user overrides if provided
        if user_config:
            default_config.update(user_config)
            
        return ResultsConfig(**default_config)


def validate_params_against_schema(config: dict, schema_path: str) -> None:
    """Validate parameters against JSON schema."""
    with open(schema_path) as f:
        schema = json.load(f)

    mode = config.get("mode", "")
    
    # Extract parameter definitions from schema based on structure
    param_defs = {}
    if mode in schema:
        # Binder schema structure: mode -> parameter definitions
        param_defs = schema[mode]
    else:
        # Nextflow schema structure: definitions -> sections -> properties
        for section in schema.get("definitions", {}).values():
            for prop_name, prop_def in section.get("properties", {}).items():
                param_defs[prop_name] = prop_def

    # Validate fixed params
    for param_name, value in config.get("fixed_params", {}).items():
        validate_param_value(param_name, value, param_defs.get(param_name))

    # Validate sweep params
    for param_name, sweep_def in config.get("sweep_params", {}).items():
        param_def = param_defs.get(param_name)
        if param_def:
            # Create sweep and validate each value
            sweep = create_sweep(sweep_def)
            for value in sweep.generate_values():
                validate_param_value(param_name, value, param_def)


def validate_param_value(
    param_name: str, value: Any, param_def: Optional[dict]
) -> None:
    """Validate a single parameter value against its schema definition."""
    if not param_def:
        return  # No schema definition, skip validation

    param_type = param_def.get("type")

    # Type validation
    if param_type == "string" and not isinstance(value, str):
        raise ValueError(f"{param_name} must be a string, got {type(value)}")
    elif param_type == "integer" and not isinstance(value, int):
        # Allow floats that are actually integers
        if isinstance(value, float) and value.is_integer():
            value = int(value)  # Convert for further validation
        else:
            raise ValueError(f"{param_name} must be an integer, got {type(value)}")
    elif param_type == "number" and not isinstance(value, (int, float)):
        raise ValueError(f"{param_name} must be a number, got {type(value)}")
    elif param_type == "boolean" and not isinstance(value, bool):
        raise ValueError(f"{param_name} must be a boolean, got {type(value)}")

    # Enum validation
    if "enum" in param_def and value not in param_def["enum"]:
        raise ValueError(
            f"{param_name} must be one of {param_def['enum']}, got {value}"
        )

    # Pattern validation
    if "pattern" in param_def and isinstance(value, str):
        if not re.match(param_def["pattern"], value):
            error_msg = f"{param_name} value '{value}' does not match pattern {param_def['pattern']}"
            
            # Add description from schema
            if "description" in param_def:
                error_msg += f"\n{param_def['description']}"
            
            # Add example from schema
            if "example" in param_def:
                error_msg += f"\nExample: {param_def['example']}"
            
            raise ValueError(error_msg)

    # Range validation
    if "minimum" in param_def and isinstance(value, (int, float)):
        if value < param_def["minimum"]:
            raise ValueError(
                f"{param_name} must be >= {param_def['minimum']}, got {value}"
            )

    if "maximum" in param_def and isinstance(value, (int, float)):
        if value > param_def["maximum"]:
            raise ValueError(
                f"{param_name} must be <= {param_def['maximum']}, got {value}"
            )




def parse_nextflow_value(value_str: str) -> Any:
    """Parse a Nextflow configuration value string."""
    value_str = value_str.strip()

    # Remove trailing comments
    value_str = re.sub(r"//.*$", "", value_str, flags=re.MULTILINE).strip()

    # Handle null
    if value_str == "null":
        return None

    # Handle booleans
    if value_str == "true":
        return True
    elif value_str == "false":
        return False

    # Handle numbers
    try:
        if "." in value_str:
            return float(value_str)
        else:
            return int(value_str)
    except ValueError:
        pass

    # Handle strings (remove quotes)
    if (value_str.startswith("'") and value_str.endswith("'")) or (
        value_str.startswith('"') and value_str.endswith('"')
    ):
        return value_str[1:-1]

    # Handle lists (basic support)
    if value_str.startswith("[") and value_str.endswith("]"):
        # Very basic list parsing - doesn't handle nested structures
        inner = value_str[1:-1].strip()
        if not inner:
            return []
        # Split by comma and clean up
        items = [item.strip().strip("'\"") for item in inner.split(",")]
        return items

    # Default to string
    return value_str
