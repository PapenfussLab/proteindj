#!/usr/bin/env python3
"""Validator for binder mode parameters using the simplified schema."""

import json
import os
import re
from typing import Any


class BinderValidator:
    """Validator for binder mode parameter sweeps."""

    def __init__(self, schema_path: str = None):
        """Initialize validator with schema."""
        if schema_path is None:
            schema_path = os.path.join(os.path.dirname(__file__), "binder_schema.json")

        with open(schema_path) as f:
            self.schema = json.load(f)

    def validate_parameter(self, mode: str, param_name: str, value: Any) -> bool:
        """Validate a single parameter value."""
        if mode not in self.schema:
            raise ValueError(
                f"Unknown mode: {mode}. Supported modes: {list(self.schema.keys())}"
            )

        if param_name not in self.schema[mode]:
            raise ValueError(
                f"Parameter '{param_name}' not supported for mode '{mode}'. "
                f"Supported parameters: {list(self.schema[mode].keys())}"
            )

        param_schema = self.schema[mode][param_name]

        # No preprocessing - let validation catch incorrect values

        # Type validation - handle both single types and arrays of types
        expected_types = param_schema["type"]
        if not isinstance(expected_types, list):
            expected_types = [expected_types]
            
        type_valid = False
        for expected_type in expected_types:
            if expected_type == "string" and isinstance(value, str):
                type_valid = True
                break
            elif expected_type == "null" and (value is None or value == "null"):
                type_valid = True
                break
            elif expected_type == "number" and isinstance(value, (int, float)):
                type_valid = True
                break
            elif expected_type == "integer" and isinstance(value, int):
                type_valid = True
                break
                
        if not type_valid:
            type_names = [t if t != "null" else "None" for t in expected_types]
            raise ValueError(
                f"{param_name} must be one of types {type_names}, got {type(value).__name__}"
            )
        
        # Early return for null values - no further validation needed
        if value is None or value == "null":
            return True

        if param_schema.get("type") == "number" or "number" in expected_types:
            if not isinstance(value, (int, float)):
                raise ValueError(
                    f"{param_name} must be a number, got {type(value).__name__}"
                )
        elif param_schema["type"] == "integer":
            if not isinstance(value, int):
                raise ValueError(
                    f"{param_name} must be an integer, got {type(value).__name__}"
                )

        # Enum validation
        if "enum" in param_schema and value not in param_schema["enum"]:
            raise ValueError(
                f"{param_name} must be one of {param_schema['enum']}, got '{value}'"
            )

        # Range validation
        if "minimum" in param_schema and isinstance(value, (int, float)):
            if value < param_schema["minimum"]:
                raise ValueError(
                    f"{param_name} must be >= {param_schema['minimum']}, got {value}"
                )

        if "maximum" in param_schema and isinstance(value, (int, float)):
            if value > param_schema["maximum"]:
                raise ValueError(
                    f"{param_name} must be <= {param_schema['maximum']}, got {value}"
                )

        # Pattern validation
        if "pattern" in param_schema and isinstance(value, str):
            if not re.match(param_schema["pattern"], value):
                error_msg = f"{param_name} value '{value}' does not match pattern {param_schema['pattern']}"
                
                # Add description from schema
                if "description" in param_schema:
                    error_msg += f"\n{param_schema['description']}"
                
                # Add example from schema
                if "example" in param_schema:
                    error_msg += f"\nExample: {param_schema['example']}"
                
                raise ValueError(error_msg)

        return True

    def validate_config(
        self,
        mode: str,
        sweep_params: dict[str, Any],
        fixed_params: dict[str, Any] = None,
    ) -> bool:
        """Validate entire configuration."""
        if mode not in self.schema:
            raise ValueError(
                f"Unknown mode: {mode}. Supported modes: {list(self.schema.keys())}"
            )

        # Validate sweep parameters
        for param_name, param_def in sweep_params.items():
            # We need to generate values from the sweep definition to validate them
            from .sweep_types import create_sweep

            try:
                sweep = create_sweep(param_def)
                values = sweep.generate_values()
                for value in values:
                    self.validate_parameter(mode, param_name, value)
            except Exception as e:
                raise ValueError(str(e))

        # Validate fixed parameters if provided
        if fixed_params:
            for param_name, value in fixed_params.items():
                if param_name in self.schema[mode]:
                    self.validate_parameter(mode, param_name, value)

        return True

    def get_supported_parameters(self, mode: str) -> list[str]:
        """Get list of supported parameters for a mode."""
        if mode not in self.schema:
            raise ValueError(f"Unknown mode: {mode}")
        return list(self.schema[mode].keys())

    def get_parameter_info(self, mode: str, param_name: str) -> dict[str, Any]:
        """Get parameter schema information."""
        if mode not in self.schema:
            raise ValueError(f"Unknown mode: {mode}")
        if param_name not in self.schema[mode]:
            raise ValueError(
                f"Parameter '{param_name}' not supported for mode '{mode}'"
            )
        return self.schema[mode][param_name]

    def get_supported_modes(self) -> list[str]:
        """Get list of supported modes."""
        return list(self.schema.keys())
