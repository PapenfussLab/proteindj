#!/usr/bin/env python3
"""Parameter converters for translating sweep parameters to nextflow configs."""

from abc import ABC, abstractmethod
from typing import Any, Optional


class ParameterConverter(ABC):
    """Abstract base class for parameter conversion."""

    @abstractmethod
    def to_profile_param(self, param_name: str, value: Any) -> dict[str, Any]:
        """Convert parameter to nextflow profile format."""
        pass

    @abstractmethod
    def to_command_arg(self, param_name: str, value: Any) -> Optional[str]:
        """Convert parameter to command line argument."""
        pass

    def format_value_for_name(self, value: Any) -> str:
        """Format value for use in profile/directory names."""
        if isinstance(value, float):
            # Use rstrip to remove trailing zeros, but preserve precision for small numbers
            formatted = f"{value:.10f}".rstrip('0').rstrip('.')
            return formatted.replace(".", "")
        elif isinstance(value, list):
            return "".join(str(v).replace(" ", "").replace(".", "").replace(",", "").replace("_", "").replace("/", "").replace("[", "").replace("]", "") for v in value)
        else:
            return str(value).replace(" ", "").replace(".", "").replace(",", "").replace("_", "").replace("/", "").replace("[", "").replace("]", "")


class DefaultParameterConverter(ParameterConverter):
    """Default converter that passes parameters through unchanged."""

    def to_profile_param(self, param_name: str, value: Any) -> dict[str, Any]:
        """Direct pass-through of parameter."""
        return {param_name: value}

    def to_command_arg(self, param_name: str, value: Any) -> Optional[str]:
        """Convert to command line argument."""
        if value is None:
            return None
        return f"--{param_name} '{value}'"

class ModelsParameterConverter(ParameterConverter):
    """Converter for model selection parameter."""

    def to_profile_param(self, param_name: str, value: str) -> dict[str, Any]:
        """Convert model name to rfd_diffusion_model parameter."""
        if value == "default":
            return {}  # Don't set parameter for default model
        return {"rfd_diffusion_model": value}

    def to_command_arg(self, param_name: str, value: str) -> Optional[str]:
        """No command line argument for models (handled by profile)."""
        return None


class NoiseScaleParameterConverter(ParameterConverter):
    """Converter for noise scale parameter."""

    def to_profile_param(self, param_name: str, value: float) -> dict[str, Any]:
        """Pass through noise scale."""
        return {"rfd_noise_scale": value}

    def to_command_arg(self, param_name: str, value: float) -> Optional[str]:
        """No command line argument for noise scale (handled by profile)."""
        return None


class HotspotsParameterConverter(ParameterConverter):
    """Converter for hotspots parameter."""

    def to_profile_param(self, param_name: str, value: Any) -> dict[str, Any]:
        """Convert hotspot list or string to formatted string."""
        # Handle empty values - convert to null
        if value in [None, "null", "", "[]"] or (isinstance(value, list) and len(value) == 0):
            return {"rfd_hotspots": None}
        elif isinstance(value, str):
            # If already a string, pass through directly
            return {"rfd_hotspots": value}
        elif isinstance(value, list):
            # If a list, join with commas
            hotspot_str = ",".join(value) if value else ""
            return {"rfd_hotspots": f"[{hotspot_str}]"}
        else:
            # Fallback for other types
            return {"rfd_hotspots": str(value)}

    def to_command_arg(self, param_name: str, value: Any) -> Optional[str]:
        """No command line argument for hotspots (handled by profile)."""
        return None

    def format_value_for_name(self, value: Any) -> str:
        """Format hotspots for directory naming."""
        if not value:
            return "nohotspots"
        return "hotspots" + "".join(h.replace(" ", "").replace(",", "").replace("_", "").replace("[", "").replace("]", "") for h in value)


class InputPdbParameterConverter(ParameterConverter):
    """Converter for input PDB parameter."""

    def to_profile_param(self, param_name: str, value: str) -> dict[str, Any]:
        """Pass through input PDB path."""
        return {"rfd_input_pdb": value}

    def to_command_arg(self, param_name: str, value: str) -> Optional[str]:
        """Convert to command line argument for partial diffusion mode."""
        # Only used in command line for partial diffusion mode
        return f"--input_pdb '{value}'"


class ContigsParameterConverter(ParameterConverter):
    """Converter for contigs parameter."""

    def to_profile_param(self, param_name: str, value: str) -> dict[str, Any]:
        """Pass through contigs."""
        return {"rfd_contigs": value}

    def to_command_arg(self, param_name: str, value: str) -> Optional[str]:
        """Convert to command line argument."""
        return f"--contigs '{value}'"


# Registry of parameter converters
PARAMETER_CONVERTERS = {
    "models": ModelsParameterConverter(),
    "model": ModelsParameterConverter(),  # Alias
    "rfd_noise_scale": NoiseScaleParameterConverter(),
    "noise_scale": NoiseScaleParameterConverter(),  # Alias
    "hotspots": HotspotsParameterConverter(),
    "rfd_hotspots": HotspotsParameterConverter(),  # Alias
    "input_pdb": InputPdbParameterConverter(),
    "rfd_input_pdb": InputPdbParameterConverter(),  # Alias
    "contigs": ContigsParameterConverter(),
    "rfd_contigs": ContigsParameterConverter(),  # Alias
}

# Default converter for fallthrough
DEFAULT_CONVERTER = DefaultParameterConverter()


def get_converter(param_name: str) -> ParameterConverter:
    """Get the appropriate converter for a parameter."""
    return PARAMETER_CONVERTERS.get(param_name, DEFAULT_CONVERTER)
