#!/usr/bin/env python3
"""Sweep type definitions for parameter sweeping."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Union


class SweepType(ABC):
    """Abstract base class for parameter sweep types."""

    @abstractmethod
    def generate_values(self) -> list[Any]:
        """Generate all values for this sweep."""
        pass

    @abstractmethod
    def to_dict(self) -> dict:
        """Convert sweep definition to dictionary."""
        pass

    @classmethod
    @abstractmethod
    def from_dict(cls, data) -> "SweepType":
        """Create sweep from data."""
        pass


@dataclass
class ListSweep(SweepType):
    """Sweep through a list of discrete values."""

    values: list[Any]

    def generate_values(self) -> list[Any]:
        """Return the list of values."""
        return self.values

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {"type": "list", "values": self.values}

    @classmethod
    def from_dict(cls, data) -> "ListSweep":
        """Create from dictionary or list."""
        if isinstance(data, list):
            return cls(values=data)
        elif isinstance(data, dict):
            return cls(values=data.get("values", data))
        else:
            raise ValueError(f"Cannot create ListSweep from {type(data)}")


@dataclass
class RangeSweep(SweepType):
    """Sweep through a range of numeric values."""

    min: float
    max: float
    step: float

    def __post_init__(self):
        """Validate range parameters."""
        if self.step <= 0:
            raise ValueError("Step must be positive")
        if self.min > self.max:
            raise ValueError("Min must be less than or equal to max")

    def generate_values(self) -> list[float]:
        """Generate range values from min to max with step increments."""
        values = []
        value = self.min
        while value <= self.max + 1e-10:  # Small epsilon for floating point
            values.append(value)
            value += self.step
        return values

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {"type": "range", "min": self.min, "max": self.max, "step": self.step}

    @classmethod
    def from_dict(cls, data) -> "RangeSweep":
        """Create from dictionary."""
        if not isinstance(data, dict):
            raise ValueError(f"RangeSweep requires dictionary data, got {type(data)}")
        return cls(
            min=float(data["min"]), max=float(data["max"]), step=float(data["step"])
        )


@dataclass
class PairedSweep(SweepType):
    """Sweep through values paired with other parameters (zipped, not Cartesian product)."""

    values: list[Any]
    paired_params: dict[str, list[Any]]

    def __post_init__(self):
        """Validate paired parameter structure and lengths."""
        if not self.paired_params:
            raise ValueError(
                "PairedSweep requires at least one parameter in 'paired_with'. "
                "If no pairing is needed, use a regular list sweep instead."
            )
        primary_len = len(self.values)
        for param_name, param_values in self.paired_params.items():
            if not isinstance(param_values, list):
                raise ValueError(
                    f"Paired parameter '{param_name}' must be a list of values, "
                    f"got {type(param_values).__name__}. "
                    f"Wrap the value in a list, e.g.: [{param_values}]"
                )
            if len(param_values) != primary_len:
                raise ValueError(
                    f"Paired parameter '{param_name}' has {len(param_values)} values, "
                    f"but primary parameter has {primary_len} values. "
                    f"All paired parameters must have the same length."
                )

    def generate_values(self) -> list[Any]:
        """Return the list of values for the primary parameter."""
        return self.values

    def get_paired_value(self, param_name: str, index: int) -> Any:
        """Get the paired value for a parameter at a given index."""
        return self.paired_params[param_name][index]

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return {
            "type": "paired",
            "values": self.values,
            "paired_with": self.paired_params,
        }

    @classmethod
    def from_dict(cls, data) -> "PairedSweep":
        """Create from dictionary."""
        if not isinstance(data, dict):
            raise ValueError(f"PairedSweep requires dictionary data, got {type(data)}")
        
        values = data.get("values", [])
        paired_with = data.get("paired_with", {})
        
        return cls(values=values, paired_params=paired_with)


def create_sweep(data: Union[dict, list]) -> SweepType:
    """Factory function to create appropriate sweep type from data."""
    if isinstance(data, list):
        return ListSweep(values=data)

    if isinstance(data, dict):
        # Check for paired_with FIRST — even if 'type' is also specified,
        # paired_with takes precedence to avoid silently dropping pairings.
        if "paired_with" in data:
            return PairedSweep.from_dict(data)

        if "type" in data:
            if data["type"] == "range":
                return RangeSweep.from_dict(data)
            elif data["type"] == "list":
                return ListSweep.from_dict(data)
            elif data["type"] == "paired":
                return PairedSweep.from_dict(data)

        # Check if it's a range by presence of min/max/step
        if all(k in data for k in ["min", "max", "step"]):
            return RangeSweep.from_dict(data)

        # Check if it has values key
        if "values" in data:
            return ListSweep.from_dict(data)

    raise ValueError(f"Cannot create sweep from data: {data}")
