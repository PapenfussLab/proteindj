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


def create_sweep(data: Union[dict, list]) -> SweepType:
    """Factory function to create appropriate sweep type from data."""
    if isinstance(data, list):
        return ListSweep(values=data)

    if isinstance(data, dict):
        if "type" in data:
            if data["type"] == "range":
                return RangeSweep.from_dict(data)
            elif data["type"] == "list":
                return ListSweep.from_dict(data)

        # Check if it's a range by presence of min/max/step
        if all(k in data for k in ["min", "max", "step"]):
            return RangeSweep.from_dict(data)

        # Check if it has values key
        if "values" in data:
            return ListSweep.from_dict(data)

    raise ValueError(f"Cannot create sweep from data: {data}")
