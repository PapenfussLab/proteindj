"""Tests for sweep type definitions - FIXED VERSION."""

import pytest

from bindsweeper.sweep_types import ListSweep, RangeSweep, create_sweep


class TestListSweep:
    """Test ListSweep functionality."""

    def test_create_list_sweep(self):
        """Test creating a list sweep."""
        values = ["default", "beta", "complex_base"]
        sweep = ListSweep(values=values)
        assert sweep.values == values

    def test_generate_values(self):
        """Test generating values from list sweep."""
        values = [1, 2, 3]
        sweep = ListSweep(values=values)
        assert sweep.generate_values() == values

    def test_to_dict(self):
        """Test converting list sweep to dictionary."""
        values = ["a", "b", "c"]
        sweep = ListSweep(values=values)
        expected = {"type": "list", "values": values}
        assert sweep.to_dict() == expected

    def test_from_dict(self):
        """Test creating list sweep from dictionary."""
        data = {"values": ["x", "y", "z"]}
        sweep = ListSweep.from_dict(data)
        assert sweep.values == ["x", "y", "z"]

    def test_from_dict_direct_values(self):
        """Test creating list sweep from direct values - FIXED."""
        # The actual implementation expects data to be the values directly
        # when passed as a list to from_dict
        data = ["a", "b", "c"]
        sweep = ListSweep.from_dict(data)
        assert sweep.values == data


class TestRangeSweep:
    """Test RangeSweep functionality."""

    def test_create_range_sweep(self):
        """Test creating a range sweep."""
        sweep = RangeSweep(min=0.0, max=1.0, step=0.5)
        assert sweep.min == 0.0
        assert sweep.max == 1.0
        assert sweep.step == 0.5

    def test_invalid_step(self):
        """Test invalid step parameter."""
        with pytest.raises(ValueError, match="Step must be positive"):
            RangeSweep(min=0.0, max=1.0, step=0.0)

        with pytest.raises(ValueError, match="Step must be positive"):
            RangeSweep(min=0.0, max=1.0, step=-0.1)

    def test_invalid_range(self):
        """Test invalid range parameters."""
        with pytest.raises(ValueError, match="Min must be less than or equal to max"):
            RangeSweep(min=1.0, max=0.0, step=0.1)

    def test_generate_values(self):
        """Test generating values from range sweep."""
        sweep = RangeSweep(min=0.0, max=1.0, step=0.5)
        values = sweep.generate_values()
        expected = [0.0, 0.5, 1.0]
        assert values == expected

    def test_generate_values_precise(self):
        """Test range generation with precise boundaries."""
        sweep = RangeSweep(min=0.1, max=0.3, step=0.1)
        values = sweep.generate_values()
        assert len(values) == 3
        assert abs(values[0] - 0.1) < 1e-10
        assert abs(values[1] - 0.2) < 1e-10
        assert abs(values[2] - 0.3) < 1e-10

    def test_single_value_range(self):
        """Test range that generates single value."""
        sweep = RangeSweep(min=1.0, max=1.0, step=0.1)
        values = sweep.generate_values()
        assert values == [1.0]

    def test_to_dict(self):
        """Test converting range sweep to dictionary."""
        sweep = RangeSweep(min=0.0, max=2.0, step=1.0)
        expected = {"type": "range", "min": 0.0, "max": 2.0, "step": 1.0}
        assert sweep.to_dict() == expected

    def test_from_dict(self):
        """Test creating range sweep from dictionary."""
        data = {"min": 0.0, "max": 1.0, "step": 0.25}
        sweep = RangeSweep.from_dict(data)
        assert sweep.min == 0.0
        assert sweep.max == 1.0
        assert sweep.step == 0.25


class TestCreateSweep:
    """Test sweep factory function."""

    def test_create_from_list(self):
        """Test creating sweep from plain list."""
        data = ["a", "b", "c"]
        sweep = create_sweep(data)
        assert isinstance(sweep, ListSweep)
        assert sweep.values == data

    def test_create_from_dict_with_type_list(self):
        """Test creating list sweep from dict with type."""
        data = {"type": "list", "values": [1, 2, 3]}
        sweep = create_sweep(data)
        assert isinstance(sweep, ListSweep)
        assert sweep.values == [1, 2, 3]

    def test_create_from_dict_with_type_range(self):
        """Test creating range sweep from dict with type."""
        data = {"type": "range", "min": 0, "max": 10, "step": 2}
        sweep = create_sweep(data)
        assert isinstance(sweep, RangeSweep)
        assert sweep.min == 0
        assert sweep.max == 10
        assert sweep.step == 2

    def test_create_from_dict_infer_range(self):
        """Test creating range sweep by inferring from keys."""
        data = {"min": 0.0, "max": 1.0, "step": 0.1}
        sweep = create_sweep(data)
        assert isinstance(sweep, RangeSweep)

    def test_create_from_dict_infer_list(self):
        """Test creating list sweep by inferring from keys."""
        data = {"values": ["x", "y", "z"]}
        sweep = create_sweep(data)
        assert isinstance(sweep, ListSweep)
        assert sweep.values == ["x", "y", "z"]

    def test_create_invalid_data(self):
        """Test creating sweep from invalid data."""
        with pytest.raises(ValueError, match="Cannot create sweep from data"):
            create_sweep({"invalid": "data"})

        with pytest.raises(ValueError, match="Cannot create sweep from data"):
            create_sweep("invalid_string")

        with pytest.raises(ValueError, match="Cannot create sweep from data"):
            create_sweep(42)
