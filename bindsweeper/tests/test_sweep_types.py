"""Tests for sweep type definitions."""

import pytest

from bindsweeper.sweep_types import ListSweep, PairedSweep, RangeSweep, create_sweep


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
        """Test creating list sweep from direct values."""
        # The from_dict method expects data.get("values", data)
        # When data is a list, it should use the list as values
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

    def test_create_paired_from_dict_with_type(self):
        """Test creating paired sweep from dict with explicit type."""
        data = {
            "type": "paired",
            "values": ["a.pdb", "b.pdb"],
            "paired_with": {"msa_path": ["a.a3m", "b.a3m"]},
        }
        sweep = create_sweep(data)
        assert isinstance(sweep, PairedSweep)

    def test_create_paired_from_dict_inferred(self):
        """Test creating paired sweep inferred from paired_with key."""
        data = {
            "values": ["a.pdb", "b.pdb"],
            "paired_with": {"msa_path": ["a.a3m", "b.a3m"]},
        }
        sweep = create_sweep(data)
        assert isinstance(sweep, PairedSweep)

    def test_paired_with_takes_precedence_over_type_list(self):
        """Test that paired_with takes precedence even if type='list' is specified."""
        data = {
            "type": "list",
            "values": ["a.pdb", "b.pdb"],
            "paired_with": {"msa_path": ["a.a3m", "b.a3m"]},
        }
        sweep = create_sweep(data)
        # Should create PairedSweep, not ListSweep, to avoid silently dropping pairing
        assert isinstance(sweep, PairedSweep)
        assert sweep.generate_values() == ["a.pdb", "b.pdb"]
        assert sweep.get_paired_value("msa_path", 0) == "a.a3m"


class TestPairedSweep:
    """Test PairedSweep functionality."""

    def test_create_paired_sweep(self):
        """Test creating a paired sweep."""
        sweep = PairedSweep(
            values=["target1.pdb", "target2.pdb", "target3.pdb"],
            paired_params={
                "msa_path": ["msa1.a3m", "msa2.a3m", "msa3.a3m"],
            },
        )
        assert sweep.values == ["target1.pdb", "target2.pdb", "target3.pdb"]
        assert len(sweep.paired_params) == 1

    def test_generate_values(self):
        """Test generating values returns primary values only."""
        sweep = PairedSweep(
            values=["a.pdb", "b.pdb"],
            paired_params={"msa": ["a.a3m", "b.a3m"]},
        )
        assert sweep.generate_values() == ["a.pdb", "b.pdb"]

    def test_get_paired_value(self):
        """Test retrieving paired values by index."""
        sweep = PairedSweep(
            values=["t1.pdb", "t2.pdb", "t3.pdb"],
            paired_params={
                "msa": ["m1.a3m", "m2.a3m", "m3.a3m"],
                "chain": ["A", "B", "C"],
            },
        )
        assert sweep.get_paired_value("msa", 0) == "m1.a3m"
        assert sweep.get_paired_value("msa", 2) == "m3.a3m"
        assert sweep.get_paired_value("chain", 1) == "B"

    def test_length_mismatch_raises(self):
        """Test that mismatched lengths raise ValueError."""
        with pytest.raises(ValueError, match="must have the same length"):
            PairedSweep(
                values=["a.pdb", "b.pdb", "c.pdb"],
                paired_params={"msa": ["a.a3m", "b.a3m"]},  # 2 != 3
            )

    def test_empty_paired_params_raises(self):
        """Test that empty paired_with raises ValueError."""
        with pytest.raises(ValueError, match="at least one parameter"):
            PairedSweep(values=["a.pdb", "b.pdb"], paired_params={})

    def test_scalar_paired_value_raises(self):
        """Test that scalar (non-list) paired value raises ValueError."""
        with pytest.raises(ValueError, match="must be a list"):
            PairedSweep(
                values=["a.pdb", "b.pdb"],
                paired_params={"msa": "single_value.a3m"},
            )

    def test_to_dict(self):
        """Test converting paired sweep to dictionary."""
        sweep = PairedSweep(
            values=["a.pdb", "b.pdb"],
            paired_params={"msa": ["a.a3m", "b.a3m"]},
        )
        expected = {
            "type": "paired",
            "values": ["a.pdb", "b.pdb"],
            "paired_with": {"msa": ["a.a3m", "b.a3m"]},
        }
        assert sweep.to_dict() == expected

    def test_from_dict(self):
        """Test creating paired sweep from dictionary."""
        data = {
            "values": ["x.pdb", "y.pdb"],
            "paired_with": {"chain": ["A", "B"]},
        }
        sweep = PairedSweep.from_dict(data)
        assert sweep.values == ["x.pdb", "y.pdb"]
        assert sweep.paired_params == {"chain": ["A", "B"]}

    def test_roundtrip(self):
        """Test to_dict/from_dict roundtrip."""
        original = PairedSweep(
            values=["a.pdb", "b.pdb", "c.pdb"],
            paired_params={
                "msa": ["a.a3m", "b.a3m", "c.a3m"],
                "chain": ["A", "B", "C"],
            },
        )
        restored = PairedSweep.from_dict(original.to_dict())
        assert restored.values == original.values
        assert restored.paired_params == original.paired_params

    def test_multiple_paired_params(self):
        """Test pairing with multiple secondary parameters."""
        sweep = PairedSweep(
            values=["t1.pdb", "t2.pdb"],
            paired_params={
                "msa_path": ["m1.a3m", "m2.a3m"],
                "target_chain": ["A", "B"],
            },
        )
        assert sweep.get_paired_value("msa_path", 0) == "m1.a3m"
        assert sweep.get_paired_value("target_chain", 1) == "B"

    def test_single_value_paired(self):
        """Test paired sweep with single value still works."""
        sweep = PairedSweep(
            values=["only.pdb"],
            paired_params={"msa": ["only.a3m"]},
        )
        assert sweep.generate_values() == ["only.pdb"]
        assert sweep.get_paired_value("msa", 0) == "only.a3m"
