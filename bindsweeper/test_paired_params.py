#!/usr/bin/env python3
"""Test script for paired parameter functionality."""

import sys
from pathlib import Path

# Add bindsweeper to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from bindsweeper.sweep_types import PairedSweep, ListSweep, create_sweep


def test_paired_sweep_creation():
    """Test creating a PairedSweep."""
    print("Test 1: Creating PairedSweep...")
    
    data = {
        "values": ["file1.pdb", "file2.pdb", "file3.pdb"],
        "paired_with": {
            "msa_path": ["msa1.a3m", "msa2.a3m", "msa3.a3m"]
        }
    }
    
    sweep = create_sweep(data)
    assert isinstance(sweep, PairedSweep)
    assert sweep.generate_values() == ["file1.pdb", "file2.pdb", "file3.pdb"]
    assert sweep.get_paired_value("msa_path", 0) == "msa1.a3m"
    assert sweep.get_paired_value("msa_path", 1) == "msa2.a3m"
    assert sweep.get_paired_value("msa_path", 2) == "msa3.a3m"
    
    print("✓ PairedSweep creation works correctly")


def test_paired_sweep_validation():
    """Test that PairedSweep validates length mismatches."""
    print("\nTest 2: Validating length mismatch detection...")
    
    data = {
        "values": ["file1.pdb", "file2.pdb", "file3.pdb"],
        "paired_with": {
            "msa_path": ["msa1.a3m", "msa2.a3m"]  # Wrong length!
        }
    }
    
    try:
        sweep = create_sweep(data)
        print("✗ Should have raised ValueError for length mismatch")
        sys.exit(1)
    except ValueError as e:
        if "must have the same length" in str(e):
            print(f"✓ Length mismatch detected correctly: {e}")
        else:
            print(f"✗ Wrong error message: {e}")
            sys.exit(1)


def test_yaml_config_parsing():
    """Test parsing a YAML config with paired parameters."""
    print("\nTest 3: Parsing YAML config with paired parameters...")
    
    import tempfile
    import yaml
    from bindsweeper.sweep_config import SweepConfig
    
    config_data = {
        "mode": "bindcraft_denovo",
        "profile": "milton",
        "fixed_params": {
            "skip_fold_seq": True,
            "pred_method": "boltz"
        },
        "sweep_params": {
            "uncropped_target_pdb": {
                "values": [
                    "input/protein1.pdb",
                    "input/protein2.pdb",
                    "input/protein3.pdb"
                ],
                "paired_with": {
                    "boltz_msa_path": [
                        "input/msas/protein1.a3m",
                        "input/msas/protein2.a3m",
                        "input/msas/protein3.a3m"
                    ]
                }
            }
        }
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(config_data, f)
        temp_path = f.name
    
    try:
        config = SweepConfig.from_yaml(temp_path)
        assert config.mode == "bindcraft_denovo"
        assert "uncropped_target_pdb" in config.sweep_params
        
        sweep = config.sweep_params["uncropped_target_pdb"]
        assert isinstance(sweep, PairedSweep)
        assert len(sweep.generate_values()) == 3
        assert "boltz_msa_path" in sweep.paired_params
        
        print("✓ YAML config parsed correctly")
    finally:
        Path(temp_path).unlink()


def test_combination_generation():
    """Test that combinations are generated correctly with paired parameters."""
    print("\nTest 4: Generating combinations with paired parameters...")
    
    import tempfile
    import yaml
    from bindsweeper.sweep_config import SweepConfig
    from bindsweeper.sweep_engine import SweepEngine
    
    config_data = {
        "mode": "bindcraft_denovo",
        "profile": "milton",
        "fixed_params": {
            "skip_fold_seq": True
        },
        "sweep_params": {
            "uncropped_target_pdb": {
                "values": [
                    "protein1.pdb",
                    "protein2.pdb",
                    "protein3.pdb"
                ],
                "paired_with": {
                    "boltz_msa_path": [
                        "protein1.a3m",
                        "protein2.a3m",
                        "protein3.a3m"
                    ]
                }
            }
        }
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(config_data, f)
        temp_path = f.name
    
    try:
        config = SweepConfig.from_yaml(temp_path)
        engine = SweepEngine(config, "./test_output", "./nextflow.config")
        
        combinations = engine.generate_combinations()
        
        # Should have 3 combinations (zipped, not Cartesian)
        assert len(combinations) == 3, f"Expected 3 combinations, got {len(combinations)}"
        
        # Verify pairings
        for i, combo in enumerate(combinations):
            pdb_file = combo.swept_params["uncropped_target_pdb"]
            msa_file = combo.swept_params["boltz_msa_path"]
            
            # Check that they're paired correctly
            assert pdb_file.replace(".pdb", "") == msa_file.replace(".a3m", "")
            
        print(f"✓ Generated {len(combinations)} combinations correctly:")
        for combo in combinations:
            print(f"  - PDB: {combo.swept_params['uncropped_target_pdb']}, "
                  f"MSA: {combo.swept_params['boltz_msa_path']}")
            
    finally:
        Path(temp_path).unlink()


def test_paired_with_unpaired_cartesian():
    """Test combining paired parameters with unpaired parameters (Cartesian)."""
    print("\nTest 5: Combining paired + unpaired parameters...")
    
    import tempfile
    import yaml
    from bindsweeper.sweep_config import SweepConfig
    from bindsweeper.sweep_engine import SweepEngine
    
    config_data = {
        "mode": "bindcraft_denovo",
        "profile": "milton",
        "fixed_params": {},
        "sweep_params": {
            "uncropped_target_pdb": {
                "values": ["protein1.pdb", "protein2.pdb"],
                "paired_with": {
                    "boltz_msa_path": ["protein1.a3m", "protein2.a3m"]
                }
            },
            "rfd_noise_scale": {
                "values": [0.0, 0.1]
            }
        }
    }
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        yaml.dump(config_data, f)
        temp_path = f.name
    
    try:
        config = SweepConfig.from_yaml(temp_path)
        engine = SweepEngine(config, "./test_output", "./nextflow.config")
        
        combinations = engine.generate_combinations()
        
        # Should have 2 (paired) × 2 (unpaired) = 4 combinations
        assert len(combinations) == 4, f"Expected 4 combinations, got {len(combinations)}"
        
        print(f"✓ Generated {len(combinations)} combinations correctly:")
        for combo in combinations:
            print(f"  - PDB: {combo.swept_params['uncropped_target_pdb']}, "
                  f"MSA: {combo.swept_params['boltz_msa_path']}, "
                  f"Noise: {combo.swept_params['rfd_noise_scale']}")
            
    finally:
        Path(temp_path).unlink()


if __name__ == "__main__":
    print("="*60)
    print("Testing Paired Parameter Functionality")
    print("="*60)
    
    try:
        test_paired_sweep_creation()
        test_paired_sweep_validation()
        test_yaml_config_parsing()
        test_combination_generation()
        test_paired_with_unpaired_cartesian()
        
        print("\n" + "="*60)
        print("✓ All tests passed!")
        print("="*60)
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
