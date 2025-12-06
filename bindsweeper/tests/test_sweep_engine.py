"""Tests for sweep execution engine."""

import datetime
import subprocess
from unittest.mock import Mock, patch

import pytest

from bindsweeper.sweep_config import SweepConfig
from bindsweeper.sweep_engine import CommandResult, SweepCombination, SweepEngine


class TestCommandResult:
    """Test CommandResult dataclass."""

    def test_command_result_creation(self):
        """Test creating a command result."""
        start = datetime.datetime.now()
        end = start + datetime.timedelta(seconds=120)

        result = CommandResult(
            success=True,
            start_time=start,
            end_time=end,
            duration=120.0,
            return_code=0,
            output_dir="/test/output",
        )

        assert result.success is True
        assert result.duration == 120.0
        assert result.return_code == 0
        assert result.output_dir == "/test/output"

    def test_duration_str_formatting(self):
        """Test duration string formatting."""
        start = datetime.datetime.now()
        end = start + datetime.timedelta(seconds=3665)  # 1 hour, 1 minute, 5 seconds

        result = CommandResult(
            success=True, start_time=start, end_time=end, duration=3665.0, return_code=0
        )

        assert result.duration_str == "1:01:05"


class TestSweepCombination:
    """Test SweepCombination dataclass."""

    def test_sweep_combination_creation(self):
        """Test creating a sweep combination."""
        combo = SweepCombination(
            mode="binder_denovo",
            all_params={"num_designs": 8, "rfd_noise_scale": 0.5},
            swept_params={"rfd_noise_scale": 0.5},
            profile_name="test_profile",
            output_dir="/test/output",
            command="nextflow run main.nf",
        )

        assert combo.mode == "binder_denovo"
        assert combo.swept_params == {"rfd_noise_scale": 0.5}
        assert combo.profile_name == "test_profile"


class TestSweepEngine:
    """Test SweepEngine functionality."""

    @pytest.fixture
    def mock_config(self):
        """Mock sweep configuration."""
        config = Mock(spec=SweepConfig)
        config.mode = "binder_denovo"
        config.fixed_params = {"design_length": "50-100", "num_designs": 4}
        config.profile = "milton"
        config.pipeline_path = "main.nf"

        # Mock sweep parameters
        noise_sweep = Mock()
        noise_sweep.generate_values.return_value = [0.0, 0.5, 1.0]

        model_sweep = Mock()
        model_sweep.generate_values.return_value = ["default", "beta"]

        config.sweep_params = {"rfd_noise_scale": noise_sweep, "models": model_sweep}

        return config

    @pytest.fixture
    def sweep_engine(self, mock_config, temp_dir, config_files):
        """Create a sweep engine with mocked dependencies."""
        engine = SweepEngine(
            config=mock_config,
            base_output_dir=temp_dir,
            nextflow_config_path=config_files["nextflow_config"],
        )
        return engine

    def test_engine_initialization(self, sweep_engine, mock_config, temp_dir):
        """Test sweep engine initialization."""
        assert sweep_engine.config == mock_config
        assert sweep_engine.base_output_dir == temp_dir
        assert sweep_engine.resume == False
        assert sweep_engine.parallel == False
        assert sweep_engine.max_parallel == 4

    def test_generate_combinations(self, sweep_engine):
        """Test generating parameter combinations."""
        combinations = sweep_engine.generate_combinations()

        # Should generate 3 * 2 = 6 combinations
        assert len(combinations) == 6

        # Check first combination
        combo = combinations[0]
        assert combo.mode == "binder_denovo"
        assert "rfd_noise_scale" in combo.swept_params
        assert "models" in combo.swept_params
        assert combo.swept_params["rfd_noise_scale"] == 0.0
        assert combo.swept_params["models"] == "default"

    def test_generate_combinations_no_sweep_params(self, temp_dir, config_files):
        """Test generating combinations with no sweep parameters."""
        config = Mock(spec=SweepConfig)
        config.mode = "denovo"
        config.fixed_params = {"num_designs": 4}
        config.sweep_params = {}
        config.profile = "milton"
        config.pipeline_path = "main.nf"

        engine = SweepEngine(config, temp_dir, config_files["nextflow_config"])

        combinations = engine.generate_combinations()
        assert len(combinations) == 0

    def test_generate_output_dir(self, sweep_engine):
        """Test output directory generation."""
        swept_params = {"rfd_noise_scale": 0.5, "models": "beta"}

        with patch("bindsweeper.sweep_engine.get_converter") as mock_converter:
            # Mock parameter converters
            mock_converter.return_value.format_value_for_name.side_effect = (
                lambda x: str(x)
            )

            output_dir = sweep_engine._generate_output_dir(swept_params)

            # Should contain formatted parameter values
            assert "models_beta" in output_dir or "beta" in output_dir
            assert "noisescale_0.5" in output_dir or "0.5" in output_dir

    def test_generate_command(self, sweep_engine):
        """Test command generation."""
        command = sweep_engine._generate_command("test_profile", "/test/output")

        assert "nextflow -c bindsweeper.config run main.nf" in command
        assert "-profile milton,test_profile" in command
        assert "--out_dir '/test/output'" in command
        assert "--zip_pdbs false" in command

    @patch("bindsweeper.sweep_engine.generate_profile_content")
    def test_generate_profiles(self, mock_generate_profile, sweep_engine):
        """Test profile generation."""
        mock_generate_profile.return_value = "mock_profile_content"

        combinations = [
            SweepCombination(
                mode="test_mode",
                all_params={"param1": "value1"},
                swept_params={"param1": "value1"},
                profile_name="test_profile",
                output_dir="/test",
                command="test_command",
            )
        ]

        profiles = sweep_engine.generate_profiles(combinations)

        assert len(profiles) == 1
        assert profiles[0] == "mock_profile_content"
        mock_generate_profile.assert_called_once()

    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_execute_combination_success(
        self, mock_makedirs, mock_subprocess, sweep_engine
    ):
        """Test successful combination execution."""
        # Mock successful subprocess
        mock_process = Mock()
        mock_process.returncode = 0
        mock_subprocess.return_value = mock_process

        combination = SweepCombination(
            mode="test_mode",
            all_params={},
            swept_params={"param": "value"},
            profile_name="test_profile",
            output_dir="/test/output",
            command="test command",
        )

        result = sweep_engine.execute_combination(combination)

        assert result.success is True
        assert result.return_code == 0
        assert result.output_dir == "/test/output"
        assert result.swept_params == {"param": "value"}
        mock_makedirs.assert_called_once_with("/test/output", exist_ok=True)

    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_execute_combination_failure(
        self, mock_makedirs, mock_subprocess, sweep_engine
    ):
        """Test failed combination execution."""
        # Mock failed subprocess
        error = subprocess.CalledProcessError(1, "test command")
        error.stderr = "Test error message"
        mock_subprocess.side_effect = error

        combination = SweepCombination(
            mode="test_mode",
            all_params={},
            swept_params={"param": "value"},
            profile_name="test_profile",
            output_dir="/test/output",
            command="test command",
        )

        result = sweep_engine.execute_combination(combination)

        assert result.success is False
        assert result.return_code == 1
        assert "Test error message" in result.error_message

    def test_execute_sweep_dry_run(self, sweep_engine):
        """Test executing sweep in dry run mode."""
        combinations = [
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": "value"},
                profile_name="test_profile",
                output_dir="/test/output",
                command="test command",
            )
        ]

        results = sweep_engine.execute_sweep(combinations, dry_run=True)

        # Dry run should return empty results
        assert len(results) == 0

    @patch.object(SweepEngine, "execute_combination")
    def test_execute_sweep_continue_on_error(self, mock_execute, sweep_engine):
        """Test executing sweep with continue on error."""
        # Mock one successful and one failed execution
        success_result = CommandResult(
            success=True,
            start_time=datetime.datetime.now(),
            end_time=datetime.datetime.now(),
            duration=10.0,
            return_code=0,
        )

        failure_result = CommandResult(
            success=False,
            start_time=datetime.datetime.now(),
            end_time=datetime.datetime.now(),
            duration=5.0,
            return_code=1,
            error_message="Test error",
        )

        mock_execute.side_effect = [success_result, failure_result]

        # Create mock combinations with required attributes
        mock_combo1 = Mock()
        mock_combo1.mode = "test_mode"
        mock_combo1.swept_params = {"param1": "value1"}
        mock_combo1.output_dir = "/test/output1"
        mock_combo1.profile_name = "test_profile1"
        mock_combo1.command = "test command 1"

        mock_combo2 = Mock()
        mock_combo2.mode = "test_mode"
        mock_combo2.swept_params = {"param2": "value2"}
        mock_combo2.output_dir = "/test/output2"
        mock_combo2.profile_name = "test_profile2"
        mock_combo2.command = "test command 2"

        combinations = [mock_combo1, mock_combo2]

        results = sweep_engine.execute_sweep(
            combinations, dry_run=False, continue_on_error=True
        )

        assert len(results) == 2
        assert results[0].success is True
        assert results[1].success is False

    @patch.object(SweepEngine, "execute_combination")
    def test_execute_sweep_stop_on_error(self, mock_execute, sweep_engine):
        """Test executing sweep that stops on error."""
        failure_result = CommandResult(
            success=False,
            start_time=datetime.datetime.now(),
            end_time=datetime.datetime.now(),
            duration=5.0,
            return_code=1,
            error_message="Test error",
        )

        mock_execute.return_value = failure_result

        # Create mock combination with required attributes
        mock_combo = Mock()
        mock_combo.mode = "test_mode"
        mock_combo.swept_params = {"param": "value"}
        mock_combo.output_dir = "/test/output"
        mock_combo.profile_name = "test_profile"
        mock_combo.command = "test command"

        combinations = [mock_combo]

        with pytest.raises(RuntimeError, match="Sweep execution failed"):
            sweep_engine.execute_sweep(
                combinations, dry_run=False, continue_on_error=False
            )

    @patch("subprocess.run")
    @patch("os.makedirs")
    def test_execute_combination_with_isolated_cache(
        self, mock_makedirs, mock_subprocess, sweep_engine
    ):
        """Test executing combination with isolated cache directory."""
        mock_subprocess.return_value = Mock(returncode=0)

        combination = SweepCombination(
            mode="test_mode",
            all_params={},
            swept_params={"param": "value"},
            profile_name="test_profile",
            output_dir="/test/output",
            command="test command",
        )

        result = sweep_engine.execute_combination(combination, use_isolated_cache=True)

        assert result.success is True
        # Should create both output dir and cache dir
        assert mock_makedirs.call_count == 2
        
        # Check that subprocess was called with env containing NXF_CACHE_DIR
        call_args = mock_subprocess.call_args
        assert "env" in call_args.kwargs
        assert "NXF_CACHE_DIR" in call_args.kwargs["env"]
        assert call_args.kwargs["env"]["NXF_CACHE_DIR"] == "/test/output/.nextflow_cache"

    @patch.object(SweepEngine, "_execute_parallel")
    def test_execute_sweep_parallel_mode(self, mock_parallel, sweep_engine):
        """Test that parallel mode triggers parallel execution."""
        sweep_engine.parallel = True
        
        combinations = [
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": "value"},
                profile_name="test_profile",
                output_dir="/test/output",
                command="test command",
            )
        ]

        mock_parallel.return_value = []

        sweep_engine.execute_sweep(combinations, dry_run=False, parallel=True)

        # Verify parallel execution was called
        mock_parallel.assert_called_once_with(combinations, False)

    def test_execute_parallel_success(self, sweep_engine, temp_dir):
        """Test parallel execution with successful combinations using real commands."""
        # Use simple shell commands that will actually succeed
        combinations = [
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": f"value{i}"},
                profile_name=f"test_profile{i}",
                output_dir=f"{temp_dir}/output{i}",
                command=f"echo 'test {i}' > {temp_dir}/test{i}.txt",
            )
            for i in range(2)  # Use fewer for faster tests
        ]

        results = sweep_engine._execute_parallel(combinations, continue_on_error=False)

        assert len(results) == 2
        assert all(r.success for r in results)
        # Verify output files were created
        import os
        assert os.path.exists(f"{temp_dir}/test0.txt")
        assert os.path.exists(f"{temp_dir}/test1.txt")

    def test_execute_parallel_with_failure(self, sweep_engine, temp_dir):
        """Test parallel execution with failure and continue_on_error."""
        # Mix of successful and failing commands
        combinations = [
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": "success"},
                profile_name="test_profile_success",
                output_dir=f"{temp_dir}/output_success",
                command=f"echo 'success' > {temp_dir}/success.txt",
            ),
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": "failure"},
                profile_name="test_profile_failure",
                output_dir=f"{temp_dir}/output_failure",
                command="exit 1",  # This will fail
            ),
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": "success2"},
                profile_name="test_profile_success2",
                output_dir=f"{temp_dir}/output_success2",
                command=f"echo 'success2' > {temp_dir}/success2.txt",
            ),
        ]

        results = sweep_engine._execute_parallel(combinations, continue_on_error=True)

        assert len(results) == 3
        # Two should succeed, one should fail
        successful = [r for r in results if r.success]
        failed = [r for r in results if not r.success]
        assert len(successful) == 2
        assert len(failed) == 1

    def test_execute_parallel_stop_on_failure(self, sweep_engine, temp_dir):
        """Test parallel execution stops on failure when continue_on_error=False."""
        # Create combinations where at least one will fail
        combinations = [
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": "failure"},
                profile_name="test_profile_failure",
                output_dir=f"{temp_dir}/output_failure",
                command="exit 1",  # This will fail
            ),
            SweepCombination(
                mode="test_mode",
                all_params={},
                swept_params={"param": "other"},
                profile_name="test_profile_other",
                output_dir=f"{temp_dir}/output_other",
                command=f"echo 'test'",
            ),
        ]

        # Should raise RuntimeError when a combination fails
        with pytest.raises(RuntimeError, match="Combination failed"):
            sweep_engine._execute_parallel(combinations, continue_on_error=False)
