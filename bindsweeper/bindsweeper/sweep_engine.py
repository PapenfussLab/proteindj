#!/usr/bin/env python3
"""Core sweep execution engine."""

import datetime
import logging
import os
import subprocess
from dataclasses import dataclass
from itertools import product
from typing import Any

from .parameter_converters import get_converter
from .profile_generator import generate_profile_content, generate_profile_name
from .sweep_config import SweepConfig

logger = logging.getLogger(__name__)


@dataclass
class CommandResult:
    """Store results of command execution."""

    success: bool
    start_time: datetime.datetime
    end_time: datetime.datetime
    duration: float
    return_code: int
    error_message: str = ""
    output_dir: str = ""
    swept_params: dict[str, Any] = None

    @property
    def duration_str(self) -> str:
        """Format duration as human-readable string."""
        return str(datetime.timedelta(seconds=int(self.duration)))


@dataclass
class SweepCombination:
    """Represents one combination of parameter values."""

    mode: str
    all_params: dict[str, Any]
    swept_params: dict[str, Any]
    profile_name: str
    output_dir: str
    command: str


class SweepEngine:
    """Engine for executing parameter sweeps."""

    def __init__(
        self, config: SweepConfig, base_output_dir: str, nextflow_config_path: str, resume: bool = False
    ):
        self.config = config
        self.base_output_dir = os.path.expandvars(base_output_dir)
        self.nextflow_config_path = nextflow_config_path
        self.resume = resume


    def generate_combinations(self) -> list[SweepCombination]:
        """Generate all parameter combinations for the sweep."""
        combinations = []

        # Get base parameters (only fixed overrides, no nextflow defaults)
        base_params = {**self.config.fixed_params}

        # Generate values for each sweep parameter
        param_names = list(self.config.sweep_params.keys())
        if not param_names:
            # No parameters to sweep
            return []

        param_values = []
        for param_name in param_names:
            sweep = self.config.sweep_params[param_name]
            param_values.append(sweep.generate_values())

        # Generate cartesian product
        for value_combo in product(*param_values):
            # Create swept params dict
            swept_params = dict(zip(param_names, value_combo))

            # Combine with base params
            all_params = {**base_params, **swept_params}

            # Generate output directory
            output_dir = self._generate_output_dir(swept_params)

            # Generate profile name
            profile_name = generate_profile_name(self.config.mode, swept_params)

            # Generate parameter combination identifier
            param_combo_id = self._generate_param_combo_id(swept_params)
            
            # Generate command
            command = self._generate_command(profile_name, output_dir, param_combo_id)

            combinations.append(
                SweepCombination(
                    mode=self.config.mode,
                    all_params=all_params,
                    swept_params=swept_params,
                    profile_name=profile_name,
                    output_dir=output_dir,
                    command=command,
                )
            )

        return combinations

    def generate_quick_test_combinations(self) -> list[SweepCombination]:
        """Generate quick test combinations with reduced designs for preflight testing."""
        combinations = []

        # Get base parameters (only fixed overrides, no nextflow defaults)
        base_params = {**self.config.fixed_params}

        # Override with quick test parameters
        quick_test_params = {
            **base_params,
            "num_designs": 2,
            "seqs_per_design": 2,
        }

        # Generate values for each sweep parameter
        param_names = list(self.config.sweep_params.keys())
        if not param_names:
            # No parameters to sweep
            return []

        param_values = []
        for param_name in param_names:
            sweep = self.config.sweep_params[param_name]
            param_values.append(sweep.generate_values())

        # Generate cartesian product
        for value_combo in product(*param_values):
            # Create swept params dict
            swept_params = dict(zip(param_names, value_combo))

            # Combine with quick test params
            all_params = {**quick_test_params, **swept_params}

            # Generate output directory with "quicktest_" prefix
            output_dir = self._generate_output_dir(swept_params, prefix="quicktest_")

            # Generate profile name with "quicktest_" prefix
            profile_name = generate_profile_name(
                self.config.mode, swept_params, prefix="quicktest_"
            )

            # Generate parameter combination identifier with quicktest prefix
            param_combo_id = f"quicktest_{self._generate_param_combo_id(swept_params)}"
            
            # Generate command
            command = self._generate_command(profile_name, output_dir, param_combo_id)

            combinations.append(
                SweepCombination(
                    mode=self.config.mode,
                    all_params=all_params,
                    swept_params=swept_params,
                    profile_name=profile_name,
                    output_dir=output_dir,
                    command=command,
                )
            )

        return combinations

    def _generate_param_combo_id(self, swept_params: dict[str, Any]) -> str:
        """Generate a parameter combination identifier."""
        from .parameter_converters import get_converter
        
        components = []
        for param_name, value in sorted(swept_params.items()):
            converter = get_converter(param_name)
            formatted_value = converter.format_value_for_name(value)
            
            # Use short parameter names for combination ID
            if param_name in ["hotspots", "hotspot_residues"] and formatted_value.startswith("hotspots"):
                components.append(formatted_value)
            else:
                param_short = param_name.replace("rfd_", "r").replace("fampnn_", "f").replace("mpnn_", "m").replace("af2_", "a").replace("boltz_", "b").replace("_", "")[:6]
                components.append(f"{param_short}{formatted_value}")
        
        return "_".join(components) if components else "default"

    def _generate_output_dir(
        self, swept_params: dict[str, Any], prefix: str = ""
    ) -> str:
        """Generate simplified output directory path."""
        components = []

        # Add prefix if provided
        if prefix:
            components.append(prefix.rstrip("_"))

        # Add swept parameter values to directory name (simplified)
        for param_name, value in sorted(swept_params.items()):
            converter = get_converter(param_name)
            formatted_value = converter.format_value_for_name(value)

            # For parameters with custom formatters that include the parameter name, use formatted value directly
            if param_name in ["hotspots", "hotspot_residues"] and formatted_value.startswith("hotspots"):
                components.append(formatted_value)
            else:
                # Use very short parameter names
                param_short = param_name.replace("rfd_", "r").replace("fampnn_", "f").replace("mpnn_", "m").replace("af2_", "a").replace("boltz_", "b").replace("_", "")[
                    :6
                ]  # Max 6 chars
                components.append(f"{param_short}{formatted_value}")

        # Create simplified path structure
        dir_name = "_".join(components)
        return os.path.join(self.base_output_dir, dir_name)

    def _generate_command(self, profile_name: str, output_dir: str, param_combo: str = None) -> str:
        """Generate nextflow command for a combination."""
        profiles = [profile_name]

        # Add additional profile if specified in config
        if self.config.profile:
            profiles.insert(0, self.config.profile)

        profile_str = ",".join(profiles)

        # Add parameter combination identifier as a custom parameter
        param_combo_arg = f"--bindsweeper_param_combo '{param_combo}'" if param_combo else ""
        
        # Add -resume flag if resume mode is enabled
        resume_arg = "-resume" if self.resume else ""

        # Use both nextflow.config and bindsweeper.config with -C flag
        return f"nextflow -c bindsweeper.config run {self.config.pipeline_path} {resume_arg} -profile {profile_str} --out_dir '{output_dir}' --zip_pdbs false {param_combo_arg}".strip()

    def generate_profiles(self, combinations: list[SweepCombination]) -> list[str]:
        """Generate profile content for all combinations."""
        profiles = []

        for combo in combinations:
            profile_content = generate_profile_content(
                combo.mode, combo.all_params, combo.swept_params, combo.profile_name
            )
            profiles.append(profile_content)

        return profiles

    def execute_combination(self, combination: SweepCombination) -> CommandResult:
        """Execute a single parameter combination."""
        start_time = datetime.datetime.now()

        try:
            os.makedirs(combination.output_dir, exist_ok=True)

            logger.info(f"Executing: {combination.command}")

            process = subprocess.run(
                combination.command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
            )

            end_time = datetime.datetime.now()
            duration = (end_time - start_time).total_seconds()

            return CommandResult(
                success=True,
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                return_code=process.returncode,
                output_dir=combination.output_dir,
                swept_params=combination.swept_params,
            )

        except subprocess.CalledProcessError as e:
            end_time = datetime.datetime.now()
            duration = (end_time - start_time).total_seconds()

            error_msg = (
                f"Command failed with return code {e.returncode}\nSTDERR: {e.stderr}"
            )
            logger.error(error_msg)

            return CommandResult(
                success=False,
                start_time=start_time,
                end_time=end_time,
                duration=duration,
                return_code=e.returncode,
                error_message=error_msg,
                output_dir=combination.output_dir,
                swept_params=combination.swept_params,
            )

    def execute_sweep(
        self,
        combinations: list[SweepCombination],
        dry_run: bool = False,
        continue_on_error: bool = False,
        resume: bool = False,
    ) -> list[CommandResult]:
        """Execute all parameter combinations.
        
        When resume=True, adds -resume flag to Nextflow commands, allowing
        Nextflow's caching mechanism to skip cached tasks. Nextflow will
        automatically detect parameter changes and re-run affected tasks.
        """
        results = []

        logger.info(f"Executing {len(combinations)} parameter combinations")
        
        if resume:
            logger.info("Resume mode: Nextflow will use cached tasks where possible")

        for i, combo in enumerate(combinations, 1):
            logger.info(f"\nProcessing combination {i}/{len(combinations)}:")
            logger.info(f"  Mode: {combo.mode}")

            for param_name, value in combo.swept_params.items():
                logger.info(f"  {param_name}: {value}")

            logger.info(f"  Output Directory: {combo.output_dir}")
            logger.info(f"  Profile: {combo.profile_name}")
            logger.info(f"  Command: {combo.command}")

            if dry_run:
                logger.info("Dry run - skipping execution")
                continue

            result = self.execute_combination(combo)
            results.append(result)

            if not result.success:
                logger.error(f"Command failed after {result.duration_str}")
                if not continue_on_error:
                    raise RuntimeError("Sweep execution failed")
            else:
                logger.info(f"Command completed successfully in {result.duration_str}")

        return results
