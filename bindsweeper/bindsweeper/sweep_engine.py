#!/usr/bin/env python3
"""Core sweep execution engine."""

import datetime
import logging
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import product
from typing import Any

from .parameter_converters import get_converter
from .profile_generator import generate_profile_content, generate_profile_name
from .sweep_config import SweepConfig
from .sweep_types import PairedSweep

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
        self, config: SweepConfig, base_output_dir: str, nextflow_config_path: str, 
        resume: bool = False, parallel: bool = False, max_parallel: int = 4
    ):
        self.config = config
        self.base_output_dir = os.path.expandvars(base_output_dir)
        self.nextflow_config_path = nextflow_config_path
        self.resume = resume
        self.parallel = parallel
        self.max_parallel = max_parallel


    def generate_combinations(self) -> list[SweepCombination]:
        """Generate all parameter combinations for the sweep."""
        # Get base parameters (only fixed overrides, no nextflow defaults)
        base_params = {**self.config.fixed_params}

        # Separate paired and unpaired sweep parameters
        paired_params = {}
        unpaired_params = {}
        
        for param_name, sweep in self.config.sweep_params.items():
            if isinstance(sweep, PairedSweep):
                paired_params[param_name] = sweep
            else:
                unpaired_params[param_name] = sweep

        # Generate combinations differently based on whether we have paired params
        if paired_params:
            combinations = self._generate_paired_combinations_internal(base_params, paired_params, unpaired_params)
        else:
            combinations = self._generate_cartesian_combinations_internal(base_params, unpaired_params)

        return combinations

    def generate_quick_test_combinations(self) -> list[SweepCombination]:
        """Generate quick test combinations with reduced designs for preflight testing."""
        # Get base parameters (only fixed overrides, no nextflow defaults)
        base_params = {**self.config.fixed_params}

        # Override with quick test parameters
        quick_test_params = {
            **base_params,
            "num_designs": 2,
            "seqs_per_design": 2,
        }

        # Separate paired and unpaired sweep parameters
        paired_params = {}
        unpaired_params = {}
        
        for param_name, sweep in self.config.sweep_params.items():
            if isinstance(sweep, PairedSweep):
                paired_params[param_name] = sweep
            else:
                unpaired_params[param_name] = sweep

        # Generate combinations differently based on whether we have paired params
        if paired_params:
            combinations = self._generate_paired_combinations_internal(
                quick_test_params, paired_params, unpaired_params, prefix="quicktest_"
            )
        else:
            combinations = self._generate_cartesian_combinations_internal(
                quick_test_params, unpaired_params, prefix="quicktest_"
            )

        return combinations

    def _generate_cartesian_combinations_internal(
        self, base_params: dict, unpaired_params: dict, prefix: str = ""
    ) -> list[SweepCombination]:
        """Internal method to generate Cartesian product with optional prefix."""
        combinations = []
        
        # Generate values for each sweep parameter
        param_names = list(unpaired_params.keys())
        if not param_names:
            # No parameters to sweep
            return []

        param_values = []
        for param_name in param_names:
            sweep = unpaired_params[param_name]
            param_values.append(sweep.generate_values())

        # Generate cartesian product
        for value_combo in product(*param_values):
            # Create swept params dict
            swept_params = dict(zip(param_names, value_combo))

            # Combine with base params
            all_params = {**base_params, **swept_params}

            # Generate output directory with optional prefix
            output_dir = self._generate_output_dir(swept_params, prefix=prefix)

            # Generate profile name with optional prefix
            profile_name = generate_profile_name(self.config.mode, swept_params, prefix=prefix)

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

    def _generate_paired_combinations_internal(
        self, base_params: dict, paired_params: dict, unpaired_params: dict, prefix: str = ""
    ) -> list[SweepCombination]:
        """Internal method to generate paired combinations with optional prefix."""
        combinations = []
        
        # Currently only support one paired parameter group
        if len(paired_params) > 1:
            raise ValueError(
                "Currently only one paired parameter group is supported. "
                "Multiple parameters can be paired together using 'paired_with', "
                "but you cannot have multiple separate paired groups."
            )
        
        # Get the paired parameter
        primary_param_name = list(paired_params.keys())[0]
        paired_sweep = paired_params[primary_param_name]
        
        # Cache generated values to avoid repeated generation
        primary_values = paired_sweep.generate_values()
        num_paired_combos = len(primary_values)
        
        # Generate unpaired cartesian product if any
        unpaired_combos = []
        if unpaired_params:
            param_names = list(unpaired_params.keys())
            param_values = [unpaired_params[name].generate_values() for name in param_names]
            for value_combo in product(*param_values):
                unpaired_combos.append(dict(zip(param_names, value_combo)))
        else:
            unpaired_combos = [{}]  # Single empty combo if no unpaired params
        
        # Generate combinations: zip paired params, then cartesian with unpaired
        for i in range(num_paired_combos):
            # Get paired parameter values for this index
            paired_swept_params = {
                primary_param_name: primary_values[i]
            }
            
            # Add all paired parameters
            for paired_param_name in paired_sweep.paired_params.keys():
                paired_swept_params[paired_param_name] = paired_sweep.get_paired_value(paired_param_name, i)
            
            # Combine with each unpaired combination
            for unpaired_combo in unpaired_combos:
                # Merge paired and unpaired swept params
                swept_params = {**paired_swept_params, **unpaired_combo}
                
                # Combine with base params
                all_params = {**base_params, **swept_params}

                # Generate output directory with optional prefix
                output_dir = self._generate_output_dir(swept_params, prefix=prefix)

                # Generate profile name with optional prefix
                profile_name = generate_profile_name(self.config.mode, swept_params, prefix=prefix)

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

    def _generate_param_combo_id(self, swept_params: dict[str, Any]) -> str:
        """Generate a parameter combination identifier."""
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

    def execute_combination(self, combination: SweepCombination, use_isolated_cache: bool = False) -> CommandResult:
        """Execute a single parameter combination.
        
        Args:
            combination: The parameter combination to execute
            use_isolated_cache: If True, use isolated cache directory per combination
        """
        start_time = datetime.datetime.now()

        try:
            os.makedirs(combination.output_dir, exist_ok=True)

            logger.info(f"Executing: {combination.command}")

            # Set isolated cache directory for parallel execution
            env = os.environ.copy()
            if use_isolated_cache:
                cache_dir = os.path.join(combination.output_dir, ".nextflow_cache")
                os.makedirs(cache_dir, exist_ok=True)
                env["NXF_CACHE_DIR"] = cache_dir
                logger.debug(f"Using isolated cache: {cache_dir}")

            process = subprocess.run(
                combination.command,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                env=env,
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
        parallel: bool = None,
    ) -> list[CommandResult]:
        """Execute all parameter combinations.
        
        When resume=True, adds -resume flag to Nextflow commands, allowing
        Nextflow's caching mechanism to skip cached tasks. Nextflow will
        automatically detect parameter changes and re-run affected tasks.
        
        When parallel=True, executes combinations concurrently using ThreadPoolExecutor
        with isolated cache directories per combination to prevent conflicts.
        
        Args:
            combinations: List of parameter combinations to execute
            dry_run: If True, print commands without executing
            continue_on_error: If True, continue even if combinations fail
            resume: If True, add -resume flag to Nextflow commands
            parallel: If True, execute combinations in parallel. If None, uses self.parallel
        """
        # Use instance parallel setting if not explicitly provided
        if parallel is None:
            parallel = self.parallel
            
        results = []

        logger.info(f"Executing {len(combinations)} parameter combinations")
        
        if resume:
            logger.info("Resume mode: Nextflow will use cached tasks where possible")
        
        if parallel and not dry_run:
            logger.info(f"Parallel mode: Running up to {self.max_parallel} combinations concurrently")
            logger.info("Each combination uses isolated Nextflow cache to prevent conflicts")
            results = self._execute_parallel(
                combinations, continue_on_error
            )
        else:
            # Sequential execution (existing logic)
            if parallel and dry_run:
                logger.info("Parallel mode enabled but running in dry-run (sequential preview)")
            
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

                result = self.execute_combination(combo, use_isolated_cache=False)
                results.append(result)

                if not result.success:
                    logger.error(f"Command failed after {result.duration_str}")
                    if not continue_on_error:
                        raise RuntimeError("Sweep execution failed")
                else:
                    logger.info(f"Command completed successfully in {result.duration_str}")

        return results
    
    def _execute_parallel(
        self,
        combinations: list[SweepCombination],
        continue_on_error: bool = False,
    ) -> list[CommandResult]:
        """Execute combinations in parallel using ThreadPoolExecutor.
        
        Each combination gets an isolated Nextflow cache directory to prevent
        conflicts between concurrent runs. Uses threads rather than processes
        since subprocess calls release the GIL, allowing true parallelism.
        """
        results = []
        failed_count = 0
        
        # Use ThreadPoolExecutor to run combinations in parallel
        # Threads work well here since subprocess.run releases GIL
        with ThreadPoolExecutor(max_workers=self.max_parallel) as executor:
            # Submit all combinations
            future_to_combo = {
                executor.submit(self.execute_combination, combo, use_isolated_cache=True): combo
                for combo in combinations
            }
            
            # Process results as they complete
            for i, future in enumerate(as_completed(future_to_combo), 1):
                combo = future_to_combo[future]
                try:
                    result = future.result()
                    results.append(result)
                    
                    logger.info(f"\nCompleted {i}/{len(combinations)}: {combo.profile_name}")
                    
                    if result.success:
                        logger.info(f"  ✓ Success in {result.duration_str}")
                    else:
                        failed_count += 1
                        logger.error(f"  ✗ Failed after {result.duration_str}")
                        if not continue_on_error:
                            # Cancel remaining futures
                            for f in future_to_combo:
                                f.cancel()
                            raise RuntimeError(f"Combination failed: {combo.profile_name}")
                    
                except Exception as e:
                    failed_count += 1
                    logger.error(f"\nException in combination {combo.profile_name}: {e}")
                    if not continue_on_error:
                        # Cancel remaining futures
                        for f in future_to_combo:
                            f.cancel()
                        raise
        
        logger.info(f"\nParallel execution complete: {len(results) - failed_count}/{len(results)} succeeded")
        return results
