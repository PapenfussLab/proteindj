#!/usr/bin/env python3
"""Command-line interface for the sweep tool."""

import logging
import os
import shutil
import subprocess
import sys

import click

from .profile_generator import write_profiles_to_bindsweeper_config
from .results_processor import ResultsProcessor
from .sweep_config import SweepConfig
from .sweep_engine import SweepEngine
from .utils import (
    check_prerequisites,
    find_config_files,
    get_schema_path,
    prompt_for_file,
    setup_logging,
    validate_output_directory,
)

CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}
COMBINATORIAL_EXPLOSION_THRESHOLD = 8
logger = logging.getLogger(__name__)






@click.command(context_settings=CONTEXT_SETTINGS)
@click.option("--debug", is_flag=True, help="Enable debug logging")
@click.version_option(version="0.1.9")
@click.option("--dry-run", is_flag=True, help="Print commands without executing them")
@click.option(
    "--skip-sweep", is_flag=True, help="Skip parameter sweep and only process results"
)
@click.option(
    "--continue-on-error",
    is_flag=True,
    help="Continue if individual parameter sweeps fail",
)
@click.option("-y", "--yes-to-all", is_flag=True, help="Skip confirming files")
@click.option(
    "--quick-test",
    is_flag=True,
    help="Run a quick test with 2 designs and 2 sequences per design before the full run",
)
@click.option(
    "--config",
    type=click.Path(),
    help="Path to sweep configuration YAML file (default: sweep.yaml)",
)
@click.option(
    "--output-dir",
    type=click.Path(),
    help="Output directory for results (default: ./results)",
)
@click.option(
    "--pipeline-path",
    type=click.Path(),
    help="Path to the ProteinDJ main.nf file (default: ./main.nf)",
)
@click.option(
    "--nextflow-config",
    type=click.Path(exists=True),
    help="Path to nextflow.config file (default: auto-discover)",
)
def cli(
    debug: bool,
    dry_run: bool,
    skip_sweep: bool,
    continue_on_error: bool,
    yes_to_all: bool,
    quick_test: bool,
    config: str,
    output_dir: str,
    pipeline_path: str,
    nextflow_config: str,
) -> None:
    """Execute parameter sweep and process results.

    This tool automates parameter sweeps for ProteinDJ/RFdiffusion workflows.
    It reads configuration from sweep.yaml and generates appropriate nextflow
    profiles and commands to explore parameter combinations.
    """

    # Find and validate config files
    skip_confirmation = yes_to_all
    current_dir = os.getcwd()

    try:
        if dry_run:
            click.echo(f"Performing dry run to validate config and preview parameter combinations\n")

        # Use provided nextflow_config or discover it
        if nextflow_config:
            nextflow_config_path = nextflow_config
        else:
            nextflow_config_path = os.path.join(current_dir, "nextflow.config")
            if not os.path.exists(nextflow_config_path):
                # Look in parent directory
                nextflow_config_path = os.path.join(
                    os.path.dirname(current_dir), "nextflow.config"
                )
                if not os.path.exists(nextflow_config_path):
                    raise FileNotFoundError(
                        "nextflow.config not found in current or parent directory"
                    )
        
        # Use provided config (YAML) or default/prompt
        if config:
            sweep_yaml = config
        else:
            sweep_yaml = os.path.join(current_dir, "sweep.yaml")
            if os.path.exists(sweep_yaml):
                click.echo(f"✓ Automatically detected sweep configuration: {sweep_yaml}")
            else:
                if skip_sweep:
                    # In skip_sweep mode, we still need a config for validation, but make it optional
                    click.echo("✓ Skip-sweep mode: config validation will be skipped")
                    sweep_yaml = None
                else:
                    sweep_yaml = prompt_for_file("sweep.yaml", current_dir)
        
        # Use provided output_dir or parse from nextflow config or default
        if output_dir:
            out_dir = output_dir
        else:
            from .utils import parse_out_dir_from_nextflow
            try:
                out_dir = parse_out_dir_from_nextflow(nextflow_config_path, dry_run)
            except:
                # Fallback to default if parsing fails
                out_dir = os.path.join(current_dir, "results")

        config_paths = {
            "nextflow_config": nextflow_config_path,
            "sweep_yaml": sweep_yaml,  # Can be None when skip_sweep=True
            "out_dir": out_dir,
        }

        # Validate output directory
        validate_output_directory(out_dir, dry_run, skip_sweep, skip_confirmation)

        # Setup logging
        global logger
        logger = setup_logging(debug, out_dir, dry_run)

        if debug:
            logger.debug("Debug logging enabled")

        # Validate sweep configuration if provided
        config = None
        if sweep_yaml:
            click.echo(f"Validating sweep configuration: {sweep_yaml}")

            # Determine appropriate schema based on mode
            import yaml
            with open(sweep_yaml) as f:
                sweep_data = yaml.safe_load(f)
            mode = sweep_data.get("mode", "")

            if mode.startswith("binder_") or mode.startswith("bindcraft_"):
                # Use binder schema for binder and bindcraft modes
                schema_path = os.path.join(os.path.dirname(__file__), "binder_schema.json")
            else:
                # Use nextflow schema for other modes
                schema_path = get_schema_path(os.path.dirname(nextflow_config_path))

            config = SweepConfig.from_yaml(sweep_yaml, schema_path)

        results = []
        if skip_sweep:
            click.echo("✓ Skip-sweep mode enabled - processing existing results only")
        else:
            # Check prerequisites
            if not dry_run and not check_prerequisites():
                sys.exit(1)

            # Override pipeline_path if provided via CLI
            if pipeline_path:
                config.pipeline_path = pipeline_path

            # Validate pipeline path exists
            if not os.path.exists(config.pipeline_path):
                raise FileNotFoundError(
                    f"Pipeline script not found: {config.pipeline_path}. "
                    "Make sure you're running bindsweeper from the proteinDJ root directory "
                    "where main.nf is located, or use --pipeline-path to specify the correct path."
                )

            logger.info("Loaded sweep configuration:")
            logger.info(f"Mode: {config.mode}")
            logger.info(f"Fixed parameters: {len(config.fixed_params)}")
            logger.info(f"Sweep parameters: {list(config.sweep_params.keys())}")

            
            # Initialize sweep engine
            engine = SweepEngine(config, out_dir, nextflow_config_path)
        
        if not skip_sweep:
            # Run quick test if requested
            if quick_test:
                logger.info("\n=== Running Quick Test ===")
                quick_combinations = engine.generate_quick_test_combinations()
                logger.info(
                    f"Generated {len(quick_combinations)} quick test combinations"
                )

                if quick_combinations:
                    # Generate quick test profiles
                    quick_profiles = engine.generate_profiles(quick_combinations)
                    write_profiles_to_bindsweeper_config(
                        quick_profiles, "bindsweeper.config", dry_run
                    )

                    if not dry_run:
                        # Execute quick test
                        quick_results = engine.execute_sweep(
                            quick_combinations, dry_run, continue_on_error
                        )

                        # Check if quick test passed
                        failed_tests = [r for r in quick_results if not r.success]
                        if failed_tests:
                            logger.error(
                                f"Quick test failed for {len(failed_tests)} combinations. "
                                "Check configuration before proceeding with full sweep."
                            )
                            if not continue_on_error:
                                logger.error(
                                    "Stopping due to quick test failures. Use --continue-on-error to proceed anyway."
                                )
                                sys.exit(1)
                        else:
                            logger.info("✓ Quick test passed!")
                        
                        # Process quick test results
                        if quick_results:
                            results = quick_results
                            
            else:
                logger.info("\n=== Running Full Sweep ===")
                # Generate parameter combinations
                combinations = engine.generate_combinations()
                logger.info(f"\nGenerated {len(combinations)} parameter combinations")

                # Check for combinatorial explosion and warn user
                if len(combinations) > COMBINATORIAL_EXPLOSION_THRESHOLD:
                    click.echo(f"\nWARNING: Large number of parameter combinations detected!", err=True)
                    click.echo(f"This sweep will generate {len(combinations)} parameter combinations.", err=True)
                    click.echo(f"Each combination runs a full ProteinDJ pipeline which can take significant time and resources.", err=True)
                    click.echo(f"Consider using --quick-test first to validate your configuration with 2 designs per combination.\n", err=True)
                    
                    if not skip_confirmation and not dry_run and not click.confirm("Do you want to proceed with the full sweep?"):
                        click.echo("Sweep cancelled by user.")
                        sys.exit(0)

                if combinations:
                    # Generate profiles
                    profiles = engine.generate_profiles(combinations)

                    # Write profiles to bindsweeper.config
                    write_profiles_to_bindsweeper_config(profiles, "bindsweeper.config", dry_run)

                    # Execute sweep
                    results = engine.execute_sweep(combinations, dry_run, continue_on_error)

        # Process results
        if not dry_run:
            # Use default results config if skip_sweep (no config loaded)
            from .sweep_config import ResultsConfig
            results_config = config.results_config if config else ResultsConfig()
            processor = ResultsProcessor(results_config, out_dir)
            processor.process_results(results, config_paths, dry_run, skip_sweep)
            click.echo("\nSweep completed successfully!")
        else:
            click.echo("\nDry run completed successfully!")

    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        if debug:
            logger.exception("Full traceback:")
        sys.exit(1)


if __name__ == "__main__":
    cli()
