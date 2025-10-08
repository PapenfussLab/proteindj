#!/usr/bin/env python3
"""Utility functions for the sweep tool."""

import logging
import os
import shutil
import sys
from logging.handlers import RotatingFileHandler
from typing import Optional


def setup_logging(debug: bool, out_dir: str, dry_run: bool = False) -> logging.Logger:
    """Setup logging configuration with log file in output directory."""
    out_dir = os.path.expandvars(out_dir)
    log_file = os.path.join(out_dir, "sweep.log")

    if dry_run:
        print(f"Would set up logging to: {log_file}")
        return logging.getLogger()

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    file_handler = RotatingFileHandler(
        log_file,
        maxBytes=10 * 1024 * 1024,  # 10MB
        backupCount=5,
    )
    file_handler.setFormatter(formatter)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    logger.handlers = []
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def check_prerequisites() -> bool:
    """Check if required tools are available."""
    if shutil.which("nextflow") is None:
        logging.error("Nextflow is not installed or not in PATH")
        return False
    return True


def find_config_files(
    current_dir: str, dry_run: bool = False, skip_confirmation: bool = False
) -> tuple[str, str, str]:
    """Find or prompt for config files and validate nextflow config."""
    # Find nextflow.config
    nextflow_config = os.path.join(current_dir, "nextflow.config")
    if not os.path.isfile(nextflow_config):
        nextflow_config = prompt_for_file("nextflow.config", current_dir)

    # Find sweep.yaml
    sweep_yaml = os.path.join(current_dir, "sweep.yaml")
    if not os.path.isfile(sweep_yaml):
        sweep_yaml = prompt_for_file("sweep.yaml", current_dir)

    # Confirm config files with user
    if not skip_confirmation and not confirm_config_files(nextflow_config, sweep_yaml):
        sys.exit("User cancelled - config files not confirmed")

    # Parse output directory from nextflow config
    out_dir = parse_out_dir_from_nextflow(nextflow_config, dry_run)

    # Confirm output directory
    if not skip_confirmation and not confirm_output_directory(out_dir):
        sys.exit("User cancelled - output directory not confirmed")

    return nextflow_config, sweep_yaml, out_dir


def prompt_for_file(filename: str, current_dir: str) -> str:
    """Prompt user for directory containing a specific file."""
    while True:
        print(f"\n{filename} not found in {current_dir}")
        dir_path = input(
            f"Please enter the directory containing {filename} (or 'q' to quit): "
        ).strip()

        if dir_path.lower() == "q":
            sys.exit("User requested quit")

        file_path = os.path.join(dir_path, filename)
        if os.path.isfile(file_path):
            return file_path
        print(f"Error: {filename} not found in specified directory")


def confirm_config_files(nextflow_config: str, sweep_yaml: str) -> bool:
    """Confirm found configuration files with user."""
    print("\nFound in current directory:")
    print(f"- nextflow.config: {nextflow_config}")
    print(f"- sweep.yaml: {sweep_yaml}")

    confirm = input("Proceed with these files? [y/N]: ")
    return confirm.lower() == "y"


def confirm_output_directory(out_dir: str) -> bool:
    """Confirm parsed output directory with user."""
    print("\nParsed output directory from nextflow.config:")
    print(f"- Output directory: {out_dir}")

    confirm = input("Proceed with this output directory? [y/N]: ")
    return confirm.lower() == "y"


def parse_out_dir_from_nextflow(config_path: str, dry_run: bool = False) -> str:
    """Parse nextflow.config file to extract out_dir parameter."""
    import re

    try:
        with open(config_path) as f:
            content = f.read()

        # Look for out_dir parameter, excluding commented lines
        match = re.search(
            r'^[^/#]*?(?:^|\s)out_dir\s*=\s*[\'"]([^\'"]+)[\'"]', content, re.MULTILINE
        )
        if not match:
            sys.exit("Error: out_dir parameter not found in nextflow.config")

        out_dir = match.group(1)

        # Expand environment variables in the path
        out_dir = os.path.expandvars(out_dir)

        # Validate absolute path
        if not os.path.isabs(out_dir):
            sys.exit("Error: out_dir in nextflow.config must be an absolute path")

        return out_dir

    except Exception as e:
        sys.exit(f"Error parsing nextflow.config: {str(e)}")


def validate_output_directory(
    out_dir: str,
    dry_run: bool = False,
    skip_sweep: bool = False,
    skip_confirmation: bool = False,
) -> None:
    """Validate and potentially create output directory."""
    out_dir = os.path.expandvars(out_dir)

    if dry_run:
        if not os.path.exists(out_dir):
            print(f"Would create output directory: {out_dir}")
        elif os.listdir(out_dir) and not skip_sweep:
            print(f"Would prompt for confirmation (directory not empty): {out_dir}")
        return

    # Create output directory if it doesn't exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
        return  # New empty directory, no need for further validation

    if os.listdir(out_dir) and not skip_sweep and not skip_confirmation:
        confirm = input(
            f"Output directory '{out_dir}' is not empty. Contents may be overwritten. Continue? [y/N]: "
        )
        if confirm.lower() != "y":
            sys.exit("User cancelled due to existing files")


def get_schema_path(nextflow_dir: str) -> Optional[str]:
    """Try to find the nextflow schema JSON file."""
    possible_paths = [
        os.path.join(nextflow_dir, "nextflow_schema.json"),
        os.path.join(nextflow_dir, "schema.json"),
        os.path.join(os.path.dirname(nextflow_dir), "nextflow_schema.json"),
        # Look in schemas subdirectory (ProteinDJ structure)
        os.path.join(nextflow_dir, "schemas", "nextflow_schema.json"),
        os.path.join(os.path.dirname(nextflow_dir), "schemas", "nextflow_schema.json"),
    ]

    for path in possible_paths:
        if os.path.exists(path):
            return path

    return None
