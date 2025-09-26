#!/usr/bin/env python3
"""Results processing and merging for parameter sweeps."""

import logging
import os
import shutil
import zipfile
from typing import Optional

import pandas as pd

from .sweep_config import ResultsConfig
from .sweep_engine import CommandResult
from .success_rate_analyzer import SuccessRateAnalyzer

logger = logging.getLogger(__name__)


class ResultsProcessor:
    """Process and merge results from parameter sweeps."""

    def __init__(self, config: ResultsConfig, base_output_dir: str):
        self.config = config
        self.base_output_dir = os.path.expandvars(base_output_dir)
        self.success_analyzer = SuccessRateAnalyzer(base_output_dir)

    def process_results(
        self,
        results: list[CommandResult],
        config_paths: dict,
        dry_run: bool = False,
        skip_sweep: bool = False,
    ) -> None:
        """Process results after parameter sweep completion."""
        logger.info("Starting results processing...")

        # Find all CSV files
        csv_info = self.find_rank_csvs(results, skip_sweep)

        # Prepare output paths
        pdb_output_dir = os.path.join(self.base_output_dir, self.config.pdb_output_dir)
        output_csv = os.path.join(self.base_output_dir, self.config.output_csv)

        if csv_info:
            # Copy and rename PDB files
            copied_files = self.copy_and_rename_pdbs(csv_info, pdb_output_dir, dry_run)

            # Merge CSVs
            merged_df = self.merge_csvs(csv_info, dry_run)

            if not dry_run and merged_df is not None:
                merged_df.to_csv(output_csv, index=False)
                logger.info(f"\nMerged {len(csv_info)} CSV files into {output_csv}")
                logger.info(f"Total rows in merged file: {len(merged_df)}")
                logger.info("\nOrigins included:")
                for origin in sorted(merged_df["origin"].unique()):
                    count = len(merged_df[merged_df["origin"] == origin])
                    logger.info(f"  {origin}: {count} rows")

                logger.info(f"\nCopied {len(copied_files)} PDB files to {pdb_output_dir}")

                if self.config.zip_results:
                    zip_path = self.zip_results_with_configs(
                        pdb_output_dir, output_csv, config_paths
                    )
                    if zip_path:
                        logger.info(f"Created zip file: {zip_path}")
        else:
            # No successful designs found - create empty CSV
            logger.info("No successful designs found to process, creating empty results CSV")
            if not dry_run:
                # Create empty DataFrame with expected structure
                empty_df = pd.DataFrame()
                empty_df.to_csv(output_csv, index=False)
                logger.info(f"Created empty results file: {output_csv}")

        # Always process success rate analysis regardless of whether CSV files exist
        self._process_success_metrics(results, dry_run, skip_sweep)

    def find_rank_csvs(
        self, results: list[CommandResult], skip_sweep: bool = False
    ) -> list[tuple[str, str, str]]:
        """Find all CSV files in rank directories and their corresponding extract/results directories."""
        csv_info = []

        if skip_sweep:
            # Look for existing results in output directory
            for dirname in os.listdir(self.base_output_dir):
                dirpath = os.path.join(self.base_output_dir, dirname)
                if not os.path.isdir(dirpath):
                    continue

                # Try the expected structure (rank/best.csv)
                rank_dir = os.path.join(dirpath, self.config.rank_dirname)
                if os.path.exists(rank_dir) and os.path.isdir(rank_dir):
                    csv_path = os.path.join(rank_dir, self.config.csv_filename)
                    if os.path.exists(csv_path):
                        pdb_dir = os.path.join(rank_dir, self.config.results_dirname)
                        if os.path.exists(pdb_dir):
                            csv_info.append((csv_path, dirname, pdb_dir))
                            continue
                        else:
                            logger.warning(
                                f"No {self.config.results_dirname} "
                                f"directory found for {dirname}"
                            )
                # If no results found, silently skip (no passing results for this directory)
        else:
            # Process results from sweep
            for result in results:
                if not result.success:
                    continue

                root_dir = os.path.expandvars(result.output_dir)

                for dirpath, _dirnames, filenames in os.walk(root_dir):
                    if (
                        os.path.basename(dirpath) == self.config.rank_dirname
                        and self.config.csv_filename in filenames
                    ):
                        csv_path = os.path.join(dirpath, self.config.csv_filename)
                        origin = os.path.basename(os.path.dirname(dirpath))
                        pdb_dir = os.path.join(
                            dirpath,
                            self.config.results_dirname,
                        )

                        if os.path.exists(pdb_dir):
                            csv_info.append((csv_path, origin, pdb_dir))
                        else:
                            logger.warning(
                                f"No {self.config.results_dirname} "
                                f"directory found for {origin}"
                            )

        return csv_info

    def copy_and_rename_pdbs(
        self,
        csv_info: list[tuple[str, str, str]],
        output_dir: str,
        dry_run: bool = False,
    ) -> list[str]:
        """Copy and rename PDB files from extract/results directories."""
        if not dry_run:
            os.makedirs(output_dir, exist_ok=True)

        copied_files = []
        for _, origin, extract_dir in csv_info:
            for filename in os.listdir(extract_dir):
                if filename.endswith(".pdb"):
                    old_path = os.path.join(extract_dir, filename)
                    new_name = f"{os.path.splitext(filename)[0]}_{origin}.pdb"
                    new_path = os.path.join(output_dir, new_name)

                    if dry_run:
                        logger.info(f"Would copy: {old_path} to: {new_path}")
                    else:
                        shutil.copy2(old_path, new_path)
                        copied_files.append(new_path)
                        logger.debug(f"Copied {old_path} to {new_path}")

        return copied_files

    def merge_csvs(
        self, csv_info: list[tuple[str, str, str]], dry_run: bool = False
    ) -> Optional[pd.DataFrame]:
        """Merge CSVs and add origin column."""
        if dry_run:
            logger.info("\nFiles that would be merged:")
            for path, origin, _ in csv_info:
                logger.info(f"File: {path}")
                logger.info(f"Origin: {origin}\n")
            return None

        dfs = []
        for path, origin, _ in csv_info:
            try:
                df = pd.read_csv(path)
                df["origin"] = origin
                dfs.append(df)
                logger.debug(f"Successfully read and processed {path}")
            except Exception as e:
                logger.error(f"Error processing CSV file {path}: {e}")

        if not dfs:
            logger.error("No CSV files found or processed!")
            return None

        return pd.concat(dfs, ignore_index=True)

    def zip_results_with_configs(
        self, pdb_dir: str, csv_file: str, config_paths: dict
    ) -> str:
        """Create organized zip file with results and configs."""
        zip_path = os.path.join(self.base_output_dir, "merged_results.zip")
        logger.info(f"Creating zip archive: {zip_path}")

        try:
            with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
                # Add config files
                if "sweep_yaml" in config_paths and config_paths["sweep_yaml"] is not None:
                    zipf.write(
                        config_paths["sweep_yaml"],
                        f"inputs/{os.path.basename(config_paths['sweep_yaml'])}",
                    )

                # Add log file
                log_file = os.path.join(self.base_output_dir, "sweep.log")
                if os.path.exists(log_file):
                    zipf.write(log_file, f"logs/{os.path.basename(log_file)}")

                # Add CSV file
                zipf.write(csv_file, os.path.basename(csv_file))

                # Add PDB files
                for root, _, files in os.walk(pdb_dir):
                    for file in files:
                        if file.endswith(".pdb"):
                            file_path = os.path.join(root, file)
                            arcname = os.path.join("pdbs", file)
                            zipf.write(file_path, arcname)

            # Clean up original files
            shutil.rmtree(pdb_dir)
            os.remove(csv_file)
            logger.info("Cleaned up temporary files")

            return zip_path
        except Exception as e:
            logger.error(f"Error creating zip archive: {e}")
            return ""
    
    def _process_success_metrics(self, results: list[CommandResult], dry_run: bool = False, skip_sweep: bool = False) -> None:
        """Process success rate metrics and generate analysis reports."""
        if dry_run:
            logger.info("Dry run - skipping success rate analysis")
            return
        
        logger.info("\n=== Processing Success Rate Metrics ===")
        
        # Collect success metrics from all results (analyzer handles both sweep and existing results)
        success_metrics = self.success_analyzer.collect_success_metrics(results)
        
        if not success_metrics:
            logger.warning("No success metrics found to analyze")
            return
        
        # Generate success summary CSV
        success_summary_path = os.path.join(self.base_output_dir, "sweep_success_summary.csv")
        self.success_analyzer.generate_success_summary_csv(success_summary_path)
        
        # Save best parameter combination
        best_combo_path = os.path.join(self.base_output_dir, "best_parameter_combination.json")
        self.success_analyzer.save_best_parameter_combination(best_combo_path)
        
        # Generate and log comprehensive report
        report = self.success_analyzer.generate_success_report()
        
        # Log summary statistics
        if "summary" in report:
            summary = report["summary"]
            logger.info(f"\nSuccess Rate Analysis Summary:")
            logger.info(f"  Total parameter combinations: {summary['total_parameter_combinations']}")
            logger.info(f"  Mean success rate: {summary['mean_success_rate']:.1%}")
            logger.info(f"  Best success rate: {summary['max_success_rate']:.1%}")
            logger.info(f"  Worst success rate: {summary['min_success_rate']:.1%}")
            logger.info(f"  Total successful designs: {summary['total_successful_designs']}")
            logger.info(f"  Total attempted designs: {summary['total_attempted_designs']}")
        
        # Log best combination details
        if "best_combination" in report and report["best_combination"]:
            best = report["best_combination"]
            logger.info(f"\nBest Parameter Combination:")
            logger.info(f"  Parameters: {best['parameter_combination']}")
            logger.info(f"  Success rate: {best['success_rate']:.1%}")
            logger.info(f"  Successful designs: {best['successful_designs']}/{best['total_designs']}")
        
        logger.info(f"\nSuccess rate analysis complete!")
        logger.info(f"  Summary CSV: {success_summary_path}")
        logger.info(f"  Best combination JSON: {best_combo_path}")
