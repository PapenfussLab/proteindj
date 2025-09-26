#!/usr/bin/env python3
"""Success rate analysis for parameter sweeps."""

import json
import logging
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class SuccessMetrics:
    """Represents success metrics for a single parameter combination."""
    
    parameter_combination: str
    total_designs: int
    successful_designs: int
    success_rate: float
    rfd_generated: int
    rfd_filtered: int
    seq_generated: int
    seq_filtered: int
    pred_filtered: int
    timestamp: str
    pipeline_metrics: Dict[str, float]
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SuccessMetrics":
        """Create SuccessMetrics from dictionary."""
        return cls(
            parameter_combination=data["parameter_combination"],
            total_designs=data["total_designs"],
            successful_designs=data["successful_designs"],
            success_rate=data["success_rate"],
            rfd_generated=data["rfd_generated"],
            rfd_filtered=data["rfd_filtered"],
            seq_generated=data["seq_generated"],
            seq_filtered=data["seq_filtered"],
            pred_filtered=data["pred_filtered"],
            timestamp=data["timestamp"],
            pipeline_metrics=data.get("pipeline_metrics", {})
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        return {
            "parameter_combination": self.parameter_combination,
            "total_designs": self.total_designs,
            "successful_designs": self.successful_designs,
            "success_rate": self.success_rate,
            "rfd_generated": self.rfd_generated,
            "rfd_filtered": self.rfd_filtered,
            "seq_generated": self.seq_generated,
            "seq_filtered": self.seq_filtered,
            "pred_filtered": self.pred_filtered,
            "timestamp": self.timestamp,
            "pipeline_metrics": self.pipeline_metrics
        }


class SuccessRateAnalyzer:
    """Analyzes success rates from parameter sweep results."""
    
    def __init__(self, base_output_dir: str):
        self.base_output_dir = os.path.expandvars(base_output_dir)
        self.success_metrics: List[SuccessMetrics] = []
    
    def collect_success_metrics(self, results: list) -> List[SuccessMetrics]:
        """Collect success metrics from all sweep results."""
        metrics = []
        
        if results:
            # Process from sweep results
            for result in results:
                # Always check for success_metrics.json regardless of result.success status
                # The success_metrics.json file is the authoritative source of truth
                metrics_file = os.path.join(result.output_dir, "results", "success_metrics.json")
                if os.path.exists(metrics_file):
                    try:
                        with open(metrics_file, 'r') as f:
                            data = json.load(f)
                        success_metrics = SuccessMetrics.from_dict(data)
                        metrics.append(success_metrics)
                        logger.debug(f"Loaded success metrics from {metrics_file}")
                    except Exception as e:
                        logger.error(f"Error loading success metrics from {metrics_file}: {e}")
                else:
                    if not result.success:
                        logger.debug(f"Skipping result with no success metrics: {result.output_dir}")
                    else:
                        logger.warning(f"Success metrics file not found: {metrics_file}")
        else:
            # Scan existing directories for success metrics
            for dirname in os.listdir(self.base_output_dir):
                dirpath = os.path.join(self.base_output_dir, dirname)
                if os.path.isdir(dirpath) and not dirname.endswith(('.zip', '.csv', '.log')):
                    metrics_file = os.path.join(dirpath, "results", "success_metrics.json")
                    if os.path.exists(metrics_file):
                        try:
                            with open(metrics_file, 'r') as f:
                                data = json.load(f)
                            success_metrics = SuccessMetrics.from_dict(data)
                            metrics.append(success_metrics)
                            logger.debug(f"Loaded success metrics from {metrics_file}")
                        except Exception as e:
                            logger.error(f"Error loading success metrics from {metrics_file}: {e}")
        
        self.success_metrics = metrics
        logger.info(f"Collected success metrics from {len(metrics)} parameter combinations")
        return metrics
    
    def find_best_parameter_combination(self) -> Optional[SuccessMetrics]:
        """Find the parameter combination with the highest success rate."""
        if not self.success_metrics:
            logger.warning("No success metrics available")
            return None
        
        # Sort by success rate (descending), then by successful designs count (descending)
        best = max(self.success_metrics, 
                  key=lambda x: (x.success_rate, x.successful_designs))
        
        logger.info(f"Best parameter combination: {best.parameter_combination} "
                   f"(success rate: {best.success_rate:.1%}, "
                   f"successful designs: {best.successful_designs})")
        
        return best
    
    def generate_success_summary_csv(self, output_path: str) -> None:
        """Generate a CSV summary of all success rates, ranked by performance."""
        if not self.success_metrics:
            logger.warning("No success metrics to summarize")
            return
        
        # Sort by success rate descending, then by successful designs descending
        sorted_metrics = sorted(self.success_metrics, 
                              key=lambda x: (x.success_rate, x.successful_designs), 
                              reverse=True)
        
        # Create DataFrame
        summary_data = []
        for i, metrics in enumerate(sorted_metrics, 1):
            row = {
                'rank': i,
                'parameter_combination': metrics.parameter_combination,
                'success_rate': metrics.success_rate,
                'successful_designs': metrics.successful_designs,
                'total_designs': metrics.total_designs,
                'rfd_retention_rate': metrics.pipeline_metrics.get('rfd_retention_rate', 0.0),
                'seq_retention_rate': metrics.pipeline_metrics.get('seq_retention_rate', 0.0),
                'pred_retention_rate': metrics.pipeline_metrics.get('pred_retention_rate', 0.0),
                'overall_retention_rate': metrics.pipeline_metrics.get('overall_retention_rate', 0.0),
                'timestamp': metrics.timestamp
            }
            summary_data.append(row)
        
        df = pd.DataFrame(summary_data)
        df.to_csv(output_path, index=False)
        
        logger.info(f"Generated success summary CSV: {output_path}")
        logger.info(f"Top 3 parameter combinations:")
        for i, metrics in enumerate(sorted_metrics[:3], 1):
            logger.info(f"  {i}. {metrics.parameter_combination}: "
                       f"{metrics.success_rate:.1%} ({metrics.successful_designs}/{metrics.total_designs})")
    
    def save_best_parameter_combination(self, output_path: str) -> None:
        """Save the best parameter combination to a JSON file."""
        best = self.find_best_parameter_combination()
        if not best:
            logger.warning("No best parameter combination to save")
            return
        
        best_data = {
            "best_parameter_combination": best.parameter_combination,
            "success_rate": best.success_rate,
            "successful_designs": best.successful_designs,
            "total_designs": best.total_designs,
            "analysis_timestamp": pd.Timestamp.now().isoformat(),
            "detailed_metrics": best.to_dict()
        }
        
        with open(output_path, 'w') as f:
            json.dump(best_data, f, indent=2)
        
        logger.info(f"Saved best parameter combination to: {output_path}")
    
    def generate_success_report(self) -> Dict[str, Any]:
        """Generate a comprehensive success rate report."""
        if not self.success_metrics:
            return {"error": "No success metrics available"}
        
        # Calculate summary statistics
        success_rates = [m.success_rate for m in self.success_metrics]
        successful_designs = [m.successful_designs for m in self.success_metrics]
        
        best = self.find_best_parameter_combination()
        worst = min(self.success_metrics, key=lambda x: (x.success_rate, x.successful_designs))
        
        report = {
            "summary": {
                "total_parameter_combinations": len(self.success_metrics),
                "mean_success_rate": sum(success_rates) / len(success_rates),
                "max_success_rate": max(success_rates),
                "min_success_rate": min(success_rates),
                "total_successful_designs": sum(successful_designs),
                "total_attempted_designs": sum(m.total_designs for m in self.success_metrics)
            },
            "best_combination": best.to_dict() if best else None,
            "worst_combination": worst.to_dict() if worst else None,
            "all_combinations": [m.to_dict() for m in sorted(
                self.success_metrics, 
                key=lambda x: (x.success_rate, x.successful_designs), 
                reverse=True
            )]
        }
        
        return report