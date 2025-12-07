#!/usr/bin/env python3
"""
generate_success_metrics.py - Generate success metrics JSON for parameter sweeps

This script calculates success rates and generates a JSON file with metrics
that can be consumed by BindSweeper for parameter optimization.
"""

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from typing import Dict, Any, Optional


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate success metrics JSON from pipeline results"
    )
    parser.add_argument(
        "--fold-count", 
        type=int,
        required=True,
        help="Number of designs generated"
    )
    parser.add_argument(
        "--filter-fold-count",
        type=int, 
        required=True,
        help="Number of designs after filtering"
    )
    parser.add_argument(
        "--seq-count",
        type=int,
        required=True, 
        help="Number of sequences generated"
    )
    parser.add_argument(
        "--filter-seq-count",
        type=int,
        required=True,
        help="Number of sequences after filtering"  
    )
    parser.add_argument(
        "--pred-count",
        type=int,
        required=True,
        help="Number of predictions generated"
    )
    parser.add_argument(
        "--filter-pred-count",
        type=int,
        required=True,
        help="Number of predictions after filtering"
    )
    parser.add_argument(
        "--filter-analysis-count",
        type=int,
        required=True,
        help="Number of designs after analysis filtering"
    )
    parser.add_argument(
        "--final-designs-count",
        type=int,
        required=True,
        help="Number of final successful designs"
    )
    parser.add_argument(
        "--parameter-combination",
        default="unknown",
        help="Parameter combination identifier"
    )
    parser.add_argument(
        "--output",
        default="success_metrics.json",
        help="Output JSON file (default: success_metrics.json)"
    )
    return parser.parse_args()


def calculate_success_rate(successful_designs: int, total_designs: int) -> float:
    """Calculate success rate, handling division by zero."""
    if total_designs == 0:
        return 0.0
    return successful_designs / total_designs


def calculate_overall_success_rate(args) -> tuple:
    """
    Calculate overall success rate using the first non-zero count as denominator.
    This handles cases where pipeline stages are skipped.
    
    Returns:
        tuple: (success_rate, total_designs_count) where total_designs_count is the entry point count
        
    Raises:
        ValueError: If the pipeline stage configuration is invalid or unexpected
    """
    # Determine the entry point based on which stages were run
    # Use the earliest unfiltered count as the denominator
    
    if args.fold_count > 0:
        # Full pipeline - start from fold stage
        total = args.fold_count
    elif args.seq_count > 0:
        # skip_fold=true - start from sequence stage
        total = args.seq_count
    elif args.pred_count > 0:
        # skip_fold_seq=true - start from prediction stage (unfiltered)
        total = args.pred_count
    elif args.filter_pred_count > 0:
        # skip_fold_seq_pred=true - start from filtered predictions (user input)
        total = args.filter_pred_count
    else:
        # No valid entry point found - this indicates a configuration error
        raise ValueError(
            "No valid pipeline entry point found. At least one of fold_count, seq_count, "
            "pred_count, or filter_pred_count must be greater than 0. "
            f"Current values: fold={args.fold_count}, seq={args.seq_count}, "
            f"pred={args.pred_count}, filter_pred={args.filter_pred_count}"
        )
    
    return calculate_success_rate(args.final_designs_count, total), total


def generate_success_metrics(args) -> Dict[str, Any]:
    """Generate success metrics dictionary."""
    
    # Calculate overall success rate based on pipeline entry point
    # This correctly handles skipped stages by using the first non-zero count
    success_rate, total_designs = calculate_overall_success_rate(args)
    
    # Build pipeline metrics, marking stages as None when they weren't run
    pipeline_metrics = {}
    
    # Fold retention rate (only if fold stage was run)
    if args.fold_count > 0:
        pipeline_metrics["fold_retention_rate"] = round(
            calculate_success_rate(args.filter_fold_count, args.fold_count), 4
        )
    else:
        pipeline_metrics["fold_retention_rate"] = None
    
    # Sequence retention rate (only if seq stage was run)
    if args.seq_count > 0:
        pipeline_metrics["seq_retention_rate"] = round(
            calculate_success_rate(args.filter_seq_count, args.seq_count), 4
        )
    else:
        pipeline_metrics["seq_retention_rate"] = None
    
    # Prediction retention rate (only if pred stage was run)
    # pred_count represents predictions before filtering
    if args.pred_count > 0:
        # Normal case: predictions were generated and counted
        pipeline_metrics["pred_retention_rate"] = round(
            calculate_success_rate(args.filter_pred_count, args.pred_count), 4
        )
    else:
        # Prediction stage was skipped
        pipeline_metrics["pred_retention_rate"] = None
    
    # Analysis retention rate (only if analysis stage was run)
    if args.filter_pred_count > 0:
        pipeline_metrics["analysis_retention_rate"] = round(
            calculate_success_rate(args.filter_analysis_count, args.filter_pred_count), 4
        )
    else:
        pipeline_metrics["analysis_retention_rate"] = None
    
    # Overall retention rate (always calculated from entry point to final)
    pipeline_metrics["overall_retention_rate"] = round(success_rate, 4)
    
    metrics = {
        "parameter_combination": args.parameter_combination,
        "total_designs": total_designs,
        "successful_designs": args.final_designs_count,
        "success_rate": round(success_rate, 4),
        "fold_generated": args.fold_count,
        "fold_filtered": args.filter_fold_count,
        "seq_generated": args.seq_count,
        "seq_filtered": args.filter_seq_count,
        "pred_generated": args.pred_count,
        "pred_filtered": args.filter_pred_count,
        "analysis_filtered": args.filter_analysis_count,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "pipeline_metrics": pipeline_metrics
    }
    
    return metrics


def main():
    """Main function to generate success metrics JSON."""
    args = parse_arguments()
    
    try:
        # Generate metrics
        metrics = generate_success_metrics(args)
        
        # Write to JSON file
        with open(args.output, 'w') as f:
            json.dump(metrics, f, indent=2)
        
        print(f"Generated success metrics: {args.output}")
        print(f"Success rate: {metrics['success_rate']:.1%} ({metrics['successful_designs']}/{metrics['total_designs']})")
        
        return 0
        
    except Exception as e:
        print(f"Error generating success metrics: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())