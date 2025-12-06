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


def generate_success_metrics(args) -> Dict[str, Any]:
    """Generate success metrics dictionary."""
    
    # Calculate overall success rate based on final designs vs total sequences generated
    # This gives the true success rate: how many of the generated sequences made it through
    success_rate = calculate_success_rate(args.final_designs_count, args.seq_count)
    
    metrics = {
        "parameter_combination": args.parameter_combination,
        "total_designs": args.seq_count,
        "successful_designs": args.final_designs_count,
        "success_rate": round(success_rate, 4),
        "fold_generated": args.fold_count,
        "fold_filtered": args.filter_fold_count,
        "seq_generated": args.seq_count,
        "seq_filtered": args.filter_seq_count,
        "pred_filtered": args.filter_pred_count,
        "analysis_filtered": args.filter_analysis_count,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "pipeline_metrics": {
            "fold_retention_rate": round(calculate_success_rate(args.filter_fold_count, args.fold_count), 4),
            "seq_retention_rate": round(calculate_success_rate(args.filter_seq_count, args.seq_count), 4),
            "pred_retention_rate": round(calculate_success_rate(args.filter_pred_count, args.filter_seq_count), 4) if args.filter_seq_count > 0 else 0.0,
            "analysis_retention_rate": round(calculate_success_rate(args.filter_analysis_count, args.filter_pred_count), 4) if args.filter_pred_count > 0 else 0.0,
            "overall_retention_rate": round(calculate_success_rate(args.final_designs_count, args.seq_count), 4)
        }
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