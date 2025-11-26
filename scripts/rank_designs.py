#!/usr/bin/env python3
"""
rank_designs.py - Rank and select top protein designs based on prediction metrics

This script ranks designs from best_designs.csv based on prediction quality metrics,
renames PDB files with ranked prefixes, and outputs the top N designs.
"""

import argparse
import os
import shutil
import pandas as pd
import sys
from pathlib import Path


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Rank protein designs and select top designs"
    )
    parser.add_argument(
        "--csv", 
        required=True, 
        help="Input CSV file with design metrics (best_designs.csv)"
    )
    parser.add_argument(
        "--pdb-dir", 
        required=True, 
        help="Directory containing PDB files"
    )
    parser.add_argument(
        "--output-csv",
        default="ranked_designs.csv",
        help="Output CSV file with rankings (default: ranked_designs.csv)"
    )
    parser.add_argument(
        "--output-dir",
        default="ranked_designs",
        help="Output directory for ranked PDB files (default: ranked_designs)"
    )
    parser.add_argument(
        "--ranking-metric",
        required=True,
        choices=["af2_pae_interaction", "boltz_ptm_interface"],
        help="Metric to use for ranking"
    )
    parser.add_argument(
        "--max-designs",
        type=int,
        default=None,
        help="Maximum number of top designs to output (default: all)"
    )
    
    return parser.parse_args()


def rank_designs(df, metric):
    """
    Rank designs based on the specified metric.
    
    Args:
        df: pandas DataFrame with design metrics
        metric: column name to rank by
    
    Returns:
        DataFrame sorted by metric with rank column added
    """
    # Check if metric exists in dataframe
    if metric not in df.columns:
        print(f"Error: Metric '{metric}' not found in CSV file.", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Remove rows with NaN values in the ranking metric
    initial_count = len(df)
    df = df.dropna(subset=[metric])
    removed_count = initial_count - len(df)
    
    if removed_count > 0:
        print(f"Removed {removed_count} designs with missing {metric} values")
    
    if len(df) == 0:
        print(f"Error: No designs have valid {metric} values", file=sys.stderr)
        sys.exit(1)
    
    # Sort by metric
    # For PAE metrics: lower is better (ascending=True)
    # For PTM metrics: higher is better (ascending=False)
    if "pae" in metric.lower():
        ascending = True
        print(f"Ranking by {metric} (lower is better)")
    else:  # ptm, plddt, confidence, etc.
        ascending = False
        print(f"Ranking by {metric} (higher is better)")
    
    df_sorted = df.sort_values(by=metric, ascending=ascending).reset_index(drop=True)
    
    # Add rank column (1-indexed)
    df_sorted.insert(0, 'rank', range(1, len(df_sorted) + 1))
    
    return df_sorted


def generate_pdb_filename_from_row(row):
    """
    Generate expected PDB filename from CSV row data.
    
    Args:
        row: pandas Series with fold_id and seq_id columns
    
    Returns:
        Expected PDB filename
    """
    fold_id = row.get('fold_id')
    seq_id = row.get('seq_id')
    
    # Handle different prediction methods
    if pd.notna(row.get('boltz_ptm')):
        # Boltz prediction
        return f"fold_{fold_id}_seq_{seq_id}_boltzpred.pdb"
    else:
        # AF2 prediction
        return f"fold_{fold_id}_seq_{seq_id}_af2pred.pdb"


def copy_and_rename_pdbs(df, pdb_dir, output_dir, max_designs=None):
    """
    Copy and rename PDB files with ranked prefixes.
    
    Args:
        df: DataFrame with rankings and design information
        pdb_dir: Directory containing original PDB files
        output_dir: Directory to output ranked PDB files
        max_designs: Maximum number of designs to copy (None for all)
    
    Returns:
        Number of successfully copied files
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Limit to top N designs if specified
    if max_designs is not None:
        df = df.head(max_designs)
        print(f"Selecting top {max_designs} designs")
    
    copied_count = 0
    missing_files = []
    
    for idx, row in df.iterrows():
        rank = row['rank']
        
        # Generate expected PDB filename
        original_filename = generate_pdb_filename_from_row(row)
        original_path = os.path.join(pdb_dir, original_filename)
        
        # Check if file exists
        if not os.path.exists(original_path):
            missing_files.append(original_filename)
            continue
        
        # Generate ranked filename with zero-padded prefix
        # Determine padding width based on total number of designs
        padding_width = len(str(len(df)))
        ranked_filename = f"{str(rank).zfill(padding_width)}_{original_filename}"
        ranked_path = os.path.join(output_dir, ranked_filename)
        
        # Copy file
        shutil.copy2(original_path, ranked_path)
        copied_count += 1
    
    if missing_files:
        print(f"\nWarning: {len(missing_files)} PDB files not found:", file=sys.stderr)
        for filename in missing_files[:5]:  # Show first 5
            print(f"  - {filename}", file=sys.stderr)
        if len(missing_files) > 5:
            print(f"  ... and {len(missing_files) - 5} more", file=sys.stderr)
    
    return copied_count


def main():
    args = parse_arguments()
    
    print(f"Reading design metrics from {args.csv}")
    
    # Read CSV
    try:
        df = pd.read_csv(args.csv)
    except Exception as e:
        print(f"Error reading CSV file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(df)} designs in CSV")
    
    # Rank designs
    df_ranked = rank_designs(df, args.ranking_metric)
    
    print(f"Ranked {len(df_ranked)} designs by {args.ranking_metric}")
    
    # Apply max_designs limit if specified
    output_count = len(df_ranked) if args.max_designs is None else min(len(df_ranked), args.max_designs)
    
    # Save ranked CSV
    output_csv_path = args.output_csv
    df_ranked.to_csv(output_csv_path, index=False)
    print(f"Saved ranked CSV to {output_csv_path}")
    
    # Copy and rename PDB files
    print(f"\nCopying and renaming PDB files to {args.output_dir}")
    copied_count = copy_and_rename_pdbs(
        df_ranked,
        args.pdb_dir,
        args.output_dir,
        args.max_designs
    )
    
    print(f"\nRanking complete:")
    print(f"  - Total designs ranked: {len(df_ranked)}")
    print(f"  - Designs output: {output_count}")
    print(f"  - PDB files copied: {copied_count}")
    
    if copied_count == 0:
        print("\nError: No PDB files were copied", file=sys.stderr)
        sys.exit(1)
    
    # Show top 5 designs
    print(f"\nTop 5 designs by {args.ranking_metric}:")
    print(df_ranked[['rank', 'fold_id', 'seq_id', args.ranking_metric]].head(5).to_string(index=False))


if __name__ == "__main__":
    main()
