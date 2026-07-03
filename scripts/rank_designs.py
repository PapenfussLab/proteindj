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
        help="Metric to use for ranking (e.g., 'af2_pae_interaction', 'boltz_ptm_interface', 'boltz_ipSAE_min')"
    )
    parser.add_argument(
        "--max-designs",
        type=int,
        default=None,
        help="Maximum number of top designs to output (default: all)"
    )
    parser.add_argument(
        "--max-seqs-per-fold",
        type=int,
        default=None,
        help="Maximum number of sequences per fold_id to keep (default: no limit)"
    )
    
    return parser.parse_args()


def rank_designs(df, metric, max_seqs_per_fold=None):
    """
    Rank designs based on the specified metric.
    
    Args:
        df: pandas DataFrame with design metrics
        metric: column name to rank by
        max_seqs_per_fold: maximum number of sequences to keep per fold_id (None for no limit)
    
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
    
    # Determine sort direction based on metric type
    # Lower is better for: PAE, RMSD, PDE
    # Higher is better for: PTM, pLDDT, confidence, ipSAE, LIS, pDockQ
    metric_lower = metric.lower()
    
    lower_is_better_keywords = ['pae', 'rmsd', 'pde']
    higher_is_better_keywords = ['ptm', 'plddt', 'conf', 'ipsae', 'lis', 'pdockq']
    
    # Check if metric contains keywords for direction
    is_lower_better = any(keyword in metric_lower for keyword in lower_is_better_keywords)
    is_higher_better = any(keyword in metric_lower for keyword in higher_is_better_keywords)
    
    # Handle special case where both might match (e.g., unsat_hbonds has 'hbond')
    # In this case, prioritize the lower_is_better keywords
    if is_lower_better:
        ascending = True
        print(f"Ranking by {metric} (lower is better)")
    elif is_higher_better:
        ascending = False
        print(f"Ranking by {metric} (higher is better)")
    else:
        # Default to higher is better if we can't determine
        ascending = False
        print(f"Warning: Could not determine direction for metric '{metric}', assuming higher is better")
        print(f"Ranking by {metric} (higher is better)")
    
    df_sorted = df.sort_values(by=metric, ascending=ascending).reset_index(drop=True)
    
    # Apply max_seqs_per_fold filter if specified
    if max_seqs_per_fold is not None:
        print(f"\nApplying max_seqs_per_fold filter: keeping top {max_seqs_per_fold} sequences per fold")
        
        # Check if fold_id column exists
        if 'fold_id' not in df_sorted.columns:
            print("Warning: fold_id column not found, skipping max_seqs_per_fold filter", file=sys.stderr)
        else:
            # Group by fold_id and keep top N sequences per fold
            df_filtered = df_sorted.groupby('fold_id', group_keys=False).apply(
                lambda x: x.head(max_seqs_per_fold)
            ).reset_index(drop=True)
            
            original_count = len(df_sorted)
            filtered_count = len(df_filtered)
            removed = original_count - filtered_count
            
            print(f"  - Original designs: {original_count}")
            print(f"  - Designs after filter: {filtered_count}")
            print(f"  - Designs removed: {removed}")
            
            # Show fold distribution
            fold_counts = df_filtered.groupby('fold_id').size()
            print(f"  - Number of unique folds: {len(fold_counts)}")
            print(f"  - Sequences per fold: min={fold_counts.min()}, max={fold_counts.max()}, mean={fold_counts.mean():.1f}")
            
            # Re-sort by metric to restore global ranking order after per-fold filtering
            df_sorted = df_filtered.sort_values(by=metric, ascending=ascending).reset_index(drop=True)
    
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
    fold_id = int(row.get('fold_id'))
    seq_id = int(row.get('seq_id'))
    
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
    df_ranked = rank_designs(df, args.ranking_metric, args.max_seqs_per_fold)
    
    print(f"\nRanked {len(df_ranked)} designs by {args.ranking_metric}")
    
    # Apply max_designs as the global cap — takes priority over max_seqs_per_fold
    if args.max_designs is not None and len(df_ranked) > args.max_designs:
        print(f"Applying max_designs limit: selecting top {args.max_designs} of {len(df_ranked)} designs")
        df_ranked = df_ranked.head(args.max_designs)

    output_count = len(df_ranked)

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
    
    # Show top designs
    display_count = min(5, len(df_ranked))
    print(f"\nTop {display_count} designs by {args.ranking_metric}:")
    print(df_ranked[['rank', 'fold_id', 'seq_id', args.ranking_metric]].head(display_count).to_string(index=False))


if __name__ == "__main__":
    main()
