#!/usr/bin/env python3
"""
filter_best_designs.py - Filter CSV based on fold_id and seq_id from PDB files

This script filters a CSV file containing protein design metrics to include only
entries that correspond to final PDB files. It matches based on fold_id and seq_id
values extracted from PDB filenames in a specified directory, or just fold_id if
seq_id column is not present in the CSV.
"""

import argparse
import os
import shutil
import re
import pandas as pd
import sys


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Filter CSV based on fold_id and seq_id from PDB files in a directory"
    )
    parser.add_argument(
        "--csv", 
        required=True, 
        help="Input CSV file with design metrics"
    )
    parser.add_argument(
        "--pdb-dir", 
        required=True, 
        help="Directory containing PDB files to filter by"
    )
    parser.add_argument(
        "--output-csv", 
        default="best_designs.csv", 
        help="Output CSV file (default: best_designs.csv)"
    )
    parser.add_argument(
        "--output-dir", 
        default="best_designs", 
        help="Directory to copy PDB files to (default: best_designs)"
    )
    return parser.parse_args()


def extract_ids_from_pdb(pdb_filename):
    """Extract fold_id and seq_id from PDB filename."""
    # First try to match fold_id and seq_id pattern
    match = re.search(r'fold_(\d+)_seq_(\d+)', os.path.basename(pdb_filename))
    if match:
        fold_id = int(match.group(1))
        seq_id = int(match.group(2))
        return fold_id, seq_id
    
    # If that fails, try to match just fold_id pattern
    match = re.search(r'fold_(\d+)(?!.*seq)', os.path.basename(pdb_filename))
    if match:
        fold_id = int(match.group(1))
        return fold_id, None
    
    return None, None


def main():
    """Main function to filter CSV and copy PDB files."""
    args = parse_arguments()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read the full CSV file
    try:
        df = pd.read_csv(args.csv)
        print(f"Read CSV with {len(df)} entries")
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return 1
    
    # Check if seq_id column exists in the CSV
    has_seq_id = 'seq_id' in df.columns
    print(f"CSV has seq_id column: {has_seq_id}")
    
    # Extract fold_id and seq_id from PDB filenames in the input directory
    pdb_info = []
    mismatch_errors = []
    
    for filename in os.listdir(args.pdb_dir):
        if filename.endswith('.pdb'):
            pdb_path = os.path.join(args.pdb_dir, filename)
            fold_id, seq_id = extract_ids_from_pdb(filename)
            
            if fold_id is not None:
                # Check for mismatch between CSV and PDB format
                if has_seq_id and seq_id is None:
                    mismatch_errors.append(f"CSV has seq_id column but PDB file '{filename}' doesn't have seq_id in filename")
                elif not has_seq_id and seq_id is not None:
                    mismatch_errors.append(f"CSV lacks seq_id column but PDB file '{filename}' has seq_id in filename")
                else:
                    pdb_info.append((fold_id, seq_id, pdb_path))
            else:
                print(f"Warning: Could not extract fold_id from {filename}")
    
    # If there are any mismatches, report them and exit with error
    if mismatch_errors:
        print("\nERROR: Mismatch between CSV format and PDB filename format:")
        for error in mismatch_errors:
            print(f"  - {error}")
        print("\nPlease ensure that:")
        print("  - If CSV has seq_id column, all PDB files should follow 'fold_X_seq_Y.pdb' format")
        print("  - If CSV lacks seq_id column, all PDB files should follow 'fold_X.pdb' format")
        return 1
    
    print(f"Found {len(pdb_info)} valid PDB files")
    
    # Filter the dataframe based on available columns and PDB file patterns
    filtered_rows = []
    
    for fold_id, seq_id, pdb_path in pdb_info:
        if has_seq_id and seq_id is not None:
            # Both CSV and PDB have seq_id, match on both
            matching_rows = df[(df['fold_id'] == fold_id) & (df['seq_id'] == seq_id)]
            if not matching_rows.empty:
                filtered_rows.append(matching_rows)
            else:
                print(f"Warning: No matching entry found for fold_id={fold_id}, seq_id={seq_id}")
        
        elif not has_seq_id and seq_id is None:
            # Neither CSV nor PDB have seq_id, match only on fold_id
            matching_rows = df[df['fold_id'] == fold_id]
            if not matching_rows.empty:
                filtered_rows.append(matching_rows)
            else:
                print(f"Warning: No matching entry found for fold_id={fold_id}")
    
    # Combine all matching rows
    filtered_df = pd.concat(filtered_rows)
    print(f"Found {len(filtered_df)} matching entries in CSV")
    # Save the filtered dataframe to a new CSV file
    filtered_df.to_csv(args.output_csv, index=False)
    print(f"Saved filtered CSV to {args.output_csv}")
    
    # Copy the PDB files to the output directory
    for fold_id, seq_id, pdb_path in pdb_info:
        dest = os.path.join(args.output_dir, os.path.basename(pdb_path))
        shutil.copy(pdb_path, dest)
        print(f"Copied {os.path.basename(pdb_path)} to {args.output_dir}")
    
    print("Done!")
    return 0


if __name__ == "__main__":
    exit(main())
