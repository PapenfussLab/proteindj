#!/usr/bin/env python3

import argparse
import sys
import os
import json
from pathlib import Path
import shutil


def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter designs based on biophysical and sequence analysis metrics.')
    
    # Biophysical Analysis Metrics
    parser.add_argument('--pr-min-helices', type=int, help='Minimum number of alpha-helices')
    parser.add_argument('--pr-max-helices', type=int, help='Maximum number of alpha-helices')
    parser.add_argument('--pr-min-strands', type=int, help='Minimum number of beta-strands')
    parser.add_argument('--pr-max-strands', type=int, help='Maximum number of beta-strands')
    parser.add_argument('--pr-min-total-ss', type=int, help='Minimum total secondary structures')
    parser.add_argument('--pr-max-total-ss', type=int, help='Maximum total secondary structures')
    parser.add_argument('--pr-min-rog', type=float, help='Minimum radius of gyration (Å)')
    parser.add_argument('--pr-max-rog', type=float, help='Maximum radius of gyration (Å)')
    parser.add_argument('--pr-min-intface-bsa', type=float, help='Minimum buried surface area (Å²)')
    parser.add_argument('--pr-min-intface-shpcomp', type=float, help='Minimum interface shape complementarity')
    parser.add_argument('--pr-min-intface-hbonds', type=int, help='Minimum interface hydrogen bonds')
    parser.add_argument('--pr-max-intface-unsat-hbonds', type=int, help='Maximum unsatisfied hydrogen bonds')
    parser.add_argument('--pr-max-intface-deltag', type=float, help='Maximum interface deltaG (REU)')
    parser.add_argument('--pr-max-intface-deltagtobsa', type=float, help='Maximum interface deltaGtoBSA ratio')
    parser.add_argument('--pr-min-intface-packstat', type=float, help='Minimum interface packing statistic')
    parser.add_argument('--pr-max-tem', type=float, help='Maximum total energy metric (REU)')
    parser.add_argument('--pr-max-surfhphobics', type=float, help='Maximum surface hydrophobics percentage')
    
    # Sequence Analysis Metrics
    parser.add_argument('--seq-min-ext-coef', type=float, help='Minimum extinction coefficient')
    parser.add_argument('--seq-max-ext-coef', type=float, help='Maximum extinction coefficient')
    parser.add_argument('--seq-min-pi', type=float, help='Minimum isoelectric point')
    parser.add_argument('--seq-max-pi', type=float, help='Maximum isoelectric point')
    
    # Input/Output
    parser.add_argument('--jsonl-file', required=True, help='Path to JSONL file containing analysis data')
    parser.add_argument('--pdb-directory', required=True, help='Directory containing PDB files')
    parser.add_argument('--output-directory', default='output', help='Directory to copy passing PDB files to')
    parser.add_argument('--output-score-file', default='filtered.jsonl', help='Output JSONL file for filtered results')
    parser.add_argument('--num-to-extract', type=int, help='Number of designs to extract (extracts all if not specified)')
    
    return parser.parse_args()


def read_jsonl(file_path):
    """
    Read JSONL file and return list of entries.
    
    Args:
        file_path: Path to JSONL file
        
    Returns:
        List of data entries
    """
    try:
        data = []
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue
                try:
                    data.append(json.loads(line))
                except json.JSONDecodeError as e:
                    print(f"Error parsing JSON on line {line_num}: {e}", file=sys.stderr)
                    
        print(f"Loaded {len(data)} entries from {file_path}")
        return data
    except Exception as e:
        print(f"Error reading JSONL file: {e}", file=sys.stderr)
        sys.exit(1)


def passes_filter(value, min_val, max_val, metric_name):
    """
    Check if value passes min/max filter criteria.
    
    Args:
        value: Value to check
        min_val: Minimum threshold (None if not set)
        max_val: Maximum threshold (None if not set)
        metric_name: Name of metric for error reporting
        
    Returns:
        Tuple of (passes: bool, failure_reason: str or None)
    """
    if min_val is not None and max_val is not None:
        if not (min_val <= value <= max_val):
            return False, f"{metric_name} {value} not in range [{min_val}, {max_val}]"
    elif min_val is not None:
        if value < min_val:
            return False, f"{metric_name} {value} < {min_val}"
    elif max_val is not None:
        if value > max_val:
            return False, f"{metric_name} {value} > {max_val}"
    
    return True, None


def filter_data(data, args):
    """
    Filter data based on provided criteria.
    
    Args:
        data: List of data entries
        args: Parsed command line arguments
        
    Returns:
        Tuple of (passed_designs, passed_entries)
    """
    passed_designs = []
    passed_entries = []
    
    # Define all filters to check
    filters = [
        # Biophysical metrics
        ('pr_helices', args.pr_min_helices, args.pr_max_helices, int),
        ('pr_strands', args.pr_min_strands, args.pr_max_strands, int),
        ('pr_total_ss', args.pr_min_total_ss, args.pr_max_total_ss, int),
        ('pr_RoG', args.pr_min_rog, args.pr_max_rog, float),
        ('pr_intface_BSA', args.pr_min_intface_bsa, None, float),
        ('pr_intface_shpcomp', args.pr_min_intface_shpcomp, None, float),
        ('pr_intface_hbonds', args.pr_min_intface_hbonds, None, int),
        ('pr_intface_unsat_hbonds', None, args.pr_max_intface_unsat_hbonds, int),
        ('pr_intface_deltaG', None, args.pr_max_intface_deltag, float),
        ('pr_intface_deltaGtoBSA', None, args.pr_max_intface_deltagtobsa, float),
        ('pr_intface_packstat', args.pr_min_intface_packstat, None, float),
        ('pr_TEM', None, args.pr_max_tem, float),
        ('pr_surfhphobics', None, args.pr_max_surfhphobics, float),
        # Sequence metrics
        ('seq_ext_coef', args.seq_min_ext_coef, args.seq_max_ext_coef, float),
        ('seq_pI', args.seq_min_pi, args.seq_max_pi, float),
    ]
    
    for entry in data:
        try:
            description = entry.get('description')
            if not description:
                print(f"design rejected: missing description")
                continue
            
            failures = []
            
            # Check each filter
            for metric_key, min_val, max_val, value_type in filters:
                # Skip if no filter is set for this metric
                if min_val is None and max_val is None:
                    continue
                
                # Check if metric exists in entry
                if metric_key not in entry:
                    failures.append(f"missing metric '{metric_key}'")
                    continue
                
                try:
                    # Convert to appropriate type
                    value = value_type(entry[metric_key])
                    
                    # Check filter
                    passes, reason = passes_filter(value, min_val, max_val, metric_key)
                    if not passes:
                        failures.append(reason)
                        
                except (ValueError, TypeError) as e:
                    failures.append(f"invalid {metric_key} value: {entry[metric_key]}")
            
            # If no failures, accept the design
            if not failures:
                print(f"design '{description}' accepted")
                passed_designs.append(description)
                passed_entries.append(entry)
            else:
                failure_reasons = ', '.join(failures)
                print(f"design '{description}' rejected: {failure_reasons}")
                
        except Exception as e:
            print(f"design rejected: error processing entry - {e}")
            continue
    
    return passed_designs, passed_entries


def copy_pdb_files(passed_designs, pdb_directory, output_dir):
    """
    Copy PDB files from the PDB directory to the output directory.
    
    Args:
        passed_designs: List of design names that passed filtering
        pdb_directory: Source directory containing PDB files
        output_dir: Destination directory to copy files to
        
    Returns:
        List of successfully copied design names
    """
    successfully_copied = []
    pdb_dir_path = Path(pdb_directory)
    
    for design_name in passed_designs:
        # Append .pdb to the description to get the filename
        source_path = pdb_dir_path / f"{design_name}.pdb"
        
        if source_path.exists():
            dest_path = Path(output_dir) / source_path.name
            shutil.copy2(source_path, dest_path)
            print(f"Copied {source_path} to {dest_path}")
            successfully_copied.append(design_name)
        else:
            print(f"Warning: Could not find PDB file for design '{design_name}' at {source_path}")
    
    return successfully_copied


def write_filtered_jsonl(passed_entries, successfully_copied, output_score_file):
    """
    Write filtered entries to JSONL file.
    
    Args:
        passed_entries: List of entries that passed filtering
        successfully_copied: List of design names that were successfully copied
        output_score_file: Path to output JSONL file
    """
    with open(output_score_file, 'w') as outfile:
        for entry in passed_entries:
            if entry['description'] in successfully_copied:
                outfile.write(json.dumps(entry) + '\n')
    
    print(f"Filtered score file created: {output_score_file}")


def main():
    args = parse_arguments()
    
    # Create output directory
    output_dir = Path(args.output_directory)
    output_dir.mkdir(exist_ok=True)
    
    # Read JSONL file
    data = read_jsonl(args.jsonl_file)
    
    if not data:
        print("No data found in the JSONL file. Exiting.")
        sys.exit(0)
    
    # Check if any filter is applied
    filter_args = [
        args.pr_min_helices, args.pr_max_helices,
        args.pr_min_strands, args.pr_max_strands,
        args.pr_min_total_ss, args.pr_max_total_ss,
        args.pr_min_rog, args.pr_max_rog,
        args.pr_min_intface_bsa,
        args.pr_min_intface_shpcomp,
        args.pr_min_intface_hbonds,
        args.pr_max_intface_unsat_hbonds,
        args.pr_max_intface_deltag,
        args.pr_max_intface_deltagtobsa,
        args.pr_min_intface_packstat,
        args.pr_max_tem,
        args.pr_max_surfhphobics,
        args.seq_min_ext_coef, args.seq_max_ext_coef,
        args.seq_min_pi, args.seq_max_pi,
    ]
    
    any_filter_applied = any(f is not None for f in filter_args)
    
    if not any_filter_applied:
        print("Warning: No filter parameters specified. All designs will pass.")
    
    # Filter data
    passed_designs, passed_entries = filter_data(data, args)
    
    if not passed_designs:
        print("No designs passed filtering.")
        sys.exit(0)
    
    # Limit number of designs if specified
    if args.num_to_extract is not None:
        if args.num_to_extract < len(passed_designs):
            print(f"Limiting to top {args.num_to_extract} designs out of {len(passed_designs)} that passed filters")
            passed_designs = passed_designs[:args.num_to_extract]
            passed_entries = passed_entries[:args.num_to_extract]
        else:
            print(f"Requested {args.num_to_extract} designs but only {len(passed_designs)} passed filters")
    
    # Copy PDB files
    print(f"\nCopying PDB files...")
    print(f"Total designs to process: {len(passed_designs)}")
    
    successfully_copied_designs = copy_pdb_files(passed_designs, args.pdb_directory, args.output_directory)
    
    # Print summary
    print(f"\nFile Copy Summary:")
    print(f"Total designs that passed filtering: {len(passed_designs)}")
    print(f"Successfully copied PDB files: {len(successfully_copied_designs)}")
    
    # Create filtered score file
    write_filtered_jsonl(passed_entries, successfully_copied_designs, args.output_score_file)


if __name__ == '__main__':
    main()
