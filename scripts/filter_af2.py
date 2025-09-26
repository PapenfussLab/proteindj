#!/usr/bin/env python3

import argparse
import sys
import os
import json
from pathlib import Path
import shutil
import glob

def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter designs and copy top hits from PDB files.')
    parser.add_argument('--json-directory', required=True, help='Path to directory containing JSON score files')
    parser.add_argument('--af2-max-pae-interaction', type=float, help='Maximum PAE interaction score')
    parser.add_argument('--af2-max-pae-overall', type=float, help='Maximum PAE overall score')
    parser.add_argument('--af2-max-pae-binder', type=float, help='Maximum PAE binder score')
    parser.add_argument('--af2-max-pae-target', type=float, help='Maximum PAE target score')
    parser.add_argument('--af2-min-plddt-overall', type=float, help='Minimum overall pLDDT score')
    parser.add_argument('--af2-min-plddt-binder', type=float, help='Minimum binder pLDDT score')
    parser.add_argument('--af2-min-plddt-target', type=float, help='Minimum target pLDDT score')
    parser.add_argument('--af2-max-rmsd-overall', type=float, help='Max overall RMSD when all chains are aligned')
    parser.add_argument('--af2-max-rmsd-binder-bndaln', type=float, help='Max binder RMSD when binder chains are aligned')
    parser.add_argument('--af2-max-rmsd-binder-tgtaln', type=float, help='Max binder RMSD when target chains are aligned')
    parser.add_argument('--af2-max-rmsd-target', type=float, help='Max target RMSD when target chains are aligned')
    parser.add_argument('--output-directory', default='output', help='Directory to copy passing PDB files to')
    parser.add_argument('--output-score-file', default='filtered.jsonl')
    parser.add_argument('--num-to-extract', type=int, help='Number of designs to extract (extracts all if not specified)')
    parser.add_argument('--json-pattern', default='*.json', help='Pattern to match JSON files in the directory')
    return parser.parse_args()

def read_data_from_directory(directory_path, pattern='*.json'):
    """
    Read all JSON files from a directory.
    
    Args:
        directory_path: Path to directory containing JSON files
        pattern: Glob pattern to match JSON files
        
    Returns:
        List of data entries from all JSON files
    """
    try:
        data = []
        json_files = glob.glob(os.path.join(directory_path, pattern))
        
        if not json_files:
            print(f"No JSON files found in {directory_path} matching pattern {pattern}", file=sys.stderr)
            sys.exit(1)
            
        print(f"Found {len(json_files)} JSON files to process")
        
        for json_file in json_files:
            try:
                with open(json_file, 'r') as f:
                    file_content = f.read().strip()
                    if not file_content:
                        print(f"Skipping empty file: {json_file}")
                        continue
                        
                    # Check if the file contains one JSON object per line (JSONL)
                    if file_content.startswith('[') or file_content.startswith('{'):
                        # Single JSON object or array
                        file_data = json.loads(file_content)
                        if isinstance(file_data, list):
                            data.extend(file_data)
                        else:
                            data.append(file_data)
                    else:
                        # JSONL format (one JSON object per line)
                        for line in file_content.splitlines():
                            if line.strip():
                                data.append(json.loads(line))
                                
                print(f"Processed {json_file}")
            except json.JSONDecodeError as e:
                print(f"Error parsing JSON in file {json_file}: {e}", file=sys.stderr)
            except Exception as e:
                print(f"Error processing file {json_file}: {e}", file=sys.stderr)
                
        return data
    except Exception as e:
        print(f"Error reading JSON files from directory: {e}", file=sys.stderr)
        sys.exit(1)

def filter_data(data, args):
    passed_designs = []
    passed_entries = []
    
    for entry in data:
        try:
            description = entry.get('description')
            if not description:
                print(f"design rejected: missing description")
                continue

            # Check each criterion individually and print specific failures
            failures = []
            
            # PAE interaction check
            if args.af2_max_pae_interaction is not None:
                pae_interaction = float(entry['af2_pae_interaction'])
                if pae_interaction > args.af2_max_pae_interaction:
                    failures.append(f"pae_interaction {pae_interaction} > {args.af2_max_pae_interaction}")
            
            # PAE overall check
            if args.af2_max_pae_overall is not None:
                pae_overall = float(entry['af2_pae_overall'])
                if pae_overall > args.af2_max_pae_overall:
                    failures.append(f"pae_overall {pae_overall} > {args.af2_max_pae_overall}")

            # PAE binder check
            if args.af2_max_pae_binder is not None:
                pae_binder = float(entry['af2_pae_binder'])
                if pae_binder > args.af2_max_pae_binder:
                    failures.append(f"pae_binder {pae_binder} > {args.af2_max_pae_binder}")
            
            # PAE target check
            if args.af2_max_pae_target is not None:
                pae_target = float(entry['af2_pae_target'])
                if pae_target > args.af2_max_pae_target:
                    failures.append(f"pae_target {pae_target} > {args.af2_max_pae_target}")

            # Overall RMSD check
            if args.af2_max_rmsd_overall is not None:
                rmsd_overall = float(entry['af2_rmsd_overall'])
                if rmsd_overall > args.af2_max_rmsd_overall:
                    failures.append(f"rmsd_overall {rmsd_overall} > {args.af2_max_rmsd_overall}")
                        
            # Binder Binder-Aligned RMSD check
            if args.af2_max_rmsd_binder_bndaln is not None:
                rmsd_binder_bndaln = float(entry['af2_rmsd_binder_bndaln'])
                if rmsd_binder_bndaln > args.af2_max_rmsd_binder_bndaln:
                    failures.append(f"rmsd_binder_bndaln {rmsd_binder_bndaln} > {args.af2_max_rmsd_binder_bndaln}")
            
            # Binder Target-Aligned RMSD check
            if args.af2_max_rmsd_binder_tgtaln is not None:
                rmsd_binder_tgtaln = float(entry['af2_rmsd_binder_tgtaln'])
                if rmsd_binder_tgtaln > args.af2_max_rmsd_binder_tgtaln:
                    failures.append(f"rmsd_binder_tgtaln {rmsd_binder_tgtaln} > {args.af2_max_rmsd_binder_tgtaln}")
            
            # Target RMSD check
            if args.af2_max_rmsd_target is not None:
                rmsd_target = float(entry['af2_rmsd_target'])
                if rmsd_target > args.af2_max_rmsd_target:
                    failures.append(f"rmsd_target {rmsd_target} > {args.af2_max_rmsd_target}")

            # pLDDT overall check
            if args.af2_min_plddt_overall is not None:
                plddt_overall = float(entry['af2_plddt_overall'])
                if plddt_overall < args.af2_min_plddt_overall:
                    failures.append(f"plddt_overall {plddt_overall} < {args.af2_min_plddt_overall}")
            
            # pLDDT binder check
            if args.af2_min_plddt_binder is not None:
                plddt_binder = float(entry['af2_plddt_binder'])
                if plddt_binder < args.af2_min_plddt_binder:
                    failures.append(f"plddt_binder {plddt_binder} < {args.af2_min_plddt_binder}")
            
            # pLDDT target check
            if args.af2_min_plddt_target is not None:
                plddt_target = float(entry['af2_plddt_target'])
                if plddt_target < args.af2_min_plddt_target:
                    failures.append(f"plddt_target {plddt_target} < {args.af2_min_plddt_target}")

            # If no failures, accept the design
            if not failures:
                print(f"design '{description}' accepted")
                passed_designs.append(description)
                passed_entries.append(entry)
            else:
                failure_reasons = ', '.join(failures)
                print(f"design '{description}' rejected: {failure_reasons}")
        except (ValueError, KeyError) as e:
            print(f"design rejected: invalid data - {e}")
            continue
    
    return passed_designs, passed_entries

def copy_pdb_files(passed_designs, output_dir):
    """
    Copy PDB files from the current directory to the output directory.
    
    Args:
        passed_designs: List of design names that passed filtering
        output_dir: Directory to copy files to
        
    Returns:
        List of successfully copied design names
    """
    successfully_copied = []
    
    for design_name in passed_designs:
        # Append .pdb to the description to get the filename
        source_path = Path(f"{design_name}.pdb")
        
        if source_path.exists():
            dest_path = Path(output_dir) / source_path.name
            shutil.copy2(source_path, dest_path)
            print(f"Copied {source_path} to {dest_path}")
            successfully_copied.append(design_name)
        else:
            print(f"Warning: Could not find PDB file for design '{design_name}'")
    
    return successfully_copied

def write_filtered_score_file(passed_entries, successfully_copied, output_score_file):
    with open(output_score_file, 'w') as outfile:
        for entry in passed_entries:
            if entry['description'] in successfully_copied:
                outfile.write(json.dumps(entry) + '\n')
    
    print(f"filtered score file created: {output_score_file}")

def main():
    args = parse_arguments()

    # Create output directory
    output_dir = Path(args.output_directory)
    output_dir.mkdir(exist_ok=True)
    
    # Read and filter data from directory
    data = read_data_from_directory(args.json_directory, args.json_pattern)
    
    if not data:
        print("No data found in the JSON files. Exiting.")
        sys.exit(0)
        
    print(f"Total entries loaded from all JSON files: {len(data)}")
    
    # Check if any filter is applied
    any_filter_applied = (
        args.af2_max_pae_interaction is not None or
        args.af2_max_pae_overall is not None or
        args.af2_max_pae_binder is not None or
        args.af2_max_pae_target is not None or
        args.af2_max_rmsd_overall is not None or
        args.af2_max_rmsd_binder_bndaln is not None or
        args.af2_max_rmsd_binder_tgtaln is not None or
        args.af2_max_rmsd_target is not None or
        args.af2_min_plddt_overall is not None or
        args.af2_min_plddt_binder is not None or
        args.af2_min_plddt_target is not None
    )
    
    if not any_filter_applied:
        print("Warning: No filter parameters specified. All designs will pass.")
    
    passed_designs, passed_entries = filter_data(data, args)

    if not passed_designs:
        print("no designs passed filtering.")
        sys.exit(0)
        
    if args.num_to_extract is not None:
        if args.num_to_extract < len(passed_designs):
            print(f"Limiting to top {args.num_to_extract} designs out of {len(passed_designs)} that passed filters")
            passed_designs = passed_designs[:args.num_to_extract]
            passed_entries = passed_entries[:args.num_to_extract]
        else:
            print(f"Requested {args.num_to_extract} designs but only {len(passed_designs)} passed filters")

    # Copy PDB files sequentially
    print(f"\nCopying PDB files...")
    print(f"Total designs to process: {len(passed_designs)}")
    
    successfully_copied_designs = copy_pdb_files(passed_designs, args.output_directory)
    
    # Print summary
    print(f"\nFile Copy Summary:")
    print(f"Total designs that passed filtering: {len(passed_designs)}")
    print(f"Successfully copied PDB files: {len(successfully_copied_designs)}")
    
    # Create filtered score file
    write_filtered_score_file(passed_entries, successfully_copied_designs, args.output_score_file)

if __name__ == '__main__':
    main()
