#!/usr/bin/env python3

import argparse
import sys
import os
import json
import re
from pathlib import Path
import shutil
import glob

def parse_arguments():
    parser = argparse.ArgumentParser(description='Filter designs and copy top hits from PDB files.')
    parser.add_argument('--json-directory', required=True, help='Path to directory containing JSON score files')
    parser.add_argument('--boltz-max-overall-rmsd', type=float,
                    help='Maximum allowed RMSD value')
    parser.add_argument('--boltz-max-binder-rmsd', type=float,
                        help='Maximum allowed binder chain RMSD value')
    parser.add_argument('--boltz-max-target-rmsd', type=float,
                        help='Maximum allowed target chain RMSD value')
    parser.add_argument('--boltz-min-conf-score', type=float, 
                    help='Minimum confidence score')
    parser.add_argument('--boltz-min-ptm', type=float, 
                    help='Minimum pTM score')
    parser.add_argument('--boltz-min-ptm-interface', type=float,
                    help='Minimum interface pTM score')
    parser.add_argument('--boltz-min-plddt', type=float,
                    help='Minimum complex pLDDT score')
    parser.add_argument('--boltz-min-plddt-interface', type=float,
                    help='Minimum interface-weighted pLDDT')
    parser.add_argument('--boltz-max-pde', type=float,
                    help='Maximum complex PDE score')
    parser.add_argument('--boltz-max-pde-interface', type=float,
                    help='Maximum interface-weighted PDE')
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
                # Extract metadata from filename
                filename = Path(json_file).stem
                match = re.match(r'fold_(\d+)_seq_(\d+)_boltzpred', filename)
                if not match:
                    print(f"Skipping invalid filename format: {filename}")
                    continue
                    
                fold_id = int(match.group(1))
                seq_id = int(match.group(2))
                
                with open(json_file, 'r') as f:
                    file_content = f.read().strip()
                    if not file_content:
                        print(f"Skipping empty file: {json_file}")
                        continue

                    # Add filename metadata to all entries
                    if file_content.startswith('['):  # JSON array
                        file_data = json.loads(file_content)
                        for entry in file_data:
                            entry.update({
                                "description": filename,
                                "fold_id": fold_id,
                                "seq_id": seq_id
                            })
                        data.extend(file_data)
                    elif file_content.startswith('{'):  # Single JSON object
                        entry = json.loads(file_content)
                        entry.update({
                            "description": filename,
                            "fold_id": fold_id,
                            "seq_id": seq_id
                        })
                        data.append(entry)
                    else:  # JSONL format
                        for line in file_content.splitlines():
                            if line.strip():
                                entry = json.loads(line)
                                entry.update({
                                    "description": filename,
                                    "fold_id": fold_id,
                                    "seq_id": seq_id
                                })
                                data.append(entry)
                                
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
            failures = []
            
            # Overall RMSD check
            if args.boltz_max_overall_rmsd:
                rmsd = entry.get('boltz_overall_rmsd', 1000)
                if rmsd > args.boltz_max_overall_rmsd:
                    failures.append(f"overall_rmsd {rmsd:.2f} > {args.boltz_max_overall_rmsd}")
            
            # Binder RMSD check
            if args.boltz_max_binder_rmsd:
                binder_rmsd = entry.get('boltz_binder_rmsd', 1000)
                if binder_rmsd > args.boltz_max_binder_rmsd:
                    failures.append(f"binder_rmsd {binder_rmsd:.2f} > {args.boltz_max_binder_rmsd}")
            
            # Target RMSD check
            if args.boltz_max_target_rmsd:
                target_rmsd = entry.get('boltz_target_rmsd', 1000)
                if target_rmsd > args.boltz_max_target_rmsd:
                    failures.append(f"target_rmsd {target_rmsd:.2f} > {args.boltz_max_target_rmsd}")

            # Confidence score check
            if args.boltz_min_conf_score:
                score = entry.get('boltz_conf_score', 0)
                if score < args.boltz_min_conf_score:
                    failures.append(f"conf_score {score:.3f} < {args.boltz_min_conf_score}")
            
            # PTM/iPTM checks
            if args.boltz_min_ptm:
                ptm = entry.get('boltz_ptm', 0)
                if ptm < args.boltz_min_ptm:
                    failures.append(f"ptm {ptm:.3f} < {args.boltz_min_ptm}")
            
            if args.boltz_min_ptm_interface:
                ptm_interface = entry.get('boltz_ptm_interface', 0)
                if ptm_interface < args.boltz_min_ptm_interface:
                    failures.append(f"ptm_interface {ptm_interface:.3f} < {args.boltz_min_ptm_interface}")
            
            # Structure quality checks
            if args.boltz_min_plddt:
                plddt = entry.get('boltz_plddt', 0)
                if plddt < args.boltz_min_plddt:
                    failures.append(f"plddt {plddt:.3f} < {args.boltz_min_plddt}")
            
            if args.boltz_min_plddt_interface:
                plddt_interface = entry.get('boltz_plddt_interface', 0)
                if plddt_interface < args.boltz_min_plddt_interface:
                    failures.append(f"plddt_interface {plddt_interface:.3f} < {args.boltz_min_plddt_interface}")
            
            # Energy metric checks
            if args.boltz_max_pde:
                pde = entry.get('boltz_pde', 0)
                if pde > args.boltz_max_pde:
                    failures.append(f"pde {pde:.2f} > {args.boltz_max_pde}")
            
            if args.boltz_max_pde_interface:
                pde_interface = entry.get('boltz_pde_interface', 0)
                if pde_interface > args.boltz_max_pde_interface:
                    failures.append(f"pde_interface {pde_interface:.2f} > {args.boltz_max_pde_interface}")

            if not failures:
                passed_designs.append(entry['description'])
                passed_entries.append(entry)
            else:
                print(f"Rejected {entry['description']}: {', '.join(failures)}")
                
        except KeyError as e:
            print(f"Invalid entry missing key {e}")
    
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
        args.boltz_max_overall_rmsd is not None or
        args.boltz_max_binder_rmsd is not None or
        args.boltz_max_target_rmsd is not None or
        args.boltz_min_conf_score is not None or
        args.boltz_min_ptm is not None or
        args.boltz_min_ptm_interface is not None or
        args.boltz_min_plddt is not None or
        args.boltz_min_plddt_interface is not None or
        args.boltz_max_pde is not None or
        args.boltz_max_pde_interface is not None
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
