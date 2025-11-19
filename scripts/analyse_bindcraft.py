#!/usr/bin/env python3
"""
BindCraft Design Processor
Processes BindCraft batch output PDB and CSV files for downstream compatibility.
- Handles multiple designs per batch CSV file
- Swaps chains: A→B, B→A in PDB files
- Ensures chain order is A then B
- Renumbers chain B residues to avoid overlap with chain A
- Extracts and transforms metadata from CSV to JSON
- Generates RFdiffusion inpaint_seq mask array
- Assigns sequential fold_id's starting from 0
"""

import argparse
import csv
import json
import re
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain


def parse_trajectory_time(time_str):
    """
    Convert TrajectoryTime string to seconds.
    Expected format: "X hours, Y minutes, Z seconds"
    
    Args:
        time_str: Time string from CSV
    
    Returns:
        int: Total time in seconds
    """
    total_seconds = 0
    
    # Extract hours
    hours_match = re.search(r'(\d+)\s*hour', time_str)
    if hours_match:
        total_seconds += int(hours_match.group(1)) * 3600
    
    # Extract minutes
    minutes_match = re.search(r'(\d+)\s*minute', time_str)
    if minutes_match:
        total_seconds += int(minutes_match.group(1)) * 60
    
    # Extract seconds
    seconds_match = re.search(r'(\d+)\s*second', time_str)
    if seconds_match:
        total_seconds += int(seconds_match.group(1))
    
    return total_seconds


def get_chain_length(chain):
    """
    Get the number of residues in a chain.
    
    Args:
        chain: BioPython Chain object
    
    Returns:
        int: Number of residues
    """
    return len([res for res in chain if res.id[0] == ' '])  # Only standard residues


def swap_and_renumber_chains(input_pdb_path, output_pdb_path):
    """
    Swap chains (A→B, B→A), ensure order is A then B, and renumber chain B
    to start after chain A ends.
    
    Args:
        input_pdb_path: Path to input PDB file
        output_pdb_path: Path to output PDB file
    
    Returns:
        tuple: (chain_a_length, chain_b_length) - lengths of the output chains
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb_path)
    
    # Create new structure with swapped and renumbered chains
    new_structure = Structure.Structure('protein')
    new_model = Model.Model(0)
    new_structure.add(new_model)
    
    original_chain_a = None
    original_chain_b = None
    
    # Get original chains
    for model in structure:
        for chain in model:
            if chain.id == 'A':
                original_chain_a = chain
            elif chain.id == 'B':
                original_chain_b = chain
    
    if not original_chain_a or not original_chain_b:
        raise ValueError("Input PDB must contain both chain A and chain B")
    
    # Original chain B becomes new chain A (binder becomes chain A)
    new_chain_a = Chain.Chain('A')
    res_num = 1
    for residue in original_chain_b:
        new_res = residue.copy()
        hetfield, _, insertion = new_res.id
        new_res.id = (hetfield, res_num, insertion)
        new_chain_a.add(new_res)
        if hetfield == ' ':  # Only increment for standard residues
            res_num += 1
    
    new_model.add(new_chain_a)
    chain_a_length = get_chain_length(new_chain_a)
    
    # Original chain A becomes new chain B (target becomes chain B)
    # Chain B starts numbering after chain A ends
    new_chain_b = Chain.Chain('B')
    res_num = chain_a_length + 1
    for residue in original_chain_a:
        new_res = residue.copy()
        hetfield, _, insertion = new_res.id
        new_res.id = (hetfield, res_num, insertion)
        new_chain_b.add(new_res)
        if hetfield == ' ':  # Only increment for standard residues
            res_num += 1
    
    new_model.add(new_chain_b)
    chain_b_length = get_chain_length(new_chain_b)
    
    # Write output
    io = PDBIO()
    io.set_structure(new_structure)
    io.save(str(output_pdb_path))
    
    return chain_a_length, chain_b_length


def transform_interface_residues(interface_res_str):
    """
    Transform InterfaceResidues from B to A notation.
    E.g., "B1,B2,B3" → "A1,A2,A3"
    
    Args:
        interface_res_str: Comma-separated interface residues string
    
    Returns:
        str: Transformed interface residues with A instead of B
    """
    residues = interface_res_str.split(',')
    transformed = []
    
    for res in residues:
        res = res.strip()
        if res.startswith('B'):
            # Replace B with A
            transformed.append('A' + res[1:])
        else:
            # Keep other chains unchanged
            transformed.append(res)
    
    return ','.join(transformed)


def parse_interface_residues(interface_res_str):
    """
    Parse interface residues string to extract residue numbers.
    E.g., "A1,A2,A3" → [1, 2, 3]
    
    Args:
        interface_res_str: Comma-separated interface residues string (e.g., "A1,A2,A3")
    
    Returns:
        set: Set of residue numbers (1-indexed)
    """
    residues = interface_res_str.split(',')
    residue_numbers = set()
    
    for res in residues:
        res = res.strip()
        # Extract the numeric part (e.g., "A123" → 123)
        match = re.search(r'[A-Z](\d+)', res)
        if match:
            residue_numbers.add(int(match.group(1)))
    
    return residue_numbers


def generate_inpaint_seq(chain_a_length, chain_b_length, interface_res_str=None, 
                         fix_interface=False):
    """
    Generate RFdiffusion inpaint_seq mask array.
    
    By default:
    - Chain A (binder) residues: False
    - Chain B (target) residues: True
    
    If fix_interface is True:
    - Interface residues in Chain A are set to True
    
    Args:
        chain_a_length: Number of residues in chain A (binder)
        chain_b_length: Number of residues in chain B (target)
        interface_res_str: Comma-separated interface residues (e.g., "A1,A2,A3")
        fix_interface: If True, set interface residues to True
    
    Returns:
        list: Boolean array for rfd_inpaint_seq
    """
    # Initialize array: False for chain A, True for chain B
    inpaint_seq = [False] * chain_a_length + [True] * chain_b_length
    
    # If fix_interface is True, set interface residues in chain A to True
    if fix_interface and interface_res_str:
        interface_residues = parse_interface_residues(interface_res_str)
        
        for res_num in interface_residues:
            # Convert to 0-indexed array position
            # res_num is 1-indexed position in chain A
            if 1 <= res_num <= chain_a_length:
                array_index = res_num - 1
                inpaint_seq[array_index] = True
    
    return inpaint_seq


def parse_design_name(design_name):
    """
    Parse the Design column to extract batch number, length, and seed.
    E.g., "Batch_0_l146_s894229" → (0, 146, 894229)
    
    Args:
        design_name: Design name from CSV (e.g., "Batch_0_l146_s894229")
    
    Returns:
        tuple: (batch_num, length, seed) or None if parsing fails
    """
    # Pattern: Batch_Z_lXXX_sYYY (case insensitive)
    match = re.search(r'batch_(\d+)_l(\d+)_s(\d+)', design_name, re.IGNORECASE)
    if match:
        batch_num = int(match.group(1))
        length = int(match.group(2))
        seed = int(match.group(3))
        return batch_num, length, seed
    return None


def construct_pdb_filename(batch_num, length, seed):
    """
    Construct PDB filename from batch number, length, and seed.
    E.g., (0, 146, 894229) → "batch_0_l146_s894229.pdb"
    
    Args:
        batch_num: Batch number
        length: Design length
        seed: Random seed
    
    Returns:
        str: PDB filename
    """
    return f"batch_{batch_num}_l{length}_s{seed}.pdb"


def create_metadata_from_row(row, fold_id, chain_a_length, chain_b_length, 
                              fix_interface=False):
    """
    Create metadata dictionary from CSV row.
    
    Args:
        row: Dictionary representing a CSV row
        fold_id: New sequential fold_id to assign
        chain_a_length: Length of chain A (binder)
        chain_b_length: Length of chain B (target)
        fix_interface: If True, set interface residues to True in inpaint_seq
    
    Returns:
        dict: Metadata dictionary
    """
    # Transform interface residues
    bc_intface_res = transform_interface_residues(row['InterfaceResidues'])
    
    # Generate inpaint_seq mask
    rfd_inpaint_seq = generate_inpaint_seq(
        chain_a_length, 
        chain_b_length, 
        bc_intface_res,
        fix_interface
    )
    
    # Extract metadata
    metadata = {
        'fold_id': fold_id,
        'bc_length': int(row['Length']),
        'bc_plddt': float(row['pLDDT']),
        'bc_target_rmsd': float(row['Target_RMSD']),
        'bc_intface_res': bc_intface_res,
        'bc_time': parse_trajectory_time(row['TrajectoryTime']),
        'rfd_inpaint_seq': rfd_inpaint_seq
    }
    
    return metadata


def collect_all_designs(input_dir):
    """
    Collect all designs from all batch CSV files.
    
    Args:
        input_dir: Directory containing batch CSV files
    
    Returns:
        list: List of tuples (csv_row, pdb_filename)
    """
    designs = []
    csv_files = sorted(input_dir.glob('batch_*.csv'))
    
    for csv_path in csv_files:
        with open(csv_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Parse the Design column
                design_info = parse_design_name(row['Design'])
                
                if design_info is None:
                    print(f"Warning: Could not parse design name '{row['Design']}', skipping...")
                    continue
                
                batch_num, length, seed = design_info
                pdb_filename = construct_pdb_filename(batch_num, length, seed)
                
                designs.append((row, pdb_filename))
    
    return designs


def process_design(row, pdb_path, fold_id, output_dir, fix_interface=False):
    """
    Process a single design.
    
    Args:
        row: CSV row dictionary
        pdb_path: Path to input PDB file
        fold_id: Sequential fold_id for this design
        output_dir: Output directory path
        fix_interface: If True, set interface residues to True in inpaint_seq
    """
    # Define output paths
    output_pdb = output_dir / f"fold_{fold_id}.pdb"
    output_json = output_dir / f"fold_{fold_id}.json"
    
    # Swap and renumber chains in PDB, get chain lengths
    chain_a_length, chain_b_length = swap_and_renumber_chains(pdb_path, output_pdb)
    
    # Create metadata from CSV row
    metadata = create_metadata_from_row(
        row, 
        fold_id, 
        chain_a_length, 
        chain_b_length,
        fix_interface
    )
    
    # Write metadata to JSON
    with open(output_json, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"Processed: {pdb_path.name} → fold_{fold_id}.pdb, fold_{fold_id}.json")


def main():
    parser = argparse.ArgumentParser(
        description='Process BindCraft batch design outputs (PDB + CSV files)'
    )
    parser.add_argument(
        '--input_dir',
        type=str,
        required=True,
        help='Input directory containing batch CSV and PDB files'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory for processed files'
    )
    parser.add_argument(
        '--fix_interface_residues',
        action='store_true',
        default=False,
        help='If set, interface residues in chain A will be marked as True in rfd_inpaint_seq'
    )
    
    args = parser.parse_args()
    
    # Convert to Path objects
    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    
    # Create output directory if it doesn't exist
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Collecting designs from batch CSV files...")
    designs = collect_all_designs(input_dir)
    
    print(f"Found {len(designs)} designs across all batches")
    print(f"Output directory: {output_dir}")
    print(f"Fix interface residues: {args.fix_interface_residues}")
    print("-" * 60)
    
    # Process each design with sequential fold_id starting from 0
    fold_id = 0
    for row, pdb_filename in designs:
        pdb_path = input_dir / pdb_filename
        
        if not pdb_path.exists():
            print(f"Warning: PDB file {pdb_filename} not found, skipping...")
            continue
        
        try:
            process_design(row, pdb_path, fold_id, output_dir, args.fix_interface_residues)
            fold_id += 1
        except Exception as e:
            print(f"Error processing {pdb_filename}: {e}")
            continue
    
    print("-" * 60)
    print(f"Processing complete! Generated {fold_id} designs.")


if __name__ == '__main__':
    main()
