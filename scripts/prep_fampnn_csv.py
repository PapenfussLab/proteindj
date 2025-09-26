import json
import argparse
from pathlib import Path
from collections import defaultdict

def parse_pdb_chains(pdb_path):
    """Extract chain IDs and residue numbers with actual PDB numbering, preserving order."""
    chains = defaultdict(list)
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                chain = line[21]
                resnum = int(line[22:26].strip())  # PDB residue number (1-based)
                res_id = (chain, resnum)
                if res_id not in seen:
                    chains[chain].append(resnum)
                    seen.add(res_id)
    return chains

def indices_to_chain_ranges(indices, chain_residues):
    """Convert indices to FA-MPNN ranges using actual PDB residue numbers."""
    if not indices:
        return ""
    
    # Create a list of (chain, resnum) in PDB order
    index_map = []
    for chain in chain_residues:
        for resnum in chain_residues[chain]:
            index_map.append((chain, resnum))
    
    # Get residues at valid indices
    valid_indices = [i for i in indices if i < len(index_map)]
    fixed_residues = [index_map[i] for i in valid_indices]
    
    # Group into contiguous ranges per chain
    ranges = []
    current_chain, current_start, current_end = None, None, None
    
    for chain, resnum in fixed_residues:
        if chain == current_chain and resnum == current_end + 1:
            current_end = resnum
        else:
            if current_chain is not None:
                ranges.append(f"{current_chain}{current_start}-{current_end}")
            current_chain, current_start, current_end = chain, resnum, resnum
    
    if current_chain is not None:
        ranges.append(f"{current_chain}{current_start}-{current_end}")
    
    return ",".join(ranges)

def process_file_pair(json_path, pdb_path):
    """Process a single PDB/JSON pair."""
    with open(json_path) as f:
        data = json.load(f)
    
    # Get indices where sequence was fixed (True)
    inpaint_seq = data.get('rfd_inpaint_seq', [])
    fixed_indices = [i for i, val in enumerate(inpaint_seq) if val]
    
    # Get chain residue mapping from PDB
    chain_residues = parse_pdb_chains(pdb_path)
    
    return indices_to_chain_ranges(fixed_indices, chain_residues)

def main():
    parser = argparse.ArgumentParser(description='Generate FAMPNN CSV from RFdiffusion outputs')
    parser.add_argument('--input_dir', required=True, help='Directory containing PDB/JSON files')
    parser.add_argument('--out_csv', required=True, help='Output CSV path')
    parser.add_argument('--fix_target_sidechains', action='store_true', help='Fix target sidechains using the same ranges as fixed sequence positions')
    args = parser.parse_args()

    with open(args.out_csv, 'w') as csv_out:
        csv_out.write("pdb,fixed_seq_positions,fixed_sidechains\n")
        
        for json_file in Path(args.input_dir).glob('*.json'):
            # Find matching PDB file
            pdb_file = json_file.with_suffix('.pdb')
            
            if not pdb_file.exists():
                print(f"Warning: No matching PDB for {json_file.name}")
                continue

            # Process the file pair
            seq_ranges = process_file_pair(json_file, pdb_file)
            
            if args.fix_target_sidechains:
                # Write CSV row (PDB name and ranges to fix sequence and side-chains)
                csv_out.write(f'"{pdb_file.stem}","{seq_ranges}","{seq_ranges}"\n')
            else:
                # Write CSV row (PDB name and ranges to fix sequence only)
                csv_out.write(f'"{pdb_file.stem}","{seq_ranges}",""\n')

if __name__ == '__main__':
    main()