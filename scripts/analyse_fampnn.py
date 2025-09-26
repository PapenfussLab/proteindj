import gemmi
import argparse
import os
import json

def get_chain_sequence(chain):
    """Extract amino acid sequence from a Gemmi chain object"""
    sequence = []
    for residue in chain:
        # Get residue info using built-in chemical data
        res_info = gemmi.find_tabulated_residue(residue.name)
        if not res_info or not res_info.is_amino_acid():
            continue
        # Use fasta_code() to get X for non-standard residues
        aa = res_info.fasta_code()
        sequence.append(aa if aa != ' ' else 'X')
    return ''.join(sequence)

def average_per_residue_bfactor(input_dir, chain_id, ignore_cbeta, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    results = {}
    
    for pdb_file in os.listdir(input_dir):
        if pdb_file.endswith('.pdb'):
            full_path = os.path.join(input_dir, pdb_file)
            try:
                structure = gemmi.read_structure(full_path)
            except FileNotFoundError:
                print(f"File {full_path} not found, skipping")
                continue
                
            if ignore_cbeta:
                backbone_atoms = {'C', 'N', 'O', 'CA', 'CB'}
            else:
                backbone_atoms = {'C', 'N', 'O', 'CA'}
                
            design_name = os.path.splitext(pdb_file)[0]  # Remove .pdb extension
            
            # If processing of all chains is requested
            if chain_id == 'all_chains':
                all_chains_residue_averages = []
                all_chains_sequences = {}
                chain_averages = {}
                
                for model in structure:
                    for chain in model:
                        chain_id_current = chain.name
                        sequence_str = get_chain_sequence(chain)
                        all_chains_sequences[chain_id_current] = sequence_str
                        
                        chain_residue_averages = []
                        for residue in chain:
                            residue_total = 0.0
                            residue_count = 0
                            
                            for atom in residue:
                                if atom.name not in backbone_atoms:
                                    residue_total += atom.b_iso
                                    residue_count += 1
                            
                            if residue_count > 0:
                                chain_residue_averages.append(residue_total / residue_count)
                        
                        if chain_residue_averages:
                            all_chains_residue_averages.extend(chain_residue_averages)
                            chain_avg = sum(chain_residue_averages) / len(chain_residue_averages)
                            chain_averages[chain_id_current] = round(chain_avg, 2)
                            print(f"Chain {chain_id_current} average pSCE: {chain_avg:.2f}")
                
                if not all_chains_residue_averages:
                    print(f"No atoms with pSCE score in {pdb_file}")
                    continue
                
                # Calculate the overall average across all chains
                average_psce = sum(all_chains_residue_averages) / len(all_chains_residue_averages)
                print(f"Overall average pSCE for all chains in {pdb_file}: {average_psce:.2f}")

                # Create JSON output with sequences from all chains
                output_data = {
                    "design": design_name,
                    "sequence": '|'.join(f"{chain_id}:{seq}" for chain_id, seq in all_chains_sequences.items()),
                    "chain_avg_psce": chain_averages,
                    "fampnn_avg_psce": round(average_psce, 2)
                }
            
            # If a specific chain is requested, process only that chain
            else:
                residue_averages = []
                sequence_str = ""
                
                for model in structure:
                    chain = model.find_chain(chain_id)
                    if not chain:
                        continue

                    # extract sequence for chain_id
                    sequence_str = get_chain_sequence(chain)    
                    for residue in chain:
                        residue_total = 0.0
                        residue_count = 0
                        
                        for atom in residue:
                            if atom.name not in backbone_atoms:
                                residue_total += atom.b_iso
                                residue_count += 1
                        
                        if residue_count > 0:
                            residue_averages.append(residue_total / residue_count)
                
                if not residue_averages:
                    print(f"No atoms with pSCE score in chain {chain_id} of {pdb_file}")
                    continue  # Skip to next file
                
                average_psce = sum(residue_averages) / len(residue_averages)
                print(f"Average pSCE (on a per-residue basis) for side chains in chain {chain_id} of {pdb_file}: {average_psce:.2f}")
                
                # Create JSON output for single chain
                output_data = {
                    "design": design_name,
                    "sequence": sequence_str,
                    "fampnn_avg_psce": round(average_psce, 2)
                }
            
            # Write JSON file
            output_filename = f"{design_name}.json"
            output_path = os.path.join(out_dir, output_filename)
            with open(output_path, 'w') as f:
                f.write(json.dumps(output_data) + '\n')
            
            results[pdb_file] = average_psce
            print(f"Created {output_filename}")
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Calculates average pSCE metric (Predicted Sidechain Confidence Error) of designs by Full-Atom MPNN')
    parser.add_argument('--input_dir', required=True, help='Input directory containing PDB files')
    parser.add_argument('--chain_id', required=False, default='all', 
                        help='Chain ID of designed chain. If not provided, will calculate average over all chains (default: all)')
    parser.add_argument('--ignore_cbeta', action='store_true', help='Ignore C-beta atom types when calculating average pSCE')
    parser.add_argument('--out_dir', default='./averagePSCE', help='Output directory for scores in JSON format (default: ./averagePSCE)')
    args = parser.parse_args()

    average_psce = average_per_residue_bfactor(args.input_dir, args.chain_id, args.ignore_cbeta, args.out_dir)

if __name__ == "__main__":
    main()