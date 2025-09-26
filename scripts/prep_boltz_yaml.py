import os
import argparse
from collections import defaultdict
from Bio import SeqIO
import yaml

def extract_chain_groups(pdb_path):
    """Group chains by identical sequences, returns {sequence: [chain_ids]}"""
    sequence_map = defaultdict(list)
    
    for record in SeqIO.parse(pdb_path, 'pdb-atom'):
        chain_id = record.annotations.get('chain', '')
        sequence = str(record.seq)
        sequence_map[sequence].append(chain_id)
    
    return sequence_map

def generate_yaml_config(chain_groups):
    """Create YAML structure with combined chain IDs for identical sequences"""
    return {
        'sequences': [
            {
                'protein': {
                    'id': sorted(chain_ids),
                    'sequence': sequence,
                    'msa': 'empty'
                }
            }
            for sequence, chain_ids in chain_groups.items()
        ]
    }

def process_pdb_files(input_dir, output_dir):
    """Process PDB files with flexible chain handling"""
    os.makedirs(output_dir, exist_ok=True)
    
    for filename in os.listdir(input_dir):
        if filename.startswith('fold_') and filename.endswith('.pdb'):
            pdb_path = os.path.join(input_dir, filename)
            
            # Group chains by identical sequences
            chain_groups = extract_chain_groups(pdb_path)
            
            # Generate YAML config
            yaml_config = generate_yaml_config(chain_groups)
            
            # Create output path
            yaml_filename = filename.replace('.pdb', '.yaml')
            yaml_path = os.path.join(output_dir, yaml_filename)
            
            # Write YAML file
            with open(yaml_path, 'w') as yaml_file:
                yaml.dump(yaml_config, yaml_file, sort_keys=False)
            
            print(f'Generated: {yaml_path}')

def main():
    """Command-line interface setup"""
    parser = argparse.ArgumentParser(
        description='Generate Boltz-2 YAML configs with automatic chain grouping',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing PDB files')
    parser.add_argument('-o', '--output', default=None,
                        help='Output directory for YAML files')
    
    args = parser.parse_args()
    output_dir = args.output if args.output else args.input
    
    process_pdb_files(args.input, output_dir)

if __name__ == '__main__':
    main()
