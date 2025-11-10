import os
import shutil
import argparse
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select
import yaml


class ChainSelect(Select):
    """Select all atoms (used for clean PDB writing)."""
    def accept_residue(self, residue):
        return 1


def add_seqres_to_pdb(input_pdb, output_pdb):
    """
    Add SEQRES records to a PDB file using BioPython.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file with SEQRES records
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb)
    
    # Extract sequences for each chain (keep 3-letter codes)
    chain_sequences = {}
    for model in structure:
        for chain in model:
            sequence = []
            for residue in chain:
                if residue.id[0] == ' ':  # Standard residue
                    resname = residue.resname
                    sequence.append(resname)
            if sequence:
                chain_sequences[chain.id] = sequence
    
    # Write PDB with SEQRES records
    with open(output_pdb, 'w') as out_file:
        # Write SEQRES records
        for chain_id, sequence in chain_sequences.items():
            # SEQRES records have max 13 residues per line
            num_residues = len(sequence)
            for line_num, i in enumerate(range(0, num_residues, 13), start=1):
                seq_chunk = sequence[i:i+13]
                # Format: SEQRES serNum chainID numRes resName1 resName2 ...
                residues = ' '.join(seq_chunk)
                out_file.write(f"SEQRES {line_num:>3} {chain_id} {num_residues:>4}  {residues}\n")
        
        # Write the structure using BioPython (this includes ATOM records)
        io = PDBIO()
        io.set_structure(structure)
        
        # Save to a temporary string buffer
        import io as io_module
        buffer = io_module.StringIO()
        io.save(buffer, ChainSelect())
        
        # Get the content and write it (skip the END record from BioPython)
        buffer.seek(0)
        for line in buffer:
            if not line.startswith('END'):
                out_file.write(line)
        
        # Write final END record
        out_file.write('END\n')


def extract_sequences(pdb_path):
    """Extract sequences from PDB file, returns list of {id, sequence} dicts."""
    sequences = []
    
    for record in SeqIO.parse(pdb_path, 'pdb-atom'):
        chain_id = record.annotations.get('chain', '')
        sequence = str(record.seq)
        sequences.append({
            'id': chain_id,
            'sequence': sequence,
            'msa': 'empty'
        })
    
    return sequences


def get_chain_ids(sequences):
    """Extract all chain IDs from sequences list."""
    chain_ids = []
    for seq in sequences:
        chain_ids.extend(seq['id'])
    return chain_ids


def generate_yaml_config(sequences, use_template=False, pdb_filename=None, 
                        template_chain='B', template_force=False, template_threshold=None):
    """
    Create YAML structure with individual chain sequences.
    
    Args:
        sequences: List of sequence dicts with id, sequence, and msa
        use_template: Whether to include template information
        pdb_filename: Name of the PDB template file
        template_chain: Chain ID to use as template (default 'B' for target in binder mode)
        template_force: Whether to enforce template with potential
        template_threshold: Distance threshold in Angstroms for template deviation
    """
    config = {
        'sequences': [
            {'protein': seq}
            for seq in sequences
        ]
    }
    
    # Add template configuration if requested
    if use_template and pdb_filename:
        template_config = {
            'pdb': f'templates/{pdb_filename}',
            'chain_id': template_chain
        }
        
        # Add optional template parameters if specified
        if template_force:
            template_config['force'] = True
        if template_threshold is not None:
            template_config['threshold'] = float(template_threshold)
        
        config['templates'] = [template_config]
    
    return config


def process_pdb_files(input_dir, output_dir, use_template=False, template_chain='B',
                      template_force=False, template_threshold=None):
    """
    Process PDB files and generate YAML configs with optional PDB templates.
    
    Args:
        input_dir: Input directory containing PDB files
        output_dir: Output directory for YAML files
        use_template: Whether to create template directory and include in YAML
        template_chain: Chain ID to use as template (default 'B' for target in binder mode)
        template_force: Whether to enforce template with potential
        template_threshold: Distance threshold in Angstroms for template deviation
    """
    # Create output directory for YAMLs
    os.makedirs(output_dir, exist_ok=True)
    
    # Create templates subdirectory if using templates
    templates_dir = None
    if use_template:
        templates_dir = os.path.join(output_dir, 'templates')
        os.makedirs(templates_dir, exist_ok=True)
    
    for filename in os.listdir(input_dir):
        if filename.startswith('fold_') and filename.endswith('.pdb'):
            pdb_path = os.path.join(input_dir, filename)
            
            # Extract sequences from PDB
            sequences = extract_sequences(pdb_path)
            
            # Check if template chain exists when using templates
            if use_template:
                chain_ids = get_chain_ids(sequences)
                if template_chain not in chain_ids:
                    print(f'ERROR: Template chain "{template_chain}" not found in {filename}')
                    print(f'       Available chains: {", ".join(chain_ids)}')
                    print(f'       Skipping this file...\n')
                    continue
                
                # Add SEQRES records and copy to templates directory
                template_path = os.path.join(templates_dir, filename)
                add_seqres_to_pdb(pdb_path, template_path)
                print(f'Copied template with SEQRES: {filename} -> templates/{filename}')
            
            # Generate YAML config
            yaml_config = generate_yaml_config(
                sequences,
                use_template=use_template,
                pdb_filename=filename if use_template else None,
                template_chain=template_chain,
                template_force=template_force,
                template_threshold=template_threshold
            )
            
            # Create output path for YAML
            yaml_filename = filename.replace('.pdb', '.yaml')
            yaml_path = os.path.join(output_dir, yaml_filename)
            
            # Write YAML file
            with open(yaml_path, 'w') as yaml_file:
                yaml.dump(yaml_config, yaml_file, sort_keys=False)
            
            print(f'Generated YAML: {yaml_filename}\n')


def main():
    """Command-line interface setup"""
    parser = argparse.ArgumentParser(
        description='Generate Boltz-2 YAML configs with optional PDB templates',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing PDB files')
    parser.add_argument('-o', '--output', default=None,
                        help='Output directory for YAML and template files')
    parser.add_argument('--use-template', action='store_true',
                        help='Create templates directory and include PDB templates in YAML configs')
    parser.add_argument('--template-chain', default='B',
                        help='Chain ID to use as template (default: B for target in binder mode)')
    parser.add_argument('--template-force', action='store_true',
                        help='Use potential to enforce template structure')
    parser.add_argument('--template-threshold', type=float, default=None,
                        help='Distance threshold (Angstroms) for template deviation')
    
    args = parser.parse_args()
    output_dir = args.output if args.output else args.input
    
    process_pdb_files(
        args.input, 
        output_dir,
        use_template=args.use_template,
        template_chain=args.template_chain,
        template_force=args.template_force,
        template_threshold=args.template_threshold
    )


if __name__ == '__main__':
    main()
