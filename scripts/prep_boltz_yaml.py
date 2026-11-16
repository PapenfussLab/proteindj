import os
import argparse
import string
import io
from copy import deepcopy
from Bio.PDB import PDBParser, PDBIO, Select, PPBuilder
from Bio.PDB.Chain import Chain as PDBChain
from Bio.PDB.Model import Model as PDBModel
from Bio.PDB.Structure import Structure as PDBStructure
import yaml


class ChainSelect(Select):
    """Select all atoms (used for clean PDB writing)."""
    def accept_residue(self, residue):
        return 1


def get_available_chain_ids(used_ids):
    """Return uppercase letters not currently used as chain IDs, in alphabetical order."""
    return [c for c in string.ascii_uppercase if c not in used_ids]


def split_structure_chains(structure, max_break_distance=3.0):
    """
    Detect chain breaks (C→N > max_break_distance Å) and split chains accordingly.

    Returns a tuple of:
        new_structure: BioPython Structure with split chains assigned new chain IDs
        chain_id_map: dict mapping original chain ID → [list of new chain IDs]
    """
    pp_builder = PPBuilder(radius=max_break_distance)
    used_ids = set(chain.id for model in structure for chain in model)
    available = iter(get_available_chain_ids(used_ids))

    new_structure = PDBStructure('structure')
    new_model = PDBModel(0)
    new_structure.add(new_model)

    chain_id_map = {}

    for model in structure:
        for chain in model:
            peptides = pp_builder.build_peptides(chain)

            if len(peptides) <= 1:
                # No breaks: copy chain unchanged
                new_chain = deepcopy(chain)
                new_model.add(new_chain)
                chain_id_map[chain.id] = [chain.id]
            else:
                # Chain has breaks: assign new IDs to each sub-chain
                new_ids = [chain.id] + [next(available) for _ in range(len(peptides) - 1)]
                chain_id_map[chain.id] = new_ids

                for peptide, new_id in zip(peptides, new_ids):
                    new_chain = PDBChain(new_id)
                    for residue in peptide:
                        new_chain.add(deepcopy(residue))
                    new_model.add(new_chain)
        break  # Only process first model

    return new_structure, chain_id_map


def add_seqres_to_pdb(input_pdb, output_pdb, max_break_distance=3.0):
    """
    Add SEQRES records to a PDB file. Chains with peptide bond breaks
    (C→N > max_break_distance Å) are split into separate chains before writing.

    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file with SEQRES records
        max_break_distance: Maximum C→N distance (Å) to consider a chain break
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_pdb)

    new_structure, _ = split_structure_chains(structure, max_break_distance)

    # Build chain sequences (3-letter codes) from the new split structure
    chain_sequences = {}
    for model in new_structure:
        for chain in model:
            seq = [r.resname for r in chain if r.id[0] == ' ']
            if seq:
                chain_sequences[chain.id] = seq
        break

    with open(output_pdb, 'w') as out_file:
        # Write SEQRES records
        for chain_id, sequence in chain_sequences.items():
            num_residues = len(sequence)
            for line_num, i in enumerate(range(0, num_residues, 13), start=1):
                seq_chunk = sequence[i:i + 13]
                residues_str = ' '.join(seq_chunk)
                out_file.write(f"SEQRES {line_num:>3} {chain_id} {num_residues:>4}  {residues_str}\n")

        # Write ATOM records via PDBIO
        pdb_io = PDBIO()
        pdb_io.set_structure(new_structure)
        buffer = io.StringIO()
        pdb_io.save(buffer, ChainSelect())
        buffer.seek(0)
        for line in buffer:
            if not line.startswith('END'):
                out_file.write(line)

        out_file.write('END\n')


def extract_sequences(pdb_path, msa_files=None, msa_chain=None, msa_exclude_chain=None, max_break_distance=3.0):
    """
    Extract sequences from a PDB file, splitting chains at peptide bond breaks
    (C→N > max_break_distance Å).

    Args:
        pdb_path: Path to PDB file
        msa_files: Optional ordered list of .a3m MSA file paths to apply to the chain(s)
                   selected by msa_chain/msa_exclude_chain, matched positionally (one file
                   per chain, in chain order). If a single file is given but multiple chains
                   are selected, it is applied to all of them (e.g. an identical-sequence
                   multi-chain target). If more/fewer files than selected chains are given,
                   they are matched positionally and the rest are left without an MSA.
        msa_chain: Chain ID to apply MSA to (default None applies to all, unless
                   msa_exclude_chain is set)
        msa_exclude_chain: If set, apply MSA to every chain except this one (e.g. binder
                           mode: exclude the binder chain so all target chain(s) get the
                           MSA, however many there are). Takes precedence over msa_chain.
        max_break_distance: Maximum C→N distance (Å) to consider a chain break

    Returns:
        sequences: list of dicts with keys id, sequence, msa — one entry per (sub-)chain.
                   Sub-chains from a split always receive msa: empty.
        chain_id_map: dict mapping original chain ID → [list of new chain IDs]
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_path)

    pp_builder = PPBuilder(radius=max_break_distance)
    used_ids = set(chain.id for model in structure for chain in model)
    available = iter(get_available_chain_ids(used_ids))

    # Resolve which original chain(s) get an MSA, and which file goes to which chain
    all_chain_ids = [chain.id for chain in structure[0]]
    if msa_exclude_chain is not None:
        eligible_ids = [cid for cid in all_chain_ids if cid != msa_exclude_chain]
    elif msa_chain is not None:
        eligible_ids = [msa_chain] if msa_chain in all_chain_ids else []
    else:
        eligible_ids = []

    chain_msa_map = {}
    if msa_files:
        if len(msa_files) == 1:
            for cid in eligible_ids:
                chain_msa_map[cid] = msa_files[0]
        else:
            if len(msa_files) != len(eligible_ids):
                print(f'WARNING: {len(msa_files)} MSA file(s) provided but {len(eligible_ids)} '
                      f'target chain(s) found in {os.path.basename(pdb_path)}; matching '
                      f'positionally, extra file(s)/chain(s) will be unmatched.')
            for cid, msa_path in zip(eligible_ids, msa_files):
                chain_msa_map[cid] = msa_path

    sequences = []
    chain_id_map = {}

    for model in structure:
        for chain in model:
            peptides = pp_builder.build_peptides(chain)

            if not peptides:
                chain_id_map[chain.id] = [chain.id]
                continue

            if len(peptides) == 1:
                # No breaks: original behaviour, MSA applies if chain was matched to a file
                msa_path = chain_msa_map.get(chain.id)
                seq = str(peptides[0].get_sequence())
                msa_value = os.path.basename(msa_path) if msa_path else 'empty'
                sequences.append({'id': chain.id, 'sequence': seq, 'msa': msa_value})
                chain_id_map[chain.id] = [chain.id]
            else:
                # Has breaks: split into sub-chains, all receive msa: empty
                new_ids = [chain.id] + [next(available) for _ in range(len(peptides) - 1)]
                chain_id_map[chain.id] = new_ids
                print(f'  Chain {chain.id}: detected {len(peptides) - 1} break(s) '
                      f'→ sub-chains: {", ".join(new_ids)}')
                for peptide, new_id in zip(peptides, new_ids):
                    seq = str(peptide.get_sequence())
                    sequences.append({'id': new_id, 'sequence': seq, 'msa': 'empty'})
        break  # Only process first model

    return sequences, chain_id_map


def get_chain_ids(sequences):
    """Extract all chain IDs from sequences list."""
    chain_ids = []
    for seq in sequences:
        chain_ids.extend(seq['id'])
    return chain_ids


def generate_yaml_config(sequences, use_template=False, pdb_filename=None,
                         template_chains=None, template_force=False, template_threshold=None):
    """
    Create YAML structure with individual chain sequences.

    Args:
        sequences: List of sequence dicts with id, sequence, and msa
        use_template: Whether to include template information
        pdb_filename: Name of the PDB template file
        template_chains: List of chain IDs to include as templates (one entry each)
        template_force: Whether to enforce template with potential
        template_threshold: Distance threshold in Angstroms for template deviation
    """
    config = {
        'sequences': [
            {'protein': seq}
            for seq in sequences
        ]
    }

    if use_template and pdb_filename and template_chains:
        template_list = []
        for tc in template_chains:
            template_config = {
                'pdb': f'templates/{pdb_filename}',
                'chain_id': tc
            }
            if template_force:
                template_config['force'] = True
            if template_threshold is not None:
                template_config['threshold'] = float(template_threshold)
            template_list.append(template_config)
        config['templates'] = template_list

    return config


def process_pdb_files(input_dir, output_dir, use_template=False, template_chain='B',
                      template_exclude_chain=None, template_force=False, template_threshold=None,
                      msa_files=None, msa_chain=None, msa_exclude_chain=None,
                      max_break_distance=3.0, binder_only_chain=None):
    """
    Process PDB files and generate YAML configs with optional PDB templates and MSA.

    Args:
        input_dir: Input directory containing PDB files
        output_dir: Output directory for YAML files
        use_template: Whether to create template directory and include in YAML
        template_chain: Original chain ID to use as the sole template (ignored if
                        template_exclude_chain is set)
        template_exclude_chain: If set, template every chain except this one (e.g. binder
                                mode: exclude the binder chain so all target chain(s) are
                                templated, however many there are). Takes precedence over
                                template_chain.
        template_force: Whether to enforce template with potential
        template_threshold: Distance threshold in Angstroms for template deviation
        msa_files: Ordered list of MSA file(s) (.a3m format) to apply to the selected chain(s)
                   (see extract_sequences for the positional matching rules)
        msa_chain: Chain ID to apply MSA to (default None)
        msa_exclude_chain: If set, apply MSA to every chain except this one (see extract_sequences)
        max_break_distance: Maximum C→N distance (Å) to consider a chain break
        binder_only_chain: If set, restrict the YAML to only this original chain's (sub-)chains
                            (e.g. for an unbound/target-free binder-only prediction), and disable
                            templates/MSA since there is no target to template/align against
    """
    if binder_only_chain:
        use_template = False
        msa_files = None
        msa_chain = None
        msa_exclude_chain = None

    os.makedirs(output_dir, exist_ok=True)

    templates_dir = None
    if use_template:
        templates_dir = os.path.join(output_dir, 'templates')
        os.makedirs(templates_dir, exist_ok=True)

    if msa_files:
        import shutil
        for msa_file in msa_files:
            msa_dest = os.path.join(output_dir, os.path.basename(msa_file))
            shutil.copy2(msa_file, msa_dest)
            print(f'Copied MSA file: {msa_file} -> {os.path.basename(msa_file)}\n')

    for filename in os.listdir(input_dir):
        if filename.startswith('fold_') and filename.endswith('.pdb'):
            pdb_path = os.path.join(input_dir, filename)

            # Extract sequences with chain-break detection
            sequences, chain_id_map = extract_sequences(
                pdb_path, msa_files=msa_files, msa_chain=msa_chain, msa_exclude_chain=msa_exclude_chain,
                max_break_distance=max_break_distance
            )

            if binder_only_chain:
                keep_ids = set(chain_id_map.get(binder_only_chain, [binder_only_chain]))
                sequences = [seq for seq in sequences if seq['id'] in keep_ids]
                if not sequences:
                    print(f'ERROR: Binder-only chain "{binder_only_chain}" not found in {filename}')
                    print(f'       Skipping this file...\n')
                    continue

            if use_template:
                if template_exclude_chain is not None:
                    target_orig_ids = [cid for cid in chain_id_map if cid != template_exclude_chain]
                    if not target_orig_ids:
                        print(f'ERROR: No target chains found (excluding "{template_exclude_chain}") in {filename}')
                        print(f'       Skipping this file...\n')
                        continue
                    template_chains = [new_id for cid in target_orig_ids for new_id in chain_id_map[cid]]
                else:
                    if template_chain not in chain_id_map:
                        all_seq_ids = [seq['id'] for seq in sequences]
                        print(f'ERROR: Template chain "{template_chain}" not found in {filename}')
                        print(f'       Available chains: {", ".join(all_seq_ids)}')
                        print(f'       Skipping this file...\n')
                        continue
                    template_chains = chain_id_map[template_chain]

                # Write template PDB with SEQRES records and split chains
                template_path = os.path.join(templates_dir, filename)
                add_seqres_to_pdb(pdb_path, template_path, max_break_distance=max_break_distance)
                print(f'Copied template with SEQRES: {filename} -> templates/{filename}')
            else:
                template_chains = chain_id_map.get(template_chain, [template_chain])

            yaml_config = generate_yaml_config(
                sequences,
                use_template=use_template,
                pdb_filename=filename if use_template else None,
                template_chains=template_chains,
                template_force=template_force,
                template_threshold=template_threshold
            )

            yaml_filename = filename.replace('.pdb', '.yaml')
            yaml_path = os.path.join(output_dir, yaml_filename)

            with open(yaml_path, 'w') as yaml_file:
                yaml.dump(yaml_config, yaml_file, sort_keys=False)

            print(f'Generated YAML: {yaml_filename}\n')


def main():
    """Command-line interface setup"""
    parser = argparse.ArgumentParser(
        description='Generate Boltz-2 YAML configs with optional PDB templates and MSA',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input directory containing PDB files')
    parser.add_argument('-o', '--output', default=None,
                        help='Output directory for YAML and template files')
    parser.add_argument('--use-template', action='store_true',
                        help='Create templates directory and include PDB templates in YAML configs')
    parser.add_argument('--template-chain', default='B',
                        help='Chain ID to use as the sole template (ignored if --template-exclude-chain is set)')
    parser.add_argument('--template-exclude-chain', default=None,
                        help='Template every chain except this one (e.g. "A" in binder mode, so all '
                             'target chain(s) are templated however many there are). Takes precedence '
                             'over --template-chain.')
    parser.add_argument('--template-force', action='store_true',
                        help='Use potential to enforce template structure')
    parser.add_argument('--template-threshold', type=float, default=None,
                        help='Distance threshold (Angstroms) for template deviation')
    parser.add_argument('--msa-files', default=None,
                        help='Comma-separated list of MSA file(s) (.a3m format) to apply to the '
                             'selected chain(s), matched positionally in chain order. A single '
                             'file is broadcast to all selected chains (e.g. for an identical-'
                             'sequence multi-chain target).')
    parser.add_argument('--msa-chain', default=None,
                        help='Chain ID to apply MSA to (default: same as template-chain)')
    parser.add_argument('--msa-exclude-chain', default=None,
                        help='Apply MSA to every chain except this one (default: same as '
                             '--template-exclude-chain). Takes precedence over --msa-chain.')
    parser.add_argument('--max-break-distance', type=float, default=3.0,
                        help='Maximum C→N distance (Å) to consider a chain break (default: 3.0)')
    parser.add_argument('--binder-only-chain', default=None,
                        help='Restrict YAML to only this original chain (e.g. "A") for an '
                             'unbound/target-free binder-only prediction. Disables templates/MSA.')

    args = parser.parse_args()
    output_dir = args.output if args.output else args.input

    msa_files = [f for f in args.msa_files.split(',') if f] if args.msa_files else None
    msa_exclude_chain = args.msa_exclude_chain if args.msa_exclude_chain else args.template_exclude_chain
    msa_chain = args.msa_chain if args.msa_chain else args.template_chain

    process_pdb_files(
        args.input,
        output_dir,
        use_template=args.use_template,
        template_chain=args.template_chain,
        template_exclude_chain=args.template_exclude_chain,
        template_force=args.template_force,
        template_threshold=args.template_threshold,
        msa_files=msa_files,
        msa_chain=msa_chain if msa_files else None,
        msa_exclude_chain=msa_exclude_chain if msa_files else None,
        max_break_distance=args.max_break_distance,
        binder_only_chain=args.binder_only_chain
    )


if __name__ == '__main__':
    main()