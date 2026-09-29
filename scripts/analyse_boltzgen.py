#!/usr/bin/env python3
"""
BoltzGen Design Processor
Processes BoltzGen `design` step output (.cif + .npz sidecar pairs) for downstream compatibility.
- Converts CIF -> PDB (BioPython)
- Relabels chains A, B, C..., placing BoltzGen's own chain 'A' (the designed binder,
  or chain A for redesign) first regardless of its physical position in the CIF -
  BoltzGen preserves original chain IDs but does not reorder chains by role, so for
  boltzgen_motifscaff the designed chain can legitimately appear anywhere in the file
- Extracts the `design_mask` array from the .npz sidecar and inverts it into
  `bg_inpaint_seq` (True = fixed/not to be redesigned by ProteinMPNN/FAMPNN,
  False = designable). This follows the same per-tool naming convention as
  RFdiffusion's `rfd_inpaint_seq` and BindCraft's `bc_inpaint_seq`; downstream
  consumer scripts check all three possible keys since only one is present per
  fold_N.json depending on which tool produced the design.
- Assigns sequential fold_id's starting from 0
"""

import argparse
import json
from pathlib import Path

import numpy as np
from Bio.PDB import Chain, MMCIFParser, Model, PDBIO, Structure


def get_protein_chains_in_order(model):
    """
    Return chains that contain at least one standard residue, with BoltzGen's own
    chain 'A' (the designed chain) placed first and all others kept in file order.

    BoltzGen preserves the original author chain ID on write but does not reorder
    chains by role, so chain 'A' can appear anywhere in the file - relabeling by raw
    file order would silently swap the designed chain with a fixed context chain.

    Args:
        model: BioPython Model object

    Returns:
        list: Chain objects, chain 'A' first, remaining chains in file order
    """
    chains = [chain for chain in model if any(res.id[0] == ' ' for res in chain)]
    if not chains:
        return chains
    chain_a = [chain for chain in chains if chain.id == 'A']
    if not chain_a:
        raise ValueError("No chain 'A' (the designed chain) found among protein chains")
    others = [chain for chain in chains if chain.id != 'A']
    return chain_a + others


def cif_to_relabelled_pdb(cif_path, pdb_path):
    """
    Convert a BoltzGen-generated CIF file to a PDB file, relabeling chains
    A, B, C... with BoltzGen's chain 'A' first (see get_protein_chains_in_order)
    and remaining chains in file order. HETATM/non-standard residues are dropped.
    Residues are renumbered sequentially and continuously across chains
    (chain A: 1..n_res_A, chain B: n_res_A+1..n_res_A+n_res_B, etc.) so residue
    numbers don't overlap between chains - this matches the convention used by
    RFdiffusion/BindCraft outputs, which downstream tools (e.g. dl_binder_design's
    FIXED-residue-label parsing) rely on to disambiguate chains by residue number.

    Args:
        cif_path: Path to input CIF file
        pdb_path: Path to output PDB file

    Returns:
        list: Number of standard residues in each output chain, in order (e.g. [134, 162])
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure('design', str(cif_path))
    model = next(structure.get_models())

    orig_chains = get_protein_chains_in_order(model)
    if not orig_chains:
        raise ValueError(f"No protein chains found in {cif_path}")

    new_structure = Structure.Structure('design')
    new_model = Model.Model(0)
    new_structure.add(new_model)

    chain_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    chain_lengths = []
    next_resnum = 1

    for letter, orig_chain in zip(chain_letters, orig_chains):
        new_chain = Chain.Chain(letter)
        n_res = 0
        for residue in orig_chain:
            if residue.id[0] != ' ':  # Skip HETATM/water
                continue
            new_residue = residue.copy()
            new_residue.id = (' ', next_resnum, ' ')
            new_chain.add(new_residue)
            n_res += 1
            next_resnum += 1
        new_model.add(new_chain)
        chain_lengths.append(n_res)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(str(pdb_path))

    return chain_lengths


def load_design_mask(npz_path):
    """
    Load the `design_mask` array from a BoltzGen .npz metadata sidecar.

    Args:
        npz_path: Path to .npz file

    Returns:
        numpy.ndarray: Boolean array, True = designable residue
    """
    data = np.load(npz_path)
    if 'design_mask' not in data:
        raise KeyError(f"No 'design_mask' array found in {npz_path}")
    return data['design_mask'].astype(bool)


def build_metadata(fold_id, design_mask, design_mode):
    """
    Build the metadata dictionary for a single design.

    BoltzGen's design_mask[i] == True means residue i was (re)designed by BoltzGen.
    ProteinDJ's bg_inpaint_seq[i] == True means residue i is FIXED (sequence held)
    during downstream ProteinMPNN/FAMPNN sequence design, so the mask is inverted.

    Args:
        fold_id: Sequential fold_id to assign
        design_mask: Boolean numpy array from the .npz sidecar
        design_mode: 'boltzgen_denovo' or 'boltzgen_motifscaff'

    Returns:
        dict: Metadata dictionary
    """
    bg_inpaint_seq = [not bool(v) for v in design_mask]

    return {
        'fold_id': fold_id,
        'bg_design_mode': design_mode,
        'bg_inpaint_seq': bg_inpaint_seq,
    }


def collect_design_pairs(input_dir):
    """
    Collect (cif_path, npz_path) pairs for all BoltzGen designs in a directory.

    BoltzGen writes a `<name>.npz` metadata sidecar only for the generated design
    (not for the accompanying `<name>_native.cif` reference structure), so using
    the .npz files as the authoritative list naturally excludes native structures.

    Args:
        input_dir: Directory containing .cif and .npz files

    Returns:
        list: Sorted list of (cif_path, npz_path) tuples
    """
    pairs = []
    for npz_path in sorted(input_dir.glob('*.npz')):
        cif_path = npz_path.with_suffix('.cif')
        if not cif_path.exists():
            print(f"Warning: no matching .cif for {npz_path.name}, skipping")
            continue
        pairs.append((cif_path, npz_path))
    return pairs


def process_design(cif_path, npz_path, fold_id, output_dir, design_mode):
    """
    Process a single BoltzGen design.

    Args:
        cif_path: Path to input CIF file
        npz_path: Path to input .npz metadata sidecar
        fold_id: Sequential fold_id for this design
        output_dir: Output directory path
        design_mode: 'boltzgen_denovo' or 'boltzgen_motifscaff'
    """
    output_pdb = output_dir / f"fold_{fold_id}.pdb"
    output_json = output_dir / f"fold_{fold_id}.json"

    chain_lengths = cif_to_relabelled_pdb(cif_path, output_pdb)
    design_mask = load_design_mask(npz_path)

    if sum(chain_lengths) != len(design_mask):
        output_pdb.unlink(missing_ok=True)
        raise ValueError(
            f"Residue count mismatch (PDB={sum(chain_lengths)}, "
            f"design_mask={len(design_mask)})"
        )

    metadata = build_metadata(fold_id, design_mask, design_mode)
    with open(output_json, 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"Processed: {cif_path.name} -> fold_{fold_id}.pdb, fold_{fold_id}.json (chains={chain_lengths})")


def main():
    parser = argparse.ArgumentParser(
        description='Process BoltzGen design step outputs (.cif + .npz pairs)'
    )
    parser.add_argument(
        '--input_dir',
        type=str,
        required=True,
        help='Input directory containing .cif and .npz files from all batches'
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory for processed files'
    )
    parser.add_argument(
        '--design_mode',
        type=str,
        default='boltzgen_denovo',
        choices=['boltzgen_denovo', 'boltzgen_motifscaff'],
        help='BoltzGen design mode used to generate these designs'
    )

    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Collecting designs from .cif/.npz pairs...")
    design_pairs = collect_design_pairs(input_dir)

    print(f"Found {len(design_pairs)} designs")
    print(f"Output directory: {output_dir}")
    print(f"Design mode: {args.design_mode}")
    print("-" * 60)

    fold_id = 0
    for cif_path, npz_path in design_pairs:
        try:
            process_design(cif_path, npz_path, fold_id, output_dir, args.design_mode)
            fold_id += 1
        except Exception as e:
            print(f"Error processing {cif_path.name}: {e}")
            continue

    print("-" * 60)
    print(f"Processing complete! Generated {fold_id} designs.")


if __name__ == '__main__':
    main()
