#!/usr/bin/env python3

import os
import sys
import logging
from pathlib import Path
import numpy as np
import argparse
from multiprocessing import Pool
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.PDB.Selection import unfold_entities
from Bio.Align import PairwiseAligner

def setup_logging():
    """Setup logging configuration"""
    log_file = "binder_merge.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def get_sequence_and_atoms(structure, chain_ids, atom_type="CA"):
    """Extract concatenated sequence and matching atoms across the given chain(s), in order.

    Supports multiple chain IDs (e.g. all of a design's own target chains: ['B', 'C']) so
    alignment can use every target chain, not just one hardcoded chain. Residues missing
    the requested atom_type are skipped from both outputs together, so sequence[i] always
    corresponds to atoms[i].
    """
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    sequence = ""
    atoms = []
    for chain_id in chain_ids:
        try:
            chain = structure[0][chain_id]
        except KeyError:
            continue
        for residue in chain:
            if residue.id[0] != ' ':  # Standard amino acid only
                continue
            atom = next((a for a in residue if a.get_name() == atom_type), None)
            if atom is None:
                continue
            sequence += aa_dict.get(residue.get_resname(), 'X')
            atoms.append(atom)
    return sequence, atoms

def get_target_chain_ids(structure, binder_chain='A'):
    """Return all chain IDs in structure except the binder chain, in structure order."""
    return [chain.id for chain in structure[0] if chain.id != binder_chain]

def find_best_alignment_region(ref_seq, target_seq, min_overlap=10):
    """Find the best alignment region between two sequences"""
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    alignments = aligner.align(ref_seq, target_seq)

    try:
        best_alignment = alignments[0]
    except IndexError:
        return None, None, None

    score = best_alignment.score

    # best_alignment.aligned gives matching (ref_ranges, target_ranges) blocks of
    # equal length without gaps (may include mismatches within a block)
    ref_indices = []
    target_indices = []
    ref_ranges, target_ranges = best_alignment.aligned
    for (r_start, r_end), (t_start, t_end) in zip(ref_ranges, target_ranges):
        for offset in range(r_end - r_start):
            ref_indices.append(r_start + offset)
            target_indices.append(t_start + offset)

    if len(ref_indices) < min_overlap:
        return None, None, None

    return ref_indices, target_indices, score

def merge_chains_to_single_chain(pdb_path, output_path, target_chain='B'):
    """Merge all chains in PDB to a single chain in alphabetical order"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", pdb_path)
    
    # Get all chains and sort alphabetically
    chains = sorted([chain.id for chain in structure[0]])
    
    with open(output_path, 'w') as f_out:
        atom_serial = 1
        res_num = 1
        
        for chain_id in chains:
            chain = structure[0][chain_id]
            
            for residue in chain:
                if residue.id[0] != ' ':  # Skip HETATM
                    continue
                    
                for atom in residue:
                    coord = atom.get_coord()
                    occupancy = atom.get_occupancy()
                    bfactor = atom.get_bfactor()
                    element = atom.element if hasattr(atom, 'element') else atom.name[0]
                    
                    atom_line = (
                        f"ATOM  {atom_serial:5d}  "
                        f"{atom.name:<4s}"
                        f"{residue.resname:>3s} "
                        f"{target_chain}"
                        f"{res_num:4d}    "
                        f"{coord[0]:8.3f}"
                        f"{coord[1]:8.3f}"
                        f"{coord[2]:8.3f}"
                        f"{occupancy:6.2f}"
                        f"{bfactor:6.2f}"
                        f"          "
                        f"{element:>2s}\n"
                    )
                    f_out.write(atom_line)
                    atom_serial += 1
                
                res_num += 1
        
        f_out.write("END\n")

def align_structures_by_sequence(reference_structure, target_structure, ref_chains, target_chain='B'):
    """Align structures using sequence-based matching.

    ref_chains: chain ID(s) in reference_structure to use, concatenated in order (e.g. all of
                the design's own target chains, ['B'] or ['B', 'C', ...] for multi-target-chain
                designs), so the alignment fit accounts for every target chain.
    """
    # Get sequences
    ref_seq, ref_atoms = get_sequence_and_atoms(reference_structure, ref_chains, "CA")
    target_seq, target_atoms = get_sequence_and_atoms(target_structure, [target_chain], "CA")
    
    logger = logging.getLogger(__name__)
    logger.info(f"Reference sequence length: {len(ref_seq)}")
    logger.info(f"Target sequence length: {len(target_seq)}")
    
    # Find best alignment
    ref_indices, target_indices, score = find_best_alignment_region(ref_seq, target_seq)
    
    if ref_indices is None:
        logger.warning("No good sequence alignment found")
        return target_structure, 999.0
    
    logger.info(f"Alignment score: {score}, aligned residues: {len(ref_indices)}")
    
    # Select atoms based on alignment (ref_atoms/target_atoms already line up with their sequences)
    aligned_ref_atoms = [ref_atoms[i] for i in ref_indices if i < len(ref_atoms)]
    aligned_target_atoms = [target_atoms[i] for i in target_indices if i < len(target_atoms)]
    
    min_atoms = min(len(aligned_ref_atoms), len(aligned_target_atoms))
    if min_atoms < 5:
        logger.warning(f"Too few atoms for alignment: {min_atoms}")
        return target_structure, 999.0
    
    aligned_ref_atoms = aligned_ref_atoms[:min_atoms]
    aligned_target_atoms = aligned_target_atoms[:min_atoms]
    
    # Perform superimposition
    superimposer = Superimposer()
    superimposer.set_atoms(aligned_ref_atoms, aligned_target_atoms)
    
    # Apply transformation to all atoms in target structure
    superimposer.apply(target_structure.get_atoms())
    
    logger.info(f"Alignment RMSD: {superimposer.rms:.3f} Å")
    
    return target_structure, superimposer.rms

def create_combined_pdb_fixed_coords(rf_pdb, aligned_uncropped_pdb, output_path):
    """Create combined PDB maintaining RFdiffusion coordinate system"""
    parser = PDBParser(QUIET=True)
    
    # Load RFdiffusion structure (this is our reference coordinate system)
    rf_structure = parser.get_structure("rf", rf_pdb)
    
    # Load aligned uncropped structure
    uncropped_structure = parser.get_structure("uncropped", aligned_uncropped_pdb)
    
    with open(output_path, 'w') as f_out:
        atom_serial = 1
        res_num = 1
        
        # First, write the binder (chain A) from RFdiffusion - keep original coordinates
        if 'A' in [chain.id for chain in rf_structure[0]]:
            chain_a = rf_structure[0]['A']
            
            for residue in chain_a:
                if residue.id[0] != ' ':
                    continue
                    
                for atom in residue:
                    coord = atom.get_coord()
                    occupancy = atom.get_occupancy()
                    bfactor = atom.get_bfactor()
                    element = atom.element if hasattr(atom, 'element') else atom.name[0]
                    
                    atom_line = (
                        f"ATOM  {atom_serial:5d}  "
                        f"{atom.name:<4s}"
                        f"{residue.resname:>3s} "
                        f"A"
                        f"{res_num:4d}    "
                        f"{coord[0]:8.3f}"
                        f"{coord[1]:8.3f}"
                        f"{coord[2]:8.3f}"
                        f"{occupancy:6.2f}"
                        f"{bfactor:6.2f}"
                        f"          "
                        f"{element:>2s}\n"
                    )
                    f_out.write(atom_line)
                    atom_serial += 1
                
                res_num += 1
        
        # Write TER record
        f_out.write(f"TER   {atom_serial:5d}      A {res_num-1:4d}\n")
        atom_serial += 1
        
        # Now write the aligned uncropped target (chain B)
        if 'B' in [chain.id for chain in uncropped_structure[0]]:
            chain_b = uncropped_structure[0]['B']
            
            for residue in chain_b:
                if residue.id[0] != ' ':
                    continue
                    
                for atom in residue:
                    coord = atom.get_coord()
                    occupancy = atom.get_occupancy()
                    bfactor = atom.get_bfactor()
                    element = atom.element if hasattr(atom, 'element') else atom.name[0]
                    
                    atom_line = (
                        f"ATOM  {atom_serial:5d}  "
                        f"{atom.name:<4s}"
                        f"{residue.resname:>3s} "
                        f"B"
                        f"{res_num:4d}    "
                        f"{coord[0]:8.3f}"
                        f"{coord[1]:8.3f}"
                        f"{coord[2]:8.3f}"
                        f"{occupancy:6.2f}"
                        f"{bfactor:6.2f}"
                        f"          "
                        f"{element:>2s}\n"
                    )
                    f_out.write(atom_line)
                    atom_serial += 1
                
                res_num += 1
        
        f_out.write(f"TER   {atom_serial:5d}      B {res_num-1:4d}\n")
        f_out.write("END\n")

def process_single_rf_file(args):
    """Process a single RFdiffusion file"""
    merged_uncropped_pdb, rf_file, output_dir, reference_rf_structure = args
    
    logger = logging.getLogger(__name__)
    
    try:
        parser = PDBParser(QUIET=True)
        
        # Load the merged uncropped structure
        uncropped_structure = parser.get_structure("uncropped", merged_uncropped_pdb)
        
        # Load current RFdiffusion file
        rf_structure = parser.get_structure("rf", rf_file)
        
        # Align uncropped structure to this design's full set of target chain(s) (everything
        # except the binder chain 'A'), not just a single hardcoded chain
        target_chain_ids = get_target_chain_ids(rf_structure, binder_chain='A')
        if not target_chain_ids:
            raise ValueError(f"No target chain(s) found in {rf_file.name} (only binder chain 'A')")
        
        aligned_uncropped, rmsd = align_structures_by_sequence(
            rf_structure, uncropped_structure, target_chain_ids, 'B'
        )
        
        # Save aligned uncropped structure temporarily
        temp_aligned = output_dir / f"temp_aligned_{rf_file.stem}.pdb"
        io = PDBIO()
        io.set_structure(aligned_uncropped)
        io.save(str(temp_aligned))
        
        # Create final combined structure
        output_path = output_dir / rf_file.name
        create_combined_pdb_fixed_coords(rf_file, temp_aligned, output_path)
        
        # Clean up temporary file
        temp_aligned.unlink()
        
        logger.info(f"Successfully processed {rf_file.name} (RMSD: {rmsd:.3f} Å)")
        return (rf_file.name, None)
        
    except Exception as e:
        logger.error(f"Failed to process {rf_file.name}: {str(e)}")
        return (rf_file.name, str(e))

def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(description="Align and merge uncropped target chain(s) with RFdiffusion binders")
    parser.add_argument("--uncropped_pdb", required=True, type=Path, 
                       help="Path to uncropped target PDB file")
    parser.add_argument("--input_dir", required=True, type=Path, 
                       help="Directory containing RFdiffusion output PDB files")
    parser.add_argument("--output_dir", type=Path, default="merged_outputs", 
                       help="Directory to save output PDB files (default: merged_outputs)")
    parser.add_argument("--ncpus", type=int, default=1, 
                       help="Number of CPUs for parallel processing (default: 1)")
    
    args = parser.parse_args()
    
    if not args.uncropped_pdb.is_file():
        raise ValueError(f"Uncropped PDB file not found: {args.uncropped_pdb}")
    if not args.input_dir.is_dir():
        raise ValueError(f"RFdiffusion directory not found: {args.input_dir}")
    
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    return args

def main():
    args = parse_arguments()
    logger = setup_logging()
    
    # Get RFdiffusion files
    rf_files = list(args.input_dir.glob("*.pdb"))
    if not rf_files:
        logger.error(f"No PDB files found in {args.input_dir}")
        sys.exit(1)
    
    logger.info(f"Found {len(rf_files)} RFdiffusion files to process")
    
    # Step 1: Merge all chains in uncropped PDB to single chain B
    logger.info("Step 1: Merging all chains in uncropped PDB to single chain B")
    merged_uncropped = args.output_dir / "merged_uncropped.pdb"
    merge_chains_to_single_chain(args.uncropped_pdb, merged_uncropped, 'B')
    
    # Step 2: Process each RFdiffusion file individually
    logger.info(f"Step 2: Processing {len(rf_files)} files using {args.ncpus} CPUs")
    
    # Load reference structure for comparison
    parser = PDBParser(QUIET=True)
    reference_rf_structure = parser.get_structure("ref", rf_files[0])
    
    process_args = [(merged_uncropped, rf_file, args.output_dir, reference_rf_structure) 
                   for rf_file in rf_files]
    
    if args.ncpus == 1:
        results = [process_single_rf_file(arg) for arg in process_args]
    else:
        with Pool(processes=args.ncpus) as pool:
            results = pool.map(process_single_rf_file, process_args)
    
    # Log summary
    successful = sum(1 for _, error in results if error is None)
    failed = sum(1 for _, error in results if error is not None)
    
    logger.info(f"\nProcessing Summary:")
    logger.info(f"Total files processed: {len(rf_files)}")
    logger.info(f"Successful: {successful}")
    logger.info(f"Failed: {failed}")
    
    if failed > 0:
        logger.info("\nFailed files:")
        for filename, error in results:
            if error:
                logger.info(f"  {filename}: {error}")
    
    # Clean up
    if merged_uncropped.exists():
        merged_uncropped.unlink()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        sys.exit(1)
