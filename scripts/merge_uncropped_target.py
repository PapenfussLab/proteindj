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
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

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

def get_sequence_from_structure(structure, chain_id):
    """Extract sequence from structure"""
    aa_dict = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    try:
        chain = structure[0][chain_id]
        sequence = ""
        for residue in chain:
            if residue.id[0] == ' ':  # Standard amino acid
                res_name = residue.get_resname()
                sequence += aa_dict.get(res_name, 'X')
        return sequence
    except KeyError:
        return ""

def get_chain_atoms(structure, chain_id, atom_type="CA"):
    """Get specified atoms from a chain"""
    try:
        chain = structure[0][chain_id]
        atoms = []
        for residue in chain:
            if residue.id[0] == ' ':  # Standard amino acid
                for atom in residue:
                    if atom.get_name() == atom_type:
                        atoms.append(atom)
                        break
        return atoms
    except KeyError:
        return []

def find_best_alignment_region(ref_seq, target_seq, min_overlap=10):
    """Find the best alignment region between two sequences"""
    alignments = pairwise2.align.localms(ref_seq, target_seq, 2, -1, -0.5, -0.1)
    
    if not alignments:
        return None, None, None
    
    best_alignment = alignments[0]
    aligned_ref, aligned_target, score, start, end = best_alignment
    
    # Find the alignment region without gaps
    ref_indices = []
    target_indices = []
    ref_pos = 0
    target_pos = 0
    
    for i in range(len(aligned_ref)):
        if aligned_ref[i] != '-' and aligned_target[i] != '-':
            ref_indices.append(ref_pos)
            target_indices.append(target_pos)
        
        if aligned_ref[i] != '-':
            ref_pos += 1
        if aligned_target[i] != '-':
            target_pos += 1
    
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

def align_structures_by_sequence(reference_structure, target_structure, ref_chain='B', target_chain='B'):
    """Align structures using sequence-based matching"""
    # Get sequences
    ref_seq = get_sequence_from_structure(reference_structure, ref_chain)
    target_seq = get_sequence_from_structure(target_structure, target_chain)
    
    logger = logging.getLogger(__name__)
    logger.info(f"Reference sequence length: {len(ref_seq)}")
    logger.info(f"Target sequence length: {len(target_seq)}")
    
    # Find best alignment
    ref_indices, target_indices, score = find_best_alignment_region(ref_seq, target_seq)
    
    if ref_indices is None:
        logger.warning("No good sequence alignment found")
        return target_structure, 999.0
    
    logger.info(f"Alignment score: {score}, aligned residues: {len(ref_indices)}")
    
    # Get CA atoms for aligned regions
    ref_atoms = get_chain_atoms(reference_structure, ref_chain, "CA")
    target_atoms = get_chain_atoms(target_structure, target_chain, "CA")
    
    # Select atoms based on alignment
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
        
        # Align uncropped structure to this specific RFdiffusion file's chain B
        aligned_uncropped, rmsd = align_structures_by_sequence(
            rf_structure, uncropped_structure, 'B', 'B'
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
