#!/usr/bin/env python3

import sys
import os
import logging
from pathlib import Path
import numpy as np
import argparse
from multiprocessing import Pool
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.PDB.Selection import unfold_entities

def setup_logging():
    """Setup logging configuration"""
    log_file = "alignment.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

# Initialize logger globally
logger = setup_logging()

def calculate_rmsd(atoms1, atoms2):
    """Calculate RMSD between two sets of atoms"""
    if len(atoms1) != len(atoms2):
        raise ValueError(f"Atom lists have different lengths: {len(atoms1)} vs {len(atoms2)}")
    
    superimposer = Superimposer()
    superimposer.set_atoms(atoms1, atoms2)
    return superimposer.rms

def get_chain_atoms(structure, chain_id, atom_type="CA"):
    """Get specified atoms from a chain"""
    try:
        chain = structure[0][chain_id]  # Model 0, Chain ID
        atoms = [atom for atom in unfold_entities(chain, 'A') if atom.get_name() == atom_type]
        if not atoms:
            raise ValueError(f"No {atom_type} atoms found in chain {chain_id}")
        return atoms
    except KeyError:
        raise ValueError(f"Chain {chain_id} not found in structure")

def align_pdb(args):
    """
    Align target PDB to reference PDB using chain B
    Returns tuple of (target_name, rmsd, error)
    """
    reference_pdb, target_pdb, output_path = args
    
    try:
        # Parse structures
        parser = PDBParser(QUIET=True)
        ref_structure = parser.get_structure("reference", reference_pdb)
        target_structure = parser.get_structure("target", target_pdb)
        
        # Get CA atoms from chain B of both structures
        ref_atoms = get_chain_atoms(ref_structure, 'B')
        target_atoms = get_chain_atoms(target_structure, 'B')
        
        # Ensure we have the same number of atoms for alignment
        min_length = min(len(ref_atoms), len(target_atoms))
        ref_atoms = ref_atoms[:min_length]
        target_atoms = target_atoms[:min_length]
        
        # Calculate initial RMSD
        initial_rmsd = calculate_rmsd(ref_atoms, target_atoms)
        
        # Perform alignment
        superimposer = Superimposer()
        superimposer.set_atoms(ref_atoms, target_atoms)
        
        # Apply rotation/translation to all atoms in the target structure
        superimposer.apply(target_structure.get_atoms())
        
        # Calculate final RMSD
        final_rmsd = superimposer.rms
        
        # Log alignment results
        logger.info(f"Initial RMSD between {Path(reference_pdb).name} and {Path(target_pdb).name}: {initial_rmsd:.3f} Å")
        logger.info(f"Final RMSD between {Path(reference_pdb).name} and {Path(target_pdb).name}: {final_rmsd:.3f} Å")

        # Save aligned structure
        io = PDBIO()
        io.set_structure(target_structure)
        io.save(str(output_path))
        
        return (Path(target_pdb).name, final_rmsd, None)
        
    except Exception as e:
        return (Path(target_pdb).name, None, str(e))

def parse_arguments():
    """Parse command-line arguments"""
    parser = argparse.ArgumentParser(description="Align PDB files using BioPython")
    parser.add_argument("--input_dir", required=True, type=Path, help="Directory containing input PDB files")
    parser.add_argument("--output_dir", type=Path, default="aligned", 
                        help="Directory to save aligned PDB files (default: aligned)")
    parser.add_argument("--ncpus", type=int, default=1, help="Number of CPUs for parallel processing (default: 1)")
    parser.add_argument("--reference", type=str, help="Reference PDB file (default: first PDB in input_dir)")
    return parser.parse_args()

def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Get all PDB files from input directory
    pdb_files = list(args.input_dir.glob("*.pdb"))
    
    if not pdb_files:
        logger.error(f"No PDB files found in {args.input_dir}")
        sys.exit(0)
    
    # Determine reference PDB
    if args.reference:
        reference_pdb = args.input_dir / args.reference
        if not reference_pdb.is_file():
            logger.error(f"Reference PDB file not found: {reference_pdb}")
            sys.exit(0)
    else:
        reference_pdb = pdb_files[0]
    
    logger.info(f"Starting alignment of {len(pdb_files)} PDB files using {args.ncpus} CPUs")
    logger.info(f"Reference structure: {reference_pdb.name}")
    
    # Create output directory if it doesn't exist
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare arguments for parallel processing
    process_args = []
    for target_pdb in pdb_files:
        output_path = args.output_dir / f"{target_pdb.stem}.pdb"
        process_args.append((reference_pdb, target_pdb, output_path))
    
    # Process PDB files in parallel
    results = []
    with Pool(processes=args.ncpus) as pool:
        results = pool.map(align_pdb, process_args)
    
    # Log summary
    successful_alignments = sum(1 for _, rmsd, error in results if rmsd is not None)
    failed_alignments = sum(1 for _, rmsd, error in results if rmsd is None)
    
    logger.info("\nAlignment Summary:")
    logger.info(f"Total structures processed: {len(pdb_files)}")
    logger.info(f"Successful alignments: {successful_alignments}")
    logger.info(f"Failed alignments: {failed_alignments}")
    
    # Log details of failed alignments
    if failed_alignments > 0:
        logger.info("\nFailed alignments:")
        for name, _, error in results:
            if error:
                logger.info(f"  {name}: {error}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.error(f"Error: {str(e)}")
        sys.exit(0)
