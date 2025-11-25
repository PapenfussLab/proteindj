#!/usr/bin/env python3
import sys
import os
import re
import json
import logging
from pathlib import Path
import argparse
from multiprocessing import Pool
from Bio.PDB import PDBParser, PDBIO, Superimposer
from Bio.PDB.Selection import unfold_entities
from copy import deepcopy

def setup_logging():
    """Configure logging"""
    logger = logging.getLogger(__name__)
    log_file = "alignment.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logger

logger = setup_logging()

def get_all_ca_atoms(structure):
    """Collect all CA atoms from all chains in structure"""
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    ca_atoms.append(residue['CA'])
    if not ca_atoms:
        raise ValueError("No CA atoms found in structure")
    return ca_atoms

def get_chain_ca_atoms(structure, chain_id):
    """Collect CA atoms from specific chain in structure"""
    ca_atoms = []
    for model in structure:
        # Get all chains in the model
        for chain in model.get_chains():
            if chain.id == chain_id:
                for residue in chain:
                    if 'CA' in residue:
                        ca_atoms.append(residue['CA'])
    if not ca_atoms:
        raise ValueError(f"No CA atoms found in chain {chain_id}")
    return ca_atoms

def align_structures(args):
    """Align Boltz structure to Design template with chain-specific handling"""
    (design_path, boltz_path, out_pdb, src_json, dst_json, 
     fold_id, seq_id, design_type) = args  # Added design_type
    
    try:
        parser = PDBParser(QUIET=True)
        ref_structure = parser.get_structure("design", design_path)
        boltz_structure = parser.get_structure("boltz", boltz_path)

        if design_type == 'binder':
            # 1. Align chain B (target) for final structure
            ref_chainB = get_chain_ca_atoms(ref_structure, 'B')
            boltz_chainB = get_chain_ca_atoms(boltz_structure, 'B')
            
            superimposer = Superimposer()
            superimposer.set_atoms(ref_chainB, boltz_chainB)
            superimposer.apply(boltz_structure.get_atoms())
            rmsd_target = superimposer.rms 

            # 2. Calculate overall RMSD (all CA atoms)
            ref_all_ca = get_all_ca_atoms(ref_structure)
            boltz_all_ca = get_all_ca_atoms(boltz_structure)
            superimposer_all = Superimposer()
            superimposer_all.set_atoms(ref_all_ca, boltz_all_ca)
            rmsd_overall = superimposer_all.rms

            # 3. Calculate binder RMSD (chain A)
            ref_chainA = get_chain_ca_atoms(ref_structure, 'A')
            boltz_chainA = get_chain_ca_atoms(boltz_structure, 'A')
            superimposer_a = Superimposer()
            superimposer_a.set_atoms(ref_chainA, boltz_chainA)
            rmsd_binder = superimposer_a.rms

            rmsd_data = {
                "boltz_overall_rmsd": round(rmsd_overall, 2),
                "boltz_target_rmsd": round(rmsd_target, 2),
                "boltz_binder_rmsd": round(rmsd_binder, 2)
            }

        else:  # Monomer design
            ref_atoms = get_all_ca_atoms(ref_structure)
            boltz_atoms = get_all_ca_atoms(boltz_structure)
            
            superimposer = Superimposer()
            superimposer.set_atoms(ref_atoms, boltz_atoms)
            superimposer.apply(boltz_structure.get_atoms())
            
            rmsd_data = {
                "boltz_overall_rmsd": round(superimposer.rms, 2)
            }

        # Save aligned structure (always chain B aligned for binder)
        io = PDBIO()
        io.set_structure(boltz_structure)
        io.save(str(out_pdb))

        # Write new JSON with only the requested fields
        if src_json.exists():
            with open(src_json, 'r') as f:
                data = json.load(f)

            # Build output dictionary
            if design_type == 'binder':
                out_json = {
                    "fold_id": fold_id,
                    "seq_id": seq_id,
                    "description": boltz_path.name,
                    "boltz_overall_rmsd": round(data.get("boltz_overall_rmsd", rmsd_data.get("boltz_overall_rmsd", 0)), 2),
                    "boltz_target_rmsd": round(rmsd_data.get("boltz_target_rmsd", 0), 2),
                    "boltz_binder_rmsd": round(rmsd_data.get("boltz_binder_rmsd", 0), 2),
                    "boltz_conf_score": round(data.get("confidence_score", 0), 3),
                    "boltz_ptm": round(data.get("ptm", 0), 3),
                    "boltz_ptm_interface": round(data.get("iptm", 0), 3),
                    "boltz_plddt": round(data.get("complex_plddt", 0), 3),
                    "boltz_plddt_interface": round(data.get("complex_iplddt", 0), 3),
                    "boltz_pde": round(data.get("complex_pde", 0), 2),
                    "boltz_pde_interface": round(data.get("complex_ipde", 0), 2)
                }
            else:
                out_json = {
                    "fold_id": fold_id,
                    "seq_id": seq_id,
                    "description": boltz_path.name,
                    "boltz_overall_rmsd": round(data.get("boltz_overall_rmsd", rmsd_data.get("boltz_overall_rmsd", 0)), 2),
                    "boltz_conf_score": round(data.get("confidence_score", 0), 3),
                    "boltz_ptm": round(data.get("ptm", 0), 3),
                    "boltz_plddt": round(data.get("complex_plddt", 0), 3),
                    "boltz_pde": round(data.get("complex_pde", 0), 2),
                }

            with open(dst_json, 'w') as f:
                json.dump(out_json, f, indent=2)

        return (boltz_path.name, rmsd_data.get('boltz_overall_rmsd'), None)

    except Exception as e:
        logger.error(f"Failed {boltz_path.name}: {str(e)}")
        return (boltz_path.name, None, str(e))

def main():
    parser = argparse.ArgumentParser(description="Align Boltz predictions to designs")
    parser.add_argument("--design_dir", type=Path, required=True, 
                      help="Directory with Design PDBs (fold_*_seq_*.pdb)")
    parser.add_argument("--boltz_dir", type=Path, required=True,
                      help="Directory with Boltz PDBs and JSONs (fold_*_seq_*_boltzpred.pdb)")
    parser.add_argument("--output_dir", type=Path, default="aligned",
                      help="Output directory for results")
    parser.add_argument("--design_type", choices=['binder', 'monomer'], required=True,
                      help="Design type: 'binder' (A/B chains) or 'monomer (A chain)'")
    parser.add_argument("--ncpus", type=int, default=1,
                      help="Number of CPUs for parallel processing")
    args = parser.parse_args()
    
    # Validate input directories
    if not args.design_dir.exists():
        logger.error(f"Design directory not found: {args.design_dir}")
        sys.exit(1)
        
    if not args.boltz_dir.exists():
        logger.error(f"Boltz directory not found: {args.boltz_dir}")
        sys.exit(1)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Map Design files
    design_files = {}
    for design_file in args.design_dir.glob("fold_*_seq_*.pdb"):
        match = re.match(r"fold_(\d+)_seq_(\d+)\.pdb", design_file.name)
        if match:
            fold_id = int(match.group(1))
            seq_id = int(match.group(2))
            design_files[(fold_id, seq_id)] = design_file
    
    # Prepare processing tasks
    tasks = []
    for boltz_file in args.boltz_dir.glob("fold_*_seq_*_boltzpred.pdb"):
        match = re.match(r"fold_(\d+)_seq_(\d+)_.*\.pdb", boltz_file.name)
        if not match:
            continue
            
        fold_id = int(match.group(1))
        seq_id = int(match.group(2))
        key = (fold_id, seq_id)
        
        if key not in design_files:
            logger.warning(f"No design file for fold {fold_id} seq {seq_id}, skipping {boltz_file.name}")
            continue
            
        # Generate paths
        base_name = boltz_file.stem  # fold_X_seq_Y_boltzpred
        src_json = args.boltz_dir / f"{base_name}.json"
        out_pdb = args.output_dir / f"{base_name}.pdb"
        dst_json = args.output_dir / f"{base_name}.json"
        
        tasks.append((
            design_files[key],
            boltz_file,
            out_pdb,
            src_json,
            dst_json,
            fold_id,
            seq_id,
            args.design_type
        ))
    
    # Log processing start
    logger.info(f"Starting alignment of {len(tasks)} Boltz structures")
    logger.info(f"Using design directory: {args.design_dir}")
    logger.info(f"Output directory: {args.output_dir}")
    
    # Process tasks in parallel
    with Pool(args.ncpus) as pool:
        results = pool.map(align_structures, tasks)
    
    # Report summary
    successes = sum(1 for r in results if r[1] is not None)
    failures = len(results) - successes
    
    logger.info("\n=== Alignment Summary ===")
    logger.info(f"Total structures processed: {len(tasks)}")
    logger.info(f"Successful alignments: {successes}")
    logger.info(f"Failed alignments: {failures}")
    
    if failures > 0:
        logger.info("\nFailed cases:")
        for name, _, error in results:
            if error:
                logger.info(f"  {name}: {error}")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error(f"Fatal error: {str(e)}")
        sys.exit(1)
