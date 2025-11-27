#!/usr/bin/env python3

import pandas as pd
import json
from Bio import PDB
from Bio.PDB import Polypeptide
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from functools import partial
from pathlib import Path
import logging
from multiprocessing import Pool
import pyrosetta as pr
from typing import List, Dict, Tuple, Optional
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.simple_filters import BuriedUnsatHbondFilter
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector
from pyrosetta.rosetta.core.scoring import dssp
import os
import traceback
import re

# setup logging - file and stream
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('analysis.log', mode='w'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def derive_ids_from_filename(filename):
    """Extract fold_id and seq_id from filename format: fold_X_seq_Y_*.pdb"""
    basename = Path(filename).stem
    fold_match = re.search(r'fold_(\d+)', basename)
    seq_match = re.search(r'seq_(\d+)', basename)
    
    fold_id = int(fold_match.group(1)) if fold_match else None
    seq_id = int(seq_match.group(1)) if seq_match else None
    
    return fold_id, seq_id

def calculate_interface_metrics(pose, chain1='A', chain2='B'):
    """Calculate interface metrics between specified chains"""
    # Analyze interface statistics
    iam = InterfaceAnalyzerMover()
    iam.set_interface(f"{chain1}_{chain2}")
    scorefxn = pr.get_fa_scorefxn()
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_pack_separated(True)
    iam.set_compute_interface_energy(True)
    iam.set_calc_hbond_sasaE(True)
    iam.set_compute_interface_sc(True)
    iam.apply(pose)

    # Retrieve statistics
    metrics = iam.get_all_data()
    
    # Add suffix for metrics if not the default A_B interface
    suffix = '' if (chain1 == 'A' and chain2 == 'B') else f'_{chain1}_{chain2}'
    
    return {
        f'pr_intface_BSA{suffix}': round(metrics.dSASA[1] / 2),
        f'pr_intface_shpcomp{suffix}': round(metrics.sc_value, 3),
        f'pr_intface_deltaG{suffix}': round(metrics.dG[1], 1),
        f'pr_intface_deltaGtoBSA{suffix}': round(metrics.dG_dSASA_ratio, 3),
        f'pr_intface_hbonds{suffix}': metrics.interface_hbonds,
        f'pr_intface_unsat_hbonds{suffix}': metrics.delta_unsat_hbonds,
        f'pr_intface_packstat{suffix}': round(metrics.packstat, 3)
    }

def get_chain_ids(pose):
    """Get unique chain IDs from a pose using PDBInfo or chain indices"""
    chain_ids = []
    if pose.pdb_info():
        # Get chain IDs from PDB information
        seen_chains = set()
        for res in range(1, pose.total_residue() + 1):
            chain = pose.pdb_info().chain(res)
            if chain not in seen_chains:
                seen_chains.add(chain)
                chain_ids.append(chain)
    else:
        # Fallback to numeric chain indices
        chain_ids = [str(i) for i in range(1, pose.num_chains() + 1)]
    return chain_ids

def calculate_chain_metrics(pose, chain_id, pdb_path):
    """Calculate all metrics for a single chain"""
    # Find chain index using PDBInfo or position
    chain_index = 1
    for i in range(1, pose.num_chains() + 1):
        start_res = pose.chain_begin(i)
        current_chain_id = str(pose.pdb_info().chain(start_res)) if pose.pdb_info() else str(i)
        if current_chain_id == chain_id:
            chain_index = i
            break
    
    # Split pose into chains
    chains = pose.split_by_chain()
    if chain_index > len(chains):
        raise ValueError(f"Chain {chain_id} not found in pose")
    
    chain_pose = chains[chain_index]
       
    # Add suffix for metrics if not chain A
    suffix = '' if chain_id == 'A' else f'_{chain_id}'
    
    # Calculate secondary structure metrics
    ss_metrics = count_secondary_structures(chain_pose)
    ss_metrics = {f'{k}{suffix}': v for k, v in ss_metrics.items()}
    
    # Calculate RoG
    rog_metrics = calculate_rog(chain_pose, chain_id)
    rog_metrics = {f'{k}{suffix}': v for k, v in rog_metrics.items()}
    
    # Calculate surface chemistry
    surface_metrics = calculate_surface_chemistry(chain_pose)
    surface_metrics = {f'{k}{suffix}': v for k, v in surface_metrics.items()}
    
    # Calculate total energy metric
    tem_metrics = calculate_tem(chain_pose)
    tem_metrics = {f'{k}{suffix}': v for k, v in tem_metrics.items()}
    
    # Get sequence metrics
    seq_metrics = {}
    sequence = get_chain_sequence(pdb_path, chain_id)
    if sequence:
        seq_metrics = calculate_seq_metrics(sequence)
        seq_metrics = {f'{k}{suffix}': v for k, v in seq_metrics.items()}
    
    # Combine all metrics
    metrics = {}
    metrics.update(ss_metrics)
    metrics.update(rog_metrics)
    metrics.update(surface_metrics)
    metrics.update(tem_metrics)
    metrics.update(seq_metrics)
    
    return metrics

def calculate_rog(pose, chain_id):
    """Calculate radius of gyration for a pose"""
    selected_chain = ChainSelector(chain_id).apply(pose)
    center_of_mass = pr.rosetta.core.pose.center_of_mass(pose, selected_chain)
    rog = pr.rosetta.core.pose.radius_of_gyration(pose, center_of_mass, selected_chain)
    
    return {
        'pr_RoG': round(rog, 2)
    }

def calculate_surface_chemistry(pose):
    """Calculate surface chemistry metrics"""
    layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    layer_sel.set_layers(pick_core=False, pick_boundary=False, pick_surface=True)
    surface_res = layer_sel.apply(pose)

    surface_hydrophobic_count = 0
    surface_polar_count = 0
    total_count = 0

    # count apolar and aromatic residues at the surface
    for i in range(1, len(surface_res) + 1):
        if surface_res[i] == True:
            res = pose.residue(i)
            if res.is_apolar() == True or res.name() in ['PHE', 'TRP', 'TYR']:
                surface_hydrophobic_count += 1
            total_count += 1

    return {
        'pr_surfhphobics': round(surface_hydrophobic_count * 100 / total_count, 1)
    }

def calculate_tem(pose):
    """Calculate total energy metric for a pose"""
    scorefxn = pr.get_fa_scorefxn()
    tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    tem.set_scorefunction(scorefxn)
    
    return {
        'pr_TEM': round(tem.calculate(pose))
    }

def count_secondary_structures(pose):
    """Count secondary structures in a pose"""
    dssp_obj = dssp.Dssp(pose)
    dssp_string = dssp_obj.get_dssp_secstruct()
    
    helix_count = 0
    strand_count = 0
    current_helix = False
    current_strand = False
    
    for ss in dssp_string:
        if ss == 'H':
            if not current_helix:
                helix_count += 1
                current_helix = True
            current_strand = False
        elif ss == 'E':
            if not current_strand:
                strand_count += 1
                current_strand = True
            current_helix = False
        else:
            current_helix = False
            current_strand = False
    
    return {
        'pr_helices': helix_count,
        'pr_strands': strand_count,
        'pr_total_ss': helix_count + strand_count
    }

def get_chain_sequence(pdb_path, chain_id):
    """Extract sequence for a specific chain from PDB file"""
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('temp', str(pdb_path))
        
        # Convert numeric chain IDs to letters if needed
        if chain_id.isdigit():
            chain_id = chr(64 + int(chain_id))  # Convert 1->A, 2->B etc
            
        for chain in structure.get_chains():
            if chain.id == chain_id:
                residues = [r for r in chain.get_residues() if PDB.is_aa(r)]
                return ''.join([Polypeptide.protein_letters_3to1[r.get_resname()] for r in residues])
    except Exception as e:
        logger.error(f"Error getting sequence: {str(e)}")
    return None

def calculate_seq_metrics(sequence):
    """Calculate metrics for a protein sequence"""
    analysis = ProteinAnalysis(sequence)
    extinction_coef = analysis.molar_extinction_coefficient()[0]
    molecular_weight = int(analysis.molecular_weight())

    return {
        'sequence': sequence,
        'seq_ext_coef': extinction_coef,
        'seq_length': len(sequence),
        'seq_MW': molecular_weight,
        'seq_pI': round(analysis.isoelectric_point(), 2)
    }

def calculate_whole_pose_metrics(pose):
    """Calculate metrics for the entire pose (all chains combined)"""
    # Create residue selector vector for all residues
    all_residues = pr.rosetta.utility.vector1_bool(pose.total_residue())
    for i in range(1, pose.total_residue() + 1):
        all_residues[i] = True  # Select all residues
        
    # Calculate center of mass
    com = pr.rosetta.core.pose.center_of_mass(pose, all_residues)

    # Radius of gyration for full structure
    rog = pr.rosetta.core.pose.radius_of_gyration(pose, com, all_residues)

    # Surface chemistry for full structure
    surface_metrics = calculate_surface_chemistry(pose)
    
    # Total energy for full structure
    tem_metrics = calculate_tem(pose)
    
    return {
        'pr_RoG_total': round(rog, 2),
        **{f'pr_surfhphobics_total': surface_metrics['pr_surfhphobics']},
        **tem_metrics
    }

def process_single_pdb(pdb_path):
    """Process a single PDB file and return all calculated metrics"""
    
    logger.info(f"Processing PDB file: {pdb_path}")
    try:
        pose = pr.pose_from_pdb(str(pdb_path))
        
        # Initialize PDBInfo if missing
        if not pose.pdb_info():
            pose.pdb_info(pr.rosetta.core.pose.PDBInfo(pose))
        
        metrics = {'description': pdb_path.stem}
        chain_ids = get_chain_ids(pose)

        # Handle metrics calculation based on design type
        if len(chain_ids) == 1:
            # Monomer design            
            # Calculate metrics for the first only chain
            primary_chain = chain_ids[0]
            chain_metrics = calculate_chain_metrics(pose, primary_chain, pdb_path)
            metrics.update(chain_metrics)

        elif len(chain_ids) == 2:
            # Calculate binder chain (A) sequence metrics
            sequence = get_chain_sequence(pdb_path, 'A')
            if sequence:
                seq_metrics = calculate_seq_metrics(sequence)
                metrics.update(seq_metrics)
            
            # Calculate interface metrics for A-B
            interface_metrics = calculate_interface_metrics(pose, 'A', 'B')
            metrics.update(interface_metrics)
            
            # Calculate chain metrics for binder chain only
            chain_metrics = calculate_chain_metrics(pose, 'A', pdb_path)
            metrics.update(chain_metrics)
            
        elif len(chain_ids) >= 3:
            # Oligomer design: Handle any number of chains, calculate all interfaces
            
            # Calculate interface metrics for all chain pairs if multi-chain
            if len(chain_ids) > 1:
                for i, chain1 in enumerate(chain_ids):
                    for chain2 in chain_ids[i+1:]:
                        interface_metrics = calculate_interface_metrics(pose, chain1, chain2)
                        metrics.update(interface_metrics)

            # Calculate per-chain metrics
            all_chain_metrics = {}
            total_helices = 0
            total_strands = 0
            total_ss = 0
            
            for chain_id in chain_ids:
                chain_metrics = calculate_chain_metrics(pose, chain_id, pdb_path)
                
                # Track secondary structure aggregates
                suffix = '' if chain_id == 'A' else f'_{chain_id}'
                total_helices += chain_metrics.get(f'pr_helices{suffix}', 0)
                total_strands += chain_metrics.get(f'pr_strands{suffix}', 0)
                total_ss += chain_metrics.get(f'pr_total_ss{suffix}', 0)
                
                all_chain_metrics.update(chain_metrics)

            # Add whole-pose metrics
            whole_pose_metrics = calculate_whole_pose_metrics(pose)
            all_chain_metrics.update(whole_pose_metrics)
            
            # Add aggregated secondary structure metrics
            all_chain_metrics.update({
                'pr_helices_allchains': total_helices,
                'pr_strands_allchains': total_strands,
                'pr_total_ss_allchains': total_ss
            })
            
            metrics.update(all_chain_metrics)
                   
        return metrics
    except Exception as e:
        logger.error(f"Error processing {pdb_path}: {str(e)}\n{traceback.format_exc()}")
        return {}

def read_jsonl(file_path):
    """Read JSONL file and return list of records"""
    records = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError as e:
                    logger.error(f"Error parsing JSONL line: {e}")
    return records

def write_jsonl(file_path, records):
    """Write records to JSONL file"""
    with open(file_path, 'w') as f:
        for record in records:
            f.write(json.dumps(record) + '\n')

def process_pdbs(pdb_dir=None, output_path=None, num_processes=4):
    """Main function to process multiple PDB files in parallel and create enriched JSONL"""
    
    logger.info(f"Starting processing of PDB files in directory: {pdb_dir}")
    try:
        # Get all PDB files in the directory
        pdb_paths = list(Path(pdb_dir).glob('*.pdb'))
        logger.info(f"Found {len(pdb_paths)} PDB files to process")
        
        if not pdb_paths:
            logger.error(f"No PDB files found in directory: {pdb_dir}")
            return []
        
        # Process PDB files in parallel
        logger.info(f"Processing PDB files with {num_processes} processes")
        process_func = partial(process_single_pdb)
        
        with Pool(processes=num_processes) as pool:
            results = pool.map(process_func, pdb_paths)
        
        # Create enriched records with fold_id and seq_id derived from filenames
        enriched_records = []
        for i, pdb_path in enumerate(pdb_paths):
            result = results[i]
            if result:  # Skip failed processing results
                # Extract fold_id and seq_id from filename
                fold_id, seq_id = derive_ids_from_filename(pdb_path.name)
                
                # Create record with derived IDs and calculated metrics
                record = {
                    'description': pdb_path.stem,
                    'fold_id': fold_id,
                    'seq_id': seq_id
                }
                
                # Add all calculated metrics (excluding duplicate description)
                for key, value in result.items():
                    if key != 'description':
                        record[key] = value
                
                enriched_records.append(record)
            else:
                logger.warning(f"Failed to process {pdb_path}, skipping")
        
        # Save enriched records to output JSONL
        logger.info(f"Saving enriched JSONL to {output_path}")
        write_jsonl(output_path, enriched_records)
        logger.info(f"Processing completed successfully. Processed {len(enriched_records)} files")
        
        return enriched_records
        
    except Exception as e:
        logger.error(f"Error in main processing: {str(e)}\n{traceback.format_exc()}")
        return []

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Process PDB files and create enriched JSONL')
    parser.add_argument('--pdb_dir', required=True, help='Directory containing PDB files')
    parser.add_argument('--output', required=True, help='Path for output enriched JSONL file')
    parser.add_argument('--num_processes', type=int, default=4, help='Number of processes to use')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Initialize PyRosetta
    logger.info("Initializing PyRosetta")
    pr.init("-out:levels all:error")
    
    # Process files and save to specified output path
    enriched_records = process_pdbs(
        pdb_dir=args.pdb_dir,
        output_path=args.output,
        num_processes=args.num_processes
    )