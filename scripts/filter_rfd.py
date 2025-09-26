import argparse
import os
from pathlib import Path
import pyrosetta as pr
from pyrosetta.rosetta.core.scoring import dssp
from multiprocessing import Pool
import logging
import sys
import uuid
import json
import shutil
import re

def setup_logger():
    """Configure logging to output to both file and stdout"""
    logger = logging.getLogger('filter_rfd')
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    
    # Create unique log filename
    unique_id = str(uuid.uuid4())[:8]
    log_filename = f'filter_rfd_{unique_id}.log'
    file_handler = logging.FileHandler(log_filename)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    return logger

def extract_fold_id(pdb_filename):
    """Extract fold ID from PDB filename"""
    match = re.search(r'fold_(\d+)', pdb_filename.name)
    if match:
        return int(match.group(1))
    else:
        return None

def analyze_structure(args):
    """Analyze structure with automatic chain detection"""
    (pdb_file, rfd_min_ss, rfd_max_ss, rfd_min_helices, rfd_max_helices, 
     rfd_min_strands, rfd_max_strands, rfd_min_rog, rfd_max_rog, 
     output_dir, json_dir) = args

    logger = logging.getLogger('filter_rfd')

    try:
        # Load structure with PyRosetta
        pose = pr.pose_from_pdb(str(pdb_file))

        # Autodetect chain structure
        chains = pose.split_by_chain()
        num_chains = len(chains)

        if num_chains == 1:
            # Monomer - analyze entire structure (single chain)
            pose_to_analyze = pose
            logger.info(f"{pdb_file.name}: Single chain detected - treating as monomer")
        elif num_chains == 2:
            # Binder - analyze first chain only
            pose_to_analyze = chains[1]
            logger.info(f"{pdb_file.name}: Two chains detected ({num_chains}) - treating as binder, analyzing first chain only")
        elif num_chains >= 3:
            # Oligomer - analyse all chains
            pose_to_analyze = pose
            logger.info(f"{pdb_file.name}: More than two chains detected - treating as oligomer, analysing all chains")
        else:
            logger.error(f"{pdb_file.name}: No chains found, skipping")
            return None

        # Calculate radius of gyration
        scorefxn = pr.ScoreFunction()
        scorefxn.set_weight(pr.rosetta.core.scoring.rg, 1.0)
        rog = round(scorefxn(pose_to_analyze), 2)

        # Count secondary structures
        dssp_obj = dssp.Dssp(pose_to_analyze)
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

        total_ss = helix_count + strand_count
        fold_id = extract_fold_id(pdb_file)

        # Helper function for cleaner filter checking
        def passes_filter(value, min_val, max_val):
            if min_val is not None and max_val is not None:
                return min_val <= value <= max_val
            if min_val is not None:
                return value >= min_val
            if max_val is not None:
                return value <= max_val
            return True

        # Apply all filters
        passes_ss_filter = passes_filter(total_ss, rfd_min_ss, rfd_max_ss)
        passes_helix_filter = passes_filter(helix_count, rfd_min_helices, rfd_max_helices)
        passes_strand_filter = passes_filter(strand_count, rfd_min_strands, rfd_max_strands)
        passes_rog_filter = passes_filter(rog, rfd_min_rog, rfd_max_rog)

        passes_all_filters = (passes_ss_filter and passes_helix_filter and 
                             passes_strand_filter and passes_rog_filter)

        # Create filter description for logging
        applied_filters = []
        if rfd_min_ss is not None or rfd_max_ss is not None:
            applied_filters.append(f"SS: {rfd_min_ss or 'None'}-{rfd_max_ss or 'None'}")
        if rfd_min_helices is not None or rfd_max_helices is not None:
            applied_filters.append(f"Helix: {rfd_min_helices or 'None'}-{rfd_max_helices or 'None'}")
        if rfd_min_strands is not None or rfd_max_strands is not None:
            applied_filters.append(f"Strand: {rfd_min_strands or 'None'}-{rfd_max_strands or 'None'}")
        if rfd_min_rog is not None or rfd_max_rog is not None:
            applied_filters.append(f"RoG: {rfd_min_rog or 'None'}-{rfd_max_rog or 'None'}")

        filters_str = ", ".join(applied_filters) if applied_filters else "No filters applied"

        logger.info(
            f"{pdb_file.name}: Analysis complete - "
            f"Helices={helix_count}, Strands={strand_count}, "
            f"Total SS={total_ss}, RoG={rog} Å, "
            f"Filters=[{filters_str}], Passed={passes_all_filters}"
        )

        # Copy files if they pass filters
        if passes_all_filters:
            # Copy PDB file
            output_path = output_dir / pdb_file.name
            shutil.copy2(pdb_file, output_path)

            # Copy corresponding JSON file if it exists
            if json_dir:
                json_file = Path(json_dir) / f"{pdb_file.stem}.json"
                if json_file.exists():
                    json_output_path = output_dir / json_file.name
                    shutil.copy2(json_file, json_output_path)
                    logger.info(f"{json_file.name}: Corresponding JSON copied to output directory")
                else:
                    logger.warning(f"{json_file.name}: Corresponding JSON file not found in {json_dir}")

            logger.info(f"{pdb_file.name}: Copied to output directory")

        return {
            "fold_id": fold_id,
            "rfd_helices": helix_count,
            "rfd_strands": strand_count,
            "rfd_total_ss": total_ss,
            "rfd_RoG": rog,
        }

    except Exception as e:
        logger.error(f"{pdb_file.name}: Failed - {str(e)}", exc_info=True)
        return None

def main():
    parser = argparse.ArgumentParser(description='Filter PDB files based on secondary structure content')
    parser.add_argument('--input-dir', type=str, help='Directory containing PDB files')
    parser.add_argument('--json-dir', type=str, default=None, 
                       help='Directory containing JSON files corresponding to PDB files (default: same as input-dir)')
    parser.add_argument('--rfd-min-ss', type=int, default=None, help='Minimum total secondary structure elements')
    parser.add_argument('--rfd-max-ss', type=int, default=None, help='Maximum total secondary structure elements')
    parser.add_argument('--rfd-min-helices', type=int, default=None, help='Minimum number of alpha helices')
    parser.add_argument('--rfd-max-helices', type=int, default=None, help='Maximum number of alpha helices')
    parser.add_argument('--rfd-min-strands', type=int, default=None, help='Minimum number of beta strands')
    parser.add_argument('--rfd-max-strands', type=int, default=None, help='Maximum number of beta strands')
    parser.add_argument('--rfd-min-rog', type=float, default=None, help='Minimum radius of gyration')
    parser.add_argument('--rfd-max-rog', type=float, default=None, help='Maximum radius of gyration')
    parser.add_argument('--ncpus', type=int, default=1, help='Number of CPU cores to use')
    parser.add_argument('--output-dir', type=Path, default="filtered", 
        help="Directory to save filtered PDB files (default: filtered)")
    args = parser.parse_args()
    
    # setup logging
    logger = setup_logger()
    logger.info(f"Starting analysis with parameters: {vars(args)}")
    
    # initialise pyrosetta
    pr.init("-out:levels all:error")
    
    # create output directory
    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True)
    
    # Set json_dir to input_dir if not specified
    json_dir = args.json_dir if args.json_dir else args.input_dir
    
    pdb_files = list(Path(args.input_dir).glob('*.pdb'))
    if not pdb_files:
        logger.error(f"No PDB files found in {args.input_dir}")
        return
    
    logger.info(f"Found {len(pdb_files)} PDB files to analyze")
    logger.info(f"Looking for JSON files in {json_dir}")
    
    process_args = [(
        pdb, 
        args.rfd_min_ss, 
        args.rfd_max_ss,
        args.rfd_min_helices, 
        args.rfd_max_helices,
        args.rfd_min_strands, 
        args.rfd_max_strands,
        args.rfd_min_rog,
        args.rfd_max_rog,
        args.output_dir,
        json_dir
    ) for pdb in pdb_files]
    
    with Pool(processes=args.ncpus) as pool:
        results = pool.map(analyze_structure, process_args)
    
    # filter out None results
    valid_results = [result for result in results if result is not None]
    
    # save analysis data to JSONL
    output_filename = f'rfd_data_{str(uuid.uuid4())[:8]}.jsonl'
    with open(output_filename, 'w') as f:
        for result in valid_results:
            json.dump(result, f)
            f.write('\n')
    
    logger.info(f"Analysis complete. Results saved to {output_filename}")
    
    # Print summary statistics
    passed_pdbs = len(list(output_dir.glob('*.pdb')))
    passed_jsons = len(list(output_dir.glob('*.json')))
    logger.info(f"Summary: {passed_pdbs} out of {len(valid_results)} structures passed all filters")
    logger.info(f"Copied {passed_jsons} corresponding JSON files to the output directory")

if __name__ == "__main__":
    main()
