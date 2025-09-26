#!/usr/bin/env python3
import os
import argparse
import json
import shutil

def load_json_scores(json_dir):
    """Load FAMPNN average PSCE scores from JSON metadata files"""
    scores_map = {}
    
    for json_file in os.listdir(json_dir):
        if json_file.endswith('.json'):
            with open(os.path.join(json_dir, json_file)) as f:
                for line in f:
                    try:
                        data = json.loads(line.strip())
                        design_name = data['design']
                        scores_map[design_name] = float(data['fampnn_avg_psce'])
                    except (json.JSONDecodeError, KeyError) as e:
                        print(f"Warning: Error parsing line in {json_file}: {e}")
    return scores_map

def filter_fampnn_scores(scores_map, max_score):
    """Filter FAMPNN designs based on fampnn_avg_psce"""
    if max_score is None:
        return scores_map
    return {design: score for design, score in scores_map.items() if score <= max_score}

def copy_filtered_designs(filtered_designs, pdb_dir, json_dir, output_dir):
    """Copy matching PDBs and JSONs to output directory"""
    os.makedirs(output_dir, exist_ok=True)
    copied_count = 0
    
    for design in filtered_designs:
        # Copy PDB file
        pdb_file = os.path.join(pdb_dir, f"{design}.pdb")
        if os.path.exists(pdb_file):
            shutil.copy2(pdb_file, output_dir)
            copied_count += 1
            
        # Copy JSON metadata
        json_file = os.path.join(json_dir, f"{design}.json")  # Now matches 'design' field
        if os.path.exists(json_file):
            shutil.copy2(json_file, output_dir)
        else:
            print(f"Warning: JSON file for {design} not found")
    
    return copied_count

def main():
    parser = argparse.ArgumentParser(description="Filter FAMPNN designs using JSON metadata")
    parser.add_argument("--jsons", required=True, help="Directory containing JSON metadata files")
    parser.add_argument("--pdbs", required=True, help="Directory containing PDB files")
    parser.add_argument("--fampnn-max-psce", type=float,
                      help="Maximum FAMPNN average PSCE score (optional; copies all designs if not provided)")
    parser.add_argument("--output-dir", default="filtered_output_fampnn",
                      help="Output directory for filtered designs (default: filtered_output_fampnn)")
    
    args = parser.parse_args()
    
    if args.fampnn_max_psce is not None:
        print(f"Filtering designs with FAMPNN average PSCE score ≤ {args.fampnn_max_psce}")
    else:
        print("No FAMPNN score filtering applied; copying all designs.")
    
    # Load and filter scores
    scores = load_json_scores(args.jsons)
    print(f'Pre-filter designs: {len(scores)}')
    
    filtered = filter_fampnn_scores(scores, args.fampnn_max_psce)
    print(f'Post-filter designs: {len(filtered)}')
  
    # Copy matching files
    copied_count = copy_filtered_designs(filtered.keys(), args.pdbs, args.jsons, args.output_dir)
    
    print(f"\nResults: {len(filtered)} designs found, {copied_count} PDB files copied to {args.output_dir}")

if __name__ == "__main__":
    main()
