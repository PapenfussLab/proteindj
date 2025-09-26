import pyrosetta
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Restores side-chains to PDB files after RFdiffusion processing')
    parser.add_argument('--input_dir', required=True, help='Input directory containing PDB files')
    parser.add_argument('--out_dir', default='./outputpdbs', help='Output directory for updated PDB files')
    args = parser.parse_args()
    
    # keep pyrosetta quiet
    pyrosetta.init("-out:levels all:error")
        
    # Make list of input PDBs
    input_dir = Path(args.input_dir)
    pdb_files = list(input_dir.glob("*.pdb"))

    # Create the ouput directory
    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True)

    for pdb_file in pdb_files:
        # import design. pyrosetta will automatically restore missing side-chains
        pose_design = pyrosetta.pose_from_pdb(str(pdb_file))
        # output designs
        output_path = out_dir / pdb_file.name
        print(f"Outputting PDB file: {output_path}")
        pose_design.dump_pdb(str(output_path))

if __name__ == "__main__":
    main()
