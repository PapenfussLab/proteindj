#!/usr/bin/env python3
"""
Generate RFdiffusion contig strings from PDB files.
Handles missing residues and chainbreaks according to RFdiffusion requirements.
Supports different design modes: binder_denovo, binder_partialdiff, monomer_partialdiff, monomer_denovo.
"""

import argparse
from Bio.PDB import PDBParser
from collections import defaultdict
import re


def get_protein_chains(structure):
    """
    Extract protein chains from structure, ignoring ligands and non-protein residues.
    
    Returns:
        dict: Dictionary with chain IDs as keys and lists of residue numbers as values
    """
    chain_residues = defaultdict(list)
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                # Filter out hetero residues (ligands, waters) and keep only amino acids
                hetfield = residue.id[0]
                if hetfield == ' ':  # Standard amino acid residue
                    resnum = residue.id[1]
                    chain_residues[chain_id].append(resnum)
    
    # Sort residue numbers for each chain
    for chain_id in chain_residues:
        chain_residues[chain_id].sort()
    
    return chain_residues


def find_continuous_ranges(residue_list):
    """
    Find continuous residue ranges, breaking at gaps (missing residues).
    
    Args:
        residue_list: Sorted list of residue numbers
        
    Returns:
        list: List of tuples (start, end) representing continuous ranges
    """
    if not residue_list:
        return []
    
    ranges = []
    start = residue_list[0]
    prev = residue_list[0]
    
    for resnum in residue_list[1:]:
        # If there's a gap of more than 1, start a new range
        if resnum != prev + 1:
            ranges.append((start, prev))
            start = resnum
        prev = resnum
    
    # Add the final range
    ranges.append((start, prev))
    
    return ranges


def format_chain_contig(chain_id, ranges, include_chain_id=True):
    """
    Format contig string for a single chain.
    
    Args:
        chain_id: Chain identifier
        ranges: List of (start, end) tuples
        include_chain_id: Whether to include chain ID in output
        
    Returns:
        str: Formatted contig string for the chain
    """
    chain_contigs = []
    for start, end in ranges:
        if include_chain_id:
            if start == end:
                chain_contigs.append(f"{chain_id}{start}")
            else:
                chain_contigs.append(f"{chain_id}{start}-{end}")
        else:
            if start == end:
                chain_contigs.append(f"{start}")
            else:
                chain_contigs.append(f"{start}-{end}")
    
    return "/".join(chain_contigs)


def validate_design_length(design_length):
    """
    Validate design_length format.
    
    Args:
        design_length: String in format 'N' or 'N-M'
        
    Returns:
        bool: True if valid
        
    Raises:
        ValueError: If format is invalid
    """
    if not design_length:
        return True
    
    # Check for single number or range
    pattern = r'^\d+(-\d+)?$'
    if not re.match(pattern, design_length):
        raise ValueError(
            f"Invalid design_length format: '{design_length}'. "
            "Must be a single number (e.g., '60') or range (e.g., '70-80')"
        )
    
    # If it's a range, validate that start <= end
    if '-' in design_length:
        start, end = map(int, design_length.split('-'))
        if start > end:
            raise ValueError(
                f"Invalid design_length range: '{design_length}'. "
                "Start must be less than or equal to end"
            )
    
    return True


def generate_contig_binder_denovo(chain_residues, design_length=None):
    """
    Generate contig string for binder_denovo mode.
    Full chain notation with chain IDs for all chains, optionally with design_length appended.
    
    Example: '[A1-77/0 B23-77/B80-105/0 60]' or '[A1-77/0 B23-77/B80-105/0 70-80]'
    """
    contig_parts = []
    
    for chain_id in sorted(chain_residues.keys()):
        residues = chain_residues[chain_id]
        ranges = find_continuous_ranges(residues)
        chain_string = format_chain_contig(chain_id, ranges, include_chain_id=True)
        contig_parts.append(chain_string)
    
    # Join chains with /0 (space before next chain)
    contig = "/0 ".join(contig_parts) + "/0"
    
    # Append design_length if provided
    if design_length:
        contig += f" {design_length}"
    
    return f"[{contig}]"


def generate_contig_monomer_denovo(design_length=None):
    """
    Generate contig string for monomer_denovo mode.
    Uses design_length exclusively if provided.
    
    Example: '[60]' or '[70-80]'
    """
    if not design_length:
        raise ValueError(
            "design_length is required for monomer_denovo mode. "
            "Provide a length (e.g., '60') or range (e.g., '70-80')"
        )
    
    return f"[{design_length}]"


def generate_contig_monomer_partialdiff(chain_residues):
    """
    Generate contig string for monomer_partialdiff mode.
    Total number of residues without chain IDs.
    
    Example: '[155-155]'
    """
    total_residues = sum(len(residues) for residues in chain_residues.values())
    return f"[{total_residues}-{total_residues}]"


def generate_contig_binder_partialdiff(chain_residues):
    """
    Generate contig string for binder_partialdiff mode.
    Chain A as residue count without chain ID, chain B with full notation including breaks.
    
    Example: "[88-88/0 B89-203]"
    """
    sorted_chains = sorted(chain_residues.keys())
    
    if len(sorted_chains) < 2:
        raise ValueError(
            "binder_partialdiff mode requires at least 2 chains (binder and target)"
        )
    
    # Chain A: just count of residues without chain ID
    chain_a = sorted_chains[0]
    chain_a_count = len(chain_residues[chain_a])
    contig_parts = [f"{chain_a_count}-{chain_a_count}"]
    
    # Chain B (and others): full notation with chain IDs and breaks
    for chain_id in sorted_chains[1:]:
        residues = chain_residues[chain_id]
        ranges = find_continuous_ranges(residues)
        chain_string = format_chain_contig(chain_id, ranges, include_chain_id=True)
        contig_parts.append(chain_string)
    
    # Join with /0 separator
    return f"[{'/0 '.join(contig_parts)}]"


def generate_contig_string(pdb_file, design_mode='binder_denovo', design_length=None):
    """
    Generate RFdiffusion contig string from PDB file.
    
    Args:
        pdb_file: Path to PDB file
        design_mode: Design mode - 'binder_denovo', 'binder_partialdiff', 
                     'monomer_partialdiff', or 'monomer_denovo'
        design_length: Optional length specification (e.g., '60' or '70-80')
        
    Returns:
        str: Contig string formatted for RFdiffusion with square brackets
    """
    # Validate design_length format if provided
    if design_length:
        validate_design_length(design_length)
    
    # For monomer_denovo, we don't need to parse the PDB if design_length is provided
    if design_mode == 'monomer_denovo':
        return generate_contig_monomer_denovo(design_length)
    
    # Parse PDB for other modes
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    # Get protein chains and their residues
    chain_residues = get_protein_chains(structure)
    
    if not chain_residues:
        raise ValueError("No protein chains found in PDB file")
    
    # Generate contig based on design mode
    if design_mode == 'binder_denovo':
        return generate_contig_binder_denovo(chain_residues, design_length)
    elif design_mode == 'monomer_partialdiff':
        return generate_contig_monomer_partialdiff(chain_residues)
    elif design_mode == 'binder_partialdiff':
        return generate_contig_binder_partialdiff(chain_residues)
    else:
        raise ValueError(
            f"Invalid design_mode: {design_mode}. "
            "Must be 'binder_denovo', 'binder_partialdiff', 'monomer_partialdiff', or 'monomer_denovo'"
        )


def main():
    parser = argparse.ArgumentParser(
        description='Generate RFdiffusion contig strings from PDB files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Design Modes:
  binder_denovo:         Full chain notation for all chains (de novo binder design)
                         Example: '[A1-77/0 B23-77/B80-105/0 60]'
  
  monomer_denovo:        Design length only (de novo monomer design)
                         Example: '[60]' or '[70-80]'
                         Requires --design_length parameter
  
  monomer_partialdiff:   Total residue count without chain IDs (partial diffusion of monomer)
                         Example: '[155-155]'
  
  binder_partialdiff:    Chain A as count, chain B+ with full notation (partial diffusion of binder)
                         Example: '[88-88/0 B89-203]'

Design Length:
  For binder_denovo:     Appends the length of the designed binder to the contig
  For monomer_denovo:    Specifies the entire contig (required)
  Format:                Single number (e.g., '60') or range (e.g., '70-80')

Examples:
  %(prog)s input.pdb --design_mode binder_denovo --design_length 60
  %(prog)s input.pdb --design_mode monomer_denovo --design_length 70-80
  %(prog)s input.pdb --design_mode monomer_partialdiff
  %(prog)s input.pdb --design_mode binder_partialdiff -o contigs.txt
        """
    )
    parser.add_argument('pdb_file', help='Input PDB file')
    parser.add_argument(
        '--design_mode',
        choices=['binder_denovo', 'binder_partialdiff', 'monomer_partialdiff', 'monomer_denovo'],
        default='binder_denovo',
        help='Design mode (default: binder_denovo)'
    )
    parser.add_argument(
        '--design_length',
        help='Design length specification (e.g., "60" or "70-80"). Required for monomer_denovo mode.'
    )
    parser.add_argument('-o', '--output', help='Output file (default: print to stdout)')
    
    args = parser.parse_args()
    
    try:
        contig_string = generate_contig_string(
            args.pdb_file, 
            args.design_mode,
            args.design_length
        )
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write(contig_string + '\n')
            print(f"Contig string written to {args.output}")
        else:
            print(contig_string)
            
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
