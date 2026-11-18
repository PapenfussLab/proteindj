#!/usr/bin/env python3
"""Prepare a one-chain PDB target as an RFD4 binder condition spec."""
import argparse
import json
import re
from pathlib import Path

from Bio.PDB import MMCIFIO, PDBParser
from Bio.SeqUtils import seq1


def resolve_residue_spec(spec, original_chain_id, renumber, description):
    """Resolve a ProteinDJ residue specification to target.cif residue IDs."""
    resolved = set()
    for token in filter(None, (item.strip() for item in spec.split(","))):
        match = re.fullmatch(r"([A-Za-z]+)(?:(\d+)(?:-(\d+))?)?", token)
        if not match or match.group(1) != original_chain_id:
            raise ValueError(f"Invalid {description} {token!r} for target chain {original_chain_id}")
        if match.group(2) is None:
            resolved.update(renumber.values())
            continue
        start, end = int(match.group(2)), int(match.group(3) or match.group(2))
        if end < start:
            raise ValueError(f"Invalid {description} range: {token}")
        for number in range(start, end + 1):
            key = (number, " ")
            if key not in renumber:
                raise ValueError(f"{description.capitalize()} {original_chain_id}{number} is absent from the target")
            resolved.add(renumber[key])
    return sorted(resolved)


def backbone_selections(residue_ids):
    """Build simple selectors for backbone atoms in the requested residues."""
    return [{"select": f"sel('A/*/{residue_id}/[N,CA,C,O]')"} for residue_id in residue_ids]


def prepare(input_pdb, length, hotspots, output_dir, flexible_residues=""):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    structure = PDBParser(QUIET=True).get_structure("target", input_pdb)
    chains = [c for c in structure[0] if any(r.id[0] == " " for r in c)]
    if len(chains) != 1:
        raise ValueError("rfd4_denovo currently requires a single protein target chain")
    chain = chains[0]
    original_id = chain.id
    residues = [r for r in chain if r.id[0] == " "]
    if not residues or any(seq1(r.resname, custom_map={"UNK": "X"}) == "X" for r in residues):
        raise ValueError("The target must contain standard amino acids")
    renumber = {(r.id[1], r.id[2]): i for i, r in enumerate(residues, 1)}
    # A continuous label_seq_id avoids dependence on the input PDB's auth numbering.
    for child in list(chain):
        if child.id[0] != " ":
            chain.detach_child(child.id)
    for index, residue in enumerate(residues, 1):
        residue.id = (" ", index, " ")
    chain.id = "A"
    cif = output_dir / "target.cif"
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(cif))

    hotspot_ids = resolve_residue_spec(hotspots, original_id, renumber, "hotspot")
    flexible_ids = resolve_residue_spec(flexible_residues, original_id, renumber, "flexible residue")
    hotspot_selections = [{"select": f"sel('A/*/{residue_id}')"} for residue_id in hotspot_ids]

    spec = {
        "binder": {
            "input": "target.cif",
            "contig": f"{length}, /0, A1-{len(residues)}",
            "conditions": {
                "C_CRDGRP": {1: [{"select": "sel('A/*/*/[N,CA,C,O]')"}]},
                "C_SEQ": {True: [{"select": "sel('A')"}]},
            },
        }
    }
    if flexible_ids:
        # Condition-spec assignments are applied in declaration order. Removing these atoms
        # from coordinate group 1 leaves their structure free while C_SEQ keeps their identity.
        spec["binder"]["conditions"]["C_CRDGRP"][None] = backbone_selections(flexible_ids)
    if hotspot_selections:
        spec["binder"]["conditions"]["C_HOT"] = {True: hotspot_selections}
    (output_dir / "binder.json").write_text(json.dumps(spec, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", required=True)
    parser.add_argument("--length", required=True)
    parser.add_argument("--hotspots", default="")
    parser.add_argument("--flexible-residues", default="")
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    prepare(args.input_pdb, args.length, args.hotspots, args.output_dir, args.flexible_residues)
