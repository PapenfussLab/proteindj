#!/usr/bin/env python3
"""Convert RFD4-generated complexes to ProteinDJ fold records."""
import argparse
import json
from pathlib import Path

from Bio.PDB import Chain, MMCIFParser, Model, PDBIO, Structure
from Bio.SeqUtils import seq1


def sequence(chain):
    return "".join(seq1(r.resname, custom_map={"UNK": "X"}) for r in chain if r.id[0] == " ")


def convert(samples_dir, target_cif, start_id, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    parser = MMCIFParser(QUIET=True)
    target = next(parser.get_structure("target", str(target_cif)).get_chains())
    target_sequence = sequence(target)
    samples = sorted(Path(samples_dir).rglob("*_generated.cif"))
    if not samples:
        raise ValueError("RFD4 did not produce any *_generated.cif files")
    metadata = []
    for index, sample in enumerate(samples):
        model = next(parser.get_structure("design", str(sample)).get_models())
        chains = [c for c in model if sequence(c)]
        matching = [c for c in chains if sequence(c) == target_sequence]
        if len(chains) != 2 or len(matching) != 1:
            raise ValueError(f"Expected one binder and one unchanged target in {sample}")
        target_chain = matching[0]
        binder = next(c for c in chains if c is not target_chain)
        if "X" in sequence(binder):
            raise ValueError(f"Non-standard binder sequence in {sample}")
        result = Structure.Structure("design")
        result_model = Model.Model(0)
        result.add(result_model)
        next_number = 1
        for name, source in (("A", binder), ("B", target_chain)):
            new_chain = Chain.Chain(name)
            for residue in source:
                if residue.id[0] != " ":
                    continue
                item = residue.copy()
                item.id = (" ", next_number, " ")
                new_chain.add(item)
                next_number += 1
            result_model.add(new_chain)
        fold_id = start_id + index
        pdb = output_dir / f"fold_{fold_id}.pdb"
        writer = PDBIO()
        writer.set_structure(result)
        writer.save(str(pdb))
        record = {"fold_id": fold_id, "rfd4_source": sample.name, "sequence": sequence(binder)}
        (output_dir / f"fold_{fold_id}.json").write_text(json.dumps(record))
        metadata.append(record)
    with (output_dir / "rfd4_metadata.jsonl").open("w") as handle:
        for record in metadata:
            handle.write(json.dumps(record) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-dir", required=True)
    parser.add_argument("--target-cif", required=True)
    parser.add_argument("--start-id", type=int, required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    convert(args.samples_dir, args.target_cif, args.start_id, args.output_dir)
