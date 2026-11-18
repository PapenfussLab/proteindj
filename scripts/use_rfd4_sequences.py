#!/usr/bin/env python3
"""Pass RFD4's designed sequence to prediction and compute sequence filters."""
import argparse
import json
import shutil
from pathlib import Path

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def run(input_dir, output_dir, min_ext, max_ext, min_pi, max_pi):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for pdb in sorted(Path(input_dir).glob("fold_*.pdb")):
        fold_id = int(pdb.stem.split("_")[1])
        chain = PDBParser(QUIET=True).get_structure("design", str(pdb))[0]["A"]
        sequence = "".join(seq1(r.resname) for r in chain if r.id[0] == " ")
        analysis = ProteinAnalysis(sequence)
        ext = analysis.molar_extinction_coefficient()[0]
        pi = round(analysis.isoelectric_point(), 2)
        name = f"fold_{fold_id}_seq_0"
        record = {"description": name, "fold_id": fold_id, "seq_id": 0,
                  "seq_length": len(sequence), "seq_ext_coef": ext,
                  "seq_MW": int(analysis.molecular_weight()), "seq_pI": pi}
        records.append(record)
        if ((min_ext is None or ext >= min_ext) and (max_ext is None or ext <= max_ext)
                and (min_pi is None or pi >= min_pi) and (max_pi is None or pi <= max_pi)):
            shutil.copyfile(pdb, output_dir / f"{name}.pdb")
    with (output_dir / "rfd4_seq_metadata.jsonl").open("w") as handle:
        for record in records:
            handle.write(json.dumps(record) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    for name in ("min-ext", "max-ext", "min-pi", "max-pi"):
        parser.add_argument("--" + name, type=float)
    args = parser.parse_args()
    run(args.input_dir, args.output_dir, args.min_ext, args.max_ext, args.min_pi, args.max_pi)
