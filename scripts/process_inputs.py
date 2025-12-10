#!/usr/bin/env python3
"""
Step-1 pipeline (per your latest rules):

Per PDB:
  1) ALWAYS merge all non-A chains into B in the copied PDBs (so A stays A; any of B/C/D/... -> B).
  2) ALWAYS rewrite PDB in place with _rewrite_pdb_continuous_and_fix_ter (no further chain-ID changes).
  3) Targets & MSA:
     - With the PDB now having only A and (maybe) B, detect chain breaks on B and create VIRTUAL segments
       named B, C, D, ... in order (PDB remains unchanged). If no breaks, the single target is B.
     - Chain A is always A:no_msa by default.
  4) target_id uniqueness:
     - Across the whole run, rows that share the same ordered target sequences (e.g., B or B+C+D)
       get the same target_id: target_1, target_2, ...

HYBRID mode overrides:
  - If the input CSV provides any of: binder_chain, target_id, target_chains, msa_info
    (and they are non-empty), use them to overwrite the inferred values for that row.
    Missing bits are inferred. If msa_info lacks an A directive, we inject "A:no_msa".

Outputs:
  - output_dir/input_pdbs/ (postprocessed PDBs in-place)
  - output_dir/run.csv
  - output_dir/Binder_seq.fasta
  - output_dir/unique_msa/ (FASTA per target chain/segment) and msa_path_* columns in run.csv
"""

import argparse
import csv
import json
import os
import re
import sys
import tempfile
import string
from pathlib import Path
from typing import Dict, List, Tuple

import shutil
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.Chain import Chain

# ---------------------------
# Small helpers
# ---------------------------

def print_info(msg: str):
    print(msg, file=sys.stderr)

def sanitize_stem(stem: str) -> str:
    # Convert to lowercase
    stem = stem.lower()
    # Replace any character not a-z, 0-9, or _ with underscore
    return re.sub(r"[^a-z0-9_]", "_", stem)

def int_str_or_blank(v) -> str:
    if v in ("", None):
        return ""
    try:
        return str(int(float(v)))
    except Exception:
        try:
            return str(int(v))
        except Exception:
            return ""

class KeepAll(Select):
    def accept_atom(self, atom):
        return 1

# ---------------------------
# PDB post-processing
# ---------------------------

def _rewrite_pdb_continuous_and_fix_ter(pdb_in: str, pdb_out: str):
    """
    Rewrite a PDB so that:
      - If any chains other than A/B exist, remap all non-A/B chains to 'B' (A stays 'A').
      - Residues are renumbered *continuously across chains* (order-of-appearance) for ATOM records.
      - Existing TER lines are removed; a single TER is inserted between chains and at the end.
      - Atom serial numbers are renumbered sequentially.

    HETATM residue numbers are kept out of the continuous scheme; ATOM-only gets continuous ids.
    """
    with open(pdb_in, "r") as fh:
        lines = fh.readlines()

    def _is_atom(ln: str) -> bool:
        rec = ln[:6]
        return rec.startswith("ATOM") or rec.startswith("HETATM")

    # Detect if we need to merge non-A/B into B
    seen_chains = set()
    for ln in lines:
        if _is_atom(ln):
            seen_chains.add(ln[21])
    merge_nonab_to_b = any(ch not in ("A", "B") for ch in seen_chains)

    out = []
    atom_serial = 1
    global_resseq = 0
    last_chain_written = None

    # IMPORTANT: key by ORIGINAL chain id to avoid collisions after merge
    # (orig_chain, orig_resseq, icode, resname) → new_resseq
    resid_map = {}

    for ln in lines:
        rec = ln[:6]

        # Drop existing TER; we will reinsert clean TERs ourselves
        if rec.startswith("TER"):
            continue

        if _is_atom(ln):
            orig_chain = ln[21]                 # chain letter in the input line
            write_chain = orig_chain
            if merge_nonab_to_b and orig_chain not in ("A", "B"):
                write_chain = "B"
                ln = f"{ln[:21]}B{ln[22:]}"     # overwrite chain ID in the output line

            # Insert TER when the (possibly merged) chain switches
            if last_chain_written is not None and write_chain != last_chain_written:
                out.append("TER\n")
            last_chain_written = write_chain

            # Parse residue identity (use ORIGINAL CHAIN for the key to avoid collisions)
            try:
                orig_resseq = int(ln[22:26])
            except ValueError:
                orig_resseq = 0
            icode = ln[26]
            resname = ln[17:20]

            if rec.startswith("ATOM"):
                resid_key = (orig_chain, orig_resseq, icode, resname)
                if resid_key not in resid_map:
                    global_resseq += 1
                    resid_map[resid_key] = global_resseq
                new_resseq = resid_map[resid_key]
                ln = f"{ln[:22]}{new_resseq:4d}{ln[26:]}"  # cols 23–26

            # Renumber atom serial (cols 7–11) for both ATOM and HETATM
            ln = f"{ln[:6]}{atom_serial:5d}{ln[11:]}"
            atom_serial += 1

            out.append(ln)
        else:
            out.append(ln)

    # Final TER if file ends with ATOM/HETATM
    if out and (out[-1][:6].startswith("ATOM") or out[-1][:6].startswith("HETATM")):
        out.append("TER\n")

    with open(pdb_out, "w") as fo:
        fo.writelines(out)




def combine_subchains(folder_path: str):
    """
    For each PDB in folder:
      1) If there is more than one chain, remove all TER lines that are NOT for chain A
         (i.e., only keep a TER that immediately follows atoms from chain 'A').
      2) Change every chain with an ID other than 'A' or 'B' to 'B'.
    Saves in-place using Biopython.
    """
    parser = PDBParser(QUIET=True)
    io = PDBIO()

    for filename in os.listdir(folder_path):
        if not filename.lower().endswith(".pdb"):
            continue
        fp = os.path.join(folder_path, filename)

        try:
            # First, see how many chains there are (model 0)
            structure = parser.get_structure(filename, fp)
            model = next(iter(structure))
            chain_ids = [ch.id for ch in model]
            multi_chain = len(chain_ids) > 1

            # If multi-chain, strip TER lines that are not for chain A (text-level pass)
            if multi_chain:
                try:
                    with open(fp, "r") as f:
                        lines = f.readlines()
                    cleaned = []
                    last_chain = None

                    def is_atom(line: str) -> bool:
                        r = line[:6]
                        return r.startswith("ATOM") or r.startswith("HETATM")

                    for ln in lines:
                        rec = ln[:6]
                        if is_atom(ln):
                            last_chain = ln[21]
                            cleaned.append(ln)
                        elif rec.startswith("TER"):
                            # Keep TER only if it closes chain A; otherwise drop it
                            if last_chain == "A":
                                cleaned.append(ln)
                            # else: skip
                        else:
                            cleaned.append(ln)
                    with open(fp, "w") as f:
                        f.writelines(cleaned)
                except Exception as e:
                    print_info(f"[subchains] WARN {filename} (TER-strip): {e}")

                # Re-parse after text edit to ensure we merge based on the updated file
                structure = parser.get_structure(filename, fp)
                model = next(iter(structure))

            # Merge all non-A chains into B
            for model in structure:
                for chain in model:
                    if chain.id not in ("A", "B"):
                        chain.id = "B"

            io.set_structure(structure)
            io.save(fp)

        except Exception as e:
            print_info(f"[subchains] WARN {filename}: {e}")


def postprocess_pdbs_merge_then_renumber(pdb_dir: str):
    """Merge non-A chains into B for all PDBs, then rewrite/renumber/TER-clean in-place."""
    if not os.path.isdir(pdb_dir):
        return
    # 1) Merge (B/C/D/... -> B)
    # 2) Renumber/TER clean
    for fn in sorted(os.listdir(pdb_dir)):
        if not fn.lower().endswith(".pdb"):
            continue
        fp = os.path.join(pdb_dir, fn)
        try:
            _rewrite_pdb_continuous_and_fix_ter(fp, fp)
        except Exception as e:
            print_info(f"[pdb-fix] WARN {fn}: {e}")

# ---------------------------
# Chain-break detection (virtual segments; PDB not modified)
# ---------------------------

def _dist(a, b) -> float:
    dx=a[0]-b[0]; dy=a[1]-b[1]; dz=a[2]-b[2]
    return (dx*dx + dy*dy + dz*dz) ** 0.5

def _break_between(prev_res, cur_res, distance_threshold: float = 4.2) -> bool:
    # residue-number gap
    if (cur_res.id[1] - prev_res.id[1]) > 1:
        return True
    # geometric hop via CA–CA
    if ("CA" in prev_res) and ("CA" in cur_res):
        d = _dist(prev_res["CA"].get_vector().get_array(),
                  cur_res["CA"].get_vector().get_array())
        if d > distance_threshold:
            return True
    return False

def segment_chain_by_breaks(model, chain_id: str, ppb: PPBuilder, distance_threshold: float = 4.2):
    """
    Return list of segments for a chain WITHOUT changing the PDB:
      [(seg_id, seq)], where seg_id is 'B', 'C', 'D' ... depending on segment order.
    Only polymer residues (hetflag == " ") are considered.
    """
    if chain_id not in model:
        return []
    chain = model[chain_id]
    residues = [r for r in chain.get_unpacked_list() if r.id[0] == " "]
    if not residues:
        return []

    # split residues into segments at breaks
    segments = []
    start = 0
    for i in range(1, len(residues)):
        if _break_between(residues[i-1], residues[i], distance_threshold):
            segments.append(residues[start:i])
            start = i
    segments.append(residues[start:])

    # assign chain IDs sequentially starting from the given chain_id (e.g., B -> B,C,D,...)
    alphabet = list(string.ascii_uppercase)
    start_idx = alphabet.index(chain_id)
    seg_out = []
    for i, seg in enumerate(segments):
        seg_id = alphabet[start_idx + i]  # B, C, D, ...
        tmp = Chain(chain_id)  # temporary holder for PPBuilder
        for res in seg:
            tmp.add(res)
        peptides = ppb.build_peptides(tmp)
        seq = "".join(str(p.get_sequence()) for p in peptides) if peptides else ""
        if seq:
            seg_out.append((seg_id, seq))
    return seg_out

# ---------------------------
# Sequence helpers
# ---------------------------

def chain_sequence(model, chain_id: str, ppb: PPBuilder) -> str:
    """Return amino-acid sequence for a chain (may be empty)."""
    if chain_id not in model:
        return ""
    chain = model[chain_id]
    peptides = ppb.build_peptides(chain)
    return "".join(str(peptide.get_sequence()) for peptide in peptides) if peptides else ""

# ---------------------------
# Unique target_id assignment (based on target sequences)
# ---------------------------

def assign_target_ids(rows: List[Dict[str, str]]) -> None:
    """
    Mutates rows in-place to assign target_id = target_1, target_2, ...
    based on the ordered list of target IDs and their sequences.
    Signature key: "B:SEQ|C:SEQ|D:SEQ" using row['segment_ids'] order and
    row['target_subchain_<ID>_seq'] values. Empty sequences are ignored.
    If a row already has a non-empty target_id, it is preserved.
    """
    signature_to_id: Dict[str, str] = {}
    next_idx = 1

    for r in rows:
        if str(r.get("target_id", "")).strip():
            continue

        seg_ids_json = r.get("segment_ids", "[]")
        try:
            seg_ids = json.loads(seg_ids_json)
        except Exception:
            seg_ids = []

        parts = []
        for sid in seg_ids:
            seq = str(r.get(f"target_subchain_{sid}_seq", "") or "")
            if seq:
                parts.append(f"{sid}:{seq}")
        signature = "|".join(parts)

        if not signature:
            r["target_id"] = f"target_{next_idx}"
            next_idx += 1
            continue

        if signature in signature_to_id:
            r["target_id"] = signature_to_id[signature]
        else:
            tid = f"target_{next_idx}"
            signature_to_id[signature] = tid
            r["target_id"] = tid
            next_idx += 1

# ---------------------------
# Builders (pdb_only + hybrid)
# ---------------------------

def build_run_csv_from_pdbs(pdb_folder: Path, out_csv: Path, distance_threshold: float = 4.2):
    """
    Build run.csv directly from PDBs after merging (non-A -> B), with A:no_msa and unique target_id.
    - PDBs now have only A and (maybe) B.
    - If B has breaks, use virtual segments B,C,D,... (PDB unchanged).
    - target_id grouped by unique ordered target sequences across the dataset.
    """
    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    out_rows: List[Dict[str, str]] = []

    for pdb in sorted(pdb_folder.glob("*.pdb")):
        stem = pdb.stem
        try:
            structure = parser.get_structure(stem, str(pdb))
            model = structure[0]

            seq_A = chain_sequence(model, "A", ppb)

            seg_columns = {}
            segment_ids: List[str] = []

            if "B" in model:
                segs = segment_chain_by_breaks(model, "B", ppb, distance_threshold=distance_threshold)
                if not segs:
                    full = chain_sequence(model, "B", ppb)
                    if full:
                        segs = [("B", full)]
                for sid, seq in segs:
                    seg_columns[f"target_subchain_{sid}_seq"] = seq
                    seg_columns[f"target_subchain_{sid}_len"] = int_str_or_blank(len(seq))
                    segment_ids.append(sid)

            msa_list = ["A:no_msa"] + [f"{sid}:run_msa" for sid in segment_ids]

            row = {
                "binder_id": stem,
                "binder_chain": "A",
                "A_seq": seq_A,
                "A_length": int_str_or_blank(len(seq_A)),
                "segment_ids": json.dumps(segment_ids),
                "target_chains": json.dumps(segment_ids),  # downstream expects targets list
                "target_chain_range": json.dumps([f"1:{len(str(seg_columns.get(f'target_subchain_{sid}_seq','')))}" for sid in segment_ids]),
                "msa_info": json.dumps(msa_list),
            }
            row.update(seg_columns)
            out_rows.append(row)

        except Exception as e:
            print_info(f"[pdb_only] FAILED parsing {pdb.name}: {e}")

    # Assign target_id globally by unique target sequences
    assign_target_ids(out_rows)

    # headers
    base = [
        "binder_id", "binder_chain", "target_id",
        "A_seq", "A_length",
        "target_chains", "target_chain_range",
        "segment_ids",
        "msa_info",
    ]
    dyn = set()
    for r in out_rows:
        for k in r.keys():
            if k not in base:
                dyn.add(k)
    fieldnames = base + sorted(dyn)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(out_csv.parent), newline="") as tmpf:
        writer = csv.DictWriter(tmpf, fieldnames=fieldnames)
        writer.writeheader()
        for r in out_rows:
            writer.writerow({fn: r.get(fn, "") for fn in fieldnames})
        Path(tmpf.name).replace(out_csv)
    print(f"[pdb_only] Wrote: {out_csv}")

def update_run_csv_hybrid(input_csv: Path, pdb_folder: Path, output_csv: Path, after_sanitize_map: Dict[str, Path], distance_threshold: float = 4.2):
    """
    HYBRID: Start from legacy CSV rows; infer from PDBs (after merging non-A->B); then, if these columns
    exist and are non-empty, use them to overwrite: binder_chain, target_id, target_chains, msa_info.
    Missing bits are inferred. If msa_info lacks an A directive, we inject "A:no_msa".
    """
    with input_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if not rows:
        raise ValueError("Input CSV is empty.")

    parser = PDBParser(QUIET=True)
    ppb = PPBuilder()
    exact_map = {p.stem: p for p in pdb_folder.glob("*.pdb")}
    out_rows: List[Dict[str, str]] = []

    for base_row in rows:
        binder_id = base_row.get("binder_id", "")
        sanitized_id = sanitize_stem(binder_id)


        binder_chain = (base_row.get("binder_chain", "") or "A").strip().upper()
        user_target_id = str(base_row.get("target_id", "") or "")
        user_target_chains_raw = base_row.get("target_chains", "")
        user_msa_info_raw = base_row.get("msa_info", "")

        # Locate PDB: exact stem or sanitized stem
        pdb_path = exact_map.get(binder_id)
        if pdb_path is None:
            pdb_path = after_sanitize_map.get(sanitized_id)

        row_out = {}
        row_out.update(base_row)  # start with user row
        # ensure binder_id is sanitized
        row_out["binder_id"] = sanitized_id


        if pdb_path is None or not pdb_path.exists():
            # No PDB available: ensure A:no_msa in msa_info if missing, then keep user data
            try:
                mi = json.loads(user_msa_info_raw) if user_msa_info_raw else []
            except Exception:
                mi = []
            if not any(str(x).startswith("A:") for x in mi):
                mi = ["A:no_msa"] + [x for x in mi if not str(x).startswith("A:")]
            row_out["msa_info"] = json.dumps(mi)
            out_rows.append(row_out)
            continue

        structure = parser.get_structure(pdb_path.stem, str(pdb_path))
        model = structure[0]
        seq_A = chain_sequence(model, binder_chain, ppb)
        row_out["binder_chain"] = binder_chain
        row_out["A_seq"] = seq_A
        row_out["A_length"] = int_str_or_blank(len(seq_A))

        # Infer segments from merged PDB (only B exists as non-A)
        seg_columns = {}
        segment_ids: List[str] = []
        if "B" in model:
            segs = segment_chain_by_breaks(model, "B", ppb, distance_threshold=distance_threshold)
            if not segs:
                full = chain_sequence(model, "B", ppb)
                if full:
                    segs = [("B", full)]
            for sid, seq in segs:
                seg_columns[f"target_subchain_{sid}_seq"] = seq
                seg_columns[f"target_subchain_{sid}_len"] = int_str_or_blank(len(seq))
                segment_ids.append(sid)

        inferred_segment_ids = segment_ids[:]
        row_out["segment_ids"] = json.dumps(segment_ids)
        row_out.update(seg_columns)

        # Default msa_info (A:no_msa + inferred segments)
        default_msa_info = ["A:no_msa"] + [f"{sid}:run_msa" for sid in segment_ids]

        # Apply user overrides when present
        used_segment_ids = segment_ids

        VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")  # canonical amino acids
        if user_target_chains_raw:
            try:
                user_target_chains = json.loads(user_target_chains_raw)
                if isinstance(user_target_chains, list) and user_target_chains:
                    used_segment_ids = [str(x) for x in user_target_chains]
                    # recompute seg_columns to align with user-provided IDs
                    seg_columns = {}
                    for sid in used_segment_ids:
                        # user-provided sequence
                        user_seq = str(base_row.get(f"target_subchain_{sid}_seq", "") or "")
                        # PDB-extracted sequence
                        if sid in model:
                            inferred_seq = chain_sequence(model, sid, ppb)
                        elif sid in inferred_segment_ids:
                            inferred_seq = row_out.get(f"target_subchain_{sid}_seq", "")
                        else:
                            inferred_seq = ""
                        
                        # store PDB sequence in _not_used column if user_seq exists
                        if user_seq:
                            seg_columns[f"pdb_extracted_trg_subch_{sid}_not_used"] = inferred_seq

                            # Validate user sequence
                            invalid_chars = [c for c in user_seq if c.upper() not in VALID_AA]
                            if invalid_chars:
                                print(f"WARNING: Invalid characters {invalid_chars} in user sequence for {sid}, binder_id {binder_id}")
                                print("Check your sequences!!")
                                # Optionally, you could remove or replace invalid chars:
                                # user_seq = ''.join(c for c in user_seq if c.upper() in VALID_AA)
                        
                        # overwrite sequence (user sequence preferred)
                        if not user_seq:
                            user_seq = inferred_seq

                        seg_columns[f"target_subchain_{sid}_seq"] = user_seq
                        seg_columns[f"target_subchain_{sid}_len"] = int_str_or_blank(len(user_seq)) if user_seq else ""

                    row_out["segment_ids"] = json.dumps(used_segment_ids)
                    row_out["target_chains"] = json.dumps(used_segment_ids)
                    row_out.update(seg_columns)
            except Exception:
                pass

       
        else:
            row_out["target_chains"] = json.dumps(used_segment_ids)

        if user_msa_info_raw:
            try:
                user_msa = json.loads(user_msa_info_raw)
                if not any(str(x).startswith("A:") for x in user_msa):
                    user_msa = ["A:no_msa"] + [x for x in user_msa if not str(x).startswith("A:")]
                row_out["msa_info"] = json.dumps(user_msa)
            except Exception:
                row_out["msa_info"] = json.dumps(default_msa_info)
        else:
            row_out["msa_info"] = json.dumps(default_msa_info)

        # target_chain_range for used_segment_ids
        ranges = []
        for sid in json.loads(row_out["segment_ids"]):
            sseq = str(row_out.get(f"target_subchain_{sid}_seq", "") or "")
            ranges.append(f"1:{len(sseq)}" if sseq else "")
        row_out["target_chain_range"] = json.dumps(ranges)

        out_rows.append(row_out)

    # Assign target_id globally (preserve user-provided non-empty values)
    assign_target_ids(out_rows)

    # headers: original + new
    fieldnames = list(rows[0].keys())
    for add in ["binder_chain", "A_seq", "A_length", "target_chains", "target_chain_range", "segment_ids", "msa_info", "target_id"]:
        if add not in fieldnames:
            fieldnames.append(add)
    dyn_keys = set()
    for r in out_rows:
        for k in r.keys():
            if k not in fieldnames:
                dyn_keys.add(k)
    fieldnames += sorted(dyn_keys)

    # Write CSV atomically
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=str(output_csv.parent), newline="") as tmpf:
        writer = csv.DictWriter(tmpf, fieldnames=fieldnames)
        writer.writeheader()
        for r in out_rows:
            writer.writerow({fn: r.get(fn, "") for fn in fieldnames})
        Path(tmpf.name).replace(output_csv)

# ---------------------------
# FASTA + MSA path generation
# ---------------------------

def generate_a_seq_fasta(csv_file: str, fasta_path: str, id_col: str = "binder_id", seq_col: str = "A_seq"):
    df = pd.read_csv(csv_file)
    n_written = 0
    with open(fasta_path, "w") as out:
        for idx, row in df.iterrows():
            seq = str(row.get(seq_col, "") or "").strip().upper()
            if not seq:
                continue
            bid = str(row.get(id_col, "") or "").strip()
            if not bid:
                bid = f"row_{idx}"
            bid = re.sub(r"\s+", "_", bid)
            out.write(f">{bid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")
            n_written += 1
    print(f"[A_seq.fasta] wrote {n_written} entries → {fasta_path}")

def process_csv_and_generate_msa(csv_file: str, outpath: str):
    """
    For each row and each msa_info directive:
      - 'no_msa'  -> msa_path_<chain> = 'no_msa'
      - 'run_msa' -> write FASTA once per unique sequence under outpath/, and set msa_path_<chain>
                     to 'unique_msa/msa/<...>.a3m'
    Updates the CSV in-place.
    """
    os.makedirs(outpath, exist_ok=True)
    os.makedirs(os.path.join(outpath, "msa"), exist_ok=True)

    df = pd.read_csv(csv_file)
    msa_lookup = {}

    for idx in df.index:
        binder_id = df.at[idx, 'binder_id']
        target_id = df.at[idx, 'target_id']
        safe_target_id = sanitize_stem(str(target_id))
        binder_chain = str(df.at[idx, 'binder_chain'])

        try:
            msa_info_list = json.loads(df.at[idx, 'msa_info'])
        except Exception as e:
            print(f"Error parsing msa_info for row {idx} ({binder_id}): {e}")
            msa_info_list = []

        for item in msa_info_list:
            try:
                chain, directive = item.split(':')
                chain = chain.strip()
                directive = directive.strip()
            except Exception as e:
                print(f"Error splitting msa_info item '{item}' for row {idx} ({binder_id}): {e}")
                continue

            col_name = f"msa_path_{chain}"
            if directive == "no_msa":
                df.at[idx, col_name] = "no_msa"
                continue

            if directive != "run_msa":
                print(f"Unrecognized directive '{directive}' for chain {chain} in row {idx} ({binder_id}).")
                df.at[idx, col_name] = ""
                continue

            if chain == binder_chain:
                seq = ""
                if f"{binder_chain}_seq" in df.columns and pd.notna(df.at[idx, f"{binder_chain}_seq"]) and df.at[idx, f"{binder_chain}_seq"] != "":
                    seq = df.at[idx, f"{binder_chain}_seq"]
                elif "A_seq" in df.columns and pd.notna(df.at[idx, "A_seq"]) and df.at[idx, "A_seq"] != "":
                    seq = df.at[idx, "A_seq"]
            else:
                seq = ""
                sub_col = f"target_subchain_{chain}_seq"
                if sub_col in df.columns and pd.notna(df.at[idx, sub_col]) and df.at[idx, sub_col] != "":
                    seq = df.at[idx, sub_col]
                else:
                    full_col = f"{chain}_seq"
                    if full_col in df.columns and pd.notna(df.at[idx, full_col]) and df.at[idx, full_col] != "":
                        seq = df.at[idx, full_col]

            if not seq or pd.isna(seq):
                print(f"Warning: No sequence found for chain {chain} in row {idx} ({binder_id}).")
                df.at[idx, col_name] = ""
                continue

            key = str(seq)
            if key in msa_lookup:
                df.at[idx, col_name] = msa_lookup[key]
                continue

            fasta_name = f"{binder_id}_{safe_target_id}_chain_{chain}.fasta"
            fasta_path = os.path.join(outpath, fasta_name)
            a3m_rel = f"msa/{binder_id}_{safe_target_id}_chain_{chain}.a3m"
            a3m_path = os.path.join(outpath, a3m_rel)

            if not os.path.exists(fasta_path):
                with open(fasta_path, "w") as f:
                    f.write(f">{binder_id}_{safe_target_id}_chain_{chain}\n")
                    f.write(str(seq))

            msa_lookup[key] = a3m_path
            df.at[idx, col_name] = a3m_path

    df.to_csv(csv_file, index=False)
    print(f"Processed CSV saved to {csv_file}")

# ---------------------------
# Utilities
# ---------------------------
def copy_and_sanitize_pdbs(input_pdbs: Path, output_pdbs: Path) -> Dict[str, Path]:
    """
    Copy PDBs from input_pdbs → output_pdbs with sanitized names (collision-safe).
    - Overwrites existing files if rerun with same input (with warning).
    - Only adds _1, _2 if *multiple inputs* collide in the same run.
    """
    output_pdbs.mkdir(parents=True, exist_ok=True)
    mapping = {}

    for pdb_file in sorted(input_pdbs.glob("*.pdb")):
        san = sanitize_stem(pdb_file.stem)
        candidate = san

        target_path = output_pdbs / f"{candidate}.pdb"
        if target_path.exists():
            print_info(f"[warning] Overwriting existing {target_path.name} in output folder")
            shutil.copy2(pdb_file, target_path)
            mapping[candidate] = target_path.resolve()
            continue

        # Handle collisions within this run
        k = 0
        while (output_pdbs / f"{candidate}.pdb").exists():
            k += 1
            candidate = f"{san}_{k}"
            target_path = output_pdbs / f"{candidate}.pdb"

        shutil.copy2(pdb_file, target_path)
        if pdb_file.name != target_path.name:
            print_info(f"[copy+sanitize] {pdb_file.name} -> {target_path.name}")
        else:
            print_info(f"[copy] {pdb_file.name}")
        mapping[candidate] = target_path.resolve()

    return mapping



def sanitize_pdb_names_inplace(folder: Path) -> Dict[str, Path]:
    """
    Rename PDB files in-place to sanitized basenames (collision-safe).
    Returns a mapping: sanitized_stem -> full Path.
    """
    pdbs = sorted(folder.glob("*.pdb"))
    mapping = {}

    for p in pdbs:
        old_stem = p.stem
        san = sanitize_stem(old_stem)

        # If already sanitized, just keep as is
        if san == old_stem:
            mapping[old_stem] = p.resolve()
            continue

        candidate = san
        k = 0
        while (folder / f"{candidate}.pdb").exists() and candidate != old_stem:
            k += 1
            candidate = f"{san}_{k}"

        new_path = folder / f"{candidate}.pdb"
        if new_path != p:
            print_info(f"[sanitize] {p.name} -> {new_path.name}")
            p.rename(new_path)
            p = new_path  # update reference

        mapping[candidate] = p.resolve()

    return mapping

# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Copy PDBs → merge non-A→B → rewrite/renumber/TER → run.csv + unique_msa/ with break-aware targets and unique target_id")
    ap.add_argument("--mode", choices=["hybrid", "pdb_only", "seq_only_csv"], default="pdb_only")
    ap.add_argument("--input_pdbs", help="Folder with input PDB files")
    ap.add_argument("--input_csv", help="Legacy CSV (hybrid) or sequences CSV (seq_only_csv)")
    ap.add_argument("--output_dir", required=True, help="Where to write run.csv and unique_msa/")
    ap.add_argument("--distance-threshold", type=float, default=4.2, help="Cα–Cα distance to flag a break (Å).")
    args = ap.parse_args()

    output_dir = Path(args.output_dir).resolve()
    output_pdbs = output_dir / "input_pdbs"
    output_dir.mkdir(parents=True, exist_ok=True)


    # check that no input PDBs collide by sanitation of the names
    if args.mode in ["hybrid", "pdb_only"]:
        if not args.input_pdbs:
            print_info(f"ERROR: --mode {args.mode} requires --input_pdbs")
            sys.exit(1)

        input_pdbs = Path(args.input_pdbs).resolve()
        if not input_pdbs.is_dir():
            print_info(f"ERROR: input_pdbs is not a directory: {input_pdbs}")
            sys.exit(1)

        # --- Check for input collisions before copying ---
        sanitized_counts = {}
        for p in input_pdbs.glob("*.pdb"):
            san = sanitize_stem(p.stem)
            sanitized_counts[san] = sanitized_counts.get(san, 0) + 1

        for san, count in sanitized_counts.items():
            if count > 1:
                print_info(f"[warning - data deletion!] Multiple input PDBs sanitize to the same name: {san}.pdb ({count} files)")

    if args.mode == "hybrid":
        if not args.input_pdbs or not args.input_csv:
            print_info("ERROR: --mode hybrid requires --input_pdbs and --input_csv")
            sys.exit(1)
        input_pdbs = Path(args.input_pdbs).resolve()
        input_csv = Path(args.input_csv).resolve()
        if not input_pdbs.is_dir():
            print_info(f"ERROR: input_pdbs is not a directory: {input_pdbs}")
            sys.exit(1)
        if not input_csv.is_file():
            print_info(f"ERROR: input_csv not found: {input_csv}")
            sys.exit(1)

        
        # Copy + sanitize → merge → renumber
        san_map = copy_and_sanitize_pdbs(input_pdbs, output_pdbs)
        postprocess_pdbs_merge_then_renumber(str(output_pdbs))

        run_csv = output_dir / "run.csv"
        update_run_csv_hybrid(input_csv, output_pdbs, run_csv, san_map, distance_threshold=args.distance_threshold)
        print(f"Wrote: {run_csv}")


    elif args.mode == "pdb_only":
        if not args.input_pdbs:
            print_info("ERROR: --mode pdb_only requires --input_pdbs")
            sys.exit(1)
        input_pdbs = Path(args.input_pdbs).resolve()
        if not input_pdbs.is_dir():
            print_info(f"ERROR: input_pdbs is not a directory: {input_pdbs}")
            sys.exit(1)

        # Copy + sanitize → merge → renumber
        _ = copy_and_sanitize_pdbs(input_pdbs, output_pdbs)
        postprocess_pdbs_merge_then_renumber(str(output_pdbs))

        run_csv = output_dir / "run.csv"
        build_run_csv_from_pdbs(output_pdbs, run_csv, distance_threshold=args.distance_threshold)
        print(f"Wrote: {run_csv}")

    elif args.mode == "seq_only_csv":
        if not args.input_csv:
            print_info("ERROR: --mode seq_only_csv requires --input_csv (with seq_* columns)")
            sys.exit(1)
        input_csv = Path(args.input_csv).resolve()
        if not input_csv.is_file():
            print_info(f"ERROR: input_csv not found: {input_csv}")
            sys.exit(1)

        run_csv = output_dir / "run.csv"
        df = pd.read_csv(input_csv)

        # --- Check required columns ---
        required_cols = {"binder_id", "target_id","A_seq"}
        missing = required_cols - set(df.columns)
        if missing:
            print_info(f"ERROR: Missing required columns in input_csv: {', '.join(missing)}")
            sys.exit(1)

        # --- Sanitize binder_id and target_id ---
        df["binder_id"] = df["binder_id"].astype(str).map(sanitize_stem)
        df["target_id"] = df["target_id"].astype(str).map(sanitize_stem)

        df.to_csv(run_csv, index=False)
        print(f"[seq_only_csv] Wrote: {run_csv}")


    # Common tail: binder FASTA + unique MSA and set msa_path_* columns
    run_csv = output_dir / "run.csv"
    a_seq_fasta = output_dir / "Binder_seq.fasta"
    generate_a_seq_fasta(str(run_csv), str(a_seq_fasta))

    unique_msa_dir = output_dir / "unique_msa"
    unique_msa_dir.mkdir(parents=True, exist_ok=True)
    (unique_msa_dir / "msa").mkdir(parents=True, exist_ok=True)

    process_csv_and_generate_msa(str(run_csv), str(unique_msa_dir))
    print(f"MSA files under: {unique_msa_dir}")

if __name__ == "__main__":
    main()

