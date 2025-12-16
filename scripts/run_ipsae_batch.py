#!/usr/bin/env python3
"""
MIT License

Copyright (c) 2025 Digital Biotechnology Lab

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

run_ipsae_batch.py — compute ipSAE/pDockQ/LIS/ipae metrics from structure files

- Auto-discovers `binder_id` from PDB/NPZ files in `--input-dir`
- Processes **in parallel** using a process pool (configurable with `--max-workers`)
- Generates/uses `*_paeXX_distYY.txt` via **ipsae_w_ipae.py**
- Outputs metrics to JSONL format with fold_id and seq_id extracted from filenames

Example:
  python run_ipsae_batch.py \
    --input-dir /path/to/boltz_predictions/ \
    --out-jsonl ipsae_metrics.jsonl \
    --ipsae-script-path ./ipsae_w_ipae.py \
    --pae-cutoff 10 --dist-cutoff 10

"""
from __future__ import annotations

import os
import argparse
import pandas as pd
import glob
import subprocess
import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# -------------------------
# CLI
# -------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute ipSAE metrics from structure files and output to JSONL"
    )
    p.add_argument("--input-dir", required=True, help="Directory containing PDB and NPZ files")
    p.add_argument("--out-jsonl", required=True, help="Path to write the calculated metrics JSONL")
    p.add_argument("--ipsae-script-path", default="ipsae_w_ipae.py", help="Path to ipsae_w_ipae.py")
    p.add_argument("--pae-cutoff", type=float, default=10.0, help="PAE cutoff")
    p.add_argument("--dist-cutoff", type=float, default=10.0, help="Distance cutoff")
    p.add_argument("--overwrite-ipsae", action="store_true", help="Recompute even if *.txt already exists")
    p.add_argument("--max-workers", type=int, default=None, help="Parallel workers (defaults to CPU count)")
    p.add_argument("--verbose", action="store_true", help="Print more info while processing")

    return p.parse_args()

# -------------------------
# File indexing (only for provided dirs)
# -------------------------

def build_file_index(input_dir: str, verbose: bool = False) -> Dict[str, Dict[str, str]]:
    """
    Build file index from a flat directory structure.
    
    Expects files like:
    - pae_fold_1_seq_19_boltzpred.npz  -> binder_id: fold_1_seq_19_boltzpred
    - fold_1_seq_19_boltzpred.pdb      -> binder_id: fold_1_seq_19_boltzpred
    
    The binder_id is extracted as the base filename without pae_ prefix and without extension.
    
    Returns: Dict[binder_id, {'structure': path, 'confidence': path}]
    """
    if not os.path.isdir(input_dir):
        raise ValueError(f"Input directory does not exist or is not a directory: {input_dir}")
    
    index: Dict[str, Dict[str, str]] = {}
    
    # Get all files directly from the input_dir (no subdirectories)
    all_files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]
    
    for filename in all_files:
        full_path = os.path.join(input_dir, filename)
        
        # Match PDB files
        if filename.endswith('.pdb'):
            # Extract binder_id: just remove .pdb extension
            bid = filename.replace('.pdb', '')
            
            index.setdefault(bid, {})['structure'] = full_path
            if verbose:
                print(f"DEBUG: Added structure for '{bid}': {full_path}")
        
        # Match PAE NPZ files (with pae_ prefix)
        elif filename.startswith('pae_') and filename.endswith('.npz'):
            # Extract binder_id: remove 'pae_' prefix and '.npz' extension
            bid = filename.replace('pae_', '').replace('.npz', '')
            
            index.setdefault(bid, {})['confidence'] = full_path
            if verbose:
                print(f"DEBUG: Added confidence for '{bid}': {full_path}")

    if verbose:
        print(f"\nDEBUG: File index summary:")
        for bid, files in index.items():
            has_struct = 'structure' in files
            has_conf = 'confidence' in files
            status = "✓ complete" if (has_struct and has_conf) else "⚠ incomplete"
            print(f"  {status} | {bid}: struct={has_struct}, conf={has_conf}")
        print()
    
    return index

# -------------------------
# Locate files for a binder
# -------------------------

def locate_files(bid: str, index: Dict[str, Dict[str, str]]) -> Tuple[Optional[Tuple[str, str]], Optional[str]]:
    """
    Locate structure and confidence files for a binder_id.
    
    Returns:
        (structure_path, confidence_path) if both found, else (None, error_message)
    """
    # try direct + case variants
    variants = [bid, bid.lower(), bid.upper()]

    for v in variants:
        files = index.get(v)
        if files and 'structure' in files and 'confidence' in files:
            return (files['structure'], files['confidence']), None

    return None, f"[{bid}] missing structure or confidence files"

# -------------------------
# IPSAE invocations & parsing
# -------------------------

def expected_txt_path(struct_path: str, pae_cutoff: float, dist_cutoff: float) -> Path:
    stem = Path(struct_path).with_suffix('').name
    return Path(struct_path).with_name(f"{stem}_pae{int(pae_cutoff):02d}_dist{int(dist_cutoff):02d}.txt")


def calculate_ipsae(conf: str, struct: str, pae_cutoff: float, dist_cutoff: float, script_path: str, overwrite: bool, verbose: bool):
    out_txt = expected_txt_path(struct, pae_cutoff, dist_cutoff)
    if out_txt.exists() and not overwrite:
        if verbose:
            print(f"  - using existing {out_txt.name}")
        return
    cmd = [sys.executable, script_path, conf, struct, str(pae_cutoff), str(dist_cutoff)]
    if verbose:
        print("  - RUN:", " ".join(map(str, cmd)))
    subprocess.run(cmd, check=True)


def get_ipsae_min_max(path: str, target_chain: str = 'A'):
    """
    Parse an ipsae_w_ipae.py TXT and aggregate metrics for the target_chain.

    Rules:
    - ipSAE_max: use Type==max rows. For each partner chain, take the maximum
      ipSAE value (after excluding zero-distance pairs), then average across partners.
    - ipSAE_min: use Type==asym rows. For each partner chain, take the lower of the
      asymmetric ipSAE values (after excluding zero-distance pairs), then average across partners.
    - Exclude zero-distance pairs (dist1==0 or dist2==0) for ipSAE-family, pDockQ, LIS.
    - Compute ipae from Type==max rows WITHOUT distance gating so it remains cutoff-independent.
    """
    if target_chain != 'A':
        raise ValueError("This function only supports target_chain='A'")
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if not lines:
        return (None,) * 8
    header = lines[0].split()
    required_cols = (
        'Chn1', 'Chn2', 'Type', 'ipSAE', 'ipSAE_avg', 'ipSAE_min_in_calculation',
        'LIS', 'ipSAE_d0chn', 'ipSAE_d0dom', 'dist1', 'dist2', 'ipae'
    )
    if not all(col in header for col in required_cols):
        missing = [col for col in required_cols if col not in header]
        raise ValueError(f"Missing columns in header: {missing} for file {path}")
    idx = {col: header.index(col) for col in required_cols}

    # Accumulate ipSAE-family metrics
    # - data_max: from Type==max rows (non-zero distance)
    # - data_asym: from Type==asym rows (non-zero distance)
    data_max: Dict[str, Dict[str, List[float]]] = {}
    data_asym: Dict[str, List[float]] = {}
    # Accumulate ipae from Type==max rows, regardless of dist1/dist2
    ipae_by_partner: Dict[str, List[float]] = {}

    for ln in lines[1:]:
        parts = ln.split()
        c1, c2 = parts[idx['Chn1']], parts[idx['Chn2']]
        if target_chain not in (c1, c2):
            continue
        t = parts[idx['Type']].lower()
        partner = c2 if c1 == target_chain else c1

        # ipae: collect from Type==max rows regardless of distance counts
        if t == 'max':
            try:
                ipae_val = float(parts[idx['ipae']])
                ipae_by_partner.setdefault(partner, []).append(ipae_val)
            except ValueError:
                pass

        # ipSAE-family (and LIS): exclude zero-distance pairs
        try:
            d1 = float(parts[idx['dist1']]); d2 = float(parts[idx['dist2']])
        except ValueError:
            continue
        if d1 == 0.0 or d2 == 0.0:
            continue

        if t == 'max':
            bucket = data_max.setdefault(partner, {
                'ipSAE': [], 'ipSAE_avg': [], 'LIS': [], 'ipSAE_min_in_calculation': [],
                'ipSAE_d0chn': [], 'ipSAE_d0dom': []
            })
            try:
                bucket['ipSAE'].append(float(parts[idx['ipSAE']]))
                bucket['ipSAE_avg'].append(float(parts[idx['ipSAE_avg']]))
                bucket['LIS'].append(float(parts[idx['LIS']]))
                bucket['ipSAE_min_in_calculation'].append(float(parts[idx['ipSAE_min_in_calculation']]))
                bucket['ipSAE_d0chn'].append(float(parts[idx['ipSAE_d0chn']]))
                bucket['ipSAE_d0dom'].append(float(parts[idx['ipSAE_d0dom']]))
            except ValueError:
                continue
        elif t == 'asym':
            try:
                val = float(parts[idx['ipSAE']])
            except ValueError:
                continue
            data_asym.setdefault(partner, []).append(val)

    # Aggregate ipSAE-family metrics
    # ipSAE_min from Type==asym (lower/asymmetric per partner), averaged across partners
    if not data_asym:
        avg_min = 0.0
    else:
        per_partner_min = [min(vs) if vs else 0.0 for vs in data_asym.values()]
        avg_min = sum(per_partner_min) / len(per_partner_min)

    # ipSAE_max and other metrics from Type==max rows
    if not data_max:
        avg_max = avg_ipsae_avg = avg_lis = avg_ipsae_min = avg_d0chn = avg_d0dom = 0.0
    else:
        maxs, means_avg, means_lis, means_min = [], [], [], []
        d0chn_list, d0dom_list = [], []
        for vals in data_max.values():
            maxs.append(max(vals['ipSAE']) if vals['ipSAE'] else 0.0)
            means_avg.append(min(vals['ipSAE_avg']) if vals['ipSAE_avg'] else 0.0)
            means_lis.append(min(vals['LIS']) if vals['LIS'] else 0.0)
            means_min.append(min(vals['ipSAE_min_in_calculation']) if vals['ipSAE_min_in_calculation'] else 0.0)
            d0chn_list.append(min(vals['ipSAE_d0chn']) if vals['ipSAE_d0chn'] else 0.0)
            d0dom_list.append(min(vals['ipSAE_d0dom']) if vals['ipSAE_d0dom'] else 0.0)
        avg_max = sum(maxs) / len(maxs)
        avg_ipsae_avg = sum(means_avg) / len(means_avg)
        avg_lis = sum(means_lis) / len(means_lis)
        avg_ipsae_min = sum(means_min) / len(means_min)
        avg_d0chn = sum(d0chn_list) / len(d0chn_list)
        avg_d0dom = sum(d0dom_list) / len(d0dom_list)

    # Aggregate ipae separately, independent of distance gating
    if not ipae_by_partner:
        avg_ipae = 0.0
    else:
        ipae_vals = [min(vs) if vs else 0.0 for vs in ipae_by_partner.values()]
        avg_ipae = sum(ipae_vals) / len(ipae_vals)

    return avg_min, avg_max, avg_ipsae_avg, avg_lis, avg_ipsae_min, avg_d0chn, avg_d0dom, avg_ipae


def get_pDockQ_min_max(path, target_chain='A'):
    with open(path) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    if not lines:
        return {"pDockQ": [None, None], "pDockQ2": [None, None]}
    header = lines[0].split()
    idx = {col: header.index(col) for col in ('Chn1','Chn2','pDockQ','pDockQ2','dist1','dist2')}
    partner_data = {}
    for ln in lines[1:]:
        parts = ln.split()
        c1, c2 = parts[idx['Chn1']], parts[idx['Chn2']]
        if target_chain not in (c1, c2):
            continue
        try:
            d1 = float(parts[idx['dist1']]); d2 = float(parts[idx['dist2']])
        except ValueError:
            continue
        if d1 == 0 or d2 == 0:
            continue
        partner = c2 if c1 == target_chain else c1
        try:
            pd1 = float(parts[idx['pDockQ']]); pd2 = float(parts[idx['pDockQ2']])
        except ValueError:
            continue
        bucket = partner_data.setdefault(partner, {'pDockQ': [], 'pDockQ2': []})
        bucket['pDockQ'].append(pd1)
        bucket['pDockQ2'].append(pd2)
    if not partner_data:
        return {"pDockQ": [0, 0], "pDockQ2": [0, 0]}
    mins_p1, maxs_p1, mins_p2, maxs_p2 = [], [], [], []
    for vals in partner_data.values():
        mins_p1.append(min(vals['pDockQ']) if vals['pDockQ'] else 0.0)
        maxs_p1.append(max(vals['pDockQ']) if vals['pDockQ'] else 0.0)
        mins_p2.append(min(vals['pDockQ2']) if vals['pDockQ2'] else 0.0)
        maxs_p2.append(max(vals['pDockQ2']) if vals['pDockQ2'] else 0.0)
    avg_min_p1 = sum(mins_p1) / len(mins_p1)
    avg_max_p1 = sum(maxs_p1) / len(maxs_p1)
    avg_min_p2 = sum(mins_p2) / len(mins_p2)
    avg_max_p2 = sum(maxs_p2) / len(maxs_p2)
    return {"pDockQ":  [avg_min_p1, avg_max_p1], "pDockQ2": [avg_min_p2, avg_max_p2]}


def min_max_pae_for_chain_contacts(json_path, threshold, target_chain='A'):
    with open(json_path) as f:
        data = json.load(f)
    chain_ids = data['token_chain_ids']
    cp = np.array(data['contact_probs'])
    pae = np.array(data['pae'])
    target = np.array([c == target_chain for c in chain_ids])
    non_target = ~target
    mask = np.outer(target, non_target) & (cp > threshold)
    if not mask.any():
        forward_probs = cp[np.outer(target, non_target)]
        return 25, 25, int(forward_probs.size > 0) and 0
    vals = pae[mask]
    return float(vals.min()), float(vals.max()), int(vals.size)

def find_ipsae_txts(struct_path, bid):
    folder = os.path.dirname(struct_path)
    txts = glob.glob(os.path.join(folder, f"{bid}*.txt"))
    if not txts:
        txts = glob.glob(os.path.join(folder, f"{bid.lower()}*.txt"))
    return [f for f in txts if not (f.endswith('byres.txt') or f.endswith('done.txt'))]

# -------------------------
# Processing one binder (sequential)
# -------------------------

def process_binder(bid: str, index: Dict[str, Dict[str, str]], pae_cutoff: float, dist_cutoff: float, ipsae_script: str, overwrite: bool, verbose: bool):
    results: Dict[str, float] = {}
    notes: List[str] = []
    
    files, error = locate_files(bid, index)
    if error:
        notes.append(error)
        return results, notes
    
    struct, conf = files
    
    if verbose:
        print(f"  Processing {bid}")
        print(f"    struct={struct}")
        print(f"    conf  ={conf}")
    
    try:
        calculate_ipsae(conf, struct, pae_cutoff, dist_cutoff, ipsae_script, overwrite, verbose)
    except subprocess.CalledProcessError as e:
        notes.append(f"[{bid}] IPSAE failed: {e}")
        return results, notes
    
    txts = find_ipsae_txts(struct, bid)
    if not txts:
        notes.append(f"[{bid}] No .txt found for {struct}")
        return results, notes
    
    txt = txts[0]
    mn, mx, avg_ipsae_avg, avg_LIS, avg_min_ipsae, avg_ipSAE_d0chn, avg_ipSAE_d0dom, avg_ipae = get_ipsae_min_max(txt)
    pqq = get_pDockQ_min_max(txt)
    pdockQ_mn, pdockQ_mx = pqq["pDockQ"][0], pqq["pDockQ"][1]
    pdockQ2_mn, pdockQ2_mx = pqq["pDockQ2"][0], pqq["pDockQ2"][1]
    
    results["pDockQ_min"] = pdockQ_mn
    results["pDockQ_max"] = pdockQ_mx
    results["pDockQ2_min"] = pdockQ2_mn
    results["pDockQ2_max"] = pdockQ2_mx
    results["ipSAE_min"] = mn
    results["ipSAE_max"] = mx
    results["ipSAE_avg"] = avg_ipsae_avg
    results["LIS"] = avg_LIS
    results["ipSAE_min_in_calculation"] = avg_min_ipsae
    results["ipSAE_d0chn"] = avg_ipSAE_d0chn
    results["ipSAE_d0dom"] = avg_ipSAE_d0dom
    results["ipae"] = avg_ipae
    
    return results, notes


def write_jsonl(data_dict: Dict[str, Dict[str, float]], output_path: str):
    """
    Write metrics to JSONL format.
    Each line is a JSON object with binder_id and all metrics.
    Extracts fold_id and seq_id from binder_id (e.g., fold_9_seq_19_boltzpred -> fold_id=9, seq_id=19)
    """
    import json
    import re
    
    records_written = 0
    with open(output_path, 'w') as f:
        for binder_id, metrics in data_dict.items():
            record = {'binder_id': binder_id}
            
            # Extract fold_id and seq_id from binder_id
            # Pattern: fold_X_seq_Y_boltzpred (or any suffix)
            match = re.match(r'fold_(\d+)_seq_(\d+)', binder_id)
            if match:
                record['fold_id'] = int(match.group(1))
                record['seq_id'] = int(match.group(2))
            else:
                # If pattern doesn't match, set to None or log warning
                record['fold_id'] = None
                record['seq_id'] = None
            
            record.update(metrics)
            f.write(json.dumps(record) + '\n')
            records_written += 1
    
    print(f"\n{'='*60}")
    print(f"✓ JSONL output written successfully")
    print(f"{'='*60}")
    print(f"  File: {output_path}")
    print(f"  Records: {records_written}")
    print(f"  Format: One JSON object per line")
    print(f"  Fields: binder_id, fold_id, seq_id, [metrics]")
    print(f"{'='*60}\n")

# -------------------------
# Main
# -------------------------

def main():
    args = parse_args()

    # Build file index from input directory
    index = build_file_index(args.input_dir, verbose=args.verbose)
    
    # Extract binder_ids that have both structure and confidence files
    binder_ids = [bid for bid, files in index.items() 
                  if 'structure' in files and 'confidence' in files]
    
    if not binder_ids:
        raise SystemExit(f"No valid binder_ids found with both structure and confidence files in {args.input_dir}")
    
    print(f"Found {len(binder_ids)} binder_ids to process")

    # Process in parallel and collect metrics per binder
    from concurrent.futures import ProcessPoolExecutor, as_completed
    collected: Dict[str, Dict[str, float]] = {}
    all_notes: List[str] = []

    tasks = []
    for bid in binder_ids:
        tasks.append((bid, index, args.pae_cutoff, args.dist_cutoff, args.ipsae_script_path, args.overwrite_ipsae, args.verbose))

    with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {pool.submit(process_binder, *t): t[0] for t in tasks}
        for fut in as_completed(futures):
            bid = futures[fut]
            try:
                res, notes = fut.result()
            except Exception as e:
                all_notes.append(f"[{bid}] worker failed: {e}")
                res = {}
                notes = []  # <-- ensure defined
            collected[bid] = res
            all_notes.extend(notes)

    # Write JSONL output
    write_jsonl(collected, args.out_jsonl)

    if all_notes:
        print(f"\n{'='*60}")
        print(f"⚠ Notes / Warnings ({len(all_notes)})")
        print(f"{'='*60}")
        for n in all_notes:
            print(f"  - {n}")
        print(f"{'='*60}\n")

if __name__ == '__main__':
    main()
