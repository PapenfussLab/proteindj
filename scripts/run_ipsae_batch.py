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

run_ipsae_batch.py — read a run CSV and append ipSAE/pDockQ/LIS/ipae metrics

- Reads **--run-csv** (must contain a `binder_id` column)
- Flexibly supports any subset of inputs: **AF3**, **Boltz1**, **ColabFold**
  - Provide any combination of: `--af3-dir`, `--boltz1-dir`, `--colab-dir`
  - Or explicitly set `--sources` among: af3 boltz colab (inferred from dirs if omitted)
- Processes **in parallel** using a process pool (configurable with `--max-workers`)
- Generates/uses `*_paeXX_distYY.txt` via your **ipsae_w_ipae.py**
- Appends **prefixed columns** per source and **overwrites the same run CSV**
  - Prefixes: `af3_*`, `boltz1_*` (note the 1), `colab_*`

Example:
  python run_ipsae_batch.py \
    --run-csv run.csv \
    --boltz-dir /path/to/boltz_results_fasta_folder \
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
        description="Read a run CSV (binder_id column), compute ipSAE metrics, and append columns in-place"
    )
    p.add_argument("--run-csv", required=True, help="CSV with a 'binder_id' column; will be overwritten with appended metrics")

    # Any subset is allowed; if none given, error
    p.add_argument("--boltz-dir", help="Path to BOLTZ outputs directory (reads PDB and NPZ files directly from this folder)")

    # Optional explicit source list; otherwise inferred from dirs provided
    p.add_argument("--sources", nargs="+", choices=["boltz"], help="Which sources to scan (defaults to dirs provided)")

    p.add_argument("--ipsae-script-path", default="ipsae_w_ipae.py", help="Path to ipsae_w_ipae.py")
    p.add_argument("--pae-cutoff", type=float, default=10.0, help="PAE cutoff")
    p.add_argument("--dist-cutoff", type=float, default=10.0, help="Distance cutoff")
    p.add_argument("--overwrite-ipsae", action="store_true", help="Recompute even if *.txt already exists")
    p.add_argument("--max-workers", type=int, default=None, help="Parallel workers (defaults to CPU count)")
    p.add_argument("--backup", action="store_true", help="Write run.csv.bak before overwrite")
    p.add_argument("--verbose", action="store_true", help="Print more info while processing")
    p.add_argument("--out-csv", required=True, help="Path to write the calculated metrics CSV")
    p.add_argument("--out-jsonl", help="Path to write the calculated metrics JSONL")

    return p.parse_args()

# -------------------------
# File indexing (only for provided dirs)
# -------------------------

def build_file_index(boltz_dir: Optional[str], verbose: bool = False):
    """
    Build file index from a flat directory structure.
    
    Expects files like:
    - pae_fold_1_seq_19_boltzpred.npz  -> binder_id: fold_1_seq_19_boltzpred
    - fold_1_seq_19_boltzpred.pdb      -> binder_id: fold_1_seq_19_boltzpred
    
    The binder_id is extracted as the base filename without pae_ prefix and without extension.
    """
    index: Dict[str, Dict[str, Dict[str, object]]] = {"boltz": {}, "af3": {}, "colab": {}}
    
    if boltz_dir:
        if not os.path.isdir(boltz_dir):
            print(f"WARNING: boltz_dir does not exist or is not a directory: {boltz_dir}")
            return index
            
        # Get all files directly from the boltz_dir (no subdirectories)
        all_files = [f for f in os.listdir(boltz_dir) if os.path.isfile(os.path.join(boltz_dir, f))]
        
        for filename in all_files:
            full_path = os.path.join(boltz_dir, filename)
            
            # Match PDB files
            if filename.endswith('.pdb'):
                # Extract binder_id: just remove .pdb extension
                bid = filename.replace('.pdb', '')
                
                index['boltz'].setdefault(bid, {})['structure'] = full_path
                if verbose:
                    print(f"DEBUG: Added structure for '{bid}': {full_path}")
            
            # Match PAE NPZ files (with pae_ prefix)
            elif filename.startswith('pae_') and filename.endswith('.npz'):
                # Extract binder_id: remove 'pae_' prefix and '.npz' extension
                bid = filename.replace('pae_', '').replace('.npz', '')
                
                index['boltz'].setdefault(bid, {})['confidence'] = full_path
                if verbose:
                    print(f"DEBUG: Added confidence for '{bid}': {full_path}")

    if verbose:
        print(f"\nDEBUG: Final index['boltz'] summary:")
        for bid, files in index['boltz'].items():
            has_struct = 'structure' in files
            has_conf = 'confidence' in files
            status = "✓ complete" if (has_struct and has_conf) else "⚠ incomplete"
            print(f"  {status} | {bid}: struct={has_struct}, conf={has_conf}")
        print()
    
    return index

# -------------------------
# Locate files for a binder
# -------------------------

def locate_files(bid: str, index) -> Tuple[List[Tuple[str,str,str,str]], List[str]]:
    valid: List[Tuple[str,str,str,str]] = []
    missing: List[str] = []

    # try direct + case variants
    variants = [bid, bid.lower(), bid.upper()]

    # Boltz1
    b_found = False
    for v in variants:
     
        b = index['boltz'].get(v)
        if b and 'structure' in b and 'confidence' in b:
            valid.append((bid, 'boltz', b['structure'], b['confidence']))
            b_found = True
            break

    if index['boltz'] and not b_found:
      
        missing.append(f"[{bid}] boltz missing structure or confidence files")

    # AF3
    a_found = False
    for v in variants:
        a = index['af3'].get(v)
        if a and 'structure' in a and 'confidence' in a:
            valid.append((bid, 'af3', a['structure'], a['confidence']))
            a_found = True
            break
    if index['af3'] and not a_found:
        missing.append(f"[{bid}] af3 missing structure or confidence files")

    # Colab
    c_found = False
    for v in variants:
        c = index['colab'].get(v, {})
        jsons = c.get('jsons', [])
        pdbs = c.get('pdbs', [])
        if jsons and pdbs:
            valid.append((bid, 'colab', pdbs[0], jsons[0]))
            c_found = True
            break
    if index['colab'] and not c_found:
        missing.append(f"[{bid}] colab missing JSON or PDB files")


    return valid, missing

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

def process_binder(bid: str, index, pae_cutoff: float, dist_cutoff: float, ipsae_script: str, overwrite: bool, verbose: bool):
    results: Dict[str, float] = {}
    notes: List[str] = []
    valid, miss = locate_files(bid, index)
    notes.extend(miss)
    for _, src, struct, conf in valid:
        if verbose:
            print(f"  [{src}] struct={struct}")
            print(f"  [{src}] conf  ={conf}")
        try:
            calculate_ipsae(conf, struct, pae_cutoff, dist_cutoff, ipsae_script, overwrite, verbose)
        except subprocess.CalledProcessError as e:
            notes.append(f"[{bid}-{src}] IPSAE failed: {e}")
            continue
        txts = find_ipsae_txts(struct, bid)
        if not txts:
            notes.append(f"[{bid}-{src}] No .txt found for {struct}")
            continue
        txt = txts[0]
        mn, mx, avg_ipsae_avg, avg_LIS, avg_min_ipsae, avg_ipSAE_d0chn, avg_ipSAE_d0dom, avg_ipae = get_ipsae_min_max(txt)
        pqq = get_pDockQ_min_max(txt)
        pdockQ_mn, pdockQ_mx = pqq["pDockQ"][0], pqq["pDockQ"][1]
        pdockQ2_mn, pdockQ2_mx = pqq["pDockQ2"][0], pqq["pDockQ2"][1]
        prefix = 'boltz' if src == 'boltz' else src
        results[f"{prefix}_pDockQ_min"] = pdockQ_mn
        results[f"{prefix}_pDockQ_max"] = pdockQ_mx
        results[f"{prefix}_pDockQ2_min"] = pdockQ2_mn
        results[f"{prefix}_pDockQ2_max"] = pdockQ2_mx
        results[f"{prefix}_ipSAE_min"] = mn
        results[f"{prefix}_ipSAE_max"] = mx
        results[f"{prefix}_ipSAE_avg"] = avg_ipsae_avg
        results[f"{prefix}_LIS"] = avg_LIS
        results[f"{prefix}_ipSAE_min_in_calculation"] = avg_min_ipsae
        results[f"{prefix}_ipSAE_d0chn"] = avg_ipSAE_d0chn
        results[f"{prefix}_ipSAE_d0dom"] = avg_ipSAE_d0dom
        results[f"{prefix}_ipae"] = avg_ipae
        if src == 'af3':
            rmin, rmax, count = min_max_pae_for_chain_contacts(conf, 0.60)
            results[f"{prefix}_min_pae_contact"] = rmin
            results[f"{prefix}_max_pae_contact"] = rmax
            results[f"{prefix}_res_above_contact_thres"] = count
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
    global args
    args = parse_args()

    run_csv_path = Path(args.run_csv)
    run_df = pd.read_csv(run_csv_path, dtype=str)
    if 'binder_id' not in run_df.columns:
        raise SystemExit("--run-csv must contain a 'binder_id' column")
    binder_ids = [str(x).strip() for x in run_df['binder_id'].tolist() if str(x).strip()]

    # Determine sources: explicit or inferred from provided dirs
    sources: List[str] = args.sources if args.sources else []
    
    if not sources:
        if args.boltz_dir: sources.append('boltz')
    if not sources:
        raise SystemExit("Provide at least one of --boltz1-dir/--af3-dir/--colab-dir or specify --sources")

    index = build_file_index(args.boltz_dir)
    

    # backup
    if args.backup:
        bak = run_csv_path.with_suffix(run_csv_path.suffix + '.bak')
        run_df.to_csv(bak, index=False)

        if args.verbose:
            print(f"Backed up to {bak}")

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

    # Write metrics CSV
    res_df = pd.DataFrame.from_dict(collected, orient='index')
    res_df.reset_index(inplace=True)
    res_df = res_df.rename(columns={'index': 'binder_id'})
    res_df.to_csv(args.out_csv, index=False)
    
    print(f"\n{'='*60}")
    print(f"✓ CSV output written successfully")
    print(f"{'='*60}")
    print(f"  File: {args.out_csv}")
    print(f"  Rows: {len(res_df)}")
    print(f"  Columns: {len(res_df.columns)}")
    print(f"{'='*60}\n")

    # Write JSONL if requested
    if args.out_jsonl:
        write_jsonl(collected, args.out_jsonl)

    # Merge with original run CSV
    merged_df = run_df.merge(res_df, on='binder_id', how='left')
    merged_df.to_csv(run_csv_path, index=False)
    
    print(f"{'='*60}")
    print(f"✓ Updated run CSV with metrics")
    print(f"{'='*60}")
    print(f"  File: {run_csv_path}")
    print(f"  Rows: {len(merged_df)}")
    print(f"  Columns: {len(merged_df.columns)}")
    print(f"{'='*60}\n")

    if all_notes:
        print(f"\n{'='*60}")
        print(f"⚠ Notes / Warnings ({len(all_notes)})")
        print(f"{'='*60}")
        for n in all_notes:
            print(f"  - {n}")
        print(f"{'='*60}\n")

if __name__ == '__main__':
    main()
