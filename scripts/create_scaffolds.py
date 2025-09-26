#!/usr/bin/env python
import os
import sys
import glob
import torch
import random
import numpy as np
import argparse
from multiprocessing import Pool, cpu_count
from functools import partial

try:
    import pyrosetta
    pyrosetta.init()
    APPROX = False
except:
    print("WARNING: pyRosetta not found, will use an approximate SSE calculation")
    APPROX = True


def main():
    args = get_args()
    assert args.input_pdb or args.pdb_dir is not None, 'Need to provide either an input pdb (--input_pdb) or a path to pdbs (--pdb_dir)'
    assert not (args.input_pdb is not None and args.pdb_dir is not None), 'Need to provide either --input_pdb or --pdb_dir, not both'

    os.makedirs(args.out_dir, exist_ok=True)
    if args.pdb_dir is not None:
        pdbs = glob.glob(f'{args.pdb_dir}/*pdb')
    else:
        pdbs = [args.input_pdb]

    print(f"Using {args.num_processes} processes")

    # Create a pool of worker processes
    with Pool(processes=args.num_processes) as pool:
        # Use partial to pass the out_dir to the process_pdb function
        process_pdb_partial = partial(process_pdb, out_dir=args.out_dir)

        # Map the pdbs to the process_pdb function
        pool.map(process_pdb_partial, pdbs)


def process_pdb(pdb, out_dir):
    """
    Process a single PDB file.  This function is designed to be run in parallel.
    """
    name = os.path.split(pdb)[1][:-4]  # Extract the base name without extension

    # Check if output files already exist; skip processing if they do
    ss_file = os.path.join(out_dir, f'{name}_ss.pt')
    adj_file = os.path.join(out_dir, f'{name}_adj.pt')

    if os.path.exists(ss_file) and os.path.exists(adj_file):
        print(f"Skipping {pdb} as output files already exist.")
        return

    print(f"Processing {pdb}...")
    secstruc_dict = extract_secstruc(pdb)
    xyz, _, _ = parse_pdb_torch(pdb)
    ss, idx = ss_to_tensor(secstruc_dict)
    block_adj = construct_block_adj_matrix(torch.FloatTensor(ss), torch.tensor(xyz)).float()
    ss_tens, mask = mask_ss(ss, idx, max_mask=0)
    ss_argmax = torch.argmax(ss_tens[:, :4], dim=1).float()

    torch.save(ss_argmax, os.path.join(out_dir, f'{name}_ss.pt'))
    torch.save(block_adj, os.path.join(out_dir, f'{name}_adj.pt'))


def get_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--pdb_dir", required=False,
                        help="path to directory of pdbs. Either pass this or the path to a specific pdb (--input_pdb)",
                        default=None)
    parser.add_argument("--input_pdb", required=False,
                        help="path to input pdb. Either provide this of path to directory of pdbs (--pdb_dir)",
                        default=None)
    parser.add_argument("--out_dir", dest="out_dir", required=True, help='need to specify an output path')
    parser.add_argument("--num_processes", type=int, required=False,
                        help="Number of processes to use for parallel processing (default = 4).",
                        default=4)
    args = parser.parse_args()
    return args


def extract_secstruc(fn):
    pdb = parse_pdb(fn)
    idx = pdb['idx']
    if APPROX:
        aa_sequence = pdb["seq"]
        secstruct = get_sse(pdb["xyz"][:, 1])
    else:
        dssp = pyrosetta.rosetta.core.scoring.dssp
        pose = pyrosetta.io.pose_from_pdb(fn)
        dssp.Dssp(pose).insert_ss_into_pose(pose, True)
        aa_sequence = pose.sequence()
        secstruct = pose.secstruct()
    secstruc_dict = {'sequence': [i for i in aa_sequence],
                     'idx': [int(i) for i in idx],
                     'ss': [i for i in secstruct]}
    return secstruc_dict


def ss_to_tensor(ss):
    """
    Function to convert ss files to indexed tensors
    0 = Helix
    1 = Strand
    2 = Loop
    3 = Mask/unknown
    4 = idx for pdb
    """
    ss_conv = {'H': 0, 'E': 1, 'L': 2}
    idx = np.array(ss['idx'])
    ss_int = np.array([int(ss_conv[i]) for i in ss['ss']])
    return ss_int, idx


def mask_ss(ss, idx, min_mask=0, max_mask=1.0):
    mask_prop = random.uniform(min_mask, max_mask)
    transitions = np.where(ss[:-1] - ss[1:] != 0)[0]  # gets last index of each block of ss
    stuck_counter = 0
    while len(ss[ss == 3]) / len(ss) < mask_prop or stuck_counter > 100:
        width = random.randint(1, 9)
        start = random.choice(transitions)
        offset = random.randint(-8, 1)
        try:

            ss[start + offset:start + offset + width] = 3
        except:
            stuck_counter += 1
            pass
    ss = torch.tensor(ss)
    ss = torch.nn.functional.one_hot(ss, num_classes=4)
    ss = torch.cat((ss, torch.tensor(idx)[..., None]), dim=-1)
    #     mask = torch.where(torch.argmax(ss[:,:-1], dim=-1) == 3, False, True)
    mask = torch.tensor(np.where(np.argmax(ss[:, :-1].numpy(), axis=-1) == 3))
    return ss, mask


def generate_Cbeta(N, Ca, C):
    # recreate Cb given N,Ca,C
    b = Ca - N
    c = C - Ca
    a = torch.cross(b, c, dim=-1)
    # Cb = -0.58273431*a + 0.56802827*b - 0.54067466*c + Ca
    # fd: below matches sidechain generator (=Rosetta params)
    Cb = -0.57910144 * a + 0.5689693 * b - 0.5441217 * c + Ca

    return Cb


def get_pair_dist(a, b):
    """calculate pair distances between two sets of points

    Parameters
    ----------
    a,b : pytorch tensors of shape [batch,nres,3]
          store Cartesian coordinates of two sets of atoms
    Returns
    -------
    dist : pytorch tensor of shape [batch,nres,nres]
           stores paitwise distances between atoms in a and b
    """

    dist = torch.cdist(a, b, p=2)
    return dist


def construct_block_adj_matrix(sstruct, xyz, cutoff=6, include_loops=False):
    '''
    Given a sstruct specification and backbone coordinates, build a block adjacency matrix.

    Input:

        sstruct (torch.FloatTensor): (L) length tensor with numeric encoding of sstruct at each position

        xyz (torch.FloatTensor): (L,3,3) tensor of Cartesian coordinates of backbone N,Ca,C atoms

        cutoff (float): The Cb distance cutoff under which residue pairs are considered adjacent
                        By eye, Nate thinks 6A is a good Cb distance cutoff

    Output:

        block_adj (torch.FloatTensor): (L,L) boolean matrix where adjacent secondary structure contacts are 1
    '''

    L = xyz.shape[0]

    # three anchor atoms
    N = xyz[:, 0]
    Ca = xyz[:, 1]
    C = xyz[:, 2]

    # recreate Cb given N,Ca,C
    Cb = generate_Cbeta(N, Ca, C)

    # May need a batch dimension - NRB
    dist = get_pair_dist(Cb, Cb)  # [L,L]
    dist[torch.isnan(dist)] = 999.9

    dist += 999.9 * torch.eye(L, device=xyz.device)
    # Now we have dist matrix and sstruct specification, turn this into a block adjacency matrix
    # There is probably a way to do this in closed-form with a beautiful einsum but I am going to do the loop approach

    # First: Construct a list of segments and the index at which they begin and end
    in_segment = True
    segments = []

    begin = -1
    end = -1

    for i in range(sstruct.shape[0]):
        # Starting edge case
        if i == 0:
            begin = 0
            continue

        if not sstruct[i] == sstruct[i - 1]:
            end = i
            segments.append((sstruct[i - 1], begin, end))

            begin = i

    # Ending edge case: last segment is length one
    if not end == sstruct.shape[0]:
        segments.append((sstruct[-1], begin, sstruct.shape[0]))
    block_adj = torch.zeros_like(dist)
    for i in range(len(segments)):
        curr_segment = segments[i]

        if curr_segment[0] == 2 and not include_loops: continue

        begin_i = curr_segment[1]
        end_i = curr_segment[2]
        for j in range(i + 1, len(segments)):
            j_segment = segments[j]

            if j_segment[0] == 2 and not include_loops: continue

            begin_j = j_segment[1]
            end_j = j_segment[2]

            if torch.any(dist[begin_i:end_i, begin_j:end_j] < cutoff):
                # Matrix is symmetic
                block_adj[begin_i:end_i, begin_j:end_j] = torch.ones(end_i - begin_i, end_j - begin_j)
                block_adj[begin_j:end_j, begin_i:end_i] = torch.ones(end_j - begin_j, end_i - begin_i)
    return block_adj


def parse_pdb_torch(filename):
    lines = open(filename, 'r').readlines()
    return parse_pdb_lines_torch(lines)


# '''
def parse_pdb_lines_torch(lines):
    # indices of residues observed in the structure
    pdb_idx = []
    for l in lines:
        if l[:4] == "ATOM" and l[12:16].strip() == "CA":
            idx = (l[21:22].strip(), int(l[22:26].strip()))
            if idx not in pdb_idx:
                pdb_idx.append(idx)

    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(pdb_idx), 27, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' ' + l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain, resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm == atom:
                xyz[idx, i_atm, :] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break
    # save atom mask
    mask = np.logical_not(np.isnan(xyz[..., 0]))
    xyz[np.isnan(xyz[..., 0])] = 0.0

    return xyz, mask, np.array(pdb_idx)


def parse_pdb(filename, **kwargs):
    '''extract xyz coords for all heavy atoms'''
    lines = open(filename, 'r').readlines()
    return parse_pdb_lines(lines, **kwargs)


def parse_pdb_lines(lines, parse_hetatom=False, ignore_het_h=True):
    # indices of residues observed in the structure
    res = [(l[22:26], l[17:20]) for l in lines if l[:4] == "ATOM" and l[12:16].strip() == "CA"]
    seq = [aa2num[r[1]] if r[1] in aa2num.keys() else 20 for r in res]
    pdb_idx = [(l[21:22].strip(), int(l[22:26].strip())) for l in lines if
               l[:4] == "ATOM" and l[12:16].strip() == "CA"]  # chain letter, res num

    # 4 BB + up to 10 SC atoms
    xyz = np.full((len(res), 27, 3), np.nan, dtype=np.float32)
    for l in lines:
        if l[:4] != "ATOM":
            continue
        chain, resNo, atom, aa = l[21:22], int(l[22:26]), ' ' + l[12:16].strip().ljust(3), l[17:20]
        idx = pdb_idx.index((chain, resNo))
        for i_atm, tgtatm in enumerate(aa2long[aa2num[aa]]):
            if tgtatm is not None and tgtatm.strip() == atom.strip():  # ignore whitespace
                xyz[idx, i_atm, :] = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
                break

    # save atom mask
    mask = np.logical_not(np.isnan(xyz[..., 0]))
    xyz[np.isnan(xyz[..., 0])] = 0.0
    # remove duplicated (chain, resi)
    new_idx = []
    i_unique = []
    for i, idx in enumerate(pdb_idx):
        if idx not in new_idx:
            new_idx.append(idx)
            i_unique.append(i)

    pdb_idx = new_idx
    xyz = xyz[i_unique]
    mask = mask[i_unique]
    seq = np.array(seq)[i_unique]

    out = {'xyz': xyz,  # cartesian coordinates, [Lx14]
           'mask': mask,  # mask showing which atoms are present in the PDB file, [Lx14]
           'idx': np.array([i[1] for i in pdb_idx]),  # residue numbers in the PDB file, [L]
           'seq': np.array(seq),  # amino acid sequence, [L]
           'pdb_idx': pdb_idx,  # list of (chain letter, residue number) in the pdb file, [L]
           }
    # heteroatoms (ligands, etc)
    if parse_hetatom:
        xyz_het, info_het = [], []
        for l in lines:
            if l[:6] == 'HETATM' and not (ignore_het_h and l[77] == 'H'):
                info_het.append(dict(
                    idx=int(l[7:11]),
                    atom_id=l[12:16],
                    atom_type=l[77],
                    name=l[16:20]
                ))
                xyz_het.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])

        out['xyz_het'] = np.array(xyz_het)
        out['info_het'] = info_het

    return out


num2aa = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
    'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO',
    'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'UNK', 'MAS',
]
aa2num = {x: i for i, x in enumerate(num2aa)}
# full sc atom representation (Nx14)
aa2long = [
    (" N  ", " CA ", " C  ", " O  ", " CB ", None, None, None, None, None, None, None, None, None, " H  ", " HA ", "1HB ",
     "2HB ", "3HB ", None, None, None, None, None, None, None, None),  # ala
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " NE ", " CZ ", " NH1", " NH2", None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", "1HG ", "2HG ", "1HD ", "2HD ", " HE ", "1HH1", "2HH1", "1HH2", "2HH2"),  # arg
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " ND2", None, None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", "1HD2", "2HD2", None, None, None, None, None, None, None),  # asn
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " OD1", " OD2", None, None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", None, None, None, None, None, None, None, None, None),  # asp
    (" N  ", " CA ", " C  ", " O  ", " CB ", " SG ", None, None, None, None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", " HG ", None, None, None, None, None, None, None, None),  # cys
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " OE1", " NE2", None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", "1HG ", "2HG ", "1HE2", "2HE2", None, None, None, None, None),  # gln
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " OE1", " OE2", None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", "1HG ", "2HG ", None, None, None, None, None, None, None),  # glu
    (" N  ", " CA ", " C  ", " O  ", None, None, None, None, None, None, None, None, None, None, " H  ", "1HA ", "2HA ",
     None, None, None, None, None, None, None, None, None, None),  # gly
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " ND1", " CD2", " CE1", " NE2", None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", " HD2", " HE1", " HE2", None, None, None, None, None, None),  # his
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2", " CD1", None, None, None, None, None, None, " H  ", " HA ", " HB ",
     "1HG2", "2HG2", "3HG2", "1HG1", "2HG1", "1HD1", "2HD1", "3HD1", None, None),  # ile
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2", None, None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", " HG ", "1HD1", "2HD1", "3HD1", "1HD2", "2HD2", "3HD2", None, None),  # leu
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", " CE ", " NZ ", None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", "1HG ", "2HG ", "1HD ", "2HD ", "1HE ", "2HE ", "1HZ ", "2HZ ", "3HZ "),  # lys
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " SD ", " CE ", None, None, None, None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", "1HG ", "2HG ", "1HE ", "2HE ", "3HE ", None, None, None, None),  # met
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", None, None, None, " H  ", " HA ",
     "1HB ", "2HB ", " HD1", " HD2", " HE1", " HE2", " HZ ", None, None, None, None),  # phe
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD ", None, None, None, None, None, None, None, " HA ", "1HB ", "2HB ",
     "1HG ", "2HG ", "1HD ", "2HD ", None, None, None, None, None, None),  # pro
    (" N  ", " CA ", " C  ", " O  ", " CB ", " OG ", None, None, None, None, None, None, None, None, " H  ", " HG ", " HA ",
     "1HB ", "2HB ", None, None, None, None, None, None, None, None),  # ser
    (" N  ", " CA ", " C  ", " O  ", " CB ", " OG1", " CG2", None, None, None, None, None, None, None, " H  ", " HG1", " HA ",
     " HB ", "1HG2", "2HG2", "3HG2", None, None, None, None, None, None),  # thr
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2", " NE1", " CE2", " CE3", " CZ2", " CZ3", " CH2", " H  ", " HA ",
     "1HB ", "2HB ", " HD1", " HE1", " HZ2", " HH2", " HZ3", " HE3", None, None, None),  # trp
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG ", " CD1", " CD2", " CE1", " CE2", " CZ ", " OH ", None, None, " H  ", " HA ",
     "1HB ", "2HB ", " HD1", " HE1", " HE2", " HD2", " HH ", None, None, None, None),  # tyr
    (" N  ", " CA ", " C  ", " O  ", " CB ", " CG1", " CG2", None, None, None, None, None, None, None, " H  ", " HA ", " HB ",
     "1HG1", "2HG1", "3HG1", "1HG2", "2HG2", "3HG2", None, None, None, None),  # val
    (" N  ", " CA ", " C  ", " O  ", " CB ", None, None, None, None, None, None, None, None, None, " H  ", " HA ", "1HB ",
     "2HB ", "3HB ", None, None, None, None, None, None, None, None),  # unk
    (" N  ", " CA ", " C  ", " O  ", " CB ", None, None, None, None, None, None, None, None, None, " H  ", " HA ", "1HB ",
     "2HB ", "3HB ", None, None, None, None, None, None, None, None),  # mask
]


def get_sse(ca_coord):
    '''
    calculates the SSE of a peptide chain based on the P-SEA algorithm (Labesse 1997)
    code borrowed from biokite: https://github.com/biokit/biokit
    '''

    def vector_dot(v1, v2):
        return (v1 * v2).sum(-1)

    def norm_vector(v):
        return v / np.linalg.norm(v, axis=-1, keepdims=True)

    def displacement(atoms1, atoms2):
        v1 = np.asarray(atoms1)
        v2 = np.asarray(atoms2)
        if len(v1.shape) <= len(v2.shape):
            diff = v2

    # You will need to complete the displacement method with the equation
    # diff = v2 - v1[:v2.shape[0]]
    # to get this working

if __name__ == "__main__":
    main()
