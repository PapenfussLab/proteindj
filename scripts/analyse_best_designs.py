#!/usr/bin/env python3

import copy
import json
import re
import time
import logging
import polars as pl
from functools import partial
from multiprocessing import Pool
from pathlib import Path

import arpeggia
import openmm
from openmm import unit
from openmm.app import ForceField, NoCutoff
from openmm.app import PDBFile as OpenMMPDBFile
from pdbfixer import PDBFixer
from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Polypeptide
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.SeqUtils import seq1
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from metrics_utils import compute_dssp_chars_by_chain, count_ss_elements, calculate_rog

def setup_logging():
    """Setup logging configuration."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('analysis.log', mode='w'),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


logger = setup_logging()

# Chothia/NACCESS-like atomic radii for SASA calculation
R_CHOTHIA = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}

# Hydrophobic amino acids (matches FreeBindCraft HYDROPHOBIC_AA_SET)
HYDROPHOBIC_AA = frozenset("ACFGILMPVWY")


# ---------------------------------------------------------------------------
# Core structure utilities
# ---------------------------------------------------------------------------

def parse_structure(pdb_path):
    """Parse PDB file and return (structure, model)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", str(pdb_path))
    return structure, structure[0]


def isolate_chain(model, chain_id):
    """Create a new structure containing only the specified chain.
    Returns (isolated_structure, isolated_model, isolated_chain).
    """
    iso_struct = Structure("iso")
    iso_model = Model(0)
    iso_chain = copy.deepcopy(model[chain_id])
    iso_model.add(iso_chain)
    iso_struct.add(iso_model)
    return iso_struct, iso_model, iso_chain


def get_chain_ids(model):
    """Get ordered list of chain IDs from a BioPython model."""
    return [chain.get_id() for chain in model.get_chains()]


def derive_ids_from_filename(filename):
    """Extract fold_id and seq_id from filename format: fold_X_seq_Y_*.pdb"""
    basename = Path(filename).stem
    fold_match = re.search(r'fold_(\d+)', basename)
    seq_match = re.search(r'seq_(\d+)', basename)

    fold_id = int(fold_match.group(1)) if fold_match else None
    seq_id = int(seq_match.group(1)) if seq_match else None

    return fold_id, seq_id


# ---------------------------------------------------------------------------
# Structure relaxation (PDBFixer + OpenMM)
# ---------------------------------------------------------------------------
#
# Boltz/AF2-predicted PDBs have no explicit hydrogens and commonly omit OXT
# atoms at true chain termini, which prevents ForceField template matching.
# Every metric in this script is computed on the relaxed (hydrogens added +
# energy-minimized) structure, not the raw prediction: a short restrained
# minimization settles newly-added hydrogens (and any restored heavy atoms)
# into a sensible local energy minimum, giving more reliable geometry for
# H-bond detection, shape complementarity, SASA, etc. PDBFixer/OpenMM reset
# per-atom B-factors to 0.00 when writing output, so the original per-residue
# B-factors (pLDDT from the prediction model) are restored afterwards.

RELAXED_PDB_DIRNAME = 'relaxed_pdbs'


def _get_residue_bfactors(pdb_path):
    """Map (chain_id, resnum) -> B-factor using each residue's CA atom.

    Boltz/AF2 write a uniform per-residue pLDDT value into the B-factor
    column, so the CA atom is a sufficient representative for the whole
    residue (same approach as scripts/prep_fampnn_designs.py's
    get_residue_metadata()).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('orig', str(pdb_path))
    bfactors = {}
    for chain in structure[0]:
        for residue in chain:
            for atom in residue:
                if atom.get_name() == 'CA':
                    resnum = residue.get_id()[1]
                    bfactors[(chain.get_id(), resnum)] = atom.get_bfactor()
                    break
    return bfactors


def _restore_bfactors(relaxed_pdb_path, residue_bfactors):
    """Rewrite the B-factor column of a relaxed PDB using residue_bfactors.

    Applies the residue's B-factor to every atom in that residue, including
    hydrogens and any other atoms PDBFixer added that have no counterpart in
    the original structure.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('relaxed', str(relaxed_pdb_path))
    for chain in structure[0]:
        for residue in chain:
            resnum = residue.get_id()[1]
            bfactor = residue_bfactors.get((chain.get_id(), resnum), 0.0)
            for atom in residue:
                atom.set_bfactor(bfactor)

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(relaxed_pdb_path))


def _reorder_chains_alphabetically(pdb_path):
    """Rewrite pdb_path with chain blocks physically reordered to match alphabetical ID order.

    OpenMM's PDBFile.writeFile() assigns output chain IDs by physical position (1st
    chain in the topology -> 'A', 2nd -> 'B', ...), ignoring the original chain.id
    string. Boltz can write its split multi-chain target sub-chains in a different
    physical order than their alphabetical labels (e.g. grouping identical-sequence
    chains together), so without this normalization, relaxation silently swaps which
    residues end up under which chain ID downstream.
    """
    pdb_path = Path(pdb_path)
    with open(pdb_path) as fin:
        lines = fin.readlines()

    chain_lines = {}
    chain_order = []
    other_lines = []
    for line in lines:
        if line.startswith(('ATOM', 'HETATM', 'TER', 'ANISOU')):
            chain_id = line[21]
            chain_lines.setdefault(chain_id, []).append(line)
            if chain_id not in chain_order:
                chain_order.append(chain_id)
        elif not line.startswith('END'):
            other_lines.append(line)

    if chain_order == sorted(chain_order):
        return pdb_path

    out_path = pdb_path.parent / f'{pdb_path.stem}_chainsorted.pdb'
    with open(out_path, 'w') as fout:
        fout.writelines(other_lines)
        for chain_id in sorted(chain_order):
            fout.writelines(chain_lines[chain_id])
        fout.write('END\n')
    return out_path


def _relax_structure(pdb_path, ph=7.0, restraint_k=1000.0):
    """Add hydrogens and energy-minimize a PDB file, preserving B-factors.

    Uses PDBFixer to restore missing heavy atoms (Boltz/AF2 predictions
    commonly omit OXT at true C-termini) and add hydrogens, then runs a
    short OpenMM minimization with heavy atoms harmonically restrained to
    their original positions so only hydrogens (and any restored atoms)
    relax into sensible geometry. The relaxed structure is written to
    RELAXED_PDB_DIRNAME alongside the input file and is used for every
    metric computed by this script.

    Returns the path to the relaxed PDB file.
    """
    pdb_path = Path(pdb_path)
    residue_bfactors = _get_residue_bfactors(pdb_path)

    fixer_input_path = _reorder_chains_alphabetically(pdb_path)
    fixer = PDBFixer(filename=str(fixer_input_path))
    if fixer_input_path != pdb_path:
        fixer_input_path.unlink()  # temporary chain-sorted copy, already loaded
    # Boltz/AF2-predicted PDBs often omit OXT atoms at genuine chain termini,
    # which otherwise makes ForceField unable to tell whether a chain's last
    # residue is an internal residue (expects a peptide bond onward) or a
    # true C-terminus (expects OXT). Restore any missing heavy atoms first.
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    # Register disulfide bonds explicitly so ForceField can unambiguously
    # match the CYX (disulfide-bonded) template instead of confusing it with
    # CYM (deprotonated thiolate), which looks identical without this bond.
    fixer.topology.createDisulfideBonds(fixer.positions)

    forcefield = ForceField('amber14-all.xml')
    system = forcefield.createSystem(fixer.topology, nonbondedMethod=NoCutoff)

    # Harmonically restrain heavy atoms so minimization only relaxes hydrogens
    restraint = openmm.CustomExternalForce('k*((x-x0)^2+(y-y0)^2+(z-z0)^2)')
    restraint.addGlobalParameter('k', restraint_k * unit.kilojoule_per_mole / unit.nanometer**2)
    for name in ('x0', 'y0', 'z0'):
        restraint.addPerParticleParameter(name)
    for atom in fixer.topology.atoms():
        if atom.element is not None and atom.element.symbol != 'H':
            pos = fixer.positions[atom.index]
            restraint.addParticle(atom.index, [pos.x, pos.y, pos.z])
    system.addForce(restraint)

    integrator = openmm.VerletIntegrator(1.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator, openmm.Platform.getPlatformByName('CPU'))
    context.setPositions(fixer.positions)

    energy_before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilocalories_per_mole)

    t0 = time.time()
    openmm.LocalEnergyMinimizer.minimize(context)
    elapsed = time.time() - t0

    state_after = context.getState(getEnergy=True, getPositions=True)
    energy_after = state_after.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole)
    positions = state_after.getPositions()

    logger.info(
        f"OpenMM minimization ({pdb_path.name}): "
        f"{energy_before:.1f} -> {energy_after:.1f} kcal/mol in {elapsed:.2f}s"
    )

    relaxed_dir = pdb_path.parent / RELAXED_PDB_DIRNAME
    relaxed_dir.mkdir(exist_ok=True)
    out_path = relaxed_dir / pdb_path.name

    with open(out_path, 'w') as f:
        OpenMMPDBFile.writeFile(fixer.topology, positions, f)

    _restore_bfactors(out_path, residue_bfactors)

    return out_path


def _validate_relaxed_structure(orig_model, relaxed_model, chain_ids):
    """Raise if relaxation changed the number of residues in any chain."""
    for chain_id in chain_ids:
        orig_count = sum(1 for _ in orig_model[chain_id])
        relaxed_count = sum(1 for _ in relaxed_model[chain_id])
        if orig_count != relaxed_count:
            raise RuntimeError(
                f"Relaxation changed residue count for chain {chain_id}: "
                f"{orig_count} -> {relaxed_count}"
            )


# ---------------------------------------------------------------------------
# SASA and surface chemistry
# ---------------------------------------------------------------------------

def compute_sasa(model):
    """Compute atom-level SASA using BioPython ShrakeRupley with Chothia radii."""
    sr = ShrakeRupley(probe_radius=1.40, n_points=960, radii_dict=R_CHOTHIA)
    sr.compute(model, level='A')


def calculate_surface_chemistry(model, chain_id):
    """Calculate surface hydrophobicity using BioPython ShrakeRupley SASA.

    Computes the percentage of SASA contributed by hydrophobic residues.
    Note: this differs from PyRosetta LayerSelector (which classifies residues
    as surface/core then counts) but captures equivalent intent.
    """
    _, iso_model, iso_chain = isolate_chain(model, chain_id)
    compute_sasa(iso_model)

    hydrophobic_sasa = 0.0
    total_sasa = 0.0

    for residue in iso_chain:
        if Polypeptide.is_aa(residue, standard=True):
            res_sasa = sum(getattr(atom, 'sasa', 0.0) for atom in residue.get_atoms())
            total_sasa += res_sasa
            try:
                aa1 = seq1(residue.get_resname()).upper()
            except Exception:
                aa1 = ''
            if aa1 in HYDROPHOBIC_AA:
                hydrophobic_sasa += res_sasa

    if total_sasa > 0:
        pct = round(hydrophobic_sasa / total_sasa * 100.0, 1)
    else:
        pct = 0.0

    return {'pr_surfhphobics': pct}


# ---------------------------------------------------------------------------
# Spatial Aggregation Propensity (SAP) via arpeggia
# ---------------------------------------------------------------------------
#
# PyRosetta's per-residue SAP calculation only counts an atom's contribution
# toward its residue's score when that contribution is positive, which is
# the same behavior as arpeggia's own level='residue' aggregation. The one
# real difference is that PyRosetta reports a value (0.0 if no atom
# qualifies) for every residue being scored, while arpeggia's residue-level
# output simply omits residues where no atom qualifies. Since an omitted
# residue would contribute 0 to the sum anyway, we sum arpeggia's reported
# per-residue sap_score column but divide by the TRUE total residue count of
# the chain (not len(df)) to keep the mean consistent with PyRosetta.

def _count_chain_residues(model, chain_id):
    """Count amino-acid residues in a chain (SAP mean denominator)."""
    return sum(1 for res in model[chain_id] if Polypeptide.is_aa(res, standard=True))


def _sap_chain_sum(pdb_path, chains, target_chain, probe_radius=1.4, n_points=960, sap_radius=5.0):
    """Sum of arpeggia's per-residue SAP scores (already filtered to positive
    atom contributions per residue) for a single chain.

    `chains` restricts both the SASA calculation and the neighbor search to
    the given comma-separated chain IDs (empty string = all chains).
    """
    df = arpeggia.sap_score(
        str(pdb_path),
        level='residue',
        probe_radius=probe_radius,
        n_points=n_points,
        sap_radius=sap_radius,
        chains=chains,
    )
    if df.height == 0:
        return 0.0
    scores = df.filter(pl.col('chain') == target_chain)['sap_score']
    return float(scores.sum()) if scores.len() else 0.0


def calculate_sap_scores(model, pdb_path, num_chains, chain_id='A', complex_sap_df=None):
    """
    Calculate Spatial Aggregation Propensity (SAP) scores using arpeggia.

    For monomer: Returns pr_SAP only
    For complex: Returns pr_SAP (alone) and pr_SAP_complex (with target shielding)

    `complex_sap_df`, if provided, is a pre-computed whole-structure
    residue-level SAP dataframe (arpeggia.sap_score(pdb_path, level='residue',
    chains='')) to avoid recomputing SASA/SAP for the full complex once per
    chain in multi-chain designs.
    """
    metrics = {}
    n_res = _count_chain_residues(model, chain_id)

    # pr_SAP: chain calculated in isolation (matches old split_by_chain behavior)
    alone_sum = _sap_chain_sum(pdb_path, chains=chain_id, target_chain=chain_id)
    metrics['pr_SAP'] = round(alone_sum / n_res, 3) if n_res else 0.0

    # pr_SAP_complex: chain's SAP within full complex context (all chains as neighbors)
    if num_chains >= 2:
        if complex_sap_df is None:
            complex_sap_df = arpeggia.sap_score(
                str(pdb_path), level='residue', probe_radius=1.4, n_points=960,
                sap_radius=5.0, chains='',
            )
        complex_scores = (
            complex_sap_df.filter(pl.col('chain') == chain_id)['sap_score']
            if complex_sap_df.height else None
        )
        complex_sum = float(complex_scores.sum()) if complex_scores is not None and complex_scores.len() else 0.0
        metrics['pr_SAP_complex'] = round(complex_sum / n_res, 3) if n_res else 0.0

    return metrics


# ---------------------------------------------------------------------------
# Interface metrics
# ---------------------------------------------------------------------------

def chain_total_sasa(chain):
    """Sum SASA over all atoms in a chain."""
    return sum(getattr(atom, "sasa", 0.0) for atom in chain.get_atoms())


def calculate_shape_complementarity(pdb_path, chain1, chain2):
    """Calculate shape complementarity using arpeggia's native `sc` function
    (Lawrence & Colman 1993 algorithm), implemented in Rust and bundled with
    the same `arpeggia` package used elsewhere in this script -- no external
    binary required.

    Returns None if the two chains have no interface close enough for SC's
    own (stricter than the BSA dSASA check) contact dot generation.
    """
    try:
        sc_val = arpeggia.sc(str(pdb_path), groups=f'{chain1}/{chain2}')
    except RuntimeError as e:
        if "No molecular dots generated" in str(e):
            return None
        raise
    if not (0.0 <= sc_val <= 1.0):
        raise ValueError(f"arpeggia.sc returned out-of-range value {sc_val} for {pdb_path}")
    return sc_val


def compute_interface_bsa(model, chain1, chain2):
    """Compute buried surface area at the interface using BioPython ShrakeRupley.

    Adapted from FreeBindCraft _compute_sasa_metrics().
    Returns dSASA / 2 (to match PyRosetta InterfaceAnalyzerMover convention).
    """
    # SASA of chains within the full complex
    compute_sasa(model)
    chain1_sasa_complex = chain_total_sasa(model[chain1])
    chain2_sasa_complex = chain_total_sasa(model[chain2])

    # SASA of chain1 as isolated monomer
    _, iso_model1, iso_chain1 = isolate_chain(model, chain1)
    compute_sasa(iso_model1)
    chain1_sasa_mono = chain_total_sasa(iso_chain1)

    # SASA of chain2 as isolated monomer
    _, iso_model2, iso_chain2 = isolate_chain(model, chain2)
    compute_sasa(iso_model2)
    chain2_sasa_mono = chain_total_sasa(iso_chain2)

    # dSASA = buried area upon complex formation
    dsasa = max(chain1_sasa_mono - chain1_sasa_complex, 0.0) + \
            max(chain2_sasa_mono - chain2_sasa_complex, 0.0)

    return round(dsasa / 2)


def compute_prodigy_scores(pdb_path, chain1, chain2):
    """Predict binding affinity using Prodigy (IC-NIS model).

    Returns dict with deltaG (kcal/mol) and Kd (M), or None if the two
    chains have no interface contacts (e.g. a failed/non-binding design)
    -- Prodigy raises ValueError("No contacts found for selection") in
    that case, which is a real, expected outcome rather than a bug.
    """
    from prodigy_prot.modules.parsers import parse_structure as prodigy_parse
    from prodigy_prot.modules.prodigy import Prodigy as ProdigyPredictor

    models, _, _ = prodigy_parse(str(pdb_path))
    prodigy = ProdigyPredictor(model=models[0], selection=[chain1, chain2])
    try:
        prodigy.predict()
    except ValueError as e:
        if "No contacts found for selection" in str(e):
            return None
        raise
    return {
        'dg': prodigy.ba_val,
        'kd': prodigy.kd_val,
    }


# ---------------------------------------------------------------------------
# H-bond metrics (arpeggia, on the relaxed/hydrogenated structure)
# ---------------------------------------------------------------------------

# H-bond donor/acceptor atom names by residue type (heavy atoms only).
# Backbone N is donor, backbone O is acceptor for all residues.
SIDECHAIN_DONORS = {
    'ARG': ('NE', 'NH1', 'NH2'), 'ASN': ('ND2',), 'CYS': ('SG',),
    'GLN': ('NE2',), 'HIS': ('ND1', 'NE2'), 'LYS': ('NZ',), 'SER': ('OG',),
    'THR': ('OG1',), 'TRP': ('NE1',), 'TYR': ('OH',),
}
SIDECHAIN_ACCEPTORS = {
    'ASN': ('OD1',), 'ASP': ('OD1', 'OD2'), 'CYS': ('SG',), 'GLN': ('OE1',),
    'GLU': ('OE1', 'OE2'), 'HIS': ('ND1', 'NE2'), 'MET': ('SD',),
    'SER': ('OG',), 'THR': ('OG1',), 'TYR': ('OH',),
}


def get_interface_polar_atoms(model, chain1, chain2, contact_dist=4.0):
    """Identify polar atoms at the interface between two chains.

    Returns (donors_c1, acceptors_c1, donors_c2, acceptors_c2) where each
    is a list of Atom objects from the respective chain that are within
    contact_dist of any atom on the other chain.
    """
    atoms_c1 = list(model[chain1].get_atoms())
    atoms_c2 = list(model[chain2].get_atoms())

    # Build neighbor search on each chain to find interface atoms
    ns_c2 = NeighborSearch(atoms_c2)
    ns_c1 = NeighborSearch(atoms_c1)

    def is_interface_atom(atom, ns_other, dist):
        """Check if atom is within dist of any atom in the other chain."""
        return len(ns_other.search(atom.get_vector().get_array(), dist, 'A')) > 0

    def collect_polar(chain_obj, ns_other, dist):
        """Collect donor and acceptor atoms at the interface for a chain."""
        donors = []
        acceptors = []
        for residue in chain_obj:
            if not Polypeptide.is_aa(residue, standard=True):
                continue
            resname = residue.get_resname()
            for atom in residue:
                aname = atom.get_name()
                if not is_interface_atom(atom, ns_other, dist):
                    continue
                # Backbone
                if aname == 'N':
                    donors.append(atom)
                elif aname == 'O':
                    acceptors.append(atom)
                # Sidechain donors
                elif aname in SIDECHAIN_DONORS.get(resname, ()):
                    donors.append(atom)
                # Sidechain acceptors
                if aname in SIDECHAIN_ACCEPTORS.get(resname, ()):
                    acceptors.append(atom)
        return donors, acceptors

    donors_c1, acceptors_c1 = collect_polar(model[chain1], ns_c2, contact_dist)
    donors_c2, acceptors_c2 = collect_polar(model[chain2], ns_c1, contact_dist)

    return donors_c1, acceptors_c1, donors_c2, acceptors_c2


def compute_interface_hbonds(pdb_path, chain1, chain2):
    """Count intermolecular hydrogen bonds using arpeggia's contacts().

    `pdb_path` must be the relaxed (hydrogens-added) structure so arpeggia
    can apply its donor-H...acceptor distance + angle (>= 90 deg) criterion,
    rather than falling back to a distance-only polar contact.
    """
    df = arpeggia.contacts(str(pdb_path), groups=f'{chain1}/{chain2}')
    if df.height == 0:
        return 0
    return df.filter(pl.col('interaction') == 'HydrogenBond').height


def _get_hbond_satisfied_atoms(pdb_path):
    """Heavy donor/acceptor atoms with a satisfied hydrogen bond anywhere in
    the structure (intra- or inter-chain), per arpeggia's angle-dependent
    detection on the relaxed (hydrogens-added) structure.

    Uses the default (whole-structure) groups so that intra-chain backbone
    H-bonds (e.g. helical i, i-4 pattern) count as satisfying -- a buried
    polar atom should not be flagged as unsatisfied just because its H-bond
    partner happens to be on the same chain.

    Returns a set of (chain, resi, atom_name) tuples.
    """
    df = arpeggia.contacts(str(pdb_path))
    if df.height == 0:
        return set()

    hb = df.filter(pl.col('interaction') == 'HydrogenBond')
    satisfied = set()
    for row in hb.iter_rows(named=True):
        satisfied.add((row['from_chain'], row['from_resi'], row['from_atomn']))
        satisfied.add((row['to_chain'], row['to_resi'], row['to_atomn']))
    return satisfied


def compute_unsat_hbonds(model, pdb_path, chain1, chain2):
    """Buried unsatisfied hydrogen bonds at the interface.

    Identifies polar atoms at the interface that are buried (atom SASA <
    2.0 A^2 in the complex) and lack a satisfied hydrogen bond, per
    arpeggia's angle-dependent detection on the relaxed structure.
    """
    compute_sasa(model)

    donors_c1, acceptors_c1, donors_c2, acceptors_c2 = get_interface_polar_atoms(
        model, chain1, chain2, contact_dist=5.0)

    satisfied = _get_hbond_satisfied_atoms(pdb_path)

    burial_cutoff = 2.0  # A^2 -- atom is considered buried below this SASA

    def is_satisfied(atom):
        chain_id = atom.get_parent().get_parent().get_id()
        resi = atom.get_parent().get_id()[1]
        return (chain_id, resi, atom.get_name()) in satisfied

    unsat_count = 0
    for donor in donors_c1 + donors_c2:
        if getattr(donor, 'sasa', 999.0) >= burial_cutoff:
            continue
        if not is_satisfied(donor):
            unsat_count += 1

    for acceptor in acceptors_c1 + acceptors_c2:
        if getattr(acceptor, 'sasa', 999.0) >= burial_cutoff:
            continue
        if not is_satisfied(acceptor):
            unsat_count += 1

    return unsat_count


def calculate_interface_metrics(model, pdb_path, chain1='A', chain2='B'):
    """Calculate interface metrics between specified chains on the relaxed
    structure.

    Computes BSA, shape complementarity (arpeggia), binding deltaG (Prodigy),
    hydrogen bonds, and unsatisfied hydrogen bonds (arpeggia, angle-dependent).

    If the two chains do not contact each other at all (e.g. a failed or
    non-binding design), BSA is 0 and shape complementarity/binding
    affinity/H-bonds are undefined -- these are reported as 0 rather than
    attempting a calculation that both arpeggia and Prodigy would reject.
    """
    # Add suffix for metrics if not the default A_B interface
    suffix = '' if (chain1 == 'A' and chain2 == 'B') else f'_{chain1}_{chain2}'

    # Geometry metrics
    bsa = compute_interface_bsa(model, chain1, chain2)

    if bsa == 0:
        logger.warning(
            f"No interface contact between chain {chain1} and chain {chain2} "
            f"in {pdb_path} (BSA=0); reporting zero-valued interface metrics."
        )
        return {
            f'pr_intface_BSA{suffix}': 0,
            f'pr_intface_shpcomp{suffix}': 0.0,
            f'pr_intface_deltaG{suffix}': 0.0,
            f'pr_intface_deltaGtoBSA{suffix}': 0.0,
            f'pr_intface_hbonds{suffix}': 0,
            f'pr_intface_unsat_hbonds{suffix}': 0,
        }

    sc_val = calculate_shape_complementarity(pdb_path, chain1, chain2)
    sc = round(sc_val, 3) if sc_val is not None else 0.0

    # Binding affinity (Prodigy IC-NIS model); None if no contacts found
    pg = compute_prodigy_scores(pdb_path, chain1, chain2)
    dg = pg['dg'] if pg is not None else 0.0
    dg_to_bsa = round(dg / bsa, 4)

    # Hydrogen bonds and buried unsatisfied hydrogen bonds (arpeggia)
    hbonds = compute_interface_hbonds(pdb_path, chain1, chain2)
    unsat = compute_unsat_hbonds(model, pdb_path, chain1, chain2)

    return {
        f'pr_intface_BSA{suffix}': bsa,
        f'pr_intface_shpcomp{suffix}': sc,
        f'pr_intface_deltaG{suffix}': round(dg, 2),
        f'pr_intface_deltaGtoBSA{suffix}': dg_to_bsa,
        f'pr_intface_hbonds{suffix}': hbonds,
        f'pr_intface_unsat_hbonds{suffix}': unsat,
    }


# ---------------------------------------------------------------------------
# Sequence metrics
# ---------------------------------------------------------------------------

def get_chain_sequence(pdb_path, chain_id):
    """Extract sequence for a specific chain from PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('temp', str(pdb_path))

    # Convert numeric chain IDs to letters if needed
    if chain_id.isdigit():
        chain_id = chr(64 + int(chain_id))  # Convert 1->A, 2->B etc

    for chain in structure.get_chains():
        if chain.id == chain_id:
            residues = [r for r in chain.get_residues() if PDB.is_aa(r)]
            return ''.join([Polypeptide.protein_letters_3to1[r.get_resname()] for r in residues])

    raise ValueError(f"Chain {chain_id} not found in {pdb_path}")


def calculate_seq_metrics(sequence):
    """Calculate metrics for a protein sequence."""
    analysis = ProteinAnalysis(sequence)
    extinction_coef = analysis.molar_extinction_coefficient()[0]
    molecular_weight = int(analysis.molecular_weight())

    return {
        'sequence': sequence,
        'seq_ext_coef': extinction_coef,
        'seq_length': len(sequence),
        'seq_MW': molecular_weight,
        'seq_pI': round(analysis.isoelectric_point(), 2)
    }


# ---------------------------------------------------------------------------
# Per-chain and whole-pose metric aggregation
# ---------------------------------------------------------------------------

def calculate_chain_metrics(model, chain_id, pdb_path, dssp_chars_by_chain):
    """Calculate all metrics for a single chain.

    `dssp_chars_by_chain` is the whole-structure DSSP result from a single
    compute_dssp_chars_by_chain() call, shared across all chains so mkdssp
    only runs once per PDB file.
    """
    # Add suffix for metrics if not chain A
    suffix = '' if chain_id == 'A' else f'_{chain_id}'

    # Calculate secondary structure metrics
    helix_count, strand_count = count_ss_elements(dssp_chars_by_chain.get(chain_id, []))
    ss_metrics = {
        f'pr_helices{suffix}': helix_count,
        f'pr_strands{suffix}': strand_count,
        f'pr_total_ss{suffix}': helix_count + strand_count,
    }

    # Calculate RoG
    rog_metrics = {f'pr_RoG{suffix}': calculate_rog(model[chain_id])}

    # Calculate surface chemistry
    surface_metrics = calculate_surface_chemistry(model, chain_id)
    surface_metrics = {f'{k}{suffix}': v for k, v in surface_metrics.items()}

    # Get sequence metrics
    sequence = get_chain_sequence(pdb_path, chain_id)
    seq_metrics = calculate_seq_metrics(sequence)
    seq_metrics = {f'{k}{suffix}': v for k, v in seq_metrics.items()}

    # Combine all metrics
    metrics = {}
    metrics.update(ss_metrics)
    metrics.update(rog_metrics)
    metrics.update(surface_metrics)
    metrics.update(seq_metrics)

    return metrics


def calculate_whole_pose_metrics(model, pdb_path):
    """Calculate metrics for the entire structure (all chains combined)."""
    rog_all = calculate_rog(model)

    compute_sasa(model)
    hydrophobic_sasa = 0.0
    total_sasa = 0.0
    for chain in model.get_chains():
        for residue in chain:
            if Polypeptide.is_aa(residue, standard=True):
                res_sasa = sum(getattr(atom, 'sasa', 0.0) for atom in residue.get_atoms())
                total_sasa += res_sasa
                try:
                    aa1 = seq1(residue.get_resname()).upper()
                except Exception:
                    aa1 = ''
                if aa1 in HYDROPHOBIC_AA:
                    hydrophobic_sasa += res_sasa

    surfhphobics = round(hydrophobic_sasa / total_sasa * 100.0, 1) if total_sasa > 0 else 0.0

    return {
        'pr_RoG_total': rog_all,
        'pr_surfhphobics_total': surfhphobics,
    }


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_single_pdb(pdb_path):
    """Process a single PDB file and return all calculated metrics.

    The structure is relaxed (hydrogens added + energy-minimized) once up
    front; every metric below is computed from that relaxed structure. The
    original input PDB is only used to determine chain IDs and to validate
    that relaxation didn't change any chain's residue count.
    """
    logger.info(f"Processing PDB file: {pdb_path}")

    _, orig_model = parse_structure(pdb_path)
    chain_ids = get_chain_ids(orig_model)

    relaxed_path = _relax_structure(pdb_path)
    _, model = parse_structure(relaxed_path)
    _validate_relaxed_structure(orig_model, model, chain_ids)

    # Run mkdssp once regardless of chain count; per-chain counts are then
    # derived from the cached result instead of re-running mkdssp per chain.
    dssp_chars_by_chain = compute_dssp_chars_by_chain(model, relaxed_path)

    metrics = {'description': pdb_path.stem}

    # Handle metrics calculation based on design type
    if len(chain_ids) == 1:
        # Monomer design
        primary_chain = chain_ids[0]
        chain_metrics = calculate_chain_metrics(model, primary_chain, relaxed_path, dssp_chars_by_chain)
        metrics.update(chain_metrics)

        # Calculate SAP scores for monomer
        sap_metrics = calculate_sap_scores(model, relaxed_path, 1, primary_chain)
        metrics.update(sap_metrics)

    elif len(chain_ids) == 2:
        # Calculate binder chain (A) sequence metrics
        sequence = get_chain_sequence(relaxed_path, 'A')
        seq_metrics = calculate_seq_metrics(sequence)
        metrics.update(seq_metrics)

        # Calculate interface metrics for A-B
        interface_metrics = calculate_interface_metrics(model, relaxed_path, 'A', 'B')
        metrics.update(interface_metrics)

        # Calculate chain metrics for binder chain only
        chain_metrics = calculate_chain_metrics(model, 'A', relaxed_path, dssp_chars_by_chain)
        metrics.update(chain_metrics)

        # Calculate SAP scores for binder (alone and complex)
        sap_metrics = calculate_sap_scores(model, relaxed_path, 2, 'A')
        metrics.update(sap_metrics)

    elif len(chain_ids) >= 3:
        # Oligomer design: calculate all pairwise interfaces
        for i, c1 in enumerate(chain_ids):
            for c2 in chain_ids[i+1:]:
                interface_metrics = calculate_interface_metrics(model, relaxed_path, c1, c2)
                metrics.update(interface_metrics)

        # Calculate per-chain metrics
        all_chain_metrics = {}
        total_helices = 0
        total_strands = 0
        total_ss = 0

        for chain_id in chain_ids:
            chain_metrics = calculate_chain_metrics(model, chain_id, relaxed_path, dssp_chars_by_chain)

            # Track secondary structure aggregates
            suffix = '' if chain_id == 'A' else f'_{chain_id}'
            total_helices += chain_metrics.get(f'pr_helices{suffix}', 0)
            total_strands += chain_metrics.get(f'pr_strands{suffix}', 0)
            total_ss += chain_metrics.get(f'pr_total_ss{suffix}', 0)

            all_chain_metrics.update(chain_metrics)

        # Add whole-pose metrics
        whole_pose_metrics = calculate_whole_pose_metrics(model, relaxed_path)
        all_chain_metrics.update(whole_pose_metrics)

        # Calculate SAP for each chain (complex-context SAP computed once and reused)
        complex_sap_df = arpeggia.sap_score(
            str(relaxed_path), level='residue', probe_radius=1.4, n_points=960,
            sap_radius=5.0, chains='',
        )
        for chain_id in chain_ids:
            sap_metrics = calculate_sap_scores(
                model, relaxed_path, len(chain_ids), chain_id, complex_sap_df=complex_sap_df
            )
            # Add suffix if not chain A
            if chain_id != 'A':
                sap_metrics = {f'{k}_{chain_id}': v for k, v in sap_metrics.items()}
            all_chain_metrics.update(sap_metrics)

        # Add aggregated secondary structure metrics
        all_chain_metrics.update({
            'pr_helices_allchains': total_helices,
            'pr_strands_allchains': total_strands,
            'pr_total_ss_allchains': total_ss
        })

        metrics.update(all_chain_metrics)

    return metrics


# ---------------------------------------------------------------------------
# JSONL I/O and batch processing
# ---------------------------------------------------------------------------

def read_jsonl(file_path):
    """Read JSONL file and return list of records."""
    records = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line:  # Skip empty lines
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError as e:
                    logger.error(f"Error parsing JSONL line: {e}")
    return records


def write_jsonl(file_path, records):
    """Write records to JSONL file."""
    with open(file_path, 'w') as f:
        for record in records:
            f.write(json.dumps(record) + '\n')


def process_pdbs(pdb_dir=None, output_path=None, num_processes=4):
    """Main function to process multiple PDB files in parallel and create enriched JSONL."""

    logger.info(f"Starting processing of PDB files in directory: {pdb_dir}")

    # Get all PDB files in the directory
    pdb_paths = list(Path(pdb_dir).glob('*.pdb'))
    logger.info(f"Found {len(pdb_paths)} PDB files to process")

    if not pdb_paths:
        raise FileNotFoundError(f"No PDB files found in directory: {pdb_dir}")

    # Process PDB files in parallel
    logger.info(f"Processing PDB files with {num_processes} processes")
    process_func = partial(process_single_pdb)

    with Pool(processes=num_processes) as pool:
        results = pool.map(process_func, pdb_paths)

    # Create enriched records with fold_id and seq_id derived from filenames
    enriched_records = []
    for pdb_path, result in zip(pdb_paths, results):
        # Extract fold_id and seq_id from filename
        fold_id, seq_id = derive_ids_from_filename(pdb_path.name)

        # Create record with derived IDs and calculated metrics
        record = {
            'description': pdb_path.stem,
            'fold_id': fold_id,
            'seq_id': seq_id
        }

        # Add all calculated metrics (excluding duplicate description)
        for key, value in result.items():
            if key != 'description':
                record[key] = value

        enriched_records.append(record)

    # Save enriched records to output JSONL
    logger.info(f"Saving enriched JSONL to {output_path}")
    write_jsonl(output_path, enriched_records)
    logger.info(f"Processing completed successfully. Processed {len(enriched_records)} files")

    return enriched_records


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Process PDB files and create enriched JSONL')
    parser.add_argument('--pdb_dir', required=True, help='Directory containing PDB files')
    parser.add_argument('--output', required=True, help='Path for output enriched JSONL file')
    parser.add_argument('--num_processes', type=int, default=4, help='Number of processes to use')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Process files and save to specified output path
    enriched_records = process_pdbs(
        pdb_dir=args.pdb_dir,
        output_path=args.output,
        num_processes=args.num_processes
    )
