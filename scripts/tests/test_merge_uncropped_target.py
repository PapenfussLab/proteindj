"""Tier A unit tests for merge_uncropped_target.py.

find_best_alignment_region/align_structures_by_sequence are heuristic
(pairwise alignment + superposition) and are intentionally out of scope for
this deterministic Tier A suite.
"""
from Bio.PDB import PDBParser

from merge_uncropped_target import (
    get_sequence_and_atoms,
    get_target_chain_ids,
    merge_chains_to_single_chain,
)


def _write_pdb(tmp_path, lines, name="input.pdb"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def _parse(path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure('s', path)


# ---------------------------------------------------------------------------
# get_sequence_and_atoms
# ---------------------------------------------------------------------------

def test_get_sequence_and_atoms_maps_known_residues(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, chain='A', resnum=1, resname='ALA'),
        pdb_atom_line(serial=2, chain='A', resnum=2, resname='GLY'),
        pdb_atom_line(serial=3, chain='A', resnum=3, resname='CYS'),
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    seq, _ = get_sequence_and_atoms(structure, ['A'])
    assert seq == "AGC"


def test_get_sequence_and_atoms_unknown_residue_is_x(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, chain='A', resnum=1, resname='ALA'),
        pdb_atom_line(serial=2, chain='A', resnum=2, resname='LIG'),
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    seq, _ = get_sequence_and_atoms(structure, ['A'])
    assert seq == "AX"


def test_get_sequence_and_atoms_skips_hetatm(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, chain='A', resnum=1, resname='ALA'),
        "HETATM    2  O   HOH A   2      10.000  10.000  10.000  1.00  0.00           O",
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    seq, _ = get_sequence_and_atoms(structure, ['A'])
    assert seq == "A"


def test_get_sequence_and_atoms_missing_chain_returns_empty(tmp_path, pdb_atom_line):
    lines = [pdb_atom_line(serial=1, chain='A', resnum=1, resname='ALA')]
    structure = _parse(_write_pdb(tmp_path, lines))
    seq, atoms = get_sequence_and_atoms(structure, ['Z'])
    assert seq == ""
    assert atoms == []


def test_get_sequence_and_atoms_concatenates_multiple_chains(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, chain='B', resnum=1, resname='ALA'),
        pdb_atom_line(serial=2, chain='C', resnum=1, resname='GLY'),
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    seq, atoms = get_sequence_and_atoms(structure, ['B', 'C'])
    assert seq == "AG"
    assert [a.get_serial_number() for a in atoms] == [1, 2]


def test_get_sequence_and_atoms_returns_one_atom_per_residue(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, name='N', chain='A', resnum=1),
        pdb_atom_line(serial=2, name='CA', chain='A', resnum=1),
        pdb_atom_line(serial=3, name='C', chain='A', resnum=1),
        pdb_atom_line(serial=4, name='N', chain='A', resnum=2),
        pdb_atom_line(serial=5, name='CA', chain='A', resnum=2),
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    _, atoms = get_sequence_and_atoms(structure, ['A'], atom_type='CA')
    assert [a.get_serial_number() for a in atoms] == [2, 5]


# ---------------------------------------------------------------------------
# get_target_chain_ids
# ---------------------------------------------------------------------------

def test_get_target_chain_ids_excludes_binder_chain(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, chain='A', resnum=1, resname='ALA'),
        pdb_atom_line(serial=2, chain='B', resnum=1, resname='GLY'),
        pdb_atom_line(serial=3, chain='C', resnum=1, resname='CYS'),
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    assert get_target_chain_ids(structure, binder_chain='A') == ['B', 'C']


# ---------------------------------------------------------------------------
# merge_chains_to_single_chain
# ---------------------------------------------------------------------------

def test_merge_chains_to_single_chain_renumbers_continuously(tmp_path, pdb_atom_line):
    lines = (
        [pdb_atom_line(serial=1, chain='B', resnum=5, resname='ALA')]
        + [pdb_atom_line(serial=2, chain='A', resnum=1, resname='GLY')]
        + [pdb_atom_line(serial=3, chain='A', resnum=2, resname='CYS')]
    )
    input_path = _write_pdb(tmp_path, lines)
    output_path = tmp_path / "merged.pdb"
    merge_chains_to_single_chain(input_path, output_path, target_chain='X')

    out_lines = output_path.read_text().splitlines()
    atom_lines = [l for l in out_lines if l.startswith('ATOM')]
    # Chains merged in alphabetical order (A before B), renumbered 1..N, all as target_chain
    assert [l[17:20] for l in atom_lines] == ['GLY', 'CYS', 'ALA']
    assert [l[21] for l in atom_lines] == ['X', 'X', 'X']
    assert [l[22:26].strip() for l in atom_lines] == ['1', '2', '3']
    assert out_lines[-1] == "END"


def test_merge_chains_to_single_chain_skips_hetatm(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, chain='A', resnum=1, resname='ALA'),
        "HETATM    2  O   HOH A   2      10.000  10.000  10.000  1.00  0.00           O",
    ]
    input_path = _write_pdb(tmp_path, lines)
    output_path = tmp_path / "merged.pdb"
    merge_chains_to_single_chain(input_path, output_path, target_chain='B')

    atom_lines = [l for l in output_path.read_text().splitlines() if l.startswith('ATOM')]
    assert len(atom_lines) == 1
    assert atom_lines[0][17:20] == 'ALA'
