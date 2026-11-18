"""Tier A unit tests for align_af2.py.

align_af2.py runs `setup_logging()` at import time, which creates an
"alignment.log" FileHandler in the current working directory as a side
effect. We import it with cwd temporarily pointed at a throwaway directory
so this doesn't litter the repo (this only matters on the first import in
the whole test session, since module imports are cached).
"""
import os
import tempfile

import pytest
from Bio.PDB import Atom, PDBParser

_prev_cwd = os.getcwd()
os.chdir(tempfile.mkdtemp())
try:
    from align_af2 import calculate_rmsd, get_chain_atoms, get_target_atoms
finally:
    os.chdir(_prev_cwd)


def _atom(name, coord, serial=1):
    return Atom.Atom(name, coord, 1.0, 1.0, ' ', name, serial, name[0])


def _write_pdb(tmp_path, lines, name="input.pdb"):
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return path


def _parse(path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure('s', path)


# ---------------------------------------------------------------------------
# calculate_rmsd
# ---------------------------------------------------------------------------

def test_calculate_rmsd_mismatched_lengths_raises():
    atoms1 = [_atom('CA', (0, 0, 0))]
    atoms2 = [_atom('CA', (0, 0, 0)), _atom('CA', (1, 1, 1), serial=2)]
    with pytest.raises(ValueError, match="different lengths: 1 vs 2"):
        calculate_rmsd(atoms1, atoms2)


def test_calculate_rmsd_identical_shapes_is_zero():
    atoms1 = [_atom('CA', (0, 0, 0)), _atom('CA', (1, 0, 0), 2), _atom('CA', (0, 1, 0), 3)]
    atoms2 = [_atom('CA', (5, 5, 5)), _atom('CA', (6, 5, 5), 2), _atom('CA', (5, 6, 5), 3)]
    assert calculate_rmsd(atoms1, atoms2) == pytest.approx(0.0, abs=1e-6)


def test_calculate_rmsd_known_value():
    atoms1 = [_atom('CA', (0, 0, 0)), _atom('CA', (1, 0, 0), 2), _atom('CA', (0, 1, 0), 3)]
    atoms2 = [_atom('CA', (0, 0, 0)), _atom('CA', (1, 0, 0), 2), _atom('CA', (0, 1, 1), 3)]
    assert calculate_rmsd(atoms1, atoms2) == pytest.approx(0.1849, abs=1e-3)


# ---------------------------------------------------------------------------
# get_chain_atoms
# ---------------------------------------------------------------------------

def test_get_chain_atoms_returns_atoms_of_type(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, name='N', chain='B', resnum=1),
        pdb_atom_line(serial=2, name='CA', chain='B', resnum=1),
        pdb_atom_line(serial=3, name='CA', chain='B', resnum=2),
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    atoms = get_chain_atoms(structure, 'B')
    assert [a.get_serial_number() for a in atoms] == [2, 3]


def test_get_chain_atoms_missing_chain_raises(tmp_path, pdb_atom_line):
    lines = [pdb_atom_line(serial=1, chain='B', resnum=1)]
    structure = _parse(_write_pdb(tmp_path, lines))
    with pytest.raises(ValueError, match="Chain Z not found"):
        get_chain_atoms(structure, 'Z')


def test_get_chain_atoms_no_matching_atom_type_raises(tmp_path, pdb_atom_line):
    lines = [pdb_atom_line(serial=1, name='N', chain='B', resnum=1)]
    structure = _parse(_write_pdb(tmp_path, lines))
    with pytest.raises(ValueError, match="No CA atoms found in chain B"):
        get_chain_atoms(structure, 'B', atom_type='CA')


# ---------------------------------------------------------------------------
# get_target_atoms
# ---------------------------------------------------------------------------

def test_get_target_atoms_excludes_binder_chain(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, name='CA', chain='A', resnum=1),
        pdb_atom_line(serial=2, name='CA', chain='B', resnum=1),
        pdb_atom_line(serial=3, name='CA', chain='C', resnum=1),
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    atoms = get_target_atoms(structure, binder_chain='A')
    assert [a.get_serial_number() for a in atoms] == [2, 3]


def test_get_target_atoms_no_target_chains_raises(tmp_path, pdb_atom_line):
    lines = [pdb_atom_line(serial=1, name='CA', chain='A', resnum=1)]
    structure = _parse(_write_pdb(tmp_path, lines))
    with pytest.raises(ValueError, match="No CA atoms found in target chains"):
        get_target_atoms(structure, binder_chain='A')


def test_get_target_atoms_ignores_hetatm_calcium_named_ca(tmp_path, pdb_atom_line):
    lines = [
        pdb_atom_line(serial=1, name='CA', chain='A', resnum=1),
        pdb_atom_line(serial=2, name='CA', chain='B', resnum=1),
        "HETATM    3 CA    CA B 101      10.000  10.000  10.000  1.00  0.00          CA",
    ]
    structure = _parse(_write_pdb(tmp_path, lines))
    atoms = get_target_atoms(structure, binder_chain='A')
    assert [a.get_serial_number() for a in atoms] == [2]
