"""Tier A unit tests for analyse_boltzgen.py."""
import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, PDBParser, Residue, Structure
from Bio.PDB.mmcifio import MMCIFIO

from analyse_boltzgen import (
    build_metadata,
    load_design_mask,
    collect_design_pairs,
    get_protein_chains_in_order,
    cif_to_relabelled_pdb,
)


def _make_residue(resnum, resname='ALA', het=' '):
    res = Residue.Residue((het, resnum, ' '), resname, '')
    atom = Atom.Atom('CA', (float(resnum), 0.0, 0.0), 1.0, 1.0, ' ', 'CA', resnum, 'C')
    res.add(atom)
    return res


def _write_cif(cif_path, chains_spec):
    """chains_spec: list of (chain_id, [resnums]) in desired file order."""
    structure = Structure.Structure('s')
    model = Model.Model(0)
    structure.add(model)
    for chain_id, resnums in chains_spec:
        chain = Chain.Chain(chain_id)
        for resnum in resnums:
            chain.add(_make_residue(resnum))
        model.add(chain)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))


class TestBuildMetadata:
    def test_inverts_design_mask_to_inpaint_seq(self):
        design_mask = np.array([True, False, True])
        metadata = build_metadata(fold_id=2, design_mask=design_mask, design_mode='boltzgen_denovo')
        assert metadata == {
            'fold_id': 2,
            'bg_design_mode': 'boltzgen_denovo',
            'bg_inpaint_seq': [False, True, False],
        }

    def test_all_designed_inverts_to_all_fixed_false(self):
        design_mask = np.array([True, True, True])
        metadata = build_metadata(fold_id=0, design_mask=design_mask, design_mode='boltzgen_motifscaff')
        assert metadata['bg_inpaint_seq'] == [False, False, False]

    def test_all_fixed_inverts_to_all_true(self):
        design_mask = np.array([False, False])
        metadata = build_metadata(fold_id=1, design_mask=design_mask, design_mode='boltzgen_denovo')
        assert metadata['bg_inpaint_seq'] == [True, True]


class TestLoadDesignMask:
    def test_loads_boolean_array(self, tmp_path):
        npz_path = tmp_path / "fold_0.npz"
        np.savez(npz_path, design_mask=np.array([1, 0, 1]))
        mask = load_design_mask(str(npz_path))
        assert mask.dtype == bool
        assert list(mask) == [True, False, True]

    def test_missing_design_mask_key_raises_keyerror(self, tmp_path):
        npz_path = tmp_path / "fold_0.npz"
        np.savez(npz_path, some_other_array=np.array([1, 2, 3]))
        with pytest.raises(KeyError):
            load_design_mask(str(npz_path))


class TestCollectDesignPairs:
    def test_pairs_cif_and_npz_files(self, tmp_path):
        (tmp_path / "fold_0.cif").write_text("dummy cif")
        np.savez(tmp_path / "fold_0.npz", design_mask=np.array([True]))
        (tmp_path / "fold_1.cif").write_text("dummy cif")
        np.savez(tmp_path / "fold_1.npz", design_mask=np.array([True]))

        pairs = collect_design_pairs(tmp_path)
        assert len(pairs) == 2
        names = sorted(p[1].stem for p in pairs)
        assert names == ["fold_0", "fold_1"]

    def test_orphan_npz_without_cif_is_skipped(self, tmp_path, capsys):
        np.savez(tmp_path / "fold_0.npz", design_mask=np.array([True]))
        # No matching fold_0.cif
        pairs = collect_design_pairs(tmp_path)
        assert pairs == []

    def test_native_cif_without_npz_is_excluded(self, tmp_path):
        # BoltzGen writes <name>_native.cif with no matching .npz - must not appear in pairs
        (tmp_path / "fold_0_native.cif").write_text("dummy cif")
        (tmp_path / "fold_0.cif").write_text("dummy cif")
        np.savez(tmp_path / "fold_0.npz", design_mask=np.array([True]))

        pairs = collect_design_pairs(tmp_path)
        assert len(pairs) == 1
        assert pairs[0][0].stem == "fold_0"


class TestGetProteinChainsInOrder:
    def test_filters_chains_with_only_hetatm(self, tmp_path):
        from Bio.PDB import Structure, Model, Chain, Residue, Atom

        structure = Structure.Structure('s')
        model = Model.Model(0)
        structure.add(model)

        chain_a = Chain.Chain('A')
        res = Residue.Residue((' ', 1, ' '), 'ALA', '')
        atom = Atom.Atom('CA', (0.0, 0.0, 0.0), 1.0, 1.0, ' ', 'CA', 1, 'C')
        res.add(atom)
        chain_a.add(res)
        model.add(chain_a)

        chain_b = Chain.Chain('B')
        het_res = Residue.Residue(('H_HOH', 1, ' '), 'HOH', '')
        het_atom = Atom.Atom('O', (0.0, 0.0, 0.0), 1.0, 1.0, ' ', 'O', 1, 'O')
        het_res.add(het_atom)
        chain_b.add(het_res)
        model.add(chain_b)

        result = get_protein_chains_in_order(model)
        assert [c.id for c in result] == ['A']


class TestCifToRelabelledPdb:
    def test_relabels_chains_in_file_order_and_renumbers_sequentially(self, tmp_path):
        cif_path = tmp_path / "fold_0.cif"
        _write_cif(cif_path, [('A', [1, 2, 3]), ('B', [10, 11])])

        pdb_path = tmp_path / "fold_0.pdb"
        chain_lengths = cif_to_relabelled_pdb(cif_path, pdb_path)

        assert chain_lengths == [3, 2]

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('p', str(pdb_path))
        model = structure[0]
        chain_ids = [c.id for c in model]
        assert chain_ids == ['A', 'B']
        # Residues renumbered continuously: chain A = 1..3, chain B = 4..5
        assert [r.id[1] for r in model['A']] == [1, 2, 3]
        assert [r.id[1] for r in model['B']] == [4, 5]

    def test_chain_a_is_placed_first_regardless_of_file_order(self, tmp_path):
        cif_path = tmp_path / "fold_0.cif"
        # Chain B's ATOM records precede chain A's - BoltzGen preserves this file order
        # but does not reorder by role, so 'A' must still map to new chain 'A'.
        _write_cif(cif_path, [('B', [10, 11]), ('A', [1, 2, 3])])

        pdb_path = tmp_path / "fold_0.pdb"
        chain_lengths = cif_to_relabelled_pdb(cif_path, pdb_path)

        assert chain_lengths == [3, 2]

        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('p', str(pdb_path))
        model = structure[0]
        assert [c.id for c in model] == ['A', 'B']
        # Original chain A's 3 residues must land in new chain A, not chain B
        assert [r.id[1] for r in model['A']] == [1, 2, 3]
        assert [r.id[1] for r in model['B']] == [4, 5]

    def test_missing_chain_a_raises_value_error(self, tmp_path):
        cif_path = tmp_path / "fold_0.cif"
        _write_cif(cif_path, [('X', [1, 2, 3]), ('Y', [10, 11])])

        pdb_path = tmp_path / "fold_0.pdb"
        with pytest.raises(ValueError):
            cif_to_relabelled_pdb(cif_path, pdb_path)

    def test_hetatm_residues_are_dropped(self, tmp_path):
        cif_path = tmp_path / "fold_0.cif"
        structure = Structure.Structure('s')
        model = Model.Model(0)
        structure.add(model)
        chain = Chain.Chain('A')
        chain.add(_make_residue(1))
        chain.add(_make_residue(2, resname='HOH', het='H_HOH'))
        model.add(chain)
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(cif_path))

        pdb_path = tmp_path / "fold_0.pdb"
        chain_lengths = cif_to_relabelled_pdb(cif_path, pdb_path)
        assert chain_lengths == [1]

    def test_no_protein_chains_raises_value_error(self, tmp_path):
        cif_path = tmp_path / "fold_0.cif"
        structure = Structure.Structure('s')
        model = Model.Model(0)
        structure.add(model)
        chain = Chain.Chain('A')
        chain.add(_make_residue(1, resname='HOH', het='H_HOH'))
        model.add(chain)
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(str(cif_path))

        pdb_path = tmp_path / "fold_0.pdb"
        with pytest.raises(ValueError):
            cif_to_relabelled_pdb(cif_path, pdb_path)
