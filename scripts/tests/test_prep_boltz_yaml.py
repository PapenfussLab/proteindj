"""Tier A unit tests for prep_boltz_yaml.py.

PPBuilder (Bio.PDB) determines peptide-bond continuity purely from the
C(i)--N(i+1) atom distance vs `radius` (aliased here as max_break_distance),
so synthetic fixtures only need N/C atoms at controlled distances - full
realistic backbone geometry isn't required.
"""
from Bio.PDB import PDBParser

from prep_boltz_yaml import (
    get_available_chain_ids,
    split_structure_chains,
    add_seqres_to_pdb,
    generate_yaml_config,
    get_chain_ids,
    extract_sequences,
)


def _res_atom_lines(pdb_atom_line, chain, resnum, n_coord, c_coord, resname='ALA', serial_start=1):
    return [
        pdb_atom_line(serial=serial_start, name='N', chain=chain, resnum=resnum,
                       resname=resname, x=n_coord[0], y=n_coord[1], z=n_coord[2]),
        pdb_atom_line(serial=serial_start + 1, name='C', chain=chain, resnum=resnum,
                       resname=resname, x=c_coord[0], y=c_coord[1], z=c_coord[2]),
    ]


def _write_broken_chain_pdb(tmp_path, pdb_atom_line):
    """Chain A: residues 1-2 bonded (C-N dist 0.3), break to 3-4 (C-N dist 8.0), bonded again."""
    lines = []
    lines += _res_atom_lines(pdb_atom_line, 'A', 1, (0, 0, 0), (1, 0, 0), serial_start=1)
    lines += _res_atom_lines(pdb_atom_line, 'A', 2, (1.3, 0, 0), (2, 0, 0), serial_start=3)
    lines += _res_atom_lines(pdb_atom_line, 'A', 3, (10, 0, 0), (11, 0, 0), serial_start=5)
    lines += _res_atom_lines(pdb_atom_line, 'A', 4, (11.3, 0, 0), (12, 0, 0), serial_start=7)
    path = tmp_path / "input.pdb"
    path.write_text("\n".join(lines) + "\n")
    return path


def _parse(path):
    parser = PDBParser(QUIET=True)
    return parser.get_structure('s', path)


def _write_three_chain_pdb(tmp_path, pdb_atom_line):
    # Two-residue chains A, B, C, each a standalone (unbroken) peptide.
    lines = []
    lines += _res_atom_lines(pdb_atom_line, 'A', 1, (0, 0, 0), (1, 0, 0), serial_start=1)
    lines += _res_atom_lines(pdb_atom_line, 'A', 2, (1.3, 0, 0), (2, 0, 0), serial_start=3)
    lines += _res_atom_lines(pdb_atom_line, 'B', 1, (20, 0, 0), (21, 0, 0), serial_start=5)
    lines += _res_atom_lines(pdb_atom_line, 'B', 2, (21.3, 0, 0), (22, 0, 0), serial_start=7)
    lines += _res_atom_lines(pdb_atom_line, 'C', 1, (40, 0, 0), (41, 0, 0), serial_start=9)
    lines += _res_atom_lines(pdb_atom_line, 'C', 2, (41.3, 0, 0), (42, 0, 0), serial_start=11)
    path = tmp_path / "input.pdb"
    path.write_text("\n".join(lines) + "\n")
    return path


# ---------------------------------------------------------------------------
# get_available_chain_ids
# ---------------------------------------------------------------------------

def test_get_available_chain_ids_excludes_used():
    result = get_available_chain_ids({'A', 'B'})
    assert result[:3] == ['C', 'D', 'E']
    assert 'A' not in result
    assert 'B' not in result


# ---------------------------------------------------------------------------
# split_structure_chains
# ---------------------------------------------------------------------------

def test_split_structure_chains_splits_on_break(tmp_path, pdb_atom_line):
    structure = _parse(_write_broken_chain_pdb(tmp_path, pdb_atom_line))
    new_structure, chain_id_map = split_structure_chains(structure, max_break_distance=3.0)

    assert chain_id_map['A'] == ['A', 'B']
    new_chain_ids = sorted(chain.id for chain in new_structure[0])
    assert new_chain_ids == ['A', 'B']

    chain_a_residues = [r.id[1] for r in new_structure[0]['A']]
    chain_b_residues = [r.id[1] for r in new_structure[0]['B']]
    assert chain_a_residues == [1, 2]
    assert chain_b_residues == [3, 4]


def test_split_structure_chains_no_break_keeps_single_chain(tmp_path, pdb_atom_line):
    lines = []
    lines += _res_atom_lines(pdb_atom_line, 'A', 1, (0, 0, 0), (1, 0, 0), serial_start=1)
    lines += _res_atom_lines(pdb_atom_line, 'A', 2, (1.3, 0, 0), (2, 0, 0), serial_start=3)
    path = tmp_path / "input.pdb"
    path.write_text("\n".join(lines) + "\n")
    structure = _parse(path)

    new_structure, chain_id_map = split_structure_chains(structure, max_break_distance=3.0)
    assert chain_id_map == {'A': ['A']}
    assert [chain.id for chain in new_structure[0]] == ['A']


# ---------------------------------------------------------------------------
# add_seqres_to_pdb
# ---------------------------------------------------------------------------

def test_add_seqres_to_pdb_writes_seqres_per_split_chain(tmp_path, pdb_atom_line):
    input_path = _write_broken_chain_pdb(tmp_path, pdb_atom_line)
    output_path = tmp_path / "with_seqres.pdb"
    add_seqres_to_pdb(input_path, output_path, max_break_distance=3.0)

    text = output_path.read_text()
    seqres_lines = [l for l in text.splitlines() if l.startswith('SEQRES')]
    assert any(l.startswith('SEQRES   1 A    2  ALA ALA') for l in seqres_lines)
    assert any(l.startswith('SEQRES   1 B    2  ALA ALA') for l in seqres_lines)
    assert text.rstrip().splitlines()[-1] == 'END'


# ---------------------------------------------------------------------------
# get_chain_ids
# ---------------------------------------------------------------------------

def test_get_chain_ids_extracts_ids():
    sequences = [{'id': 'A', 'sequence': 'AG', 'msa': 'empty'},
                 {'id': 'C', 'sequence': 'GG', 'msa': 'empty'}]
    assert get_chain_ids(sequences) == ['A', 'C']


# ---------------------------------------------------------------------------
# generate_yaml_config
# ---------------------------------------------------------------------------

def test_generate_yaml_config_without_template():
    sequences = [{'id': 'A', 'sequence': 'AG', 'msa': 'empty'}]
    config = generate_yaml_config(sequences)
    assert config == {'sequences': [{'protein': {'id': 'A', 'sequence': 'AG', 'msa': 'empty'}}]}
    assert 'templates' not in config


def test_generate_yaml_config_with_template():
    sequences = [{'id': 'A', 'sequence': 'AG', 'msa': 'empty'}]
    config = generate_yaml_config(
        sequences, use_template=True, pdb_filename='fold_1.pdb',
        template_chains=['B'], template_force=True, template_threshold=2.5,
    )
    assert config['templates'] == [
        {'pdb': 'templates/fold_1.pdb', 'chain_id': 'B', 'force': True, 'threshold': 2.5}
    ]


def test_generate_yaml_config_template_requires_all_fields():
    sequences = [{'id': 'A', 'sequence': 'AG', 'msa': 'empty'}]
    # use_template True but no pdb_filename/template_chains -> no templates block
    config = generate_yaml_config(sequences, use_template=True)
    assert 'templates' not in config


# ---------------------------------------------------------------------------
# extract_sequences - MSA file matching
# ---------------------------------------------------------------------------

def _msa_by_id(sequences):
    return {seq['id']: seq['msa'] for seq in sequences}


def test_extract_sequences_single_msa_broadcasts_to_all_target_chains(tmp_path, pdb_atom_line):
    pdb_path = _write_three_chain_pdb(tmp_path, pdb_atom_line)
    sequences, _ = extract_sequences(pdb_path, msa_files=['target.a3m'], msa_exclude_chain='A')
    assert _msa_by_id(sequences) == {'A': 'empty', 'B': 'target.a3m', 'C': 'target.a3m'}


def test_extract_sequences_multiple_msa_matched_positionally(tmp_path, pdb_atom_line):
    pdb_path = _write_three_chain_pdb(tmp_path, pdb_atom_line)
    sequences, _ = extract_sequences(
        pdb_path, msa_files=['b_target.a3m', 'c_target.a3m'], msa_exclude_chain='A'
    )
    assert _msa_by_id(sequences) == {'A': 'empty', 'B': 'b_target.a3m', 'C': 'c_target.a3m'}


def test_extract_sequences_mismatched_msa_count_matches_best_effort(tmp_path, pdb_atom_line, capsys):
    pdb_path = _write_three_chain_pdb(tmp_path, pdb_atom_line)
    # Three files for two eligible target chains (B, C) -> positional match, extra file unused
    sequences, _ = extract_sequences(
        pdb_path, msa_files=['b_target.a3m', 'extra.a3m', 'unused.a3m'], msa_exclude_chain='A'
    )
    assert _msa_by_id(sequences) == {'A': 'empty', 'B': 'b_target.a3m', 'C': 'extra.a3m'}
    assert 'WARNING' in capsys.readouterr().out


def test_extract_sequences_no_msa_files_all_empty(tmp_path, pdb_atom_line):
    pdb_path = _write_three_chain_pdb(tmp_path, pdb_atom_line)
    sequences, _ = extract_sequences(pdb_path, msa_files=None, msa_exclude_chain='A')
    assert _msa_by_id(sequences) == {'A': 'empty', 'B': 'empty', 'C': 'empty'}
