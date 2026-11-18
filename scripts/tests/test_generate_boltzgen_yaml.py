"""Tier A unit tests for generate_boltzgen_yaml.py."""
from types import SimpleNamespace

import pytest

from generate_boltzgen_yaml import (
    get_protein_chains,
    get_chain_length,
    get_residue_rank_map,
    _resolve_rank,
    clean_pdb_for_boltzgen,
    parse_design_length,
    parse_residue_ranges,
    parse_flexible_spec,
    _res_index_to_ranks,
    parse_architecture_spec,
    build_structure_groups,
    build_binding_types,
    apply_binding_types,
    apply_structure_groups,
    build_denovo_spec,
    build_motifscaff_spec,
)


@pytest.fixture
def two_chain_pdb(tmp_path, pdb_atom_line):
    """Chain A: author resnums 10,11,12,15,16,17 (gap) -> ranks 1..6.
    Chain B: author resnums 100,101,102 -> ranks 1..3.
    """
    lines = [pdb_atom_line(serial=i + 1, chain='A', resnum=r)
             for i, r in enumerate([10, 11, 12, 15, 16, 17])]
    lines += [pdb_atom_line(serial=i + 7, chain='B', resnum=r, x=20.0)
              for i, r in enumerate([100, 101, 102])]
    path = tmp_path / "input.pdb"
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def rank_map(two_chain_pdb):
    return get_residue_rank_map(two_chain_pdb)


# ---------------------------------------------------------------------------
# get_protein_chains / get_chain_length / get_residue_rank_map / _resolve_rank
# ---------------------------------------------------------------------------

def test_get_protein_chains(two_chain_pdb):
    assert get_protein_chains(two_chain_pdb) == ['A', 'B']


def test_get_chain_length(two_chain_pdb):
    assert get_chain_length(two_chain_pdb, 'A') == 6
    assert get_chain_length(two_chain_pdb, 'B') == 3
    assert get_chain_length(two_chain_pdb, 'Z') == 0


def test_get_residue_rank_map(two_chain_pdb):
    rm = get_residue_rank_map(two_chain_pdb)
    assert rm['A'] == {10: 1, 11: 2, 12: 3, 15: 4, 16: 5, 17: 6}
    assert rm['B'] == {100: 1, 101: 2, 102: 3}


def test_resolve_rank_success(rank_map):
    assert _resolve_rank('A', 15, rank_map) == 4


def test_resolve_rank_unknown_chain_raises(rank_map):
    with pytest.raises(ValueError, match="Chain 'Z' not found"):
        _resolve_rank('Z', 1, rank_map)


def test_resolve_rank_unknown_residue_raises(rank_map):
    with pytest.raises(ValueError, match="Residue A99 not found"):
        _resolve_rank('A', 99, rank_map)


# ---------------------------------------------------------------------------
# clean_pdb_for_boltzgen
# ---------------------------------------------------------------------------

def test_clean_pdb_for_boltzgen_strips_header_metadata(tmp_path, pdb_atom_line):
    atom_line = pdb_atom_line(chain='A', resnum=1)
    input_path = tmp_path / "in.pdb"
    input_path.write_text(
        "SEQRES   1 A    1  ALA\n"
        "COMPND    MOL_ID: 1;\n"
        f"{atom_line}\n"
    )
    output_path = tmp_path / "out.pdb"
    clean_pdb_for_boltzgen(input_path, output_path)
    out_lines = output_path.read_text().splitlines()
    assert out_lines == [atom_line]


def test_clean_pdb_for_boltzgen_no_reorder_by_default(tmp_path, pdb_atom_line):
    """Chain B physically precedes chain A - without reordering, order is preserved."""
    b_line = pdb_atom_line(chain='B', resnum=1)
    a_line = pdb_atom_line(chain='A', resnum=1)
    input_path = tmp_path / "in.pdb"
    input_path.write_text(f"{b_line}\n{a_line}\n")
    output_path = tmp_path / "out.pdb"
    clean_pdb_for_boltzgen(input_path, output_path)
    assert output_path.read_text().splitlines() == [b_line, a_line]


def test_clean_pdb_for_boltzgen_reorders_chain_a_first(tmp_path, pdb_atom_line):
    """Chain B physically precedes chain A in the input - reorder_chain_a_first must
    move chain A's records to the front, keeping other chains in their relative order."""
    b_line = pdb_atom_line(serial=1, chain='B', resnum=1)
    b_ter = "TER     796      GLY B   1"
    a_line = pdb_atom_line(serial=2, chain='A', resnum=1)
    a_ter = "TER     797      GLY A   1"
    input_path = tmp_path / "in.pdb"
    input_path.write_text(f"{b_line}\n{b_ter}\n{a_line}\n{a_ter}\n")
    output_path = tmp_path / "out.pdb"
    clean_pdb_for_boltzgen(input_path, output_path, reorder_chain_a_first=True)
    assert output_path.read_text().splitlines() == [a_line, a_ter, b_line, b_ter]


def test_clean_pdb_for_boltzgen_reorder_without_chain_a_raises(tmp_path, pdb_atom_line):
    input_path = tmp_path / "in.pdb"
    input_path.write_text(f"{pdb_atom_line(chain='B', resnum=1)}\n")
    with pytest.raises(ValueError, match="No chain 'A'"):
        clean_pdb_for_boltzgen(input_path, tmp_path / "out.pdb", reorder_chain_a_first=True)


# ---------------------------------------------------------------------------
# parse_design_length
# ---------------------------------------------------------------------------

def test_parse_design_length_single():
    assert parse_design_length('80') == '80'


def test_parse_design_length_range():
    assert parse_design_length('60-100') == '60..100'


# ---------------------------------------------------------------------------
# parse_residue_ranges / parse_flexible_spec
# ---------------------------------------------------------------------------

def test_parse_residue_ranges_mixed_tokens(rank_map):
    result = parse_residue_ranges('A10,A15-17,B100', rank_map)
    assert result == {'A': '1,4,5,6', 'B': '1'}


def test_parse_residue_ranges_whole_chain_token(rank_map):
    result = parse_residue_ranges('A', rank_map)
    assert result == {'A': 'all'}


def test_parse_residue_ranges_unresolvable_residue_raises(rank_map):
    with pytest.raises(ValueError, match="Residue A99 not found"):
        parse_residue_ranges('A99', rank_map)


def test_parse_residue_ranges_bad_token_raises(rank_map):
    with pytest.raises(ValueError, match="Could not parse residue token"):
        parse_residue_ranges('!!!', rank_map)


def test_parse_flexible_spec_whole_chain_is_none(rank_map):
    result = parse_flexible_spec('A10-12,B', rank_map)
    assert result == {'A': '1,2,3', 'B': None}


# ---------------------------------------------------------------------------
# _res_index_to_ranks
# ---------------------------------------------------------------------------

def test_res_index_to_ranks_all():
    assert _res_index_to_ranks('all', 5) == {1, 2, 3, 4, 5}


def test_res_index_to_ranks_none():
    assert _res_index_to_ranks(None, 3) == {1, 2, 3}


def test_res_index_to_ranks_explicit_list():
    assert _res_index_to_ranks('1,4,5', 6) == {1, 4, 5}


# ---------------------------------------------------------------------------
# parse_architecture_spec
# ---------------------------------------------------------------------------

def test_parse_architecture_spec_keep_and_insert(rank_map):
    exclude_ranks, insertions, kept_ranks = parse_architecture_spec(
        'A10-12,5,A15-17,10', rank_map, chain_a_length=6
    )
    assert exclude_ranks == []
    assert insertions == [(4, '5'), (7, '10')]
    assert kept_ranks == [1, 2, 3, 4, 5, 6]


def test_parse_architecture_spec_leading_gap_excluded(rank_map):
    exclude_ranks, insertions, kept_ranks = parse_architecture_spec(
        'A15-17', rank_map, chain_a_length=6
    )
    assert exclude_ranks == [(1, 3)]
    assert insertions == []
    assert kept_ranks == [4, 5, 6]


def test_parse_architecture_spec_trailing_gap_excluded(rank_map):
    exclude_ranks, insertions, kept_ranks = parse_architecture_spec(
        'A10-12', rank_map, chain_a_length=6
    )
    assert exclude_ranks == [(4, 6)]
    assert kept_ranks == [1, 2, 3]


def test_parse_architecture_spec_insert_range_token(rank_map):
    _, insertions, _ = parse_architecture_spec('A10-12,5-10', rank_map, chain_a_length=6)
    assert insertions == [(4, '5..10')]


def test_parse_architecture_spec_out_of_order_raises(rank_map):
    with pytest.raises(ValueError, match="strictly ascending"):
        parse_architecture_spec('A15-17,A10-12', rank_map, chain_a_length=6)


def test_parse_architecture_spec_bad_token_raises(rank_map):
    with pytest.raises(ValueError, match="Could not parse motifscaff_spec token"):
        parse_architecture_spec('!!!', rank_map, chain_a_length=6)


# ---------------------------------------------------------------------------
# build_structure_groups / build_binding_types
# ---------------------------------------------------------------------------

def test_build_structure_groups():
    groups = build_structure_groups({'A': None, 'B': '1,2'})
    assert groups == [
        {'group': {'visibility': 1, 'id': 'A'}},
        {'group': {'visibility': 0, 'id': 'A'}},
        {'group': {'visibility': 1, 'id': 'B'}},
        {'group': {'visibility': 0, 'id': 'B', 'res_index': '1,2'}},
    ]


def test_build_binding_types():
    result = build_binding_types({'A': '1,2'}, {'B': '3'})
    assert result == [
        {'chain': {'id': 'A', 'binding': '1,2'}},
        {'chain': {'id': 'B', 'not_binding': '3'}},
    ]


# ---------------------------------------------------------------------------
# apply_binding_types / apply_structure_groups
# ---------------------------------------------------------------------------

def test_apply_binding_types_noop_when_unset(rank_map):
    entity = {}
    apply_binding_types(entity, '', '', ['B'], rank_map)
    assert entity == {}


def test_apply_binding_types_sets_block(rank_map):
    entity = {}
    apply_binding_types(entity, 'B100', '', ['B'], rank_map)
    assert entity['binding_types'] == [{'chain': {'id': 'B', 'binding': '1'}}]


def test_apply_binding_types_invalid_chain_raises(rank_map):
    entity = {}
    with pytest.raises(ValueError, match="reference chain"):
        apply_binding_types(entity, 'A10', '', ['B'], rank_map)


def test_apply_structure_groups_noop_when_unset(rank_map):
    entity = {}
    apply_structure_groups(entity, '', ['B'], rank_map)
    assert entity == {}


def test_apply_structure_groups_invalid_chain_raises(rank_map):
    entity = {}
    with pytest.raises(ValueError, match="not in target"):
        apply_structure_groups(entity, 'A10', ['B'], rank_map)


# ---------------------------------------------------------------------------
# build_denovo_spec
# ---------------------------------------------------------------------------

def _args(**overrides):
    defaults = dict(
        input_pdb='', design_length='80', hotspot_residues='', bg_not_binding_residues='',
        flexible_residues='', motifscaff_spec='', motifscaff_inpaint_seq='',
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def test_build_denovo_spec_monomer_no_target():
    args = _args()
    spec = build_denovo_spec(args, target_chains=[], rank_map={})
    assert spec == {'entities': [{'protein': {'id': 'A', 'sequence': '80'}}]}


def test_build_denovo_spec_monomer_with_hotspots_raises():
    args = _args(hotspot_residues='A10')
    with pytest.raises(ValueError, match="require a target"):
        build_denovo_spec(args, target_chains=[], rank_map={})


def test_build_denovo_spec_with_target_and_hotspots(two_chain_pdb, rank_map):
    args = _args(input_pdb=str(two_chain_pdb), hotspot_residues='B100')
    spec = build_denovo_spec(args, target_chains=['B'], rank_map=rank_map)
    entities = spec['entities']
    assert entities[0] == {'protein': {'id': 'A', 'sequence': '80'}}
    file_entity = entities[1]['file']
    assert file_entity['include'] == [{'chain': {'id': 'B'}}]
    assert file_entity['binding_types'] == [{'chain': {'id': 'B', 'binding': '1'}}]


# ---------------------------------------------------------------------------
# build_motifscaff_spec
# ---------------------------------------------------------------------------

def test_build_motifscaff_spec_requires_chain_a(rank_map):
    args = _args(input_pdb='dummy', motifscaff_spec='5')
    with pytest.raises(ValueError, match="requires input_pdb to contain a chain A"):
        build_motifscaff_spec(args, all_chains=['B'], rank_map=rank_map)


def test_build_motifscaff_spec_requires_at_least_one_option(two_chain_pdb, rank_map):
    args = _args(input_pdb=str(two_chain_pdb))
    with pytest.raises(ValueError, match="requires at least one of"):
        build_motifscaff_spec(args, all_chains=['A', 'B'], rank_map=rank_map)


def test_build_motifscaff_spec_with_architecture_change(two_chain_pdb, rank_map):
    args = _args(input_pdb=str(two_chain_pdb), motifscaff_spec='A15-17')
    spec = build_motifscaff_spec(args, all_chains=['A', 'B'], rank_map=rank_map)
    file_entity = spec['entities'][0]['file']
    assert file_entity['exclude'] == [{'chain': {'id': 'A', 'res_index': '1..3'}}]
    assert file_entity['reset_res_index'] == [{'chain': {'id': 'A'}}]


def test_build_motifscaff_spec_inpaint_seq_on_removed_residue_raises(two_chain_pdb, rank_map):
    args = _args(input_pdb=str(two_chain_pdb), motifscaff_spec='A15-17',
                 motifscaff_inpaint_seq='A10')
    with pytest.raises(ValueError, match="removed by motifscaff_spec"):
        build_motifscaff_spec(args, all_chains=['A', 'B'], rank_map=rank_map)


def test_build_motifscaff_spec_inpaint_seq_valid(two_chain_pdb, rank_map):
    args = _args(input_pdb=str(two_chain_pdb), motifscaff_spec='A15-17',
                 motifscaff_inpaint_seq='A15')
    spec = build_motifscaff_spec(args, all_chains=['A', 'B'], rank_map=rank_map)
    file_entity = spec['entities'][0]['file']
    # motifscaff_inpaint_seq residue tokens resolve to original chain A ranks
    # (rank_map-based), not to re-indexed post-kept positions.
    assert file_entity['design'] == [{'chain': {'id': 'A', 'res_index': '4'}}]


def test_build_motifscaff_spec_flexible_residues(two_chain_pdb, rank_map):
    args = _args(input_pdb=str(two_chain_pdb), flexible_residues='B')
    spec = build_motifscaff_spec(args, all_chains=['A', 'B'], rank_map=rank_map)
    file_entity = spec['entities'][0]['file']
    assert {'group': {'visibility': 0, 'id': 'B'}} in file_entity['structure_groups']
