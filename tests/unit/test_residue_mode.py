"""Module for testing residue mode"""
from cool_seq_tool.data_sources.residue_mode import get_inter_residue_pos
from cool_seq_tool.schemas import ResidueMode


def test_get_inter_residue_pos():
    """Test that get_inter_residue_pos method works correctly"""
    expected = (599, 600)
    resp = get_inter_residue_pos(600, 600, ResidueMode.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(599, 600, ResidueMode.INTER_RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(599, 599, ResidueMode.INTER_RESIDUE)
    assert resp == (599, 599)

    expected = (746, 751)
    resp = get_inter_residue_pos(747, 751, ResidueMode.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(746, 751, ResidueMode.INTER_RESIDUE)
    assert resp == expected
