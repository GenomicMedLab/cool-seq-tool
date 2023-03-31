"""Module for testing position functions"""
from cool_seq_tool.utils.positions import get_inter_residue_pos
from cool_seq_tool.schemas import ResidueMode


def test_get_inter_residue_pos():
    """Test that get_inter_residue_pos function works correctly"""
    expected = 599, 600
    resp = get_inter_residue_pos(600, 600, ResidueMode.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(600, 600, ResidueMode.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(599, 600, ResidueMode.INTER_RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(599, 599, ResidueMode.INTER_RESIDUE)
    assert resp == (599, 599)
