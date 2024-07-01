"""Module for testing utility functions"""

from cool_seq_tool.schemas import ResidueMode
from cool_seq_tool.utils import get_inter_residue_pos


def test_get_inter_residue_pos():
    """Test that get_inter_residue_pos function works correctly"""
    expected = (599, 600)
    resp = get_inter_residue_pos(599, 599, ResidueMode.ZERO)
    assert resp == expected

    resp = get_inter_residue_pos(600, 600, ResidueMode.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(599, 600, ResidueMode.INTER_RESIDUE)
    assert resp == expected

    expected = (746, 751)
    resp = get_inter_residue_pos(746, 750, ResidueMode.ZERO)
    assert resp == expected

    resp = get_inter_residue_pos(747, 751, ResidueMode.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(746, 751, ResidueMode.INTER_RESIDUE)
    assert resp == expected
