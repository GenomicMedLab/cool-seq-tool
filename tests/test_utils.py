"""Module for testing residue mode"""
from cool_seq_tool.schemas import ResidueMode
from cool_seq_tool.utils import get_inter_residue_pos


def test_get_inter_residue_pos():
    """Test that get_inter_residue_pos method works correctly"""
    expected = ((599, 599), None)
    resp = get_inter_residue_pos(600, ResidueMode.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(600, ResidueMode.RESIDUE, end_pos=600)
    assert resp == expected

    resp = get_inter_residue_pos(599, ResidueMode.INTER_RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(599, ResidueMode.INTER_RESIDUE, end_pos=599)
    assert resp == expected

    resp = get_inter_residue_pos(600, "test")
    assert resp == (None, "residue_mode must be either `residue` "
                          "or `inter-residue`, not `test`")
