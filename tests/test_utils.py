"""Module for testing utility functions"""

from cool_seq_tool.schemas import CoordinateType
from cool_seq_tool.utils import get_inter_residue_pos


def test_get_inter_residue_pos():
    """Test that get_inter_residue_pos function works correctly"""
    expected = (599, 600)
    resp = get_inter_residue_pos(600, 600, CoordinateType.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(599, 600, CoordinateType.INTER_RESIDUE)
    assert resp == expected

    expected = (746, 751)
    resp = get_inter_residue_pos(747, 751, CoordinateType.RESIDUE)
    assert resp == expected

    resp = get_inter_residue_pos(746, 751, CoordinateType.INTER_RESIDUE)
    assert resp == expected
