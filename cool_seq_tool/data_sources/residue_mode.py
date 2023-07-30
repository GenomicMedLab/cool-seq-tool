"""Module for converting positions to inter-residue coordinates"""
import logging
from typing import Tuple

from cool_seq_tool.schemas import ResidueMode


logger = logging.getLogger("cool_seq_tool")


def get_inter_residue_pos(
    start_pos: int, end_pos: int, residue_mode: ResidueMode
) -> Tuple[int, int]:
    """Return inter-residue position

    :param start_pos: Start position
    :param end_pos: End position
    :param residue_mode: Residue mode for `start_pos` and `end_pos`
    :return: Inter-residue coordinates for `start_pos` and `end_pos`
    """
    if residue_mode == ResidueMode.RESIDUE:
        start_pos -= 1
    return (start_pos, end_pos)
