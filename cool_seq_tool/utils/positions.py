"""Module for converting positions"""
from typing import Optional, Tuple

from cool_seq_tool.schemas import ResidueMode
from cool_seq_tool import logger


def get_inter_residue_pos(
    start_pos: int, end_pos: int, residue_mode: ResidueMode
) -> Tuple[int, int]:
    """Return inter-residue position

    :param int start_pos: Start position
    :param int end_pos: End position
    :param ResidueMode residue_mode: Residue mode for `start_pos` and `end_pos`
    :return: Inter-residue coordinates
    """
    if residue_mode == ResidueMode.RESIDUE:
        start_pos -= 1
    return (start_pos, end_pos)


def to_zero_based(
    start_pos: int, end_pos: int, residue_mode: ResidueMode = ResidueMode.RESIDUE
) -> Tuple[int, int]:
    if residue_mode == ResidueMode.RESIDUE:
        start_pos -= 1
        end_pos -= 1
    else:
        if start_pos == end_pos:
            start_pos -= 1
        end_pos -= 1
    return start_pos, end_pos


def zero_based_to_inter_residue(
    start_pos: int, end_pos: int, begin_start_eq_end: bool,
    begin_residue_mode: ResidueMode = ResidueMode.RESIDUE
) -> Tuple[int, int]:
    if begin_residue_mode == ResidueMode.RESIDUE:
        end_pos += 1
    else:
        if begin_start_eq_end:
            start_pos += 1
        end_pos += 1

    return start_pos, end_pos
