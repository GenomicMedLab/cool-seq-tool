"""Module for converting positions to inter-residue coordinates"""
from typing import Optional, Tuple

from uta_tools.schemas import ResidueMode
from uta_tools import logger


def get_inter_residue_pos(
        start_pos: int, end_pos: Optional[int], residue_mode: str
) -> Tuple[Optional[Tuple[int, int]], Optional[str]]:
    """Return inter-residue position

    :param int start_pos: Start position
    :param Optional[int] end_pos: End position
    :param str residue_mode: `inter-residue` if start/end are 0 based coords.
        `residue` if start/end are 1 based coords
    :return: Inter-residue coordinates, warning
    """
    residue_mode = residue_mode.lower()
    if residue_mode == ResidueMode.RESIDUE:
        start_pos -= 1
        if end_pos is None:
            end_pos = start_pos
        else:
            end_pos -= 1
    elif residue_mode == ResidueMode.INTER_RESIDUE:
        if end_pos is None:
            end_pos = start_pos
    else:
        msg = f"residue_mode must be either `residue` or `inter-residue`," \
              f" not `{residue_mode}`"
        logger.warning(msg)
        return None, msg
    return (start_pos, end_pos), None
