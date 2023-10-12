"""Module for common utilities used throughout the app"""
import logging
from datetime import datetime
from typing import Optional, Tuple

from cool_seq_tool.schemas import ResidueMode, ServiceMeta
from cool_seq_tool.version import __version__

logger = logging.getLogger(__name__)


def get_inter_residue_pos(
    start_pos: int, residue_mode: ResidueMode, end_pos: Optional[int] = None
) -> Tuple[Optional[Tuple[int, int]], Optional[str]]:
    """Return inter-residue position

    :param start_pos: Start position
    :param residue_mode: Residue mode for `start_pos` and `end_pos`
    :param end_pos: End position. If `None` assumes both `start` and `end` have same
        values.
    :return: Inter-residue coordinates, warning
    """
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
        msg = (
            f"residue_mode must be either `residue` or `inter-residue`,"
            f" not `{residue_mode}`"
        )
        logger.warning(msg)
        return None, msg
    return (start_pos, end_pos), None


@staticmethod
def service_meta() -> ServiceMeta:
    """Return ServiceMeta for cool_seq_tool

    :return: ServiceMeta object
    """
    return ServiceMeta(version=__version__, response_datetime=datetime.now())
