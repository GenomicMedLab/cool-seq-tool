"""Module for common utilities used throughout the app"""
import logging
from datetime import datetime
from typing import Tuple

from cool_seq_tool.schemas import ResidueMode, ServiceMeta
from cool_seq_tool.version import __version__

logger = logging.getLogger(__name__)


def get_inter_residue_pos(
    start_pos: int, end_pos: int, residue_mode: ResidueMode
) -> Tuple[int, int]:
    """Return inter-residue position

    :param start_pos: Start position
    :param end_pos: End position
    :param residue_mode: Residue mode for `start_pos` and `end_pos`
    :return: Positions represented as inter-residue coordinates
    """
    if residue_mode == ResidueMode.RESIDUE:
        start_pos -= 1
    elif residue_mode == ResidueMode.ZERO:
        end_pos += 1
    return start_pos, end_pos


@staticmethod
def service_meta() -> ServiceMeta:
    """Return ServiceMeta for cool_seq_tool

    :return: ServiceMeta object
    """
    return ServiceMeta(version=__version__, response_datetime=datetime.now())
