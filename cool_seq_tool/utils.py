"""Provide a small set of general helper functions."""
import logging
from datetime import datetime
from typing import Optional, Tuple

from cool_seq_tool.schemas import ResidueMode, ServiceMeta
from cool_seq_tool.version import __version__

logger = logging.getLogger(__name__)


def get_inter_residue_pos(
    start_pos: int, residue_mode: ResidueMode, end_pos: Optional[int] = None
) -> Tuple[Optional[Tuple[int, int]], Optional[str]]:
    """Return equivalent inter-residue position.

    Generally, we prefer to work with inter-residue coordinates where possible. Our
    rationale is detailed in an appendix to the
    `VRS docs <https://vrs.ga4gh.org/en/stable/appendices/design_decisions.html#inter-residue-coordinates>`_.
    This function is used internally to shift user-provided coordinates accordingly.

    >>> from cool_seq_tool.utils import get_inter_residue_pos
    >>> from cool_seq_tool.schemas import ResidueMode
    >>> get_inter_residue_pos(10, ResidueMode.RESIDUE)
    ((9, 9), None)

    :param start_pos: Start position
    :param residue_mode: Residue mode for ``start_pos`` and ``end_pos``
    :param end_pos: End position. If ``None`` assumes both ``start`` and ``end`` have
        same values.
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
