"""Provide a small set of general helper functions."""

import datetime
import logging

from cool_seq_tool import __version__
from cool_seq_tool.schemas import ResidueMode, ServiceMeta

_logger = logging.getLogger(__name__)


def get_inter_residue_pos(
    start_pos: int, end_pos: int, residue_mode: ResidueMode
) -> tuple[int, int]:
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
    :param end_pos: End position
    :param residue_mode: Residue mode for `start_pos` and `end_pos`
    :return: Inter-residue coordinates
    """
    if residue_mode == ResidueMode.RESIDUE:
        start_pos -= 1
    elif residue_mode == ResidueMode.ZERO:
        end_pos += 1
    return start_pos, end_pos


@staticmethod
def service_meta() -> ServiceMeta:
    """Return description of request and service, including parameters like software
    version for reproducibility.

    :return: ServiceMeta object
    """
    return ServiceMeta(
        version=__version__,
        response_datetime=datetime.datetime.now(tz=datetime.timezone.utc),
    )
