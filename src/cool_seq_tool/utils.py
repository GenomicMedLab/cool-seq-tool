"""Provide a small set of general helper functions."""

import datetime
import logging

from bioutils.accessions import chr22XY

from cool_seq_tool import __version__
from cool_seq_tool.schemas import CoordinateType, ServiceMeta

_logger = logging.getLogger(__name__)


def get_inter_residue_pos(
    start_pos: int, end_pos: int, coordinate_type: CoordinateType
) -> tuple[int, int]:
    """Return equivalent inter-residue position.

    Residue coordinates start with 1, whereas inter-residue coordinates start with 0.

    It is preferred to work with inter-residue coordinates where possible. Our
    rationale is detailed in an appendix to the
    `VRS docs <https://vrs.ga4gh.org/en/stable/appendices/design_decisions.html#inter-residue-coordinates>`_.
    This function is used internally to shift user-provided coordinates accordingly.

    >>> from cool_seq_tool.utils import get_inter_residue_pos
    >>> from cool_seq_tool.schemas import CoordinateType
    >>> get_inter_residue_pos(10, CoordinateType.RESIDUE)
    ((9, 9), None)

    :param start_pos: Start position
    :param end_pos: End position
    :param coordinate_type: Coordinate type for `start_pos` and `end_pos`
    :return: Inter-residue coordinates
    """
    if coordinate_type == CoordinateType.RESIDUE:
        start_pos -= 1
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


def process_chromosome_input(chromosome: str, context: str = "") -> str:
    """Perform processing on a chromosome arg.

    E.g.

    >>> from cool_seq_tool.utils import process_chromosome_input
    >>> process_chromosome_input("7")
    'chr7'
    >>> process_chromosome_input("x")
    'chrX'
    >>> process_chromosome_input("chr7")
    'chr7'

    In the future, we could also use this method to be more opinionated about legal
    chromosome values, or throw exceptions in the event of invalid or unrecognized
    terms.

    :param chromosome: user-provided chromosome input
    :param context: calling context to provide in log
    :return: processed chromosome value. Idempotent -- returns original value if no
        changes needed.
    """
    original_chromosome_value = chromosome
    if chromosome.lower().startswith("chr"):
        chromosome = f"chr{chromosome[3:].upper()}"
    else:
        chromosome = chromosome.upper()
    chromosome = chr22XY(chromosome)
    if original_chromosome_value != chromosome:
        _logger.warning(
            "Transformed provided chromosome value from `%s` to `%s` in `%s`",
            original_chromosome_value,
            chromosome,
            context if context else "cool_seq_tool",
        )
    return chromosome
