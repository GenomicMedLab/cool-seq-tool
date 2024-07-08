"""Module for mapping to/from human genome assemblies.

Currently only supports GRCh37 <-> GRCh38
"""

import logging
from os import environ

from agct import Converter, Genome

from cool_seq_tool.schemas import Assembly
from cool_seq_tool.utils import process_chromosome_input

# Environment variables for paths to chain files for agct
LIFTOVER_CHAIN_37_TO_38 = environ.get("LIFTOVER_CHAIN_37_TO_38")
LIFTOVER_CHAIN_38_TO_37 = environ.get("LIFTOVER_CHAIN_38_TO_37")


_logger = logging.getLogger(__name__)


class LiftOver:
    """Class for mapping to/from human genome assemblies

    Currently only supports GRCh37 <-> GRCh38
    """

    def __init__(
        self,
        chain_file_37_to_38: str | None = None,
        chain_file_38_to_37: str | None = None,
    ) -> None:
        """Initialize liftover class

        :param chain_file_37_to_38: Optional path to chain file for 37 to 38 assembly.
            This is used for ``agct``. If this is not provided, will check to see
            if ``LIFTOVER_CHAIN_37_TO_38`` env var is set. If neither is provided, will
            allow ``agct`` to download a chain file from UCSC
        :param chain_file_38_to_37: Optional path to chain file for 38 to 37 assembly.
            This is used for ``agct``. If this is not provided, will check to see
            if ``LIFTOVER_CHAIN_38_TO_37`` env var is set. If neither is provided, will
            allow ``agct`` to download a chain file from UCSC
        """
        self.from_37_to_38 = Converter(
            chainfile=chain_file_37_to_38 or LIFTOVER_CHAIN_37_TO_38,
            from_db=Genome.HG19,
            to_db=Genome.HG38,
        )
        self.from_38_to_37 = Converter(
            chainfile=chain_file_38_to_37 or LIFTOVER_CHAIN_38_TO_37,
            from_db=Genome.HG38,
            to_db=Genome.HG19,
        )

    def get_liftover(
        self, chromosome: str, pos: int, liftover_to_assembly: Assembly
    ) -> tuple[str, int] | None:
        """Get new genome assembly data for a position on a chromosome.

        Use a UCSC-style chromosome name:

        >>> from cool_seq_tool.mappers import LiftOver
        >>> from cool_seq_tool.schemas import Assembly
        >>> lo = LiftOver()
        >>> lo.get_liftover("chr7", 140453136, Assembly.GRCH38)
        ('chr7', 140753336)

        Chromosome names can also be NCBI-style, without prefixes:

        >>> lo.get_liftover("7", 140453136, Assembly.GRCH38)
        ('chr7', 140753336)

        :param chromosome: The chromosome number, e.g. ``"chr7"``, ``"chrX"``, ``"5"``.
        :param pos: Position on the chromosome
        :param liftover_to_assembly: Assembly to liftover to
        :return: Target chromosome and target position for assembly
        """
        chromosome = process_chromosome_input(chromosome, "LiftOver.get_liftover()")
        if liftover_to_assembly == Assembly.GRCH38:
            liftover = self.from_37_to_38.convert_coordinate(chromosome, pos)
        elif liftover_to_assembly == Assembly.GRCH37:
            liftover = self.from_38_to_37.convert_coordinate(chromosome, pos)
        else:
            _logger.warning("%s assembly not supported", liftover_to_assembly)
            liftover = None

        if not liftover:
            _logger.warning("%s does not exist on %s", pos, chromosome)
            return None
        return liftover[0][:2]
