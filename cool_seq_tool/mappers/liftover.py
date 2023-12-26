"""Provide position conversion between GRCh37 and GRCh38 assemblies via PyLiftover"""
import logging
from os import environ
from typing import Optional

from pyliftover import LiftOver as PyLiftOver

from cool_seq_tool.schemas import Assembly

# Environment variables for paths to chain files for pyliftover
LIFTOVER_CHAIN_37_TO_38 = environ.get("LIFTOVER_CHAIN_37_TO_38")
LIFTOVER_CHAIN_38_TO_37 = environ.get("LIFTOVER_CHAIN_38_TO_37")


logger = logging.getLogger(__name__)


class Liftover:
    """Provide position conversion between GRCh37 and GRCh38 via PyLiftover"""

    def __init__(
        self,
        chain_file_37_to_38: Optional[str] = None,
        chain_file_38_to_37: Optional[str] = None,
    ) -> None:
        """Initialize liftover class.

        :param chain_file_37_to_38: Optional path to chain file for 37 to 38 assembly.
            This is used for ``pyliftover``. If this is not provided, will check to see
            if ``LIFTOVER_CHAIN_37_TO_38`` env var is set. If neither is provided, will
            allow ``pyliftover`` to download a chain file from UCSC
        :param chain_file_38_to_37: Optional path to chain file for 38 to 37 assembly.
            This is used for ``pyliftover``. If this is not provided, will check to see
            if ``LIFTOVER_CHAIN_38_TO_37`` env var is set. If neither is provided, will
            allow ``pyliftover`` to download a chain file from UCSC
        """
        chain_file_37_to_38 = chain_file_37_to_38 or LIFTOVER_CHAIN_37_TO_38
        if chain_file_37_to_38:
            self.grch37_to_grch38 = PyLiftOver(chain_file_37_to_38)
        else:
            self.grch37_to_grch38 = PyLiftOver("hg19", "hg38")

        chain_file_38_to_37 = chain_file_38_to_37 or LIFTOVER_CHAIN_38_TO_37
        if chain_file_38_to_37:
            self.grch38_to_grch37 = PyLiftOver(chain_file_38_to_37)
        else:
            self.grch38_to_grch37 = PyLiftOver("hg38", "hg19")

    def convert_pos(
        self, chromosome: str, pos: int, target_assembly: Assembly = Assembly.GRCH38
    ) -> Optional[int]:
        """Get new genome assembly data for a position on a chromosome.

        :param chromosome: The chromosome number. Must be prefixed with ``chr``
        :param pos: Position on the chromosome
        :param target_assembly: Assembly to liftover to
        :return: Position on ``target_assembly``
        """
        if not chromosome.startswith("chr"):
            logger.warning("`chromosome` must be prefixed with chr")
            return None

        if target_assembly == Assembly.GRCH38:
            liftover = self.grch37_to_grch38.convert_coordinate(chromosome, pos)
        elif target_assembly == Assembly.GRCH37:
            liftover = self.grch38_to_grch37.convert_coordinate(chromosome, pos)
        else:
            logger.warning(f"{target_assembly} assembly not supported")
            liftover = None

        if not liftover:
            logger.warning(f"{pos} does not exist on {chromosome}")
            return None
        else:
            return liftover[0][1]
