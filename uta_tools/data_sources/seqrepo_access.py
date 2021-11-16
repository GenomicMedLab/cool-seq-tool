"""A module for accessing SeqRepo."""
from typing import Optional, List, Tuple
from biocommons.seqrepo import SeqRepo

from uta_tools import SEQREPO_DATA_PATH, logger
from os import environ
from uta_tools.schemas import ResidueMode


class SeqRepoAccess:
    """The SeqRepoAccess class."""

    def __init__(self, seqrepo_data_path: str = SEQREPO_DATA_PATH) -> None:
        """Initialize the SeqRepoAccess class.
        :param str seqrepo_data_path: The path to the seqrepo directory.
        """
        environ['SEQREPO_LRU_CACHE_MAXSIZE'] = "none"
        self.seqrepo_client = SeqRepo(seqrepo_data_path)

    def _get_start_end(
            self, start: int, end: Optional[int] = None,
            residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Tuple[int, int]:
        """Get start and end in inter-residue (0-based) coords.

        :param int start: Start pos change
        :param Optional[int] end: End pos change
        :param ResidueMode residue_mode: Residue mode for start/end positions
            Must be either `inter-residue` or `residue`
        :return: start pos, end pos
        """
        if residue_mode == ResidueMode.RESIDUE:
            if end is None:
                end = start
            else:
                end -= 1
            start -= 1
        return start, end

    def get_reference_sequence(
            self, ac: str, start: int, end: Optional[int] = None,
            residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Tuple[Optional[str], Optional[str]]:
        """Get reference sequence for transcript at given positions

        :param str ac: Accession
        :param int start: Start pos change
        :param Optional[int] end: End pos change
        :param ResidueMode residue_mode: Residue mode for start/end positions
            Must be either `inter-residue` or `residue`
        :return: Sequence at position, warning
        """
        sequence, warnings = self.check_sequence(
            ac, start, end, residue_mode)
        if warnings:
            return None, warnings
        else:
            if sequence:
                return sequence, None
            else:
                start, end = self._get_start_end(start, end, residue_mode)
                if start != end:
                    msg = f"SeqRepo returned no sequence for inter-residue " \
                          f"coordinates {start}, {end} on {ac}"
                else:
                    msg = f"Inter-residue start ({start}) and end ({end}) " \
                          f"coordinates must have different values for " \
                          f"SeqRepo to return a reference sequence"
                logger.warning(msg)
                return None, msg

    def check_sequence(
            self, ac: str, start: int = None,
            end: Optional[int] = None,
            residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Tuple[Optional[str], Optional[str]]:
        """Check that accession and positions actually exist

        :param str ac: Accession
        :param int start: Start pos change
        :param Optional[int] end: End pos change
        :param ResidueMode residue_mode: Residue mode for start/end positions
            Must be either `inter-residue` or `residue`
        :return: Sequence at position (if accession and positions actually
            exist), warning
        """
        start, end = self._get_start_end(start, end, residue_mode)
        try:
            sequence = self.seqrepo_client.fetch(ac, start=start, end=end)
        except KeyError:
            msg = f"Accession, {ac}, not found in SeqRepo"
            logger.warning(msg)
            return None, msg
        except ValueError as e:
            error = str(e)
            if error.startswith("start out of range"):
                msg = f"Start inter-residue coordinate {start} is out of " \
                      f"range on {ac}"
            elif error.startswith("stop out of range"):
                msg = f"End inter-residue coordinate {end} is out of " \
                      f"range on {ac}"
            elif error.startswith("invalid coordinates") and ">" in error:
                msg = f"Invalid inter-residue coordinates: start ({start}) " \
                      f"cannot be greater than end ({end})"
            else:
                msg = f"{e}"
            logger.warning(msg)
            return None, msg
        else:
            return sequence, None

    def is_valid_input_sequence(
            self, ac: str, start: int = None,
            end: Optional[int] = None,
            residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Tuple[bool, Optional[str]]:
        """Determine whether or not input sequence is valid.

        :param str ac: Accession
        :param int start: Start pos change
        :param Optional[int] end: End pos change
        :param ResidueMode residue_mode: Residue mode for start/end positions
            Must be either `inter-residue` or `residue`
        :return: Bool on whether or not input is valid, warning
        """
        seq, warning = self.check_sequence(
            ac, start, end, residue_mode)
        if warning:
            return False, warning
        else:
            return True, None

    def translate_identifier(
            self, ac: str, target_namespace: str = None
    ) -> Tuple[List[Optional[str]], Optional[str]]:
        """Return list of identifiers for accession.

        :param str ac: Identifier accession
        :param str target_namespace: The namespace of identifiers to return
        :return: List of identifiers, warning
        """
        try:
            ga4gh_identifiers = self.seqrepo_client.translate_identifier(
                ac, target_namespaces=target_namespace)
        except KeyError:
            msg = f"SeqRepo unable to get translated identifiers for {ac}"
            logger.warning(msg)
            return [], msg
        else:
            return ga4gh_identifiers, None

    def aliases(self,
                input_str: str) -> Tuple[List[Optional[str]], Optional[str]]:
        """Get aliases for a given input.

        :param str input_str: Input to get aliases for
        :return: List of aliases, warning
        """
        try:
            return self.seqrepo_client.translate_alias(input_str), None
        except KeyError:
            msg = f"SeqRepo could not translate alias {input_str}"
            logger.warning(msg)
            return [], msg
