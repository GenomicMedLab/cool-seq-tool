"""A module for accessing SeqRepo."""
from typing import Optional, List, Tuple, Union
from os import environ

from biocommons.seqrepo import SeqRepo

from cool_seq_tool.schemas import ResidueMode
from cool_seq_tool import SEQREPO_DATA_PATH, logger
from cool_seq_tool.data_sources.residue_mode import get_inter_residue_pos


class SeqRepoAccess:
    """The SeqRepoAccess class."""

    def __init__(self, seqrepo_data_path: str = SEQREPO_DATA_PATH) -> None:
        """Initialize the SeqRepoAccess class.
        :param str seqrepo_data_path: The path to the seqrepo directory.
        """
        environ["SEQREPO_LRU_CACHE_MAXSIZE"] = "none"
        self.seqrepo_client = SeqRepo(seqrepo_data_path)

    def get_reference_sequence(
            self, ac: str, start: Optional[int] = None, end: Optional[int] = None,
            residue_mode: str = ResidueMode.RESIDUE
    ) -> Tuple[str, Optional[str]]:
        """Get reference sequence for an accession given a start and end position.
        If `start` and `end` are not given, it will return the entire reference sequence

        :param str ac: Accession
        :param Optional[int] start: Start pos change
        :param Optional[int] end: End pos change. If `None` assumes both
            `start` and `end` have same values, if `start` exists.
        :param str residue_mode: Residue mode for start/end positions
            Must be either `inter-residue` or `residue`
        :return: Sequence at position (if accession and positions actually
            exist, else return empty string), warning if any
        """
        if start or end:
            pos, warning = get_inter_residue_pos(start, residue_mode, end_pos=end)
            if pos is None:
                return "", warning
            else:
                start, end = pos
                if start == end:
                    end += 1
        try:
            sequence = self.seqrepo_client.fetch(ac, start=start, end=end)
        except KeyError:
            msg = f"Accession, {ac}, not found in SeqRepo"
            logger.warning(msg)
            return "", msg
        except ValueError as e:
            error = str(e)
            if error.startswith("start out of range"):
                msg = f"Start inter-residue coordinate ({start}) is out of " \
                      f"index on {ac}"
            elif error.startswith("stop out of range"):
                msg = f"End inter-residue coordinate ({end}) is out of " \
                      f"index on {ac}"
            elif error.startswith("invalid coordinates") and ">" in error:
                msg = f"Invalid inter-residue coordinates: start ({start}) " \
                      f"cannot be greater than end ({end})"
            else:
                msg = f"{e}"
            logger.warning(msg)
            return "", msg
        else:
            # If start is valid, but end is invalid, SeqRepo still returns
            # the sequence from valid positions. So we want to make sure
            # that both start and end positions are valid
            if start and end:
                expected_len_of_seq = end - start
                if len(sequence) != expected_len_of_seq:
                    return "", f"End inter-residue coordinate ({end})" \
                               f" is out of index on {ac}"
            return sequence, None

    def translate_identifier(
            self, ac: str, target_namespace: Optional[Union[str, List[str]]] = None
    ) -> Tuple[List[str], Optional[str]]:
        """Return list of identifiers for accession.

        :param ac: Identifier accession
        :param target_namespace: The namespace(s) of identifier to return
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

    def chromosome_to_acs(
            self, chromosome: str
    ) -> Tuple[Optional[List[str]], Optional[str]]:
        """Get accessions for a chromosome

        :param str chromosome: Chromosome number. Must be either 1-22, X, or Y
        :return: Accessions for chromosome (ordered by latest assembly)
        """
        acs = []
        for assembly in ["GRCh38", "GRCh37"]:
            tmp_acs = self.translate_identifier(f"{assembly}:chr{chromosome}",
                                                target_namespace="refseq")[0]
            for ac in tmp_acs:
                acs.append(ac.split("refseq:")[-1])
        if acs:
            return acs, None
        else:
            return None, f"{chromosome} is not a valid chromosome"

    def ac_to_chromosome(self, ac: str) -> Tuple[Optional[str], Optional[str]]:
        """Get chromosome for accession.

        :param str ac: Accession
        :return: Chromosome, warning
        """
        aliases, warning = self.aliases(ac)
        aliases = ([a.split(":")[-1] for a in aliases
                    if a.startswith("GRCh") and "." not in a and "chr" not in a] or [None])[0]  # noqa: E501
        if aliases is None:
            return None, f"Unable to get chromosome for {ac}"
        else:
            return aliases, None
