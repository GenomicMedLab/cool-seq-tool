"""A module for accessing SeqRepo."""
import logging
from os import environ
from pathlib import Path
from typing import List, Optional, Tuple, Union

from ga4gh.vrs.dataproxy import SeqRepoDataProxy

from cool_seq_tool.schemas import ResidueMode
from cool_seq_tool.utils import get_inter_residue_pos

logger = logging.getLogger(__name__)


class SeqRepoAccess(SeqRepoDataProxy):
    """The SeqRepoAccess class."""

    environ["SEQREPO_LRU_CACHE_MAXSIZE"] = "none"

    def get_reference_sequence(
        self,
        ac: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> Tuple[str, Optional[str]]:
        """Get reference sequence for an accession given a start and end position.
        If `start` and `end` are not given, it will return the entire reference sequence

        :param ac: Accession
        :param start: Start pos change
        :param end: End pos change. If `None` assumes both `start` and `end` have same
            values, if `start` exists.
        :param residue_mode: Residue mode for `start` and `end`
        :return: Sequence at position (if accession and positions actually
            exist, else return empty string), warning if any
        """
        if start and end:
            if start > end:
                msg = f"start ({start}) cannot be greater than end ({end})"
                return "", msg

            start, end = get_inter_residue_pos(start, end, residue_mode)
            if start == end:
                end += 1
        else:
            if start is not None and residue_mode == ResidueMode.RESIDUE:
                start -= 1

        try:
            sequence = self.sr.fetch(ac, start=start, end=end)
        except KeyError:
            msg = f"Accession, {ac}, not found in SeqRepo"
            logger.warning(msg)
            return "", msg
        except ValueError as e:
            error = str(e)
            if error.startswith("start out of range"):
                msg = (
                    f"Start inter-residue coordinate ({start}) is out of index on {ac}"
                )
            elif error.startswith("stop out of range"):
                msg = (
                    f"End inter-residue coordinate ({end}) is out of " f"index on {ac}"
                )
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
                    return (
                        "",
                        f"End inter-residue coordinate ({end}) is out of index on {ac}",
                    )
            return sequence, None

    def translate_identifier(
        self, ac: str, target_namespaces: Optional[Union[str, List[str]]] = None
    ) -> Tuple[List[str], Optional[str]]:
        """Return list of identifiers for accession.

        :param ac: Identifier accession
        :param target_namespace: The namespace(s) of identifier to return
        :return: List of identifiers, warning
        """
        try:
            ga4gh_identifiers = self.sr.translate_identifier(
                ac, target_namespaces=target_namespaces
            )
        except KeyError:
            msg = f"SeqRepo unable to get translated identifiers for {ac}"
            logger.warning(msg)
            return [], msg
        else:
            return ga4gh_identifiers, None

    def translate_alias(
        self, input_str: str
    ) -> Tuple[List[Optional[str]], Optional[str]]:
        """Get aliases for a given input.

        :param str input_str: Input to get aliases for
        :return: List of aliases, warning
        """
        try:
            return self.sr.translate_alias(input_str), None
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
            tmp_acs, _ = self.translate_identifier(
                f"{assembly}:chr{chromosome}", target_namespaces="refseq"
            )
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
        aliases, _ = self.translate_alias(ac)
        aliases = (
            [
                a.split(":")[-1]
                for a in aliases
                if a.startswith("GRCh") and "." not in a and "chr" not in a
            ]
            or [None]
        )[0]
        if aliases is None:
            return None, f"Unable to get chromosome for {ac}"
        else:
            return aliases, None

    def get_fasta_file(self, sequence_id: str, outfile_path: Path) -> None:
        """Retrieve FASTA file containing sequence for requested sequence ID.
        :param sequence_id: accession ID, sans namespace, eg `NM_152263.3`
        :param outfile_path: path to save file to
        :return: None, but saves sequence data to `outfile_path` if successful
        :raise: KeyError if SeqRepo doesn't have sequence data for the given ID
        """
        sequence = self.get_reference_sequence(sequence_id)[0]
        if not sequence:
            raise KeyError

        refseq_prefixes = [
            "NC_",
            "AC_",
            "NZ_",
            "NT_",
            "NW_",
            "NG_",
            "NM_",
            "XM_",
            "NR_",
            "XR_",
            "NP_",
            "AP_",
            "XP_",
            "YP_",
            "WP_",
        ]
        ensembl_prefixes = ["ENSE", "ENSFM", "ENSG", "ENSGT", "ENSP", "ENSR", "ENST"]

        if sequence_id[:3] in refseq_prefixes:
            aliases = self.translate_identifier(sequence_id, ["ensembl", "ga4gh"])
            header = f">refseq:{sequence_id}|{'|'.join(aliases[0])}"
        elif sequence_id[:4] in ensembl_prefixes:
            aliases = self.translate_identifier(sequence_id, ["refseq", "ga4gh"])
            header = f">ensembl:{sequence_id}|{'|'.join(aliases[0])}"
        else:
            aliases = self.translate_identifier(
                sequence_id, ["ensembl", "refseq", "ga4gh"]
            )
            header = f">gnl|ID|{sequence_id}|{'|'.join(aliases[0])}"

        line_length = 60
        file_data = [header] + [
            sequence[i : i + line_length] for i in range(0, len(sequence), line_length)
        ]
        text = "\n".join(file_data)
        outfile_path.write_text(text)
