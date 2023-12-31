"""Wrap SeqRepo to provide additional lookup and identification methods on top of basic
dereferencing functions.
"""
import logging
from os import environ
from pathlib import Path
from typing import List, Optional, Union

from ga4gh.vrs.dataproxy import SeqRepoDataProxy

from cool_seq_tool.exceptions import SeqRepoAccessError
from cool_seq_tool.schemas import ResidueMode
from cool_seq_tool.utils import get_inter_residue_pos

logger = logging.getLogger(__name__)


class SeqRepoAccess(SeqRepoDataProxy):
    """Provide a wrapper around the base SeqRepoDataProxy class from ``VRS-Python`` to
    provide additional lookup and identification methods.
    """

    environ["SEQREPO_LRU_CACHE_MAXSIZE"] = "none"

    def get_reference_sequence(
        self,
        ac: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> str:
        """Get reference sequence for an accession given a start and end position. If
        ``start`` and ``end`` are not given, returns the entire reference sequence.

        >>> from cool_seq_tool.handlers import SeqRepoAccess
        >>> from biocommons.seqrepo import SeqRepo
        >>> sr = SeqRepoAccess(SeqRepo("/usr/local/share/seqrepo/latest"))
        >>> sr.get_reference_sequence("NM_002529.3", 1, 10)[0]
        'TGCAGCTGG'
        >>> sr.get_reference_sequence("NP_001341538.1", 1, 10)[0]
        'MAALSGGGG'

        :param ac: Accession
        :param start: Start pos change
        :param end: End pos change. If ``None`` assumes both ``start`` and ``end`` have
            same values, if ``start`` exists.
        :param residue_mode: Residue mode for ``start`` and ``end``
        :return: Sequence at position(s) on an accession
        :raises SeqRepoAccessError: If ``start`` is greater than than ``end``, accession
            does not exist in SeqRepo, or coordinates are out of index on an accession
        """
        if start and end:
            if start > end:
                raise SeqRepoAccessError(
                    f"start ({start}) cannot be greater than end ({end})"
                )

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
            raise SeqRepoAccessError(msg)
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
            raise SeqRepoAccessError(msg)
        else:
            # If start is valid, but end is invalid, SeqRepo still returns
            # the sequence from valid positions. So we want to make sure
            # that both start and end positions are valid
            if start and end:
                expected_len_of_seq = end - start
                if len(sequence) != expected_len_of_seq:
                    raise SeqRepoAccessError(
                        f"End inter-residue coordinate ({end}) is out of index on {ac}"
                    )

            return sequence

    def translate_identifier(
        self, ac: str, target_namespaces: Optional[Union[str, List[str]]] = None
    ) -> List[str]:
        """Return list of identifiers for accession.

        >>> from cool_seq_tool.handlers import SeqRepoAccess
        >>> from biocommons.seqrepo import SeqRepo
        >>> sr = SeqRepoAccess(SeqRepo("/usr/local/share/seqrepo/latest"))
        >>> sr.translate_identifier("NM_002529.3")[0]
        ['MD5:18f0a6e3af9e1bbd8fef1948c7156012', 'NCBI:NM_002529.3', 'refseq:NM_002529.3', 'SEGUID:dEJQBkga9d9VeBHTyTbg6JEtTGQ', 'SHA1:74425006481af5df557811d3c936e0e8912d4c64', 'VMC:GS_RSkww1aYmsMiWbNdNnOTnVDAM3ZWp1uA', 'sha512t24u:RSkww1aYmsMiWbNdNnOTnVDAM3ZWp1uA', 'ga4gh:SQ.RSkww1aYmsMiWbNdNnOTnVDAM3ZWp1uA']
        >>> sr.translate_identifier("NM_002529.3", "ga4gh")[0]
        ['ga4gh:SQ.RSkww1aYmsMiWbNdNnOTnVDAM3ZWp1uA']

        :param ac: Identifier accession
        :param target_namespace: The namespace(s) of identifier to return
        :return: List of identifiers
        :raises SeqRepoAccessError: If accession does not exist in SeqRepo
        """
        try:
            ga4gh_identifiers = self.sr.translate_identifier(
                ac, target_namespaces=target_namespaces
            )
        except KeyError:
            msg = f"SeqRepo unable to get translated identifiers for {ac}"
            logger.warning(msg)
            raise SeqRepoAccessError(msg)
        else:
            return ga4gh_identifiers

    def translate_alias(self, input_str: str) -> List[str]:
        """Get aliases for a given input.

        :param input_str: Input to get aliases for
        :return: List of aliases
        :raises: SeqRepoAccessError if ``input_str`` does not exist in SeqRepo
        """
        try:
            return self.sr.translate_alias(input_str)
        except KeyError:
            msg = f"SeqRepo could not translate alias {input_str}"
            logger.warning(msg)
            raise SeqRepoAccessError(msg)

    def get_fasta_file(self, sequence_id: str, outfile_path: Path) -> None:
        """Retrieve FASTA file containing sequence for requested sequence ID.

        >>> from pathlib import Path
        >>> from cool_seq_tool.handlers import SeqRepoAccess
        >>> from biocommons.seqrepo import SeqRepo
        >>> sr = SeqRepoAccess(SeqRepo("/usr/local/share/seqrepo/latest"))
        >>> # write to local file tpm3.fasta:
        >>> sr.get_fasta_file("NM_002529.3", Path("tpm3.fasta"))

        FASTA file headers will include GA4GH sequence digest, Ensembl accession ID,
        and RefSeq accession ID.

        :param sequence_id: accession ID, sans namespace, eg ``NM_152263.3``
        :param outfile_path: path to save file to
        :return: None, but saves sequence data to ``outfile_path`` if successful
        """
        sequence = self.get_reference_sequence(sequence_id)

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
            header = f">refseq:{sequence_id}|{'|'.join(aliases)}"
        elif sequence_id[:4] in ensembl_prefixes:
            aliases = self.translate_identifier(sequence_id, ["refseq", "ga4gh"])
            header = f">ensembl:{sequence_id}|{'|'.join(aliases)}"
        else:
            aliases = self.translate_identifier(
                sequence_id, ["ensembl", "refseq", "ga4gh"]
            )
            header = f">gnl|ID|{sequence_id}|{'|'.join(aliases)}"

        line_length = 60
        file_data = [header] + [
            sequence[i : i + line_length] for i in range(0, len(sequence), line_length)
        ]
        text = "\n".join(file_data)
        outfile_path.write_text(text)
