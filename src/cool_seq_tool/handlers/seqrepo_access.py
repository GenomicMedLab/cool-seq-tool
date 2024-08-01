"""Wrap SeqRepo to provide additional lookup and identification methods on top of basic
dereferencing functions.
"""

import logging
from os import environ
from pathlib import Path

from ga4gh.vrs.dataproxy import SeqRepoDataProxy

from cool_seq_tool.schemas import Assembly, CoordinateType
from cool_seq_tool.utils import get_inter_residue_pos, process_chromosome_input

_logger = logging.getLogger(__name__)


SEQREPO_ROOT_DIR = environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest")


class SeqRepoAccess(SeqRepoDataProxy):
    """Provide a wrapper around the base SeqRepoDataProxy class from ``VRS-Python`` to
    provide additional lookup and identification methods.
    """

    environ["SEQREPO_LRU_CACHE_MAXSIZE"] = "none"

    def get_reference_sequence(
        self,
        ac: str,
        start: int | None = None,
        end: int | None = None,
        coordinate_type: CoordinateType = CoordinateType.RESIDUE,
    ) -> tuple[str, str | None]:
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
        :param coordinate_type: Coordinate type for ``start`` and ``end``
        :return: Sequence at position (if accession and positions actually
            exist, else return empty string), warning if any
        """
        if start and end:
            if start > end:
                msg = f"start ({start}) cannot be greater than end ({end})"
                return "", msg

            start, end = get_inter_residue_pos(start, end, coordinate_type)
            if start == end:
                end += 1
        else:
            if start is not None and coordinate_type == CoordinateType.RESIDUE:
                start -= 1

        try:
            sequence = self.sr.fetch(ac, start=start, end=end)
        except KeyError:
            msg = f"Accession, {ac}, not found in SeqRepo"
            _logger.warning(msg)
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
            _logger.warning(msg)
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
        self, ac: str, target_namespaces: str | list[str] | None = None
    ) -> tuple[list[str], str | None]:
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
        :return: List of identifiers, warning
        """
        try:
            ga4gh_identifiers = self.sr.translate_identifier(
                ac, target_namespaces=target_namespaces
            )
        except KeyError:
            msg = f"SeqRepo unable to get translated identifiers for {ac}"
            _logger.warning(msg)
            return [], msg
        else:
            return ga4gh_identifiers, None

    def translate_alias(self, input_str: str) -> tuple[list[str | None], str | None]:
        """Get aliases for a given input.

        :param str input_str: Input to get aliases for
        :return: List of aliases, warning
        """
        try:
            return self.sr.translate_alias(input_str), None
        except KeyError:
            msg = f"SeqRepo could not translate alias {input_str}"
            _logger.warning(msg)
            return [], msg

    def chromosome_to_acs(self, chromosome: str) -> tuple[list[str] | None, str | None]:
        """Get accessions for a chromosome

        :param chromosome: Chromosome number. Must be either 1-22, X, or Y
        :return: Accessions for chromosome (ordered by latest assembly)
        """
        acs = []
        for assembly in reversed(Assembly.values()):
            tmp_acs, _ = self.translate_identifier(
                f"{assembly}:{process_chromosome_input(chromosome)}",
                target_namespaces="refseq",
            )
            acs += [ac.split("refseq:")[-1] for ac in tmp_acs]
        if acs:
            return acs, None
        return (
            None,
            f'Unable to find matching accessions for "{chromosome}" in SeqRepo.',
        )

    def ac_to_chromosome(self, ac: str) -> tuple[str | None, str | None]:
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
        return aliases, None

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
