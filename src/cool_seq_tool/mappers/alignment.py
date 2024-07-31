"""Module containing alignment methods for translating to and from different
reference sequences.
"""

from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.schemas import AnnotationLayer, Assembly, CoordinateType
from cool_seq_tool.sources import TranscriptMappings, UtaDatabase


class AlignmentMapper:
    """Class for translating between p --> c --> g reference sequences."""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        transcript_mappings: TranscriptMappings,
        uta_db: UtaDatabase,
    ) -> None:
        """Initialize the AlignmentMapper class.

        :param seqrepo_access: Access to seqrepo queries
        :param transcript_mappings: Access to transcript accession mappings and
            conversions
        :param uta_db: UtaDatabase instance to give access to query UTA database
        """
        self.seqrepo_access = seqrepo_access
        self.transcript_mappings = transcript_mappings
        self.uta_db = uta_db

    async def p_to_c(
        self,
        p_ac: str,
        p_start_pos: int,
        p_end_pos: int,
        coordinate_type: CoordinateType = CoordinateType.RESIDUE,
    ) -> tuple[dict | None, str | None]:
        """Translate protein representation to cDNA representation.

        :param p_ac: Protein RefSeq accession
        :param p_start_pos: Protein start position
        :param p_end_pos: Protein end position
        :param coordinate_type: Coordinate type for ``p_start_pos`` and ``p_end_pos``
        :return: Tuple containing:

        * cDNA representation (accession, codon range positions for corresponding
          change, cds start site) if able to translate. Will return positions as
          inter-residue coordinates. If unable to translate, returns ``None``.
        * Warning, if unable to translate to cDNA representation. Else ``None``
        """
        # Get cDNA accession
        temp_c_ac = await self.uta_db.p_to_c_ac(p_ac)
        if temp_c_ac:
            c_ac = temp_c_ac[-1]
        else:
            try:
                c_ac = self.transcript_mappings.np_to_nm[p_ac]
            except KeyError:
                return None, f"{p_ac} not found in transcript mappings"

        # Get coding start site
        cds_start, warning = await self._get_cds_start(c_ac)
        if not cds_start:
            return None, warning

        # For 1-based, protein pos * 3 = end of codon
        # 1 amino acid maps to 3 nucleotides in the codon
        # Since we have the end of the codon, we will subtract 2 to get the start of the
        # codon. We want to return inter-residue (0-based), so we subtract 1 from this.
        if coordinate_type == CoordinateType.RESIDUE:
            c_pos = (p_start_pos * 3) - 3, p_end_pos * 3
        else:
            if p_start_pos == p_end_pos:
                c_pos = ((p_start_pos + 1) * 3) - 3, (p_end_pos + 1) * 3
            else:
                c_pos = ((p_start_pos + 1) * 3) - 3, p_end_pos * 3

        return {
            "c_ac": c_ac,
            "c_start_pos": c_pos[0],
            "c_end_pos": c_pos[1],
            "cds_start": cds_start,
            "coordinate_type": CoordinateType.INTER_RESIDUE.value,
        }, None

    async def _get_cds_start(self, c_ac: str) -> tuple[int | None, str | None]:
        """Get CDS start for a given cDNA RefSeq accession

        :param c_ac: cDNA RefSeq accession
        :return: Tuple containing:
            - CDS start site if found. Else ``None``
            - Warning, if unable to get CDS start. Else ``None``
        """
        cds_start_end = await self.uta_db.get_cds_start_end(c_ac)
        if not cds_start_end:
            cds_start = None
            warning = f"Accession {c_ac} not found in UTA db"
        else:
            cds_start = cds_start_end[0]
            warning = None
        return cds_start, warning

    async def c_to_g(
        self,
        c_ac: str,
        c_start_pos: int,
        c_end_pos: int,
        cds_start: int | None = None,
        coordinate_type: CoordinateType = CoordinateType.RESIDUE,
        target_genome_assembly: bool = Assembly.GRCH38,
    ) -> tuple[dict | None, str | None]:
        """Translate cDNA representation to genomic representation

        :param c_ac: cDNA RefSeq accession
        :param c_start_pos: cDNA start position for codon
        :param c_end_pos: cDNA end position for codon
        :param coding_start_site: Coding start site. If not provided, this will be
            computed.
        :param target_genome_assembly: Genome assembly to get genomic data for
        :return: Tuple containing:

        * Genomic representation (ac, positions) if able to translate. Will return
          positions as inter-residue coordinates. Else ``None``.
        * Warning, if unable to translate to genomic representation. Else ``None``
        """
        if any(
            (
                c_start_pos == c_end_pos,
                (coordinate_type == CoordinateType.INTER_RESIDUE)
                and ((c_end_pos - c_start_pos) % 3 != 0),
                (coordinate_type == CoordinateType.RESIDUE)
                and ((c_end_pos - (c_start_pos - 1)) % 3 != 0),
            )
        ):
            return (
                None,
                "c_start_pos and c_end_pos are not a valid range for the codon(s)",
            )

        warning = None
        g_coords_data = None

        # Get CDS Start if it is not provided
        if not cds_start:
            cds_start, warning = await self._get_cds_start(c_ac)
            if not cds_start:
                return None, warning

        # Change to inter-residue
        if coordinate_type == CoordinateType.RESIDUE:
            c_start_pos -= 1

        # Get aligned genomic and transcript data
        genomic_tx_data = await self.uta_db.get_genomic_tx_data(
            c_ac,
            (c_start_pos + cds_start, c_end_pos + cds_start),
            AnnotationLayer.CDNA,
            target_genome_assembly=target_genome_assembly,
        )

        if not genomic_tx_data:
            warning = (
                f"Unable to find genomic and transcript data for {c_ac} at "
                f"position ({c_start_pos}, {c_end_pos})"
            )
        else:
            alt_ac = genomic_tx_data.alt_ac

            # Validate that genomic accession assembly == target_genome_assembly
            aliases, _ = self.seqrepo_access.translate_identifier(alt_ac)
            if aliases:
                grch_aliases = [a for a in aliases if a.startswith("GRCh")]
                if not grch_aliases:
                    warning = f"Unable to find associated assembly for {alt_ac}"
                else:
                    found_assembly = grch_aliases[0].split(":")[0]
                    if found_assembly != target_genome_assembly:
                        warning = (
                            f"{alt_ac} uses {found_assembly} assembly which "
                            f"does not not match the target assembly, "
                            f"{target_genome_assembly}"
                        )
                    else:
                        g_pos = genomic_tx_data.alt_pos_change_range

                        # start pos should be less than end pos in response
                        if g_pos[0] > g_pos[1]:
                            g_start_pos = g_pos[1]
                            g_end_pos = g_pos[0]
                        else:
                            g_start_pos = g_pos[0]
                            g_end_pos = g_pos[1]

                        g_coords_data = {
                            "g_ac": alt_ac,
                            "g_start_pos": g_start_pos,
                            "g_end_pos": g_end_pos,
                            "coordinate_type": CoordinateType.INTER_RESIDUE.value,
                        }
            else:
                warning = (
                    f"Unable to validate {alt_ac} matches the target assembly,"
                    f" {target_genome_assembly}"
                )

        return g_coords_data, warning

    async def p_to_g(
        self,
        p_ac: str,
        p_start_pos: int,
        p_end_pos: int,
        coordinate_type: CoordinateType = CoordinateType.INTER_RESIDUE,
        target_genome_assembly: Assembly = Assembly.GRCH38,
    ) -> tuple[dict | None, str | None]:
        """Translate protein representation to genomic representation, by way of
        intermediary conversion into cDNA coordinates.

        :param p_ac: Protein RefSeq accession
        :param p_start_pos: Protein start position
        :param p_end_pos: Protein end position
        :param coordinate_type: Coordinate type for ``p_start_pos`` and ``p_end_pos``.
        :param target_genome_assembly: Genome assembly to get genomic data for
        :return: Tuple containing:

        * Genomic representation (ac, positions) if able to translate. Will return
          positions as inter-residue coordinates. Else ``None``.
        * Warnings, if conversion to cDNA or genomic coordinates fails.
        """
        c_data, warning = await self.p_to_c(
            p_ac, p_start_pos, p_end_pos, coordinate_type=coordinate_type
        )
        if not c_data:
            return None, warning

        # p_to_c returns c_data as inter-residue
        g_data, warning = await self.c_to_g(
            c_data["c_ac"],
            c_data["c_start_pos"],
            c_data["c_end_pos"],
            c_data["cds_start"],
            coordinate_type=CoordinateType.INTER_RESIDUE,
            target_genome_assembly=target_genome_assembly,
        )
        return g_data, warning
