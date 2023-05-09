"""Module containing alignment methods for translating to and from different
reference sequences.
"""
from typing import Optional, Tuple, Dict

from cool_seq_tool.schemas import AnnotationLayer, Assembly, ResidueMode
from cool_seq_tool.data_sources import SeqRepoAccess, TranscriptMappings, \
    UTADatabase
from cool_seq_tool.utils.positions import zero_based_to_inter_residue, to_zero_based
from cool_seq_tool.utils.validation import validate_index


class AlignmentMapperError(Exception):
    """Custom exception for AlignmentMapper class"""


class AlignmentMapper:
    """Class for translating between p -> c -> g reference sequences."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 uta_db: UTADatabase) -> None:
        """Initialize the AlignmentMapper class.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo queries
        :param TranscriptMappings transcript_mappings: Access to transcript
            accession mappings and conversions
        :param UTADatabase uta_db: UTADatabase instance to give access to query
            UTA database
        """
        self.seqrepo_access = seqrepo_access
        self.transcript_mappings = transcript_mappings
        self.uta_db = uta_db

    @staticmethod
    def p_pos_to_c_pos(
        p_start_pos: int, p_end_pos: int,
        residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Tuple[int, int]:
        """Convert protein position to cDNA position on a codon(s)

        :param int p_start_pos: Protein start position
        :param int p_end_pos: Protein end position
        :param ResidueMode residue_mode: Residue mode for `p_start_pos` and `p_end_pos`
        :return: Tuple containing:
            - cDNA start position on codon
            - cDNA end position on codon
            Inter-residue coordinates will be returned
        """
        # For 1-based, protein pos * 3 = end of codon
        # 1 amino acid maps to 3 nucleotides in the codon
        # Since we have the end of the codon, we will subtract 2 to get the start of the
        # codon. We want to return inter-residue (0-based), so we subtract 1 from this.
        if residue_mode == ResidueMode.RESIDUE:
            c_pos = (p_start_pos * 3) - 3, p_end_pos * 3
        else:
            if p_start_pos == p_end_pos:
                c_pos = ((p_start_pos + 1) * 3) - 3, (p_end_pos + 1) * 3
            else:
                c_pos = ((p_start_pos + 1) * 3) - 3, p_end_pos * 3
        return c_pos

    async def p_to_c(
        self, p_ac: str, p_start_pos: int, p_end_pos: int,
        residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Tuple[Optional[Dict], Optional[str]]:
        """Translate protein representation to cDNA representation.

        :param str p_ac: Protein accession
        :param int p_start_pos: Protein start position
        :param int p_end_pos: Protein end position
        :param ResidueMode residue_mode: Residue mode for `p_start_pos` and `p_end_pos`
        :return: Tuple containing:
            - cDNA representation (accession, codon range positions for corresponding
              change, cds start site) if able to translate. Will return positions as
              inter-residue coordinates. If unable to translate, returns `None`.
            - Warning, if unable to translate to cDNA representation. Else `None`
        """
        # Get cDNA accession
        temp_c_ac = await self.uta_db.p_to_c_ac(p_ac)
        if temp_c_ac:
            c_ac = temp_c_ac[-1]
        else:
            try:
                if p_ac.startswith("NP_"):
                    c_ac = self.transcript_mappings.np_to_nm[p_ac]
                elif p_ac.startswith("ENSP"):
                    c_ac = self.transcript_mappings.ensp_to_enst[p_ac]
                else:
                    return None, f"{p_ac} is not supported"
            except KeyError:
                return None, f"{p_ac} not found in transcript mappings"

        # Get coding start site
        cds_start, warning = await self._get_cds_start(c_ac)
        if not cds_start:
            return None, warning

        c_pos = self.p_pos_to_c_pos(p_start_pos, p_end_pos, residue_mode)

        return {
            "c_ac": c_ac,
            "c_start_pos": c_pos[0],
            "c_end_pos": c_pos[1],
            "cds_start": cds_start,
            "residue_mode": ResidueMode.INTER_RESIDUE.value
        }, None

    async def _get_cds_start(self, c_ac: str) -> Tuple[Optional[int], Optional[str]]:
        """Get CDS start for a given cDNA RefSeq accession

        :param str c_ac: cDNA RefSeq accession
        :return: Tuple containing:
            - CDS start site if found. Else `None`
            - Warning, if unable to get CDS start. Else `None`
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
        self, c_ac: str, c_start_pos: int, c_end_pos: int,
        cds_start: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
        target_genome_assembly: bool = Assembly.GRCH38,
    ) -> Tuple[Optional[Dict], Optional[str]]:
        """Translate cDNA representation to genomic representation

        :param str c_ac: cDNA RefSeq accession
        :param int c_start_pos: cDNA start position
        :param int c_end_pos: cDNA end position
        :param Optional[int] coding_start_site: Coding start site. If not provided,
            this will be computed.
        :param Assembly target_genome_assembly: Genome assembly to get genomic data for
        :return: Tuple containing:
            - Genomic representation (ac, positions) if able to translate. Will return
              positions as inter-residue coordinates. Else `None`.
            - Warning, if unable to translate to genomic representation. Else `None`
        """
        warning = None
        g_coords_data = None

        # Validate ENST accessions (since uta does not have versions)
        if c_ac.startswith("ENST"):
            if not any((
                self.transcript_mappings.ensembl_transcript_version_to_gene_symbol.get(c_ac),  # noqa: E501
                self.seqrepo_access.get_reference_sequence(c_ac, 1, 2)[0]
            )):
                return None, f"Ensembl transcript not found: {c_ac}"

        # Get CDS Start if it is not provided
        if not cds_start:
            cds_start, warning = await self._get_cds_start(c_ac)
            if not cds_start:
                return None, warning

        # Change to inter-residue
        if residue_mode == ResidueMode.RESIDUE:
            c_start_pos -= 1

        # Get aligned genomic and transcript data
        genomic_tx_data = await self.uta_db.get_genomic_tx_data(
            c_ac, (c_start_pos + cds_start, c_end_pos + cds_start),
            AnnotationLayer.CDNA, target_genome_assembly=target_genome_assembly)

        if not genomic_tx_data:
            warning = f"Unable to find genomic and transcript data for {c_ac} at "\
                      f"position ({c_start_pos}, {c_end_pos})"
        else:
            alt_ac = genomic_tx_data["alt_ac"]

            # Validate that genomic accession assembly == target_genome_assembly
            aliases, _ = self.seqrepo_access.translate_identifier(alt_ac)
            if aliases:
                grch_aliases = [a for a in aliases if a.startswith("GRCh")]
                if not grch_aliases:
                    warning = f"Unable to find associated assembly for {alt_ac}"
                else:
                    found_assembly = grch_aliases[0].split(":")[0]
                    g_pos = genomic_tx_data["alt_pos_change_range"]
                    # start pos should be less than end pos in response
                    if g_pos[0] > g_pos[1]:
                        g_start_pos = g_pos[1]
                        g_end_pos = g_pos[0]
                    else:
                        g_start_pos = g_pos[0]
                        g_end_pos = g_pos[1]

                    if found_assembly != target_genome_assembly:
                        if target_genome_assembly == Assembly.GRCH38:
                            # Need to liftover
                            grch38 = await self.g_to_grch38(alt_ac, g_start_pos, g_end_pos)
                            g_coords_data = {
                                "g_ac": grch38["ac"],
                                "g_start_pos": grch38["pos"][0],
                                "g_end_pos": grch38["pos"][1],
                                "gene": genomic_tx_data["gene"],
                                "residue_mode": ResidueMode.INTER_RESIDUE.value,
                                "strand": genomic_tx_data["strand"]
                            }
                    else:
                        g_coords_data = {
                            "g_ac": alt_ac,
                            "g_start_pos": g_start_pos,
                            "g_end_pos": g_end_pos,
                            "gene": genomic_tx_data["gene"],
                            "residue_mode": ResidueMode.INTER_RESIDUE.value,
                            "strand": genomic_tx_data["strand"]
                        }
            else:
                warning = f"Unable to validate {alt_ac} matches the target assembly,"\
                          f" {target_genome_assembly}"

        return g_coords_data, warning

    async def p_to_g(
        self, p_ac: str, p_start_pos: int, p_end_pos: int,
        residue_mode: ResidueMode = ResidueMode.INTER_RESIDUE,
        target_genome_assembly: Assembly = Assembly.GRCH38
    ) -> Tuple[Optional[Dict], Optional[str]]:
        """Translate protein representation to genomic representation

        :param str p_ac: Protein RefSeq accession
        :param int p_start_pos: Protein start position
        :param int p_end_pos: Protein end position
        :param ResidueMode residue_mode: Residue mode for `p_start_pos` and `p_end_pos`.
        :param Assembly target_genome_assembly: Genome assembly to get genomic data for
        :return: Tuple containing:
            - Genomic representation (ac, positions) if able to translate. Will return
              positions as inter-residue coordinates. Else `None`.
            and warnings. The genomic data will always return inter-residue coordinates
        """
        c_data, warning = await self.p_to_c(p_ac, p_start_pos, p_end_pos,
                                            residue_mode=residue_mode)
        if not c_data:
            return None, warning

        # p_to_c returns c_data as inter-residue
        g_data, warning = await self.c_to_g(
            c_data["c_ac"], c_data["c_start_pos"], c_data["c_end_pos"],
            c_data["cds_start"], residue_mode=ResidueMode.INTER_RESIDUE,
            target_genome_assembly=target_genome_assembly)
        return g_data, warning

    async def g_to_c(
        self, g_ac: str, g_start_pos: int, g_end_pos: int, tx_ac: str,
        cds_start: Optional[int] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
        is_zero_based: bool = False
    ) -> Tuple[Optional[Dict], Optional[str]]:
        begin_residue_mode = residue_mode
        begin_start_eq_end = g_start_pos == g_end_pos

        if not is_zero_based:
            g_start_pos, g_end_pos = to_zero_based(g_start_pos, g_end_pos)

        genomic_tx_data = await self.uta_db.get_genomic_tx_data(
            tx_ac, (g_start_pos, g_end_pos), AnnotationLayer.GENOMIC, g_ac)
        if not genomic_tx_data:
            return None, "Unable to get aligned transcript and genomic data"
        strand = genomic_tx_data["strand"]

        # transcript and genomic ranges for an exon
        tx_pos_range = genomic_tx_data["tx_pos_range"]
        alt_pos_range = genomic_tx_data["alt_pos_range"]

        if strand == "-":
            # TODO: Why is this necessary?
            pos_change = ((g_start_pos + 1) - alt_pos_range[0],
                          alt_pos_range[1] - (g_end_pos + 1))
        else:
            pos_change = g_start_pos - alt_pos_range[0], alt_pos_range[1] - g_end_pos

        if not cds_start:
            cds_start, warning = await self._get_cds_start(genomic_tx_data["tx_ac"])
            if not cds_start:
                return None, warning

        if strand == "-":
            tx_pos = tx_pos_range[0] + pos_change[1], tx_pos_range[1] - pos_change[0]
        else:
            tx_pos = tx_pos_range[0] + pos_change[0], tx_pos_range[1] - pos_change[1]
        tx_pos = tx_pos[0] - cds_start, tx_pos[1] - cds_start
        if tx_pos[0] > tx_pos[1]:
            tx_pos = tx_pos[1], tx_pos[0]

        # change back to inter-residue coords
        tx_pos = zero_based_to_inter_residue(
            tx_pos[0], tx_pos[1], begin_start_eq_end, begin_residue_mode)

        return {
            "c_ac": tx_ac,
            "c_start_pos": tx_pos[0],
            "c_end_pos": tx_pos[1],
            "cds_start": cds_start,
            "residue_mode": ResidueMode.INTER_RESIDUE.value,
            "strand": strand
        }, None

    async def g_to_grch38(
        self, ac: str, start_pos: int, end_pos: int
    ) -> Optional[Dict]:
        """Liftover genomic representation to GRCh38 assembly. If GRCh38 representation
        is provided, will perform validation checks. Only supports lifting over from
        GRCh37 assembly.

        :param str ac: Genomic accession
        :param int start_pos: Genomic start position change
        :param int end_pos: Genomic end position change
        :raise ValidationError: If validation checks do not pass
        :raise AlignmentMapperError: If error happens during alignment
        :return: RefSeq genomic NC accession, start and end pos on GRCh38 assembly
        """
        # Checking to see what chromosome and assembly we're on
        descr = await self.uta_db.get_chr_assembly(ac)
        if not descr:
            # Already GRCh38 assembly
            # Will raise ValidationError if checks do not pass
            validate_index(self.seqrepo_access, ac, (start_pos, end_pos), 0)
            return {
                "ac": ac,
                "pos": (start_pos, end_pos)
            }

        chromosome, assembly = descr
        is_same_pos = start_pos == end_pos

        # Coordinate liftover
        if assembly < "GRCh37":
            raise AlignmentMapperError("Liftover from GRCh37 is only supported")

        # TODO: We should probably raise an exception in UTA DB?
        liftover_start_i = self.uta_db.get_liftover(chromosome, start_pos,
                                                    Assembly.GRCH38)
        if liftover_start_i is None:
            return None
        else:
            start_pos = liftover_start_i[1]

        if not is_same_pos:
            liftover_end_i = self.uta_db.get_liftover(chromosome, end_pos,
                                                      Assembly.GRCH38)
            if liftover_end_i is None:
                return None
            else:
                end_pos = liftover_end_i[1]
        else:
            end_pos = start_pos

        newest_ac = await self.uta_db.get_newest_assembly_ac(ac)
        if newest_ac:
            ac = newest_ac[0]
            # Will raise ValidationError if checks do not pass
            validate_index(self.seqrepo_access, ac, (start_pos, end_pos), 0)
            return {
                "ac": ac,
                "pos": (start_pos, end_pos)
            }

        return None
