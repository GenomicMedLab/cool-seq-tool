"""Module for retrieving MANE Transcript from variation on p/c/g coordinate.
Steps:
1. Map annotation layer to genome
2. Liftover to preferred genome
    We want to liftover to GRCh38. We do not support getting MANE transcripts
    for GRCh36 and earlier assemblies.
3. Select preferred compatible annotation
4. Map back to correct annotation layer
"""
import math
from typing import Optional, Set, Tuple, Dict, List, Union

import hgvs.parser
import pandas as pd

from cool_seq_tool.schemas import AnnotationLayer, Assembly, MappedManeData, \
    ResidueMode, TranscriptPriorityLabel
from cool_seq_tool.data_sources import SeqRepoAccess, TranscriptMappings, \
    MANETranscriptMappings, UTADatabase, GeneNormalizer
from cool_seq_tool.data_sources.residue_mode import get_inter_residue_pos
from cool_seq_tool import logger


class MANETranscriptError(Exception):
    """Custom exception for MANETranscript class"""

    pass


class MANETranscript:
    """Class for retrieving MANE transcripts."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 mane_transcript_mappings: MANETranscriptMappings,
                 uta_db: UTADatabase,
                 gene_normalizer: GeneNormalizer) -> None:
        """Initialize the MANETranscript class.

        :param SeqRepoAccess seqrepo_access: Access to seqrepo queries
        :param TranscriptMappings transcript_mappings: Access to transcript
            accession mappings and conversions
        :param MANETranscriptMappings mane_transcript_mappings: Access to
            MANE Transcript accession mapping data
        :param UTADatabase uta_db: UTADatabase instance to give access to query
            UTA database
        :param GeneNormalizer gene_normalizer: Access to Gene Normalizer
        """
        self.seqrepo_access = seqrepo_access
        self.hgvs_parser = hgvs.parser.Parser()
        self.transcript_mappings = transcript_mappings
        self.mane_transcript_mappings = mane_transcript_mappings
        self.uta_db = uta_db
        self.gene_normalizer = gene_normalizer

    @staticmethod
    def _get_reading_frame(pos: int) -> int:
        """Return reading frame number.
        Only used on c. coordinate

        :param int pos: cDNA position
        :return: Reading frame
        """
        pos_mod_3 = pos % 3
        if pos_mod_3 == 0:
            pos_mod_3 = 3
        return pos_mod_3

    @staticmethod
    def _p_to_c_pos(start: int, end: int) -> Tuple[int, int]:
        """Return cDNA position given a protein position.

        :param int start: Start protein position
        :param int end: End protein position
        :return: cDNA position start, cDNA position end
        """
        start_pos = start * 3 - 1
        if end != start:
            end_pos = end * 3 - 1
        else:
            end_pos = start_pos

        return start_pos - 1, end_pos + 1

    async def _p_to_c(self, ac: str, start_pos: int,
                      end_pos: int) -> Optional[Tuple[str, Tuple[int, int]]]:
        """Convert protein (p.) annotation to cDNA (c.) annotation.

        :param str ac: Protein accession
        :param int start_pos: Protein start position
        :param int end_pos: Protein end position
        :return: [cDNA transcript accession, [cDNA pos start, cDNA pos end]]
        """
        # TODO: Check version mappings 1 to 1 relationship
        temp_ac = await self.uta_db.p_to_c_ac(ac)
        if temp_ac:
            ac = temp_ac[-1]
        else:
            try:
                if ac.startswith("NP_"):
                    ac = self.transcript_mappings.np_to_nm[ac]
                elif ac.startswith("ENSP"):
                    ac = \
                        self.transcript_mappings.ensp_to_enst[ac]
                else:
                    logger.warning(f"Unable to find accession: {ac}")
                    return None
            except KeyError:
                logger.warning(f"{ac} not found in transcript_mappings")
                return None

        pos = self._p_to_c_pos(start_pos, end_pos)
        return ac, pos

    async def _c_to_g(self, ac: str, pos: Tuple[int, int]) -> Optional[Dict]:
        """Get g. annotation from c. annotation.

        :param str ac: cDNA accession
        :param Tuple[int, int] pos: [cDNA pos start, cDNA pos end]
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
        """
        # UTA does not store ENST versions
        # So we want to make sure version is valid
        if ac.startswith("ENST"):
            if not self.transcript_mappings.ensembl_transcript_version_to_gene_symbol.get(ac):  # noqa: E501
                if not self.seqrepo_access.get_reference_sequence(ac, 1)[0]:
                    logger.warning(f"Ensembl transcript not found: {ac}")
                    return None

            temp_ac = ac.split(".")[0]
        else:
            temp_ac = ac

        # c. coordinate does not contain cds start, so we need to add it
        cds_start_end = await self.uta_db.get_cds_start_end(temp_ac)
        if not cds_start_end:
            logger.warning(f"Accession {temp_ac} not found in UTA")
            return None
        coding_start_site = cds_start_end[0]
        pos = pos[0] + coding_start_site, pos[1] + coding_start_site

        genomic_tx_data = await self._get_and_validate_genomic_tx_data(
            ac, pos, AnnotationLayer.CDNA, coding_start_site=coding_start_site)
        return genomic_tx_data

    async def _get_and_validate_genomic_tx_data(
        self, tx_ac: str, pos: Tuple[int, int],
        annotation_layer: Union[AnnotationLayer.CDNA, AnnotationLayer.GENOMIC] = AnnotationLayer.CDNA,  # noqa: E501
        coding_start_site: Optional[int] = None,
        alt_ac: Optional[str] = None
    ) -> Optional[Dict]:
        """Get and validate genomic_tx_data

        :param str tx_ac: Accession on c. coordinate
        :param Tuple[int, int] pos: (start pos, end pos)
        :param Union[AnnotationLayer.CDNA, AnnotationLayer.GENOMIC] annotation_layer:
            Annotation layer for `ac` and `pos`
        :param Optional[int] coding_start_site: Coding start site
        :param Optional[str] alt_ac: Accession on g. coordinate
        :return: genomic_tx_data if found and validated, else None
        """
        genomic_tx_data = await self.uta_db.get_genomic_tx_data(
            tx_ac, pos, annotation_layer, alt_ac=alt_ac)
        if not genomic_tx_data:
            logger.warning(f"Unable to find genomic_tx_data for {alt_ac} at position"
                           f" {pos} on annotation layer {annotation_layer}")
            return None
        genomic_tx_data["coding_start_site"] = coding_start_site

        if not alt_ac:
            # Only want to liftover if alt_ac not provided. If alt_ac is provided,
            # it means user wants to stick with the queried assembly
            og_alt_exon_id = genomic_tx_data["alt_exon_id"]
            await self.uta_db.liftover_to_38(genomic_tx_data)
            liftover_alt_exon_id = genomic_tx_data["alt_exon_id"]

            # Validation check: Exon structure
            if og_alt_exon_id != liftover_alt_exon_id:
                logger.warning(f"Original alt_exon_id {og_alt_exon_id} "
                               f"does not match liftover alt_exon_id "
                               f"{liftover_alt_exon_id}")
                return None

        return genomic_tx_data

    @staticmethod
    def _get_c_data(
        gene: str, cds_start_end: Tuple[int, int], c_pos_change: Tuple[int, int],
        strand: str, status: TranscriptPriorityLabel, refseq_c_ac: str,
        ensembl_c_ac: Optional[str] = None, alt_ac: Optional[str] = None
    ) -> Dict:
        """Return transcript data on c. coordinate.

        :param str gene: Gene symbol
        :param Tuple[int, int] cds_start_end: Coding start and end site
            for transcript
        :param Tuple[int, int] c_pos_change: Start and end positions
            for change on c. coordinate
        :param str strand: Strand
        :param TranscriptPriorityLabel status: Status of transcript
        :param str refseq_c_ac: Refseq transcript
        :param Optional[str] ensembl_c_ac: Ensembl transcript
        :param Optional[str] alt_ac: Genomic accession
        :return: Transcript data on c. coord
        """
        cds_start = cds_start_end[0]
        cds_end = cds_start_end[1]
        lt_cds_start = (c_pos_change[0] < cds_start and c_pos_change[1] < cds_start)
        gt_cds_end = (c_pos_change[1] > cds_end and c_pos_change[1] > cds_end)

        if lt_cds_start or gt_cds_end:
            logger.info(f"{refseq_c_ac} with position"
                        f" {c_pos_change} is not within CDS start/end")
        return dict(
            gene=gene,
            refseq=refseq_c_ac,
            ensembl=ensembl_c_ac,
            coding_start_site=cds_start,
            coding_end_site=cds_end,
            pos=c_pos_change,
            strand=strand,
            status=status,
            alt_ac=alt_ac
        )

    @staticmethod
    def _get_mane_p(mane_data: Dict,
                    mane_c_pos_range: Tuple[int, int]) -> Dict:
        """Translate MANE Transcript c. annotation to p. annotation

        :param Dict mane_data: MANE Transcript data
        :param Tuple[int, int] mane_c_pos_range: Position change range
            on MANE Transcript c. coordinate
        :return: MANE transcripts accessions and position change on
            p. coordinate
        """
        return dict(
            gene=mane_data["symbol"],
            refseq=mane_data["RefSeq_prot"],
            ensembl=mane_data["Ensembl_prot"],
            pos=(math.ceil(mane_c_pos_range[0] / 3),
                 math.floor(mane_c_pos_range[1] / 3)),
            strand=mane_data["chr_strand"],
            status="_".join(mane_data["MANE_status"].split()).lower()
        )

    async def _g_to_c(
        self, g: Dict, refseq_c_ac: str, status: TranscriptPriorityLabel,
        ensembl_c_ac: Optional[str] = None, alt_ac: Optional[str] = None,
        found_result: bool = False
    ) -> Optional[Dict]:
        """Get transcript c. annotation data from g. annotation.

        :param Dict g: Genomic data
        :param str refseq_c_ac: Refseq transcript accession
        :param TranscriptPriorityLabel status: Status of transcript
        :param Optional[str] ensembl_c_ac: Ensembl transcript accession
        :param Optional[str] alt_ac: Genomic accession
        :param bool found_result: `True` if found result, so do not need to query
            tx_exon_aln_v table. This is because the user did not need to liftover.
            `False` if need to get result from tx_exon_aln_v table.
        :return: Transcript data
        """
        if found_result:
            tx_g_pos = g["alt_pos_range"]
            tx_pos_range = g["tx_pos_range"]
        else:
            result = await self.uta_db.get_tx_exon_aln_v_data(
                refseq_c_ac, g["alt_pos_change_range"][0],
                g["alt_pos_change_range"][1], alt_ac=alt_ac if alt_ac else g["alt_ac"],
                use_tx_pos=False)

            if not result:
                logger.warning(f"Unable to find transcript, {refseq_c_ac}, "
                               f"position change")
                return None
            else:
                result = result[-1]
                tx_g_pos = result[5], result[6]  # alt_start_i, alt_end_i
                tx_pos_range = result[2], result[3]  # tx_start_i, tx_end_i

        cds_start_end = await self.uta_db.get_cds_start_end(refseq_c_ac)
        if not cds_start_end:
            return None
        coding_start_site = cds_start_end[0]

        g_pos = g["alt_pos_change_range"]  # start/end genomic change
        g_pos_change = g_pos[0] - tx_g_pos[0], tx_g_pos[1] - g_pos[1]

        if g["strand"] == "-":
            g_pos_change = (
                tx_g_pos[1] - g_pos[0], g_pos[1] - tx_g_pos[0]
            )

        c_pos_change = (
            tx_pos_range[0] + g_pos_change[0] - coding_start_site,
            tx_pos_range[1] - g_pos_change[1] - coding_start_site
        )

        if c_pos_change[0] > c_pos_change[1]:
            c_pos_change = c_pos_change[1], c_pos_change[0]

        return self._get_c_data(
            gene=g["gene"],
            cds_start_end=cds_start_end,
            c_pos_change=c_pos_change,
            strand=g["strand"],
            alt_ac=alt_ac,
            status=status,
            refseq_c_ac=refseq_c_ac,
            ensembl_c_ac=ensembl_c_ac)

    def _validate_reading_frames(self, ac: str, start_pos: int, end_pos: int,
                                 transcript_data: Dict) -> bool:
        """Return whether reading frames are the same after translation.

        :param str ac: Query accession
        :param int start_pos: Original start cDNA position change
        :param int end_pos: Original end cDNA position change
        :param Dict transcript_data: Ensembl and RefSeq transcripts with
            corresponding position change
        :return: `True` if reading frames are the same after translation.
            `False` otherwise
        """
        for pos, pos_index in [(start_pos, 0), (end_pos, 1)]:
            if pos is not None:
                og_rf = self._get_reading_frame(pos)
                new_rf = self._get_reading_frame(transcript_data["pos"][pos_index])

                if og_rf != new_rf:
                    logger.warning(f"{ac} original reading frame ({og_rf}) "
                                   f"does not match new "
                                   f"{transcript_data['ensembl']}, "
                                   f"{transcript_data['refseq']} reading "
                                   f"frame ({new_rf})")
                    return False
            else:
                if pos_index == 0:
                    logger.warning(f"{ac} must having start position")
                    return False
        return True

    def _validate_references(self, ac: str, coding_start_site: int,
                             start_pos: int, end_pos: int,
                             mane_transcript: Dict, expected_ref: str,
                             anno: AnnotationLayer, residue_mode: str) -> bool:
        """Return whether or not reference changes are the same.

        :param str ac: Query accession
        :param int coding_start_site: ac's coding start site
        :param int start_pos: Original start position change
        :param int end_pos: Origin end position change
        :param Dict mane_transcript: Ensembl and RefSeq transcripts with
            corresponding position change
        :param str expected_ref: Reference at position given during input
        :param AnnotationLayer anno: Annotation layer we are starting from
        :param ResidueMode residue_mode: Residue mode
        :return: `True` if reference check passes. `False` otherwise.
        """
        if anno == AnnotationLayer.CDNA:
            start_pos += coding_start_site
            end_pos += coding_start_site

        ref, warnings = self.seqrepo_access.get_reference_sequence(
            ac, start_pos, end=end_pos, residue_mode=residue_mode
        )
        if ref is None:
            return False

        if mane_transcript:
            mane_start_pos = mane_transcript["pos"][0]
            mane_end_pos = mane_transcript["pos"][1]
            if anno == "c":
                mane_cds = mane_transcript["coding_start_site"]
                mane_start_pos += mane_cds
                mane_end_pos += mane_cds
            mane_ref, warnings = self.seqrepo_access.get_reference_sequence(
                mane_transcript["refseq"],
                mane_start_pos,
                end=mane_end_pos if mane_start_pos != mane_end_pos else None,
                residue_mode=residue_mode
            )
            if not mane_ref:
                logger.info("Unable to validate reference for MANE Transcript")

            if expected_ref != mane_ref:
                logger.info(f"Expected ref, {expected_ref}, but got {mane_ref}"
                            f" on MANE accession, {mane_transcript['refseq']}")

        if expected_ref != ref:
            logger.warning(f"Expected ref, {expected_ref}, but got {ref} "
                           f"on accession, {ac}")
            return False

        return True

    def _validate_index(self, ac: str, pos: Tuple[int, int],
                        coding_start_site: int) -> bool:
        """Validate that positions actually exist on accession

        :param str ac: Accession
        :param Tuple[int, int] pos: Start position change, End position change
        :param int coding_start_site: coding start site for accession
        :return: `True` if positions exist on accession. `False` otherwise
        """
        start_pos = pos[0] + coding_start_site
        end_pos = pos[1] + coding_start_site
        if self.seqrepo_access.get_reference_sequence(ac, start_pos, end_pos,
                                                      residue_mode=ResidueMode.INTER_RESIDUE)[0]:  # noqa E501
            return True
        else:
            return False

    def _get_prioritized_transcripts_from_gene(self,
                                               df: pd.core.frame.DataFrame) -> List:
        """Sort and filter transcripts from gene to get priority list

        :param pd.core.frame.DataFrame df: Data frame containing transcripts from gene
            data
        :return: List of prioritized transcripts for a given gene. Sort by latest
            assembly, longest length of transcript, with first-published transcripts
            breaking ties. If there are multiple transcripts for a given accession, the
            most recent version of a transcript associated with an assembly will be kept
        """
        copy_df = df.copy(deep=True)
        copy_df = copy_df.drop(columns="alt_ac").drop_duplicates()
        copy_df["ac_no_version_as_int"] = copy_df["tx_ac"].apply(lambda x: int(x.split(".")[0].split("NM_00")[1]))  # noqa: E501
        copy_df["ac_version"] = copy_df["tx_ac"].apply(lambda x: x.split(".")[1])
        copy_df = copy_df.sort_values(["ac_no_version_as_int", "ac_version"],
                                      ascending=[False, False])
        copy_df = copy_df.drop_duplicates(["ac_no_version_as_int"], keep="first")
        copy_df.loc[:, "len_of_tx"] = copy_df.loc[:, "tx_ac"].apply(lambda ac: len(self.seqrepo_access.get_reference_sequence(ac)[0]))  # noqa: E501
        copy_df = copy_df.sort_values(
            ["len_of_tx", "ac_no_version_as_int"], ascending=[False, True])
        return list(copy_df["tx_ac"])

    async def get_longest_compatible_transcript(
            self, gene: str, start_pos: int, end_pos: int,
            start_annotation_layer: AnnotationLayer, ref: Optional[str] = None,
            residue_mode: str = ResidueMode.RESIDUE,
            mane_transcripts: Optional[Set] = None,
            alt_ac: Optional[str] = None
    ) -> Optional[Dict]:
        """Get longest compatible transcript from a gene.
        Try GRCh38 first, then GRCh37.
        Transcript is compatible if it passes validation checks.

        :param str gene: Gene symbol
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param  AnnotationLayer start_annotation_layer: Starting annotation layer.
        :param str ref: Reference at position given during input
        :param str residue_mode: Residue mode
        :param Optional[Set] mane_transcripts: Attempted mane transcripts that were not
            compatible
        :param Optional[str] alt_ac: Genomic accession
        :return: Data for longest compatible transcript
        """
        inter_residue_pos, _ = get_inter_residue_pos(
            start_pos, residue_mode, end_pos=end_pos)
        if not inter_residue_pos:
            return None
        residue_mode = ResidueMode.INTER_RESIDUE
        start_pos, end_pos = inter_residue_pos

        is_p_or_c_start_anno = True
        if start_annotation_layer == AnnotationLayer.PROTEIN:
            c_start_pos, c_end_pos = self._p_to_c_pos(start_pos, end_pos)
        elif start_annotation_layer == AnnotationLayer.CDNA:
            c_start_pos, c_end_pos = start_pos, end_pos
        else:
            is_p_or_c_start_anno = False

        # Data Frame that contains transcripts associated to a gene
        if is_p_or_c_start_anno:
            df = await self.uta_db.get_transcripts_from_gene(
                gene, c_start_pos, c_end_pos, use_tx_pos=True, alt_ac=alt_ac)
        else:
            df = await self.uta_db.get_transcripts_from_gene(
                gene, start_pos, end_pos, use_tx_pos=False, alt_ac=alt_ac)
        if df.empty:
            logger.warning(f"Unable to get transcripts from gene {gene}")
            return None

        prioritized_tx_acs = self._get_prioritized_transcripts_from_gene(df)

        if mane_transcripts:
            # Dont check MANE transcripts since we know that are not compatible
            prioritized_tx_acs = [el for el in prioritized_tx_acs
                                  if el not in mane_transcripts]

        for tx_ac in prioritized_tx_acs:
            # Only need to check the one row since we do liftover in _c_to_g
            tmp_df = df.loc[df["tx_ac"] == tx_ac].sort_values("alt_ac", ascending=False)
            row = tmp_df.iloc[0]

            if alt_ac is None:
                alt_ac = row["alt_ac"]

            found_tx_exon_aln_v_result = False
            if is_p_or_c_start_anno:
                # Go from c -> g annotation (liftover as well)
                g = await self._c_to_g(tx_ac, (c_start_pos, c_end_pos))
            else:
                # g -> GRCh38 (if alt_ac not provided. if it is, will use that assembly)
                g = await self._get_and_validate_genomic_tx_data(
                    tx_ac, (start_pos, end_pos),
                    annotation_layer=AnnotationLayer.GENOMIC, alt_ac=alt_ac)
                found_tx_exon_aln_v_result = True
            if not g:
                continue

            # Get prioritized transcript data for gene
            # grch38 -> c
            lcr_c_data = await self._g_to_c(
                g=g, refseq_c_ac=tx_ac,
                status=TranscriptPriorityLabel.LongestCompatibleRemaining.value,
                found_result=found_tx_exon_aln_v_result)

            if not lcr_c_data:
                continue

            # Validation checks
            if is_p_or_c_start_anno:
                validate_reading_frame = self._validate_reading_frames(
                    tx_ac, c_start_pos, c_end_pos, lcr_c_data)
                if not validate_reading_frame:
                    continue

            if ref:
                if start_annotation_layer == AnnotationLayer.PROTEIN:
                    valid_references = self._validate_references(
                        row["pro_ac"], row["cds_start_i"], start_pos,
                        end_pos, {}, ref, AnnotationLayer.PROTEIN, residue_mode)
                elif start_annotation_layer == AnnotationLayer.CDNA:
                    valid_references = self._validate_references(
                        row["tx_ac"], row["cds_start_i"], c_start_pos,
                        c_end_pos, {}, ref, AnnotationLayer.CDNA, residue_mode)
                else:
                    valid_references = self._validate_references(
                        alt_ac, 0, start_pos, end_pos, {}, ref,
                        AnnotationLayer.GENOMIC, residue_mode)

                if not valid_references:
                    continue

            if start_annotation_layer == AnnotationLayer.PROTEIN:
                pos = (math.ceil(lcr_c_data["pos"][0] / 3),
                       math.floor(lcr_c_data["pos"][1] / 3))
                ac = row["pro_ac"]
                coding_start_site = 0
            else:
                # cDNA and Genomic annotations will return c. data
                pos = lcr_c_data["pos"]
                ac = tx_ac
                coding_start_site = lcr_c_data["coding_start_site"]

            if not self._validate_index(ac, pos, coding_start_site):
                logger.warning(f"{pos} are not valid positions on {ac}"
                               f"with coding start site "
                               f"{coding_start_site}")
                continue

            return dict(
                refseq=ac if ac.startswith("N") else None,
                ensembl=ac if ac.startswith("E") else None,  # TODO: issues 87, 4
                pos=pos,
                strand=g["strand"],
                status=lcr_c_data["status"]
            )
        return None

    async def get_mane_transcript(
            self, ac: str, start_pos: int, start_annotation_layer: str,
            end_pos: Optional[int] = None, gene: Optional[str] = None,
            ref: Optional[str] = None, try_longest_compatible: bool = False,
            residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Optional[Dict]:
        """Return mane transcript.

        :param str ac: Accession
        :param int start_pos: Start position change
        :param str start_annotation_layer: Starting annotation layer.
            Must be either `p`, `c`, or `g`.
        :param Optional[int] end_pos: End position change. If `None` assumes
            both  `start_pos` and `end_pos` have same values.
        :param str gene: Gene symbol
        :param str ref: Reference at position given during input
        :param bool try_longest_compatible: `True` if should try longest
            compatible remaining if mane transcript was not compatible.
            `False` otherwise.
        :param ResidueMode residue_mode: Starting residue mode for `start_pos`
            and `end_pos`. Will always return coordinates in inter-residue
        :return: MANE data or longest transcript compatible data if validation
            checks are correct. Will return inter-residue coordinates.
            Else, `None`
        """
        inter_residue_pos, warning = get_inter_residue_pos(
            start_pos, residue_mode, end_pos=end_pos)
        if not inter_residue_pos:
            return None
        start_pos, end_pos = inter_residue_pos
        residue_mode = ResidueMode.INTER_RESIDUE
        if ref:
            ref = ref[:end_pos - start_pos]

        anno = start_annotation_layer.lower()
        if anno in ["p", "c"]:
            # Get accession and position on c. coordinate
            if anno == "p":
                c = await self._p_to_c(ac, start_pos, end_pos)
                if not c:
                    return None
                c_ac, c_pos = c
            else:
                c_ac = ac
                c_pos = start_pos, end_pos
            # Go from c -> g annotation (liftover as well)
            g = await self._c_to_g(c_ac, c_pos)
            if g is None:
                return None
            # Get mane data for gene
            mane_data = self.mane_transcript_mappings.get_gene_mane_data(g["gene"])
            if not mane_data:
                return None
            mane_data_len = len(mane_data)

            # Transcript Priority (Must pass validation checks):
            #  1. MANE Select
            #  2. MANE Plus Clinical
            #  3. Longest Compatible Remaining
            #     a. If there is a tie, choose the first-published transcript among
            #        those transcripts meeting criterion
            mane_transcripts = set()
            for i in range(mane_data_len):
                index = mane_data_len - i - 1
                current_mane_data = mane_data[index]
                mane_transcripts |= set((current_mane_data["RefSeq_nuc"],
                                         current_mane_data["Ensembl_nuc"]))
                mane = await self._g_to_c(
                    g=g, refseq_c_ac=current_mane_data["RefSeq_nuc"],
                    status="_".join(current_mane_data["MANE_status"].split()).lower(),
                    ensembl_c_ac=current_mane_data["Ensembl_nuc"])
                if not mane:
                    continue

                if not mane["alt_ac"]:
                    g_alt_ac = g.get("alt_ac")
                    if g_alt_ac:
                        mane["alt_ac"] = g_alt_ac

                valid_reading_frame = self._validate_reading_frames(
                    c_ac, c_pos[0], c_pos[1], mane
                )
                if not valid_reading_frame:
                    continue

                if anno == "p":
                    mane = self._get_mane_p(current_mane_data, mane["pos"])

                if ref:
                    valid_references = self._validate_references(
                        ac, g["coding_start_site"], start_pos, end_pos,
                        mane, ref, anno, residue_mode
                    )
                    if not valid_references:
                        continue

                return mane

            if try_longest_compatible:
                if anno == "p":
                    return await self.get_longest_compatible_transcript(
                        g["gene"], start_pos, end_pos, "p", ref,
                        residue_mode=residue_mode, mane_transcripts=mane_transcripts)
                else:
                    return await self.get_longest_compatible_transcript(
                        g["gene"], c_pos[0], c_pos[1], "c", ref,
                        residue_mode=residue_mode, mane_transcripts=mane_transcripts)
            else:
                return None
        elif anno == "g":
            return await self.g_to_mane_c(ac, start_pos, end_pos, gene=gene,
                                          residue_mode=residue_mode)
        else:
            logger.warning(f"Annotation layer not supported: {anno}")

    async def g_to_grch38(self, ac: str, start_pos: int,
                          end_pos: int) -> Optional[Dict]:
        """Return genomic coordinate on GRCh38 when not given gene context.

        :param str ac: Genomic accession
        :param int start_pos: Genomic start position change
        :param int end_pos: Genomic end position change
        :return: NC accession, start and end pos on GRCh38 assembly
        """
        if end_pos is None:
            end_pos = start_pos

        # Checking to see what chromosome and assembly we're on
        descr = await self.uta_db.get_chr_assembly(ac)
        if not descr:
            # Already GRCh38 assembly
            if self._validate_index(ac, (start_pos, end_pos), 0):
                return dict(
                    ac=ac,
                    pos=(start_pos, end_pos)
                )
            else:
                return None
        chromosome, assembly = descr
        is_same_pos = start_pos == end_pos

        # Coordinate liftover
        if assembly < "GRCh37":
            logger.warning("Liftover only supported for GRCh37")
            return None

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
            ac = newest_ac[0][0]
            if self._validate_index(ac, (start_pos, end_pos), 0):
                return dict(
                    ac=ac,
                    pos=(start_pos, end_pos)
                )

        return None

    @staticmethod
    def get_mane_c_pos_change(mane_tx_genomic_data: Dict,
                              coding_start_site: int) -> Tuple[int, int]:
        """Get mane c position change

        :param Dict mane_tx_genomic_data: MANE transcript and genomic data
        :param int coding_start_site: Coding start site
        :return: cDNA pos start, cDNA pos end
        """
        tx_pos_range = mane_tx_genomic_data["tx_pos_range"]
        alt_pos_change = mane_tx_genomic_data["alt_pos_change"]

        mane_c_pos_change = (
            tx_pos_range[0] + alt_pos_change[0] - coding_start_site,
            tx_pos_range[1] - alt_pos_change[1] - coding_start_site
        )

        if mane_c_pos_change[0] > mane_c_pos_change[1]:
            mane_c_pos_change = mane_c_pos_change[1], mane_c_pos_change[0]
        return mane_c_pos_change

    async def g_to_mane_c(
        self, ac: str, start_pos: int, end_pos: int,
        gene: Optional[str] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Optional[Dict]:
        """Return MANE Transcript on the c. coordinate.
        If gene is provided, g->GRCh38->MANE c.
            If MANE c. cannot be found, we return the genomic coordinate on
                GRCh38
        If gene is not provided, g -> GRCh38

        :param str ac: Transcript accession on g. coordinate
        :param int start_pos: genomic change start position
        :param int end_pos: genomic change end position
        :param str gene: Gene symbol
        :param ResidueMode residue_mode: Starting residue mode for `start_pos`
            and `end_pos`. Will always return coordinates in inter-residue
        :return: MANE Transcripts with cDNA change on c. coordinate if gene
            is provided. Else, GRCh38 data
        """
        inter_residue_pos, _ = get_inter_residue_pos(
            start_pos, residue_mode, end_pos=end_pos)
        if not inter_residue_pos:
            return None
        start_pos, end_pos = inter_residue_pos
        residue_mode = ResidueMode.INTER_RESIDUE

        # If gene not provided, return GRCh38
        if not gene:
            grch38 = await self.g_to_grch38(ac, start_pos, end_pos)
            if not grch38:
                return None

            return dict(
                gene=None,
                refseq=grch38["ac"],
                ensembl=None,
                coding_start_site=None,
                coding_end_site=None,
                pos=grch38["pos"],
                strand=None,
                status="GRCh38",
                alt_ac=grch38["ac"]
            )

        if not await self.uta_db.validate_genomic_ac(ac):
            logger.warning(f"Genomic accession does not exist: {ac}")
            return None

        mane_data = self.mane_transcript_mappings.get_gene_mane_data(gene)
        if not mane_data:
            return None
        mane_data_len = len(mane_data)

        for i in range(mane_data_len):
            index = mane_data_len - i - 1
            current_mane_data = mane_data[index]
            mane_c_ac = current_mane_data["RefSeq_nuc"]

            # Liftover to GRCh38
            grch38 = await self.g_to_grch38(ac, start_pos, end_pos)
            mane_tx_genomic_data = None
            if grch38:
                # GRCh38 -> MANE C
                mane_tx_genomic_data = await self.uta_db.get_mane_c_genomic_data(  # noqa: E501
                    mane_c_ac, None, grch38["pos"][0], grch38["pos"][1]
                )

            if not grch38 or not mane_tx_genomic_data:
                # GRCh38 did not work, so let's try original assembly (37)
                mane_tx_genomic_data = await self.uta_db.get_mane_c_genomic_data(  # noqa: E501
                    mane_c_ac, ac, start_pos, end_pos
                )
                if not mane_tx_genomic_data:
                    continue
                else:
                    logger.info("Not using most recent assembly")

            coding_start_site = mane_tx_genomic_data["coding_start_site"]
            coding_end_site = mane_tx_genomic_data["coding_end_site"]
            mane_c_pos_change = self.get_mane_c_pos_change(
                mane_tx_genomic_data, coding_start_site)

            if not self._validate_index(mane_c_ac, mane_c_pos_change,
                                        coding_start_site):
                logger.warning(f"{mane_c_pos_change} are not valid positions"
                               f" on {mane_c_ac}with coding start site "
                               f"{coding_start_site}")
                continue

            return self._get_c_data(
                gene=current_mane_data["symbol"],
                cds_start_end=(coding_start_site, coding_end_site),
                c_pos_change=mane_c_pos_change,
                strand=current_mane_data["chr_strand"],
                status="_".join(current_mane_data["MANE_status"].split()).lower(),
                refseq_c_ac=current_mane_data["RefSeq_nuc"],
                ensembl_c_ac=current_mane_data["Ensembl_nuc"],
                alt_ac=grch38["ac"] if grch38 else None)

    async def get_mapped_mane_data(
        self, gene: str, assembly: Assembly, genomic_position: int,
        residue_mode: ResidueMode = ResidueMode.INTER_RESIDUE
    ) -> Optional[MappedManeData]:
        """Get MANE data for gene, assembly, and position. If GRCh37 assembly is given,
        will return mapped MANE data.

        :param str gene: Gene symbol or identifier
        :param Assembly assembly: Assembly for the provided genomic position
        :param int genomic_position: Position on the genomic reference sequence to find
            MANE data for
        :param ResidueMode residue_mode: Starting residue mode for `start_pos`
            and `end_pos`. Will always return coordinates in inter-residue
        :return: Mapped MANE or Longest Compatible Remaining data if found/compatible.
            MANETranscriptError will be raised if unable to get required data for
            retrieving mapped MANE data.
        """
        hgnc_gene_data = self.gene_normalizer.get_hgnc_data(gene)
        if not hgnc_gene_data:
            raise MANETranscriptError(f"Unable to get HGNC data for gene: {gene}")

        gene = hgnc_gene_data["symbol"]

        mane_data = self.mane_transcript_mappings.get_gene_mane_data(gene)
        if not mane_data:
            raise MANETranscriptError(f"Unable to get MANE data for gene: {gene}")

        mane_data_len = len(mane_data)

        alt_ac = None
        if hgnc_gene_data["locations"]:
            chr = hgnc_gene_data["locations"][0].get("chr") or ""
            alt_acs, _ = self.seqrepo_access.translate_identifier(f"{assembly}:{chr}",
                                                                  "refseq")
            if alt_acs:
                alt_ac = alt_acs[0].split(":")[1]
            else:
                raise MANETranscriptError(f"Unable to translate identifier for: "
                                          f"{assembly}:{chr}")

        inter_residue_pos, _ = get_inter_residue_pos(genomic_position, residue_mode)
        g_pos = inter_residue_pos[0]

        mane_transcripts = set()
        for i in range(mane_data_len):
            index = mane_data_len - i - 1
            current_mane_data = mane_data[index]
            mane_transcripts |= set((current_mane_data["RefSeq_nuc"],
                                     current_mane_data["Ensembl_nuc"]))
            mane_c_ac = current_mane_data["RefSeq_nuc"]

            ac_query = mane_c_ac.split(".")[0]
            tx_exon_aln_v_data = await self.uta_db.get_tx_exon_aln_v_data(
                ac_query, g_pos, g_pos, alt_ac, False, True)

            if not tx_exon_aln_v_data:
                continue
            else:
                len_of_aligned_data = len(tx_exon_aln_v_data)
                if len_of_aligned_data == 1:
                    tx_exon_aln_v_data = tx_exon_aln_v_data[0]
                else:
                    logger.debug(f"Found {len_of_aligned_data} records for aligned "
                                 f"mapped MANE data for {ac_query}, {g_pos}, {alt_ac}")

                    # Try checking for MANE match
                    filter_data = list(filter(lambda x: x[1] == mane_c_ac,
                                              tx_exon_aln_v_data))
                    if filter_data:
                        tx_exon_aln_v_data = filter_data[0]
                    else:
                        # Try checking for older versions of MANE
                        filter_data = list(filter(lambda x: x[1].startswith(
                            mane_c_ac.split(".")[0]), tx_exon_aln_v_data))
                        if filter_data:
                            filter_data.sort(key=lambda x: x[1], reverse=True)
                            tx_exon_aln_v_data = filter_data[0]
            return MappedManeData(
                gene=gene,
                refseq=current_mane_data["RefSeq_nuc"],
                ensembl=current_mane_data["Ensembl_nuc"],
                strand="-" if tx_exon_aln_v_data[7] == -1 else "+",
                status="_".join(current_mane_data["MANE_status"].split()).lower(),
                alt_ac=alt_ac,
                assembly=assembly.value
            )

        lcr_data = await self.get_longest_compatible_transcript(
            gene, g_pos, g_pos, AnnotationLayer.GENOMIC,
            residue_mode=ResidueMode.INTER_RESIDUE, mane_transcripts=mane_transcripts,
            alt_ac=alt_ac)
        if lcr_data:
            return MappedManeData(
                gene=gene,
                refseq=lcr_data["refseq"],
                ensembl=lcr_data["ensembl"],
                strand=lcr_data["strand"],
                status=lcr_data["status"],
                alt_ac=alt_ac,
                assembly=assembly.value
            )

        return None
