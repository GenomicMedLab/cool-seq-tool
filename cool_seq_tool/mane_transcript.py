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

import pandas as pd

from cool_seq_tool.schemas import AnnotationLayer, Assembly, MappedManeData, \
    ResidueMode, TranscriptPriorityLabel
from cool_seq_tool.alignment_mapper import AlignmentMapper
from cool_seq_tool.data_sources import SeqRepoAccess, TranscriptMappings, \
    MANETranscriptMappings, UTADatabase, GeneNormalizer
from cool_seq_tool.utils.positions import get_inter_residue_pos, \
    to_zero_based, zero_based_to_inter_residue
from cool_seq_tool import logger
from cool_seq_tool.utils.validation import ValidationError, validate_index, \
    validate_reading_frames, validate_references


class MANETranscriptError(Exception):
    """Custom exception for MANETranscript class"""


class MANETranscript(AlignmentMapper):
    """Class for retrieving MANE transcripts."""

    def __init__(self, seqrepo_access: SeqRepoAccess,
                 transcript_mappings: TranscriptMappings,
                 uta_db: UTADatabase,
                 mane_transcript_mappings: MANETranscriptMappings,
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
        super().__init__(seqrepo_access, transcript_mappings, uta_db)
        self.mane_transcript_mappings = mane_transcript_mappings
        self.gene_normalizer = gene_normalizer

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
        copy_df["ac_no_version_as_int"] = copy_df["tx_ac"].apply(lambda x: int(x.split(".")[0].split("NM_")[1]))  # noqa: E501
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
        begin_start_eq_end = start_pos == end_pos
        begin_residue_mode = residue_mode

        is_p_or_c_start_anno = True
        if start_annotation_layer == AnnotationLayer.PROTEIN:
            c_start_pos, c_end_pos = self.p_pos_to_c_pos(
                start_pos, end_pos, residue_mode)
        elif start_annotation_layer == AnnotationLayer.CDNA:
            c_start_pos, c_end_pos = get_inter_residue_pos(
                start_pos, end_pos, residue_mode)
        else:
            start_pos, end_pos = get_inter_residue_pos(
                start_pos, end_pos, residue_mode)
            is_p_or_c_start_anno = False

        residue_mode = ResidueMode.INTER_RESIDUE

        # Data Frame that contains transcripts associated to a gene
        if is_p_or_c_start_anno:
            c_start_pos, c_end_pos = to_zero_based(c_start_pos, c_end_pos, residue_mode)
            df = await self.uta_db.get_transcripts_from_gene(
                gene, c_start_pos, c_end_pos, use_tx_pos=True, alt_ac=alt_ac)
        else:
            start_pos, end_pos = to_zero_based(start_pos, end_pos, residue_mode)
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
            # Only need to check the one row since we do liftover in c_to_g
            tmp_df = df.loc[df["tx_ac"] == tx_ac].sort_values("alt_ac", ascending=False)
            row = tmp_df.iloc[0]

            if alt_ac is None:
                alt_ac = row["alt_ac"]
            cds_start = row["cds_start_i"]

            found_tx_exon_aln_v_result = False
            if is_p_or_c_start_anno:
                # Go from c -> g annotation (liftover as well)
                g, _ = await self.c_to_g(
                    tx_ac, c_start_pos, c_end_pos, cds_start, residue_mode=residue_mode,
                    target_genome_assembly=Assembly.GRCH38
                )
            else:
                # g -> GRCh38 (if alt_ac not provided. if it is, will use that assembly)
                grch38 = await self.g_to_grch38(alt_ac, start_pos, end_pos)
                g = {
                    "g_ac": grch38["ac"],
                    "g_start_pos": grch38["pos"][0],
                    "g_end_pos": grch38["pos"][1]
                }
                found_tx_exon_aln_v_result = True
            if not g:
                continue

            # Get prioritized transcript data for gene
            # grch38 -> c
            # TODO: FIx
            c_data, _ = await self.g_to_c(
                g["g_ac"], g["g_start_pos"], g["g_end_pos"], tx_ac, cds_start=cds_start,
                residue_mode=ResidueMode.INTER_RESIDUE, is_zero_based=True)

            if not c_data:
                continue

            # Validation checks
            if is_p_or_c_start_anno:
                try:
                    validate_reading_frames(tx_ac, c_start_pos, c_end_pos, c_data)
                except ValidationError:
                    continue

            if ref:
                try:
                    if start_annotation_layer == AnnotationLayer.PROTEIN:
                        validate_references(
                            self.seqrepo_access, row["pro_ac"], row["cds_start_i"],
                            start_pos + 1, end_pos + 1, {}, ref,
                            AnnotationLayer.PROTEIN, ResidueMode.RESIDUE
                        )
                    elif start_annotation_layer == AnnotationLayer.CDNA:
                        validate_references(
                            self.seqrepo_access, row["tx_ac"], row["cds_start_i"],
                            c_start_pos + 1, c_end_pos + 1, {}, ref,
                            AnnotationLayer.CDNA, ResidueMode.RESIDUE
                        )
                    else:
                        validate_references(
                            self.seqrepo_access, alt_ac, 0, start_pos + 1, end_pos + 1,
                            {}, ref, AnnotationLayer.GENOMIC, ResidueMode.RESIDUE
                        )
                except ValidationError:
                    continue

            c_data["c_start_pos"], c_data["c_end_pos"] = zero_based_to_inter_residue(
                c_data["c_start_pos"], c_data["c_end_pos"], begin_start_eq_end,
                begin_residue_mode)

            if start_annotation_layer == AnnotationLayer.PROTEIN:
                pos = (math.ceil(c_data["pos"][0] / 3),
                       math.floor(c_data["pos"][1] / 3))
                ac = row["pro_ac"]
                coding_start_site = 0
            else:
                # cDNA and Genomic annotations will return c. data
                pos = c_data["c_start_pos"], c_data["c_end_pos"]
                ac = tx_ac
                coding_start_site = c_data["cds_start"]

            try:
                validate_index(self.seqrepo_access, ac, pos, coding_start_site)
            except ValidationError:
                continue

            return dict(
                refseq=ac if ac.startswith("N") else None,
                ensembl=ac if ac.startswith("E") else None,  # TODO: issues 87, 4
                pos=pos,
                strand=c_data["strand"],
                status=TranscriptPriorityLabel.LongestCompatibleRemaining.value
            )
        return None

    async def get_mane_transcript(
        self, ac: str, start_pos: int, end_pos: int,
        start_annotation_layer: AnnotationLayer,
        gene: Optional[str] = None, ref: Optional[str] = None,
        try_longest_compatible: bool = False,
        residue_mode: ResidueMode = ResidueMode.RESIDUE
    ) -> Optional[Dict]:
        """Return mane transcript.

        :param str ac: Accession
        :param int start_pos: Start position change on `ac`
        :param int end_pos: End position change on `ac`
        :param AnnotationLayer start_annotation_layer: Starting annotation layer. This
            is the annotation layer for `ac`, `start_pos`, `end_pos`
        :param Optional[str] gene: Gene symbol. Only used with AnnotationLayer.GENOMIC
        :param Optional[str] ref: Reference sequence on `ac` at positions
            (`start_pos`, `end_pos`)
        :param bool try_longest_compatible: `True` if should try longest
            compatible remaining if mane transcript was not compatible.
            `False` otherwise.
        :param ResidueMode residue_mode: Starting residue mode for `start_pos`
            and `end_pos`
        :return: MANE data or longest transcript compatible data if validation
            checks are correct. Will return inter-residue coordinates.
            Else, `None`
        """
        # TODO: Still need to raise Exceptions / warnings
        inter_residue_pos = get_inter_residue_pos(
            start_pos, end_pos, residue_mode)
        start_pos, end_pos = inter_residue_pos
        residue_mode = ResidueMode.INTER_RESIDUE
        begin_start_eq_end = start_pos == end_pos
        begin_residue_mode = residue_mode

        anno = start_annotation_layer.lower()
        if anno in {AnnotationLayer.PROTEIN, AnnotationLayer.CDNA}:
            # Get accession and positions on c. coordinate
            if anno == AnnotationLayer.PROTEIN:
                c_data, _ = await self.p_to_c(ac, start_pos, end_pos,
                                              residue_mode=residue_mode)
                if not c_data:
                    return None
                # residue position
                c_ac = c_data["c_ac"]
                c_pos = c_data["c_start_pos"], c_data["c_end_pos"]
                cds_start = c_data["cds_start"]
            else:
                c_ac = ac
                c_pos = start_pos, end_pos
                cds_start, _ = await self._get_cds_start(c_ac)
                if not cds_start:
                    return None

            # Go from c -> g annotation (liftover as well)
            g_data, w = await self.c_to_g(
                c_ac, c_pos[0], c_pos[1], cds_start, residue_mode=residue_mode,
                target_genome_assembly=Assembly.GRCH38)
            if not g_data:
                return None

            # Get mane data for gene
            mane_data = self.mane_transcript_mappings.get_gene_mane_data(g_data["gene"])
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
                current_mane_data = mane_data[i]
                mane_transcripts |= set((current_mane_data["RefSeq_nuc"],
                                         current_mane_data["Ensembl_nuc"]))

                # g_data is inter-residue
                mane, _ = await self.g_to_c(
                    g_data["g_ac"], g_data["g_start_pos"], g_data["g_end_pos"],
                    current_mane_data["RefSeq_nuc"], residue_mode=residue_mode)
                if not mane:
                    continue
                cds_start = mane["cds_start"]

                # Will raise ValidationError if checks do not pass
                try:
                    validate_reading_frames(c_ac, c_pos[0], c_pos[1], mane)
                except ValidationError:
                    continue

                if anno == AnnotationLayer.PROTEIN:
                    mane = {
                        "gene": current_mane_data["symbol"],
                        "refseq": current_mane_data["RefSeq_prot"],
                        "ensembl": current_mane_data["Ensembl_prot"],
                        "pos": (math.ceil(mane["c_start_pos"] / 3),
                                math.floor(mane["c_end_pos"] / 3)),
                        "strand": current_mane_data["chr_strand"],
                        "status": "_".join(current_mane_data["MANE_status"].split()).lower()
                    }
                else:
                    mane = {
                        "alt_ac": g_data["g_ac"],
                        "refseq": current_mane_data["RefSeq_nuc"],
                        "ensembl": current_mane_data["Ensembl_nuc"],
                        "pos": (mane["c_start_pos"], mane["c_end_pos"]),
                        "status": "_".join(current_mane_data["MANE_status"].split()).lower(),
                        "strand": current_mane_data["chr_strand"],
                        "coding_start_site": mane["cds_start"],
                        "gene": current_mane_data["symbol"]
                    }

                if ref:
                    # Will raise ValidationError if checks do not pass
                    try:
                        validate_references(
                            self.seqrepo_access, ac, cds_start, start_pos, end_pos, mane,
                            ref, anno, residue_mode
                        )
                    except ValidationError:
                        continue

                return mane

            if try_longest_compatible:
                if anno == AnnotationLayer.PROTEIN:
                    return await self.get_longest_compatible_transcript(
                        g_data["gene"], start_pos, end_pos, AnnotationLayer.PROTEIN, ref,
                        residue_mode=residue_mode, mane_transcripts=mane_transcripts)
                else:
                    return await self.get_longest_compatible_transcript(
                        g_data["gene"], c_pos[0], c_pos[1], AnnotationLayer.CDNA, ref,
                        residue_mode=residue_mode, mane_transcripts=mane_transcripts)
            else:
                return None
        elif anno == AnnotationLayer.GENOMIC:
            return await self.g_to_mane_c(ac, start_pos, end_pos, gene=gene,
                                          residue_mode=residue_mode)
        else:
            logger.warning(f"Annotation layer not supported: {anno}")

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
        begin_start_eq_end = start_pos == end_pos
        start_pos, end_pos = to_zero_based(start_pos, end_pos, residue_mode)

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
            current_mane_data = mane_data[i]
            mane_c_ac = current_mane_data["RefSeq_nuc"]

            # Liftover to GRCh38
            grch38 = await self.g_to_grch38(ac, start_pos, end_pos)
            mane_tx_genomic_data = None
            if grch38:
                # GRCh38 -> MANE C
                mane_tx_genomic_data = await self.uta_db.get_genomic_tx_data(
                    mane_c_ac, grch38["pos"], AnnotationLayer.GENOMIC)

            if not grch38 or not mane_tx_genomic_data:
                # GRCh38 did not work, so let's try original assembly (37)
                mane_tx_genomic_data = await self.uta_db.get_genomic_tx_data(
                    mane_c_ac, grch38["pos"], AnnotationLayer.GENOMIC,
                    target_genome_assembly=Assembly.GRCH37
                )
                if not mane_tx_genomic_data:
                    continue
                else:
                    logger.info("Not using most recent assembly")

            coding_start_site = mane_tx_genomic_data["coding_start_site"]
            coding_end_site = mane_tx_genomic_data["coding_end_site"]
            mane_c_pos_change = self.get_mane_c_pos_change(
                mane_tx_genomic_data, coding_start_site)

            # Will raise ValidationError if validation checks dont pass
            try:
                validate_index(
                    self.seqrepo_access, mane_c_ac, mane_c_pos_change, coding_end_site)
            except ValidationError:
                continue

            mane_c_pos_change = zero_based_to_inter_residue(
                mane_c_pos_change[0], mane_c_pos_change[1], begin_start_eq_end,
                residue_mode)

            return {
                "gene": current_mane_data["symbol"],
                "refseq": current_mane_data["RefSeq_nuc"],
                "ensembl": current_mane_data["Ensembl_nuc"],
                "coding_start_site": coding_start_site,
                "coding_end_site": coding_end_site,
                "pos": mane_c_pos_change,
                "strand": current_mane_data["chr_strand"],
                "status": "_".join(current_mane_data["MANE_status"].split()).lower(),
                "alt_ac": grch38["ac"] if grch38 else None
            }

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

        inter_residue_pos = get_inter_residue_pos(
            genomic_position, genomic_position, residue_mode)
        g_pos = inter_residue_pos[0]

        mane_transcripts = set()
        for i in range(mane_data_len):
            current_mane_data = mane_data[i]
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
