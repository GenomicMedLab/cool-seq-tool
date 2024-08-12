"""Provide mapping capabilities between transcript exon and genomic coordinates."""

import logging
from typing import TypeVar

from pydantic import Field, StrictInt

from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.mappers.liftover import LiftOver
from cool_seq_tool.mappers.mane_transcript import CdnaRepresentation, ManeTranscript
from cool_seq_tool.schemas import (
    AnnotationLayer,
    Assembly,
    BaseModelForbidExtra,
    CoordinateType,
    GenomicData,
    GenomicDataResponse,
    Strand,
    TranscriptExonData,
    TranscriptExonDataResponse,
)
from cool_seq_tool.sources.mane_transcript_mappings import ManeTranscriptMappings
from cool_seq_tool.sources.uta_database import GenesGenomicAcs, UtaDatabase
from cool_seq_tool.utils import service_meta

CoordinatesResponseType = TypeVar(
    "CoordinatesResponseType", GenomicDataResponse, TranscriptExonDataResponse
)

_logger = logging.getLogger(__name__)


class ExonCoord(BaseModelForbidExtra):
    """Model for representing exon coordinate data"""

    ord: StrictInt = Field(..., description="Exon number. 0-based.")
    tx_start_i: StrictInt = Field(
        ...,
        description="Transcript start index of the exon. Inter-residue coordinates.",
    )
    tx_end_i: StrictInt = Field(
        ..., description="Transcript end index of the exon. Inter-residue coordinates."
    )
    alt_start_i: StrictInt = Field(
        ..., description="Genomic start index of the exon. Inter-residue coordinates."
    )
    alt_end_i: StrictInt = Field(
        ..., description="Genomic end index of the exon. Inter-residue coordinates."
    )
    alt_strand: Strand = Field(..., description="Strand.")


class ExonGenomicCoordsMapper:
    """Provide capabilities for mapping transcript exon representation to/from genomic
    coordinate representation.
    """

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        uta_db: UtaDatabase,
        mane_transcript: ManeTranscript,
        mane_transcript_mappings: ManeTranscriptMappings,
        liftover: LiftOver,
    ) -> None:
        """Initialize ExonGenomicCoordsMapper class.

        A lot of resources are required for initialization, so when defaults are enough,
        it's easiest to let the core CoolSeqTool class handle it for you:

        >>> from cool_seq_tool import CoolSeqTool
        >>> egc = CoolSeqTool().ex_g_coords_mapper

        Note that this class's public methods are all defined as ``async``, so they will
        need to be called with ``await`` when called from a function, or run from an
        event loop. See the :ref:`Usage section <async_note>` for more information.

        >>> import asyncio
        >>> result = asyncio.run(
        ...     egc.tx_segment_to_genomic("NM_002529.3", exon_start=2, exon_end=17)
        ... )
        >>> result.genomic_data.start, result.genomic_data.end
        (156864428, 156881456)

        :param seqrepo_access: SeqRepo instance to give access to query SeqRepo database
        :param uta_db: UtaDatabase instance to give access to query UTA database
        :param mane_transcript: Instance to align to MANE or compatible representation
        :param mane_transcript_mappings: Instance to provide access to ManeTranscriptMappings class
        :param liftover: Instance to provide mapping between human genome assemblies
        """
        self.seqrepo_access = seqrepo_access
        self.uta_db = uta_db
        self.mane_transcript = mane_transcript
        self.mane_transcript_mappings = mane_transcript_mappings
        self.liftover = liftover

    @staticmethod
    def _return_warnings(
        resp: CoordinatesResponseType, warning_msg: list[str]
    ) -> CoordinatesResponseType:
        """Add warnings to response object

        :param resp: Response object
        :param warning_msg: Warning message(s) on why ``transcript_exon_data`` or
            ``genomic_data`` field is ``None``
        :return: Response object with warning message
        """
        for msg in warning_msg:
            _logger.warning(msg)
            resp.warnings.append(msg)
        return resp

    async def tx_segment_to_genomic(
        self,
        transcript: str,
        gene: str | None = None,
        exon_start: int | None = None,
        exon_start_offset: int = 0,
        exon_end: int | None = None,
        exon_end_offset: int = 0,
    ) -> GenomicDataResponse:
        """Get genomic data given transcript segment data.

        By default, transcript data is aligned to the GRCh38 assembly.

        >>> import asyncio
        >>> from cool_seq_tool import CoolSeqTool
        >>> egc = CoolSeqTool().ex_g_coords_mapper
        >>> tpm3 = asyncio.run(
        ...     egc.tx_segment_to_genomic(
        ...         "NM_152263.3",
        ...         gene="TPM3",
        ...         exon_start=1,
        ...         exon_end=8,
        ...     )
        ... )
        >>> tpm3.genomic_data.chr, tpm3.genomic_data.start, tpm3.genomic_data.end
        ('NC_000001.11', 154192135, 154170399)

        :param transcript: Transcript accession
        :param gene: HGNC gene symbol
        :param exon_start: Starting transcript exon number (1-based). If not provided,
            must provide ``exon_end``
        :param exon_start_offset: Starting exon offset
        :param exon_end: Ending transcript exon number (1-based). If not provided, must
            provide ``exon_start``
        :param exon_end_offset: Ending exon offset
        :return: GRCh38 genomic data (inter-residue coordinates)
        """
        resp = GenomicDataResponse(
            genomic_data=None, warnings=[], service_meta=service_meta()
        )

        # Ensure valid inputs
        warnings = []
        if not transcript:
            warnings.append("Must provide `transcript`")
        else:
            transcript = transcript.strip()

        exon_start_exists, exon_end_exists = False, False
        if exon_start is not None:
            if exon_start < 1:
                warnings.append("`exon_start` cannot be less than 1")
            exon_start_exists = True

        if exon_end is not None:
            if exon_end < 1:
                warnings.append("`exon_end` cannot be less than 1")
            exon_end_exists = True

        if not exon_start_exists and not exon_end_exists:
            warnings.append("Must provide either `exon_start` or `exon_end`")
        if exon_start_exists and exon_end_exists and (exon_start > exon_end):
            warnings.append(
                f"Start exon {exon_start} is greater than end exon {exon_end}"
            )

        if warnings:
            return self._return_warnings(resp, warnings)

        # Get exon start and exon end coordinates
        (
            tx_exon_start_coords,
            tx_exon_end_coords,
            errors,
        ) = await self._get_start_end_exon_coords(
            transcript, exon_start=exon_start, exon_end=exon_end
        )
        if errors:
            return self._return_warnings(resp, errors)

        if gene:
            gene = gene.upper().strip()

        # Get aligned genomic data (hgnc gene, alt_ac, alt_start_i, alt_end_i, strand)
        # for exon(s)
        alt_ac_start_end, warning = await self._get_alt_ac_start_and_end(
            transcript, tx_exon_start_coords, tx_exon_end_coords, gene=gene
        )
        if not alt_ac_start_end:
            return self._return_warnings(resp, [warning] if warning else [])
        alt_ac_start_data, alt_ac_end_data = alt_ac_start_end

        # Get gene and chromosome data, check that at least one was retrieved
        gene = alt_ac_start_data.hgnc if alt_ac_start_data else alt_ac_end_data.hgnc
        chromosome = (
            alt_ac_start_data.alt_ac if alt_ac_start_data else alt_ac_end_data.alt_ac
        )
        if gene is None or chromosome is None:
            return self._return_warnings(
                resp,
                [
                    "Unable to retrieve `gene` or `chromosome` from genomic start and genomic end data"
                ],
            )

        g_start = alt_ac_start_data.alt_end_i - 1 if alt_ac_start_data else None
        g_end = alt_ac_end_data.alt_start_i + 1 if alt_ac_end_data else None
        strand = (
            Strand(alt_ac_start_data.alt_strand)
            if alt_ac_start_data
            else Strand(alt_ac_end_data.alt_strand)
        )

        # Using none since could set to 0
        start_exits = g_start is not None
        end_exists = g_end is not None

        # Calculate offsets
        if strand == Strand.NEGATIVE:
            start_offset = exon_start_offset * -1 if start_exits else None
            end_offset = exon_end_offset * -1 if end_exists else 0
        else:
            start_offset = exon_start_offset if start_exits else 0
            end_offset = exon_end_offset if end_exists else 0

        # Get genomic coordinates with offsets included
        g_start = g_start + start_offset if start_exits else None
        g_end = g_end + end_offset if end_exists else None

        resp.genomic_data = GenomicData(
            gene=gene,
            chr=chromosome,
            start=g_start,
            end=g_end,
            exon_start=exon_start if start_exits else None,
            exon_start_offset=exon_start_offset,
            exon_end=exon_end if end_exists else None,
            exon_end_offset=exon_end_offset,
            transcript=transcript,
            strand=strand,
        )

        return resp

    async def genomic_to_tx_segment(
        self,
        chromosome: str | None = None,
        alt_ac: str | None = None,
        genomic_start: int | None = None,
        genomic_end: int | None = None,
        transcript: str | None = None,
        get_nearest_transcript_junction: bool = False,
        gene: str | None = None,
    ) -> GenomicDataResponse:
        """Get transcript segment data for genomic data, lifted over to GRCh38.

        Must provide inter-residue coordinates.

        MANE Transcript data will be returned if and only if ``transcript`` is not
        supplied. ``gene`` must be given in order to retrieve MANE Transcript data.

        >>> import asyncio
        >>> from cool_seq_tool import CoolSeqTool
        >>> from cool_seq_tool.schemas import Strand
        >>> egc = CoolSeqTool().ex_g_coords_mapper
        >>> result = asyncio.run(
        ...     egc.genomic_to_tx_segment(
        ...         alt_ac="NC_000001.11",
        ...         genomic_start=154192136,
        ...         genomic_end=154170400,
        ...         transcript="NM_152263.3",
        ...     )
        ... )
        >>> result.genomic_data.exon_start, result.genomic_data.exon_end
        (1, 8)

        :param chromosome: e.g. ``"1"`` or ``"chr1"``. If not provided, must provide
            ``alt_ac``. If ``alt_ac`` is also provided, ``alt_ac`` will be used.
        :param alt_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``alt_ac``
            will be used.
        :param genomic_start: Genomic position where the transcript segment starts
        :param genomic_end: Genomic position where the transcript segment ends
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript. See the :ref:`Transcript Selection policy <transcript_selection_policy>`
            page.
        :param get_nearest_transcript_junction: If ``True``, this will return the
            adjacent exon if the position specified by``genomic_start`` or
            ``genomic_end`` does not occur on an exon. For the positive strand, adjacent
            is defined as the exon preceding the breakpoint for the 5' end and the exon
            following the breakpoint for the 3' end. For the negative strand, adjacent
            is defined as the exon following the breakpoint for the 5' end and the exon
            preceding the breakpoint for the 3' end.
        :param gene: gene name. Ideally, HGNC symbol. Must be given if no ``transcript``
            value is provided.
        :param coordinate_type: Coordinate type for ``genomic_start`` and
            ``genomic_end``
        :return: Genomic data (inter-residue coordinates)
        """
        resp = GenomicDataResponse(
            genomic_data=None, warnings=[], service_meta=service_meta()
        )
        warnings = []
        if genomic_start is None and genomic_end is None:
            warnings.append("Must provide either `genomic_start` or `genomic_end`")
        if chromosome is None and alt_ac is None:
            warnings.append("Must provide either `chromosome` or `alt_ac`")
        if transcript is None and gene is None:
            warnings.append("Must provide either `gene` or `transcript`")
        if warnings:
            return self._return_warnings(resp, warnings)

        params = {key: None for key in GenomicData.model_fields}
        if gene is not None:
            gene = gene.upper().strip()

        if genomic_start:
            start_data = await self._genomic_to_transcript_exon_coordinate(
                genomic_start,
                chromosome=chromosome,
                alt_ac=alt_ac,
                transcript=transcript,
                gene=gene,
                get_nearest_transcript_junction=get_nearest_transcript_junction,
                is_start=True,
            )
            if start_data.transcript_exon_data:
                start_data = start_data.transcript_exon_data.model_dump()
            else:
                return self._return_warnings(resp, [start_data.warnings[0]])
        else:
            start_data = None

        if genomic_end:
            genomic_end -= 1
            end_data = await self._genomic_to_transcript_exon_coordinate(
                genomic_end,
                chromosome=chromosome,
                alt_ac=alt_ac,
                transcript=transcript,
                gene=gene,
                get_nearest_transcript_junction=get_nearest_transcript_junction,
                is_start=False,
            )
            if end_data.transcript_exon_data:
                end_data = end_data.transcript_exon_data.model_dump()
            else:
                return self._return_warnings(resp, [end_data.warnings[0]])
        else:
            end_data = None

        for field in ["transcript", "gene", "chr", "strand"]:
            if start_data:
                if end_data and (start_data[field] != end_data[field]):
                    msg = (
                        f"Start `{field}`, {start_data[field]}, does "
                        f"not match End `{field}`, {end_data[field]}"
                    )
                    return self._return_warnings(resp, [msg])
                params[field] = start_data[field]
            else:
                params[field] = end_data[field]

        if gene and gene != params["gene"]:
            msg = (
                f"Input gene, {gene}, does not match expected output"
                f"gene, {params['gene']}"
            )
            return self._return_warnings(resp, [msg])

        for label, data in [("start", start_data), ("end", end_data)]:
            if data:
                params[label] = data["pos"]
                params[f"exon_{label}"] = data["exon"]
                params[f"exon_{label}_offset"] = data["exon_offset"]
        resp.genomic_data = GenomicData(**params)
        return resp

    async def _get_all_exon_coords(
        self, tx_ac: str, genomic_ac: str | None = None
    ) -> list[ExonCoord]:
        """Get all exon coordinate data for a transcript.

        If ``genomic_ac`` is NOT provided, this method will use the GRCh38 accession
        associated to ``tx_ac``.

        :param tx_ac: The RefSeq transcript accession to get exon data for.
        :param genomic_ac: The RefSeq genomic accession to get exon data for.
        :return: List of all exon coordinate data for ``tx_ac`` and ``genomic_ac``.
            The exon coordinate data will include the exon number, transcript and
            genomic positions for the start and end of the exon, and strand.
            The list will be ordered by ascending exon number.
        """
        if genomic_ac:
            query = f"""
                SELECT DISTINCT ord, tx_start_i, tx_end_i, alt_start_i, alt_end_i, alt_strand
                FROM {self.uta_db.schema}.tx_exon_aln_v
                WHERE tx_ac = '{tx_ac}'
                AND alt_aln_method = 'splign'
                AND alt_ac = '{genomic_ac}'
                ORDER BY ord ASC
                """  # noqa: S608
        else:
            query = f"""
                SELECT DISTINCT ord, tx_start_i, tx_end_i, alt_start_i, alt_end_i, alt_strand
                FROM {self.uta_db.schema}.tx_exon_aln_v as t
                INNER JOIN {self.uta_db.schema}._seq_anno_most_recent as s
                ON t.alt_ac = s.ac
                WHERE s.descr = ''
                AND t.tx_ac = '{tx_ac}'
                AND t.alt_aln_method = 'splign'
                AND t.alt_ac like 'NC_000%'
                ORDER BY ord ASC
                """  # noqa: S608

        results = await self.uta_db.execute_query(query)
        return [ExonCoord(**r) for r in results]

    async def _get_start_end_exon_coords(
        self,
        tx_ac: str,
        exon_start: int | None = None,
        exon_end: int | None = None,
        genomic_ac: str | None = None,
    ) -> tuple[ExonCoord | None, ExonCoord | None, list[str]]:
        """Get exon coordinates for a transcript given exon start and exon end.

        If ``genomic_ac`` is NOT provided, this method will use the GRCh38 accession
        associated to ``tx_ac``.

        :param tx_ac: The RefSeq transcript accession to get exon data for.
        :param exon_start: Start exon number to get coordinate data for. 1-based.
        :param exon_end: End exon number to get coordinate data for. 1-based.
        :param genomic_ac: The RefSeq genomic accession to get exon data for.
        :return: Tuple containing start exon coordinate data, end exon coordinate data,
            and list of errors. The exon coordinate data will include the exon number,
            transcript and genomic positions for the start and end of the exon, and
            strand.
        """
        tx_exons = await self._get_all_exon_coords(tx_ac, genomic_ac=genomic_ac)
        if not tx_exons:
            return None, None, [f"No exons found given {tx_ac}"]

        errors = []
        start_end_exons = []
        for exon_num in [exon_start, exon_end]:
            if exon_num is not None:
                try:
                    start_end_exons.append(tx_exons[exon_num - 1])
                    continue
                except IndexError:
                    errors.append(f"Exon {exon_num} does not exist on {tx_ac}")
            start_end_exons.append(None)

        if errors:
            start_end_exons = [None, None]

        return *start_end_exons, errors

    async def _get_alt_ac_start_and_end(
        self,
        tx_ac: str,
        tx_exon_start: ExonCoord | None = None,
        tx_exon_end: ExonCoord | None = None,
        gene: str | None = None,
    ) -> tuple[tuple[tuple[int, int], tuple[int, int]] | None, str | None]:
        """Get aligned genomic coordinates for transcript exon start and end.

        :param tx_ac: Transcript accession
        :param tx_exon_start: Transcript's exon start coordinates. If not provided,
            must provide ``tx_exon_end``
        :param tx_exon_end: Transcript's exon end coordinates. If not provided, must
            provide ``tx_exon_start``
        :param gene: HGNC gene symbol
        :return: Aligned genomic data, and warnings if found
        """
        if tx_exon_start is None and tx_exon_end is None:
            msg = "Must provide either `tx_exon_start` or `tx_exon_end` or both"
            _logger.warning(msg)
            return None, msg

        alt_ac_data = {"start": None, "end": None}
        for exon, key in [(tx_exon_start, "start"), (tx_exon_end, "end")]:
            if exon:
                alt_ac_val, warning = await self.uta_db.get_alt_ac_start_or_end(
                    tx_ac, exon.tx_start_i, exon.tx_end_i, gene=gene
                )
                if alt_ac_val:
                    alt_ac_data[key] = alt_ac_val
                else:
                    return None, warning

        alt_ac_data_values = alt_ac_data.values()
        # Validate that start and end alignments have matching gene, genomic accession,
        # and strand
        if all(alt_ac_data_values):
            for attr in ["hgnc", "alt_ac", "alt_strand"]:
                start_attr = getattr(alt_ac_data["start"], attr)
                end_attr = getattr(alt_ac_data["end"], attr)
                if start_attr != end_attr:
                    error = f"{attr} mismatch. {start_attr} != {end_attr}."
                    _logger.warning(
                        "%s: %s != %s",
                        error,
                        start_attr,
                        end_attr,
                    )
                    return None, error
        return tuple(alt_ac_data_values), None

    async def _genomic_to_transcript_exon_coordinate(
        self,
        genomic_pos: int,
        chromosome: str | None = None,
        alt_ac: str | None = None,
        transcript: str | None = None,
        gene: str | None = None,
        get_nearest_transcript_junction: bool = False,
        is_start: bool = True,
    ) -> TranscriptExonDataResponse:
        """Convert individual genomic data to transcript data

        :param genomic_pos: Genomic position where the transcript segment starts or ends
            (inter-residue based)
        :param chromosome: Chromosome. Must give chromosome without a prefix
            (i.e. ``1`` or ``X``). If not provided, must provide ``alt_ac``.
            If ``alt_ac`` is also provided, ``alt_ac`` will be used.
        :param alt_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``alt_ac``
            will be used.
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript
        :param gene: HGNC gene symbol
        :param get_nearest_transcript_junction: If ``True``, this will return the
            adjacent exon if the position specified by``genomic_start`` or
            ``genomic_end`` does not occur on an exon. For the positive strand, adjacent
            is defined as the exon preceding the breakpoint for the 5' end and the exon
            following the breakpoint for the 3' end. For the negative strand, adjacent
            is defined as the exon following the breakpoint for the 5' end and the exon
            preceding the breakpoint for the 3' end.
        :param is_start: ``True`` if ``genomic_pos`` is where the transcript segment starts.
            ``False`` if ``genomic_pos`` is where the transcript segment ends.
        :return: Transcript data (inter-residue coordinates)
        """
        resp = TranscriptExonDataResponse(
            transcript_exon_data=None, warnings=[], service_meta=service_meta()
        )
        params = {key: None for key in TranscriptExonData.model_fields}

        if get_nearest_transcript_junction:
            if not gene:
                return self._return_warnings(
                    resp,
                    [
                        "Gene must be provided to select the adjacent transcript junction"
                    ],
                )
            if not alt_ac:
                alt_acs, w = self.seqrepo_access.chromosome_to_acs(chromosome)

                if not alt_acs:
                    return self._return_warnings(resp, [w])
                alt_ac = alt_acs[0]

            if not transcript:
                # Select a transcript if not provided
                mane_transcripts = self.mane_transcript_mappings.get_gene_mane_data(
                    gene
                )

                if mane_transcripts:
                    transcript = mane_transcripts[0]["RefSeq_nuc"]
                else:
                    # Attempt to find a coding transcript if a MANE transcript
                    # cannot be found
                    results = await self.uta_db.get_transcripts(
                        gene=gene, alt_ac=alt_ac
                    )

                    if not results.is_empty():
                        transcript = results[0]["tx_ac"][0]
                    else:
                        # Run if gene is for a noncoding transcript
                        query = f"""
                            SELECT DISTINCT tx_ac
                            FROM {self.uta_db.schema}.tx_exon_aln_v
                            WHERE hgnc = '{gene}'
                            AND alt_ac = '{alt_ac}'
                            """  # noqa: S608
                        result = await self.uta_db.execute_query(query)

                        if result:
                            transcript = result[0]["tx_ac"]
                        else:
                            return self._return_warnings(
                                resp,
                                [f"Could not find a transcript for {gene} on {alt_ac}"],
                            )

            tx_genomic_coords = await self._get_all_exon_coords(
                tx_ac=transcript, genomic_ac=alt_ac
            )
            if not tx_genomic_coords:
                return self._return_warnings(
                    resp, [f"No exons found given {transcript} and {alt_ac}"]
                )

            strand = Strand(tx_genomic_coords[0].alt_strand)
            params["strand"] = strand

            # Check if breakpoint occurs on an exon.
            # If not, determine the adjacent exon given the selected transcript
            if not self._is_exonic_breakpoint(genomic_pos, tx_genomic_coords):
                exon = self._get_adjacent_exon(
                    tx_exons_genomic_coords=tx_genomic_coords,
                    strand=strand,
                    start=genomic_pos if is_start else None,
                    end=genomic_pos if not is_start else None,
                )

                params["exon"] = exon
                params["transcript"] = transcript
                params["gene"] = gene
                params["pos"] = genomic_pos
                params["chr"] = alt_ac

                self._set_exon_offset(
                    params=params,
                    start=tx_genomic_coords[exon - 1].alt_start_i,
                    end=tx_genomic_coords[exon - 1].alt_end_i,
                    pos=genomic_pos,
                    is_start=is_start,
                    strand=strand,
                )
                resp.transcript_exon_data = TranscriptExonData(**params)
                return resp

        if alt_ac:
            # Check if valid accession is given
            if not await self.uta_db.validate_genomic_ac(alt_ac):
                return self._return_warnings(
                    resp, [f"Invalid genomic accession: {alt_ac}"]
                )

            genes_alt_acs, warning = await self.uta_db.get_genes_and_alt_acs(
                genomic_pos, alt_ac=alt_ac, gene=gene
            )
        elif chromosome:
            # Check if just chromosome is given. If it is, we should
            # convert this to the correct accession version
            if chromosome == "X":
                chromosome = 23
            elif chromosome == "Y":
                chromosome = 24
            else:
                chromosome = int(chromosome)

            genes_alt_acs, warning = await self.uta_db.get_genes_and_alt_acs(
                genomic_pos, chromosome=chromosome, gene=gene
            )
        else:
            genes_alt_acs = None

        if not genes_alt_acs:
            return self._return_warnings(resp, [warning])

        gene_alt_ac, warning = self._get_gene_and_alt_ac(genes_alt_acs, gene)
        if not gene_alt_ac:
            return self._return_warnings(resp, [warning])
        gene, alt_ac = gene_alt_ac

        if transcript is None:
            warnings = await self._set_mane_genomic_data(
                params, gene, alt_ac, genomic_pos, is_start
            )
            if warnings:
                return self._return_warnings(resp, [warnings])
        else:
            params["transcript"] = transcript
            params["gene"] = gene
            params["pos"] = genomic_pos
            params["chr"] = alt_ac
            warning = await self._set_genomic_data(params, is_start)
            if warning:
                return self._return_warnings(resp, [warning])

        resp.transcript_exon_data = TranscriptExonData(**params)
        return resp

    @staticmethod
    def _get_gene_and_alt_ac(
        genes_alt_acs: GenesGenomicAcs, gene: str | None
    ) -> tuple[tuple[str, str] | None, str | None]:
        """Return gene and genomic accession

        :param genes_alt_acs: Genes and genomic accessions
        :param gene: Gene symbol
        :return: Gene and genomic accession if both exist and warning if found
        """
        alt_acs = genes_alt_acs.alt_acs
        len_alt_acs = len(alt_acs)
        if len_alt_acs > 1:
            return None, f"Found more than one accessions: {alt_acs}"
        if len_alt_acs == 0:
            return None, "No genomic accessions found"
        alt_ac = next(iter(alt_acs))

        genes = genes_alt_acs.genes
        len_genes = len(genes)
        input_gene = gene
        output_gene = None
        if len_genes == 1:
            output_gene = next(iter(genes))
        elif len_genes > 1:
            return None, f"Found more than one gene: {genes}"
        elif len_genes == 0:
            return None, "No genes found"

        if input_gene is not None and output_gene != input_gene.upper():
            return (
                None,
                f"Input gene, {input_gene}, does not match "
                f"expected output gene, {output_gene}",
            )

        gene = output_gene if output_gene else input_gene
        return (gene, alt_ac), None

    async def _set_mane_genomic_data(
        self,
        params: dict,
        gene: str,
        alt_ac: str,
        pos: int,
        is_start: bool,
    ) -> str | None:
        """Set genomic data in `params` found from MANE.

        :param params: Parameters for response
        :param gene: Gene symbol
        :param alt_ac: Genomic accession
        :param pos: Genomic position
        :param is_start: `True` if `pos` is start position. `False` if `pos` is end
            position.
        :return: Warnings if found
        """
        mane_data: (
            CdnaRepresentation | None
        ) = await self.mane_transcript.get_mane_transcript(
            alt_ac,
            pos,
            pos + 1,
            AnnotationLayer.GENOMIC,
            gene=gene,
            try_longest_compatible=True,
            coordinate_type=CoordinateType.INTER_RESIDUE,
        )
        if not mane_data:
            msg = f"Unable to find mane data for {alt_ac} with position {pos}"
            if gene:
                msg += f" on gene {gene}"
            _logger.warning(msg)
            return msg

        params["gene"] = mane_data.gene
        params["transcript"] = (
            mane_data.refseq
            if mane_data.refseq
            else mane_data.ensembl
            if mane_data.ensembl
            else None
        )

        tx_exons = await self._get_all_exon_coords(
            params["transcript"], genomic_ac=alt_ac
        )
        if not tx_exons:
            return f"No exons found given {params['transcript']}"
        tx_pos = mane_data.pos[0] + mane_data.coding_start_site
        params["exon"] = self._get_exon_number(tx_exons, tx_pos)

        try:
            tx_exon: ExonCoord = tx_exons[params["exon"] - 1]
        except IndexError:
            msg = (
                f"{params['transcript']} with position {tx_pos} "
                f"does not exist on exons: {tx_exons}"
            )
            _logger.warning(msg)
            return msg

        params["strand"] = Strand(mane_data.strand)
        self._set_exon_offset(
            params,
            tx_exon.tx_start_i,
            tx_exon.tx_end_i,
            tx_pos,
            is_start=is_start,
            strand=params["strand"],
        )

        # Need to check if we need to change pos for liftover
        genomic_data, warnings = await self.uta_db.get_alt_ac_start_or_end(
            params["transcript"], tx_pos, tx_pos, gene
        )
        if genomic_data is None:
            return warnings

        params["chr"] = genomic_data.alt_ac
        genomic_pos = (
            genomic_data.alt_end_i - 1 if is_start else genomic_data.alt_start_i + 1
        )
        params["pos"] = (
            genomic_pos - params["exon_offset"]
            if params["strand"] == Strand.NEGATIVE
            else genomic_pos + params["exon_offset"]
        )
        return None

    async def _set_genomic_data(self, params: dict, is_start: bool) -> str | None:
        """Set genomic data in ``params``

        :param params: Parameters for response
        :param is_start: ``True`` if ``pos`` is start position. ``False`` if ``pos`` is
            end position.
        :return: Warnings if found
        """
        # We should always try to liftover
        grch38_ac = await self.uta_db.get_newest_assembly_ac(params["chr"])
        if not grch38_ac:
            return f"Invalid genomic accession: {params['chr']}"

        grch38_ac = grch38_ac[0]
        if grch38_ac != params["chr"]:  # params["chr"] is genomic accession
            # Liftover to 38
            descr = await self.uta_db.get_chr_assembly(params["chr"])
            if descr is None:
                return f"Unable to get chromosome and assembly for " f"{params['chr']}"

            chromosome_number, assembly = descr
            liftover_data = self.liftover.get_liftover(
                chromosome_number, params["pos"], Assembly.GRCH38
            )
            if liftover_data is None:
                return (
                    f"Position {params['pos']} does not exist on "
                    f"chromosome {chromosome_number}"
                )

            params["pos"] = liftover_data[1]
            params["chr"] = grch38_ac

        tx_exons = await self._get_all_exon_coords(
            params["transcript"], genomic_ac=grch38_ac
        )
        if not tx_exons:
            return f"No exons found given {params['transcript']}"

        data = await self.uta_db.get_tx_exon_aln_v_data(
            params["transcript"],
            params["pos"],
            params["pos"],
            alt_ac=params["chr"],
            use_tx_pos=False,
        )
        if len(data) != 1:
            return (
                f"Must find exactly one row for genomic data, "
                f"but found: {len(data)}"
            )

        # Find exon number
        data = data[0]
        data_exons = data.tx_start_i, data.tx_end_i
        i = 1
        found_tx_exon = False
        for exon in tx_exons:
            if data_exons == (exon.tx_start_i, exon.tx_end_i):
                found_tx_exon = True
                break
            i += 1
        if not found_tx_exon:
            # Either first or last
            i = 1 if data_exons == (0, tx_exons[0].tx_end_i) else i - 1
        params["exon"] = i
        params["strand"] = Strand(data.alt_strand)
        if not is_start:
            # convert back to inter-residue for end position
            params["pos"] += 1
        self._set_exon_offset(
            params,
            data.alt_start_i
            if is_start
            else data.alt_start_i + 1,  # need to convert to inter-residue
            data.alt_end_i - 1
            if is_start
            else data.alt_end_i,  # need to convert to inter-residue
            params["pos"],
            is_start=is_start,
            strand=params["strand"],
        )
        return None

    @staticmethod
    def _set_exon_offset(
        params: dict, start: int, end: int, pos: int, is_start: bool, strand: Strand
    ) -> None:
        """Set value for ``exon_offset`` in ``params``.

        :param params: Parameters for response
        :param start: Start exon coord (can be transcript or aligned genomic)
        :param end: End exon coord (can be transcript or aligned genomic)
        :param pos: Position change (can be transcript or genomic)
        :param is_start: ``True`` if ``pos`` is start position. ``False`` if ``pos`` is
            end position
        :param strand: Strand
        """
        if is_start:
            if strand == Strand.NEGATIVE:
                params["exon_offset"] = end - pos
            else:
                params["exon_offset"] = pos - end
        else:
            if strand == Strand.NEGATIVE:
                params["exon_offset"] = start - pos
            else:
                params["exon_offset"] = pos - start

    @staticmethod
    def _get_exon_number(tx_exons: list[ExonCoord], tx_pos: int) -> int:
        """Find related exon number for a position

        :param tx_exons: List of exon coordinates for a transcript
        :param tx_pos: Transcript position change
        :return: Exon number associated to transcript position change. Will be 1-based
        """
        i = 1
        for coords in tx_exons:
            if coords.tx_start_i <= tx_pos <= coords.tx_end_i:
                break
            i += 1
        return i

    @staticmethod
    def _get_adjacent_exon(
        tx_exons_genomic_coords: list[ExonCoord],
        strand: Strand,
        start: int | None = None,
        end: int | None = None,
    ) -> int:
        """Return the adjacent exon given a non-exonic breakpoint. For the positive
        strand, adjacent is defined as the exon preceding the breakpoint for the 5' end
        and the exon following the breakpoint for the 3' end. For the negative strand,
        adjacent is defined as the exon following the breakpoint for the 5' end and the
        exon preceding the breakpoint for the 3' end.

        :param tx_exons_genomic_coords: Transcript exon coordinate data
        :param strand: Strand
        :param start: Genomic coordinate of breakpoint
        :param end: Genomic coordinate of breakpoint
        :return: Exon number corresponding to adjacent exon. Will be 1-based
        """
        for i in range(len(tx_exons_genomic_coords) - 1):
            exon = tx_exons_genomic_coords[i]
            next_exon = tx_exons_genomic_coords[i + 1]
            bp = start if start else end
            if strand == Strand.POSITIVE:
                lte_exon = exon
                gte_exon = next_exon
            else:
                lte_exon = next_exon
                gte_exon = exon
            if bp >= lte_exon.alt_end_i and bp <= gte_exon.alt_start_i:
                break
        # Return current exon if end position is provided, next exon if start position
        # is provided. exon[0] needs to be incremented by 1 in both cases as exons are
        # 0-based in UTA
        return exon.ord + 1 if end else exon.ord + 2

    @staticmethod
    def _is_exonic_breakpoint(pos: int, tx_genomic_coords: list[ExonCoord]) -> bool:
        """Check if a breakpoint occurs on an exon

        :param pos: Genomic breakpoint
        :param tx_genomic_coords: A list of transcript exon coordinate data
        :return: True if the breakpoint occurs on an exon
        """
        return any(
            exon.alt_start_i <= pos <= exon.alt_end_i for exon in tx_genomic_coords
        )
