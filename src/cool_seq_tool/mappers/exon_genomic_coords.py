"""Provide mapping capabilities between transcript exon and genomic coordinates."""

import logging

from ga4gh.vrs import models
from pydantic import ConfigDict, Field, StrictInt, StrictStr, model_validator

from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.mappers.liftover import LiftOver
from cool_seq_tool.mappers.mane_transcript import CdnaRepresentation, ManeTranscript
from cool_seq_tool.schemas import (
    AnnotationLayer,
    Assembly,
    BaseModelForbidExtra,
    CoordinateType,
    ServiceMeta,
    Strand,
)
from cool_seq_tool.sources.mane_transcript_mappings import ManeTranscriptMappings
from cool_seq_tool.sources.uta_database import GenesGenomicAcs, UtaDatabase
from cool_seq_tool.utils import service_meta

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


class TxSegment(BaseModelForbidExtra):
    """Model for representing transcript segment data."""

    exon_ord: StrictInt = Field(..., description="Exon number. 0-based.")
    offset: StrictInt = Field(
        0,
        description="The value added to or subtracted from the `genomic_location` to find the start or end of an exon.",
    )
    genomic_location: models.SequenceLocation = Field(
        ..., description="The genomic position of a transcript segment."
    )

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "exon_ord": 1,
                "exon_offset": 0,
                "genomic_location": {
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1",
                    },
                    "start": 9575887,
                },
            }
        }
    )


class _GenomicTxSeg(BaseModelForbidExtra):
    """Model for representing genomic and transcript segment data."""

    seg: TxSegment | None = Field(None, description="Transcript segment.")
    gene: StrictStr | None = Field(None, description="HGNC gene symbol.")
    genomic_ac: StrictStr | None = Field(None, description="RefSeq genomic accession.")
    tx_ac: StrictStr | None = Field(None, description="RefSeq transcript accession.")
    errors: list[StrictStr] = Field([], description="Error messages.")

    @model_validator(mode="before")
    def check_errors(cls, values: dict) -> dict:  # noqa: N805
        """Ensure that fields are (un)set depending on errors

        :param values: Values in model
        :raises ValueError: If `seg`, `gene`, `genomic_ac` and `tx_ac` are not
        provided when there are no errors
        :return: Values in model
        """
        if not values.get("errors") and not all(
            (
                values.get("seg"),
                values.get("gene"),
                values.get("genomic_ac"),
                values.get("tx_ac"),
            )
        ):
            err_msg = "`seg`, `gene`, `genomic_ac` and `tx_ac` must be provided"
            raise ValueError(err_msg)
        return values


class GenomicTxSegService(BaseModelForbidExtra):
    """Service model for genomic and transcript data."""

    gene: StrictStr | None = Field(None, description="HGNC gene symbol.")
    genomic_ac: StrictStr | None = Field(None, description="RefSeq genomic accession.")
    tx_ac: StrictStr | None = Field(None, description="RefSeq transcript accession.")
    seg_start: TxSegment | None = Field(None, description="Start transcript segment.")
    seg_end: TxSegment | None = Field(None, description="End transcript segment.")
    errors: list[StrictStr] = Field([], description="Error messages.")
    service_meta: ServiceMeta = Field(..., description="Service metadata.")

    @model_validator(mode="before")
    def add_meta_check_errors(cls, values: dict) -> dict:  # noqa: N805
        """Add service metadata to model and ensure that fields are (un)set depending
        on errors

        :param values: Values in model
        :raises ValueError: If `gene`, `genomic_ac`, `tx_ac` and `seg_start` or `seg_end`
            not provided when there are no errors
        :return: Values in model, including service metadata
        """
        values["service_meta"] = service_meta()
        if not values.get("errors") and not all(
            (
                values.get("gene"),
                values.get("genomic_ac"),
                values.get("tx_ac"),
                values.get("seg_start") or values.get("seg_end"),
            )
        ):
            err_msg = "`gene`, `genomic_ac`, `tx_ac` and `seg_start` or `seg_end` must be provided"
            raise ValueError(err_msg)

        return values

    # TODO:
    model_config = ConfigDict(json_schema_extra={"example": {}})


def _return_service_errors(errors: list[str]) -> GenomicTxSegService:
    """Log errors and return service object with errors.

    :param errors: Error message(s)
    :return: Service object with error messages.
    """
    for error in errors:
        _logger.warning(error)

    return GenomicTxSegService(errors=errors)


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

        >>> from cool_seq_tool.app import CoolSeqTool
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

    async def tx_segment_to_genomic(
        self,
        transcript: str,
        gene: str | None = None,
        exon_start: int | None = None,
        exon_start_offset: int = 0,
        exon_end: int | None = None,
        exon_end_offset: int = 0,
    ) -> GenomicTxSegService:
        """Get genomic data given transcript segment data.

        By default, transcript data is aligned to the GRCh38 assembly.

        >>> import asyncio
        >>> from cool_seq_tool.app import CoolSeqTool
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
        # Ensure valid inputs
        errors = []
        if not transcript:
            errors.append("Must provide `transcript`")
        else:
            transcript = transcript

        exon_start_exists, exon_end_exists = False, False
        if exon_start is not None:
            if exon_start < 1:
                errors.append("`exon_start` cannot be less than 1")
            exon_start_exists = True

        if exon_end is not None:
            if exon_end < 1:
                errors.append("`exon_end` cannot be less than 1")
            exon_end_exists = True

        if not exon_start_exists and not exon_end_exists:
            errors.append("Must provide either `exon_start` or `exon_end`")
        if exon_start_exists and exon_end_exists and (exon_start > exon_end):
            errors.append(
                f"Start exon {exon_start} is greater than end exon {exon_end}"
            )

        if errors:
            return _return_service_errors(errors)

        # Get exon start and exon end coordinates
        (
            tx_exon_start_coords,
            tx_exon_end_coords,
            errors,
        ) = await self._get_start_end_exon_coords(
            transcript, exon_start=exon_start, exon_end=exon_end
        )
        if errors:
            return _return_service_errors(errors)

        if gene:
            gene = gene.upper()

        # Get aligned genomic data (hgnc gene, alt_ac, alt_start_i, alt_end_i, strand)
        # for exon(s)
        alt_ac_start_end, err_msg = await self._get_alt_ac_start_and_end(
            transcript, tx_exon_start_coords, tx_exon_end_coords, gene=gene
        )
        if not alt_ac_start_end:
            return _return_service_errors([err_msg] if err_msg else [])
        alt_ac_start_data, alt_ac_end_data = alt_ac_start_end

        # Get gene and chromosome data, check that at least one was retrieved
        gene = alt_ac_start_data.hgnc if alt_ac_start_data else alt_ac_end_data.hgnc
        genomic_ac = (
            alt_ac_start_data.alt_ac if alt_ac_start_data else alt_ac_end_data.alt_ac
        )
        if gene is None or genomic_ac is None:
            return _return_service_errors(
                [
                    "Unable to retrieve `gene` or `genomic_ac` from genomic start and genomic end data"
                ],
            )

        strand = (
            Strand(alt_ac_start_data.alt_strand)
            if alt_ac_start_data
            else Strand(alt_ac_end_data.alt_strand)
        )

        if exon_start_exists:
            if strand == Strand.POSITIVE:
                genomic_pos = exon_start_offset + alt_ac_start_data.alt_start_i
            else:
                genomic_pos = alt_ac_start_data.alt_start_i + exon_start_offset
            start_genomic_loc, err_msg = self._get_vrs_seq_loc(
                genomic_ac,
                genomic_pos,
                is_start=True,
                strand=strand,
            )
            if err_msg:
                return _return_service_errors([err_msg])

            seg_start = TxSegment(
                exon_ord=alt_ac_start_data.ord,
                genomic_location=start_genomic_loc,
                offset=exon_start_offset,
            )
        else:
            seg_start = None

        if exon_end_exists:
            if strand == Strand.POSITIVE:
                genomic_pos = exon_end_offset + alt_ac_end_data.alt_end_i
            else:
                genomic_pos = alt_ac_end_data.alt_start_i + exon_end_offset
            end_genomic_loc, err_msg = self._get_vrs_seq_loc(
                genomic_ac,
                genomic_pos,
                is_start=False,
                strand=strand,
            )
            if err_msg:
                return _return_service_errors([err_msg])

            seg_end = TxSegment(
                exon_ord=alt_ac_end_data.ord,
                genomic_location=end_genomic_loc,
                offset=exon_end_offset,
            )
        else:
            seg_end = None

        return GenomicTxSegService(
            gene=gene,
            genomic_ac=genomic_ac,
            tx_ac=transcript,
            seg_start=seg_start,
            seg_end=seg_end,
        )

    async def genomic_to_tx_segment(
        self,
        chromosome: str | None = None,
        genomic_ac: str | None = None,
        genomic_start: int | None = None,
        genomic_end: int | None = None,
        transcript: str | None = None,
        get_nearest_transcript_junction: bool = False,
        gene: str | None = None,
    ) -> GenomicTxSegService:
        """Get transcript segment data for genomic data, lifted over to GRCh38.

        Must provide inter-residue coordinates.

        MANE Transcript data will be returned if and only if ``transcript`` is not
        supplied. ``gene`` must be given in order to retrieve MANE Transcript data.

        # TODO: update example
        >>> import asyncio
        >>> from cool_seq_tool.app import CoolSeqTool
        >>> from cool_seq_tool.schemas import Strand
        >>> egc = CoolSeqTool().ex_g_coords_mapper
        >>> result = asyncio.run(
        ...     egc.genomic_to_tx_segment(
        ...         genomic_ac="NC_000001.11",
        ...         genomic_start=154192136,
        ...         genomic_end=154170400,
        ...         transcript="NM_152263.3",
        ...     )
        ... )
        >>> result.seg_start.exon_ord, result.seg_end.exon_ord
        (0, 7)

        :param chromosome: e.g. ``"1"`` or ``"chr1"``. If not provided, must provide
            ``genomic_ac``. If ``genomic_ac`` is also provided, ``genomic_ac`` will be
            used.
        :param genomic_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided,
            ``genomic_ac`` will be used.
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
        errors = []
        if genomic_start is None and genomic_end is None:
            errors.append("Must provide either `genomic_start` or `genomic_end`")
        if chromosome is None and genomic_ac is None:
            errors.append("Must provide either `chromosome` or `alt_ac`")
        if transcript is None and gene is None:
            errors.append("Must provide either `gene` or `transcript`")
        if errors:
            return _return_service_errors(errors)

        if gene is not None:
            gene = gene.upper()

        params = {}

        if genomic_start:
            start_tx_seg_data = await self._genomic_to_tx_segment(
                genomic_start,
                chromosome=chromosome,
                genomic_ac=genomic_ac,
                transcript=transcript,
                gene=gene,
                get_nearest_transcript_junction=get_nearest_transcript_junction,
                is_start=True,
            )
            if start_tx_seg_data.errors:
                return _return_service_errors(start_tx_seg_data.errors)

            params["gene"] = start_tx_seg_data.gene
            params["genomic_ac"] = start_tx_seg_data.genomic_ac
            params["tx_ac"] = start_tx_seg_data.tx_ac
            params["seg_start"] = start_tx_seg_data.seg
        else:
            start_tx_seg_data = None

        if genomic_end:
            end_tx_seg_data = await self._genomic_to_tx_segment(
                genomic_end,
                chromosome=chromosome,
                genomic_ac=genomic_ac,
                transcript=transcript,
                gene=gene,
                get_nearest_transcript_junction=get_nearest_transcript_junction,
                is_start=False,
            )
            if end_tx_seg_data.errors:
                return _return_service_errors(end_tx_seg_data.errors)

            if start_tx_seg_data:
                # Need to check that gene, genomic_ac, tx_ac all match
                errors = []
                for attr in ["gene", "genomic_ac", "tx_ac"]:
                    start_seg_attr = params[attr]
                    end_seg_attr = getattr(end_tx_seg_data, attr)
                    if start_seg_attr != end_seg_attr:
                        errors.append(
                            f"Start end end segment mismatch for `{attr}`. {start_seg_attr} != {end_seg_attr}."
                        )
                if errors:
                    return _return_service_errors(errors)
            else:
                params["gene"] = end_tx_seg_data.gene
                params["genomic_ac"] = end_tx_seg_data.genomic_ac
                params["tx_ac"] = end_tx_seg_data.tx_ac

            params["seg_end"] = end_tx_seg_data.seg

        return GenomicTxSegService(**params)

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
        """
        if genomic_ac:
            query = f"""
                SELECT DISTINCT ord, tx_start_i, tx_end_i, alt_start_i, alt_end_i, alt_strand
                FROM {self.uta_db.schema}.tx_exon_aln_v
                WHERE tx_ac = '{tx_ac}'
                AND alt_aln_method = 'splign'
                AND alt_ac = '{genomic_ac}'
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
                    if attr == "hgnc":
                        error = "HGNC gene symbol does not match"
                    elif attr == "alt_ac":
                        error = "Genomic accession does not match"
                    else:
                        error = "Strand does not match"
                    _logger.warning(
                        "%s: %s != %s",
                        error,
                        start_attr,
                        end_attr,
                    )
                    return None, error
        return tuple(alt_ac_data_values), None

    async def _genomic_to_tx_segment(
        self,
        genomic_pos: int,
        chromosome: str | None = None,
        genomic_ac: str | None = None,
        transcript: str | None = None,
        gene: str | None = None,
        get_nearest_transcript_junction: bool = False,
        is_start: bool = True,
    ) -> _GenomicTxSeg:
        """Convert individual genomic data to transcript data

        :param genomic_pos: Genomic position where the transcript segment starts or ends
            (inter-residue based)
        :param chromosome: Chromosome. Must give chromosome without a prefix
            (i.e. ``1`` or ``X``). If not provided, must provide ``genomic_ac``.
            If ``genomic_ac`` is also provided, ``genomic_ac`` will be used.
        :param genomic_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``genomic_ac``
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
        params = {key: None for key in _GenomicTxSeg.model_fields}

        if get_nearest_transcript_junction:
            if not gene:
                return _GenomicTxSeg(
                    errors=[
                        "`gene` must be provided to select the adjacent transcript junction"
                    ]
                )

            if not genomic_ac:
                genomic_acs, err_msg = self.seqrepo_access.chromosome_to_acs(chromosome)

                if not genomic_acs:
                    return _GenomicTxSeg(
                        errors=[err_msg],
                    )
                genomic_ac = genomic_acs[0]

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
                        gene=gene, alt_ac=genomic_ac
                    )

                    if not results.is_empty():
                        transcript = results[0]["tx_ac"][0]
                    else:
                        # Run if gene is for a noncoding transcript
                        query = f"""
                            SELECT DISTINCT tx_ac
                            FROM {self.uta_db.schema}.tx_exon_aln_v
                            WHERE hgnc = '{gene}'
                            AND alt_ac = '{genomic_ac}'
                            """  # noqa: S608
                        result = await self.uta_db.execute_query(query)

                        if result:
                            transcript = result[0]["tx_ac"]
                        else:
                            return _GenomicTxSeg(
                                errors=[
                                    f"Could not find a transcript for {gene} on {genomic_ac}"
                                ]
                            )

            tx_exons = await self._get_all_exon_coords(
                tx_ac=transcript, genomic_ac=genomic_ac
            )
            if not tx_exons:
                return _GenomicTxSeg(errors=[f"No exons found given {transcript}"])

            strand = Strand(tx_exons[0].alt_strand)
            params["strand"] = strand

            # Check if breakpoint occurs on an exon.
            # If not, determine the adjacent exon given the selected transcript
            if not self._is_exonic_breakpoint(genomic_pos, tx_exons):
                exon_num = self._get_adjacent_exon(
                    tx_exons_genomic_coords=tx_exons,
                    strand=strand,
                    start=genomic_pos if is_start else None,
                    end=genomic_pos if not is_start else None,
                )

                offset = self._get_exon_offset(
                    start_i=tx_exons[exon_num].alt_start_i,
                    end_i=tx_exons[exon_num].alt_end_i,
                    strand=strand,
                    use_start_i=strand == Strand.POSITIVE
                    if is_start
                    else strand != Strand.POSITIVE,
                    is_in_exon=False,
                    start=genomic_pos if is_start else None,
                    end=genomic_pos if not is_start else None,
                )

                genomic_location, err_msg = self._get_vrs_seq_loc(
                    genomic_ac, genomic_pos, is_start, strand
                )
                if err_msg:
                    return _GenomicTxSeg(errors=[err_msg])

                return _GenomicTxSeg(
                    gene=gene,
                    genomic_ac=genomic_ac,
                    tx_ac=transcript,
                    seg=TxSegment(
                        exon_ord=exon_num,
                        offset=offset,
                        genomic_location=genomic_location,
                    ),
                )

        if genomic_ac:
            # Check if valid accession is given
            if not await self.uta_db.validate_genomic_ac(genomic_ac):
                return _GenomicTxSeg(
                    errors=[f"Invalid genomic accession: {genomic_ac}"]
                )

            genes_alt_acs, warning = await self.uta_db.get_genes_and_alt_acs(
                genomic_pos, alt_ac=genomic_ac, gene=gene
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
            return _GenomicTxSeg(errors=[warning])

        gene_genomic_ac, warning = self._get_gene_and_alt_ac(genes_alt_acs, gene)
        if not gene_genomic_ac:
            return _GenomicTxSeg(errors=[warning])
        gene, genomic_ac = gene_genomic_ac

        if transcript is None:
            return await self._get_mane_genomic_tx_seg(
                gene, genomic_ac, genomic_pos, is_start
            )

        return await self._get_genomic_tx_seg(
            genomic_ac, genomic_pos, transcript, is_start, gene
        )

    def _get_vrs_seq_loc(
        self, genomic_ac: str, genomic_pos: int, is_start: bool, strand: Strand
    ) -> tuple[models.SequenceLocation | None, str | None]:
        """Create VRS Sequence Location for genomic position where transcript segment
        occurs

        :param genomic_ac: RefSeq genomic accession
        :param genomic_pos: Genomic position where the transcript segment occurs
        :param is_start: ``True`` if ``genomic_pos`` is where the transcript segment
            starts. ``False`` if ``genomic_pos`` is where the transcript segment ends.
        :param strand: Strand
        :return: Tuple containing VRS location (if successful) and error message (if
            unable to get GA4GH identifier for ``genomic_ac``).
        """
        ga4gh_seq_id, err_msg = self.seqrepo_access.translate_identifier(
            genomic_ac, "ga4gh"
        )
        if err_msg:
            return None, err_msg

        use_start = strand == Strand.POSITIVE if is_start else strand != Strand.POSITIVE

        return models.SequenceLocation(
            sequenceReference=models.SequenceReference(
                refgetAccession=ga4gh_seq_id[0].split("ga4gh:")[-1]
            ),
            start=genomic_pos if use_start else None,
            end=genomic_pos if not use_start else None,
        ), None

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

    async def _get_mane_genomic_tx_seg(
        self,
        gene: str,
        genomic_ac: str,
        pos: int,
        is_start: bool,
    ) -> _GenomicTxSeg:
        """Set genomic data in `params` found from MANE.

        :param gene: Gene symbol
        :param genomic_ac: Genomic accession
        :param pos: Genomic position
        :param is_start: `True` if `pos` is start position. `False` if `pos` is end
            position.
        :return: Warnings if found
        """
        mane_data: (
            CdnaRepresentation | None
        ) = await self.mane_transcript.get_mane_transcript(
            genomic_ac,
            pos,
            pos,
            AnnotationLayer.GENOMIC,
            gene=gene,
            try_longest_compatible=True,
            coordinate_type=CoordinateType.INTER_RESIDUE,
        )
        if not mane_data:
            err_msg = f"Unable to find mane data for {genomic_ac} with position {pos}"
            if gene:
                err_msg += f" on gene {gene}"
            _logger.warning(err_msg)
            return _GenomicTxSeg(errors=[err_msg])

        mane_gene = mane_data.gene
        tx_ac = (
            mane_data.refseq
            if mane_data.refseq
            else mane_data.ensembl
            if mane_data.ensembl
            else None
        )
        tx_exons = await self._get_all_exon_coords(tx_ac, genomic_ac)
        if not tx_exons:
            return _GenomicTxSeg(errors=[f"No exons found given {tx_ac}"])

        tx_pos = mane_data.pos[0] + mane_data.coding_start_site
        exon_num = self._get_exon_number(tx_exons, tx_pos)
        exon = tx_exons[exon_num]

        if exon_num == -1:
            err_msg = (
                f"{tx_ac} with position {tx_pos} does not exist on exons: {tx_exons}"
            )
            _logger.warning(err_msg)
            return _GenomicTxSeg(errors=[err_msg])

        # Need to check if we need to change pos for liftover
        genomic_data, err_msg = await self.uta_db.get_alt_ac_start_or_end(
            tx_ac, tx_pos, tx_pos, gene
        )
        if genomic_data is None:
            return _GenomicTxSeg(errors=[err_msg])

        strand = Strand(mane_data.strand)
        if strand == Strand.POSITIVE:
            offset = tx_pos - exon.tx_start_i if is_start else tx_pos - exon.tx_end_i
        else:
            offset = tx_pos - exon.tx_end_i if is_start else tx_pos - exon.tx_start_i

        if strand == Strand.POSITIVE:
            genomic_pos = (
                genomic_data.alt_start_i + offset
                if is_start
                else genomic_data.alt_end_i + offset
            )
        else:
            genomic_pos = (
                genomic_data.alt_start_i - offset
                if is_start
                else genomic_data.alt_end_i - offset
            )

        genomic_ac = genomic_data.alt_ac
        genomic_location, err_msg = self._get_vrs_seq_loc(
            genomic_ac, genomic_pos, is_start, strand
        )

        if err_msg:
            return _GenomicTxSeg(errors=[err_msg])

        return _GenomicTxSeg(
            gene=mane_gene,
            genomic_ac=genomic_ac,
            tx_ac=tx_ac,
            seg=TxSegment(
                exon_ord=exon_num,
                offset=offset,
                genomic_location=genomic_location,
            ),
        )

    async def _get_genomic_tx_seg(
        self,
        genomic_ac: str,
        genomic_pos: int,
        tx_ac: str,
        is_start: bool,
        gene: str,
    ) -> _GenomicTxSeg:
        # We should always try to liftover
        grch38_ac = await self.uta_db.get_newest_assembly_ac(genomic_ac)
        if not grch38_ac:
            return _GenomicTxSeg(errors=[f"Invalid genomic accession: {genomic_ac}"])

        grch38_ac = grch38_ac[0]
        if grch38_ac != genomic_ac:  # alt_ac is genomic accession
            # Liftover to 38
            descr = await self.uta_db.get_chr_assembly(genomic_ac)
            if descr is None:
                return _GenomicTxSeg(
                    errors=[f"Unable to get chromosome and assembly for {genomic_ac}"]
                )

            chromosome_number, assembly = descr
            liftover_data = self.liftover.get_liftover(
                chromosome_number, genomic_pos, Assembly.GRCH38
            )
            if liftover_data is None:
                return _GenomicTxSeg(
                    errors=[
                        f"Position {genomic_pos} does not exist on chromosome {chromosome_number}"
                    ]
                )

            genomic_pos = liftover_data[1]
            genomic_ac = grch38_ac

        tx_exons = await self._get_all_exon_coords(tx_ac, genomic_ac=grch38_ac)
        if not tx_exons:
            return _GenomicTxSeg(errors=[f"No exons found given {tx_ac}"])

        tx_exon_aln_data = await self.uta_db.get_tx_exon_aln_v_data(
            tx_ac,
            genomic_pos,
            genomic_pos,
            alt_ac=genomic_ac,
            use_tx_pos=False,
        )
        if len(tx_exon_aln_data) != 1:
            return _GenomicTxSeg(
                errors=[
                    f"Must find exactly one row for genomic data, but found: {len(tx_exon_aln_data)}"
                ]
            )

        tx_exon_aln_data = tx_exon_aln_data[0]
        if tx_exon_aln_data.hgnc != gene:
            return _GenomicTxSeg(
                errors=[f"Expected gene, {gene}, but found {tx_exon_aln_data.hgnc}"]
            )

        offset = self._get_exon_offset(
            start_i=tx_exon_aln_data.alt_start_i,
            end_i=tx_exon_aln_data.alt_end_i,
            strand=Strand(tx_exon_aln_data.alt_strand),
            use_start_i=False,  # This doesn't impact anything since we're on the exon
            is_in_exon=True,
            start=genomic_pos if is_start else None,
            end=genomic_pos if not is_start else None,
        )

        genomic_location, err_msg = self._get_vrs_seq_loc(
            genomic_ac, genomic_pos, is_start, tx_exon_aln_data.alt_strand
        )
        if err_msg:
            return _GenomicTxSeg(errors=[err_msg])

        return _GenomicTxSeg(
            gene=tx_exon_aln_data.hgnc,
            genomic_ac=genomic_ac,
            tx_ac=tx_exon_aln_data.tx_ac,
            seg=TxSegment(
                exon_ord=tx_exon_aln_data.ord,
                offset=offset,
                genomic_location=genomic_location,
            ),
        )

    @staticmethod
    def _get_exon_offset(
        start_i: int,
        end_i: int,
        strand: Strand,
        use_start_i: bool = True,
        is_in_exon: bool = True,
        start: int | None = None,
        end: int | None = None,
    ) -> None:
        if is_in_exon:
            if start is not None:
                offset = (
                    start - start_i if strand == Strand.POSITIVE else start_i - start
                )
            else:
                offset = end - end_i if strand == Strand.POSITIVE else end_i - end
        else:
            if strand == Strand.POSITIVE:
                offset = start - start_i if use_start_i else end - end_i
            else:
                offset = start_i - end if use_start_i else end_i - start
        return offset

    @staticmethod
    def _get_exon_number(tx_exons: list[ExonCoord], tx_pos: int) -> int:
        """Find related exon number for a position

        :param tx_exons: List of exon coordinates for a transcript
        :param tx_pos: Transcript position change
        :return: Exon number associated to transcript position change. Will be 0-based.
            If there is no exon associated to the transcript position change, -1 will
            be returned.
        """
        for coords in tx_exons:
            if coords.tx_start_i <= tx_pos < coords.tx_end_i:
                return coords.ord
        return -1

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
        :return: Exon number corresponding to adjacent exon. Will be 0-based
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
        # is provided.
        return exon.ord if end else exon.ord + 1

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
