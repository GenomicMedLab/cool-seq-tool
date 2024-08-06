"""Provide mapping capabilities between transcript exon and genomic coordinates."""

import logging
from typing import ClassVar

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
from cool_seq_tool.sources.uta_database import (
    GenesGenomicAcs,
    GenomicAlnData,
    UtaDatabase,
)
from cool_seq_tool.utils import service_meta

_logger = logging.getLogger(__name__)


class ExonCoord(BaseModelForbidExtra):
    """Model for representing exon data"""

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
    """Model for representing transcript segment data"""

    exon_ord: StrictInt = Field(..., description="Exon number. 0-based.")
    offset: StrictInt = Field(0, description="TODO")
    genomic_location: models.SequenceLocation = Field(..., description="TODO")

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


class _SeqLocAndError(BaseModelForbidExtra):
    """Model for VRS Sequence Location and errors"""

    errors: ClassVar[list[StrictStr]] = []
    sequence_location: models.SequenceLocation | None = None


class _GenomicTxSeg(BaseModelForbidExtra):
    """Model for representing genomic and transcript segment data"""

    seg: TxSegment | None = Field(None, description="Transcript segment.")
    gene: StrictStr | None = Field(None, description="HGNC gene symbol.")
    genomic_accession: StrictStr | None = Field(
        None, description="RefSeq genomic accession."
    )
    tx_ac: StrictStr | None = Field(None, description="RefSeq transcript accession.")
    errors: list[StrictStr] = Field([], description="Error messages.")

    @model_validator(mode="after")
    def check_errors(cls, values: "_GenomicTxSeg") -> "_GenomicTxSeg":  # noqa: N805
        if not values.errors and not all(
            (values.seg, values.gene, values.genomic_accession, values.tx_ac)
        ):
            err_msg = "`seg`, `gene`, `genomic_accession` and `tx_ac` must be provided"
            raise ValueError(err_msg)
        return values


class GenomicTxSegService(BaseModelForbidExtra):
    """Model for Genomic and Transcript Segment Data"""

    gene: StrictStr | None = Field(None, description="HGNC gene symbol.")
    genomic_accession: StrictStr | None = Field(
        None, description="RefSeq genomic accession."
    )
    tx_ac: StrictStr | None = Field(None, description="RefSeq transcript accession.")
    seg_start: TxSegment | None = Field(None, description="Start transcript segment.")
    seg_end: TxSegment | None = Field(None, description="End transcript segment.")
    errors: list[StrictStr] = Field([], description="Error messages.")
    service_meta: ServiceMeta = Field(..., description="Service metadata.")

    @model_validator(mode="before")
    def add_service_meta(cls, values: dict) -> dict:  # noqa: N805
        """Add service metadata to model"""
        values["service_meta"] = service_meta()
        return values

    @model_validator(mode="after")
    def check_errors(cls, values: "GenomicTxSegService") -> "GenomicTxSegService":  # noqa: N805
        """Ensure that fields are (un)set depending on errors

        :param values: Values in model
        :raises ValueError: If `gene`, `genomic_accession`, `tx_ac` and `seg_start` or `seg_end`
            not provided when there are no errors
        :return: Values in model
        """
        if not values.errors and not all(
            (
                values.gene,
                values.genomic_accession,
                values.tx_ac,
                values.seg_start or values.seg_end,
            )
        ):
            err_msg = "`gene`, `genomic_accession`, `tx_ac` and `seg_start` or `seg_end` must be provided"
            raise ValueError(err_msg)

        return values

    model_config = ConfigDict(json_schema_extra={"example": {}})


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

    @staticmethod
    def _return_errors(messages: list[str]) -> GenomicTxSegService:
        """Add warnings to response object

        :param messages: Error essage(s) on why ``transcript_exon_data`` or
            ``genomic_data`` field is ``None``
        :return: Response object with warning message
        """
        errors = []
        for msg in messages:
            _logger.warning(msg)
            errors.append(msg)
        return GenomicTxSegService(errors=errors)

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
            return self._return_errors(errors)

        # Get all exons and associated start/end coordinates for transcript
        tx_exons = await self._get_exons_coords(transcript)
        if not tx_exons:
            return self._return_errors([f"Unable to get exons for {transcript}"])

        # Get exon start and exon end coordinates
        tx_exon_coords, err_msg = self._get_tx_exon_coords(
            transcript, tx_exons, exon_start, exon_end
        )
        if not tx_exon_coords:
            return self._return_errors([err_msg] if err_msg else [])
        tx_exon_start_coords, tx_exon_end_coords = tx_exon_coords

        if gene:
            gene = gene.upper().strip()

        # Get aligned genomic data (hgnc gene, alt_ac, alt_start_i, alt_end_i, strand)
        # for exon(s)
        alt_ac_start_end, err_msg = await self._get_alt_ac_start_and_end(
            transcript, tx_exon_start_coords, tx_exon_end_coords, gene=gene
        )
        if not alt_ac_start_end:
            return self._return_errors([err_msg] if err_msg else [])
        alt_ac_start_data, alt_ac_end_data = alt_ac_start_end

        # Get gene and chromosome data, check that at least one was retrieved
        gene = alt_ac_start_data.hgnc if alt_ac_start_data else alt_ac_end_data.hgnc
        genomic_accession = (
            alt_ac_start_data.alt_ac if alt_ac_start_data else alt_ac_end_data.alt_ac
        )
        if gene is None or genomic_accession is None:
            return self._return_errors(
                [
                    "Unable to retrieve `hgnc` or `genomic_accession` from genomic start and genomic end data"
                ],
            )

        strand = (
            Strand(alt_ac_start_data.alt_strand)
            if alt_ac_start_data
            else Strand(alt_ac_end_data.alt_strand)
        )

        if exon_start is not None:
            start_genomic_loc = self._get_seq_loc(
                genomic_accession,
                exon_start_offset + alt_ac_start_data.alt_start_i,
                True,
                strand,
            )
            if not start_genomic_loc:
                return self._return_errors(start_genomic_loc.errors)

            seg_start = TxSegment(
                exon_ord=alt_ac_start_data.ord,
                genomic_location=start_genomic_loc.sequence_location,
                offset=exon_start_offset,
            )
        else:
            seg_start = None

        if exon_end is not None:
            end_genomic_loc = self._get_seq_loc(
                genomic_accession,
                exon_end_offset + alt_ac_end_data.alt_end_i,
                False,
                strand,
            )
            if not end_genomic_loc:
                return self._return_errors(end_genomic_loc.errors)

            seg_end = TxSegment(
                exon_ord=alt_ac_end_data.ord,
                genomic_location=end_genomic_loc.sequence_location,
                offset=exon_end_offset,
            )
        else:
            seg_end = None

        return GenomicTxSegService(
            gene=gene,
            genomic_accession=genomic_accession,
            tx_ac=transcript,
            seg_start=seg_start,
            seg_end=seg_end,
        )

    async def genomic_to_tx_segment(
        self,
        chromosome: str | None = None,
        genomic_accession: str | None = None,
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

        >>> import asyncio
        >>> from cool_seq_tool.app import CoolSeqTool
        >>> from cool_seq_tool.schemas import Strand
        >>> egc = CoolSeqTool().ex_g_coords_mapper
        >>> result = asyncio.run(
        ...     egc.genomic_to_tx_segment(
        ...         genomic_accession="NC_000001.11",
        ...         genomic_start=154192136,
        ...         genomic_end=154170400,
        ...         transcript="NM_152263.3",
        ...     )
        ... )
        >>> result.genomic_data.exon_start, result.genomic_data.exon_end
        (1, 8)

        :param chromosome: e.g. ``"1"`` or ``"chr1"``. If not provided, must provide
            ``genomic_accession``. If ``genomic_accession`` is also provided, ``genomic_accession`` will be used.
        :param genomic_accession: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``genomic_accession``
            will be used.
        :param genomic_start: Genomic position where the transcript segment starts.
            Must provide inter-residue coordinates.
        :param genomic_end: Genomic position where the transcript segment ends.
            Must provide inter-residue coordinates.
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript. See the :ref:`Transcript Selection policy <transcript_selection_policy>`
            page.
        :param get_nearest_transcript_junction: If ``True``, this will return the
            adjacent exon if the position specified by``genomic_start`` or ``genomic_end`` does
            not occur on an exon. For the positive strand, adjacent is defined as the
            exon preceding the breakpoint for the 5' end and the exon following the
            breakpoint for the 3' end. For the negative strand, adjacent is defined as
            the exon following the breakpoint for the 5' end and the exon preceding the
            breakpoint for the 3' end.
        :param gene: gene name. Ideally, HGNC symbol. Must be given if no ``transcript``
            value is provided.
        :param coordinate_type: Coordinate type for ``genomic_start`` and ``genomic_end``
        :return: Genomic data (inter-residue coordinates)
        """
        errors = []
        if genomic_start is None and genomic_end is None:
            errors.append("Must provide either `genomic_start` or `genomic_end`")
        if chromosome is None and genomic_accession is None:
            errors.append("Must provide either `chromosome` or `genomic_accession`")
        if transcript is None and gene is None:
            errors.append("Must provide either `gene` or `transcript`")
        if errors:
            return self._return_errors(errors)

        if gene is not None:
            gene = gene.upper().strip()

        params = {}

        if genomic_start:
            start_tx_seg_data = await self._genomic_to_tx_segment(
                genomic_start,
                chromosome=chromosome,
                genomic_accession=genomic_accession,
                transcript=transcript,
                gene=gene,
                get_nearest_transcript_junction=get_nearest_transcript_junction,
                is_start=True,
            )
            if start_tx_seg_data.errors:
                return self._return_errors(start_tx_seg_data.errors)

            params["gene"] = start_tx_seg_data.gene
            params["genomic_accession"] = start_tx_seg_data.genomic_accession
            params["tx_ac"] = start_tx_seg_data.tx_ac
            params["seg_start"] = start_tx_seg_data.seg
        else:
            start_tx_seg_data = None

        if genomic_end:
            end_tx_seg_data = await self._genomic_to_tx_segment(
                genomic_end,
                chromosome=chromosome,
                genomic_accession=genomic_accession,
                transcript=transcript,
                gene=gene,
                get_nearest_transcript_junction=get_nearest_transcript_junction,
                is_start=False,
            )
            if end_tx_seg_data.errors:
                return self._return_errors(end_tx_seg_data.errors)

            if start_tx_seg_data:
                # Need to check that gene, genomic_accession, tx_ac all match
                errors = []
                for attr in ["gene", "genomic_accession", "tx_ac"]:
                    start_seg_attr = params[attr]
                    end_seg_attr = getattr(end_tx_seg_data, attr)
                    if start_seg_attr != end_seg_attr:
                        errors.append(
                            f"Start end end segment mismatch for `{attr}`. {start_seg_attr} != {end_seg_attr}."
                        )
                if errors:
                    return self._return_errors(errors)
            else:
                params["gene"] = end_tx_seg_data.gene
                params["genomic_accession"] = end_tx_seg_data.genomic_accession
                params["tx_ac"] = end_tx_seg_data.tx_ac

            params["seg_end"] = end_tx_seg_data.seg

        return GenomicTxSegService(**params)

    @staticmethod
    def _validate_exon(
        transcript: str, tx_exons: list[ExonCoord], exon_number: int
    ) -> tuple[tuple[int, int] | None, str | None]:
        """Validate that exon number exists on a given transcript

        :param transcript: Transcript accession
        :param tx_exons: List of transcript's exons and associated coordinates
        :param exon_number: Exon number to validate. 1-based.
        :return: Exon coordinates for a given exon number and warnings if found
        """
        msg = f"Exon {exon_number} does not exist on {transcript}"
        try:
            if exon_number < 1:
                return None, msg
            exon = tx_exons[exon_number - 1]
        except IndexError:
            return None, msg
        return exon, None

    def _get_tx_exon_coords(
        self,
        transcript: str,
        tx_exons: list[ExonCoord],
        exon_start: int | None = None,
        exon_end: int | None = None,
    ) -> tuple[
        tuple[ExonCoord | None, ExonCoord | None] | None,
        str | None,
    ]:
        """Get exon coordinates for ``exon_start`` and ``exon_end``

        :param transcript: Transcript accession
        :param tx_exons: List of all transcript exons and coordinates
        :param exon_start: Start exon number. 1-based.
        :param exon_end: End exon number. 1-based.
        :return: [Transcript start exon coords, Transcript end exon coords],
            and warnings if found
        """
        if exon_start is not None:
            tx_exon_start, warning = self._validate_exon(
                transcript, tx_exons, exon_start
            )
            if not tx_exon_start:
                return None, warning
        else:
            tx_exon_start = None

        if exon_end is not None:
            tx_exon_end, warning = self._validate_exon(transcript, tx_exons, exon_end)
            if not tx_exon_end:
                return None, warning
        else:
            tx_exon_end = None
        return (tx_exon_start, tx_exon_end), None

    async def _get_alt_ac_start_and_end(
        self,
        tx_ac: str,
        tx_exon_start: tuple[int, int] | None = None,
        tx_exon_end: tuple[int, int] | None = None,
        gene: str | None = None,
    ) -> tuple[GenomicAlnData | None, str | None]:
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
        genomic_accession: str | None = None,
        transcript: str | None = None,
        gene: str | None = None,
        get_nearest_transcript_junction: bool = False,
        is_start: bool = True,
    ) -> _GenomicTxSeg:
        """Convert individual genomic data to transcript data

        :param genomic_pos: Genomic position where the transcript segment starts or ends
            (inter-residue based)
        :param chromosome: Chromosome. Must give chromosome without a prefix
            (i.e. ``1`` or ``X``). If not provided, must provide ``genomic_accession``.
            If ``genomic_accession`` is also provided, ``genomic_accession`` will be used.
        :param genomic_accession: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``genomic_accession``
            will be used.
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript
        :param gene: HGNC gene symbol
        :param get_nearest_transcript_junction: If ``True``, this will return the
            adjacent exon if the position specified by``genomic_start`` or ``genomic_end`` does
            not occur on an exon. For the positive strand, adjacent is defined as the
            exon preceding the breakpoint for the 5' end and the exon following the
            breakpoint for the 3' end. For the negative strand, adjacent is defined as
            the exon following the breakpoint for the 5' end and the exon preceding the
            breakpoint for the 3' end.
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
            if not genomic_accession:
                alt_acs, err_msg = self.seqrepo_access.chromosome_to_acs(chromosome)

                if not alt_acs:
                    return _GenomicTxSeg(
                        errors=[err_msg],
                    )
                genomic_accession = alt_acs[0]

            if not transcript:
                # Select a transcript if not provided
                mane_transcripts = self.mane_transcript_mappings.get_gene_mane_data(
                    gene
                )

                if mane_transcripts:
                    transcript = mane_transcripts[0]["RefSeq_nuc"]
                    genomic_accession = mane_transcripts[0]["GRCh38_chr"]
                else:
                    # Attempt to find a coding transcript if a MANE transcript
                    # cannot be found
                    results = await self.uta_db.get_transcripts(
                        gene=gene, alt_ac=genomic_accession
                    )

                    if not results.is_empty():
                        transcript = results[0]["tx_ac"][0]
                    else:
                        # Run if gene is for a noncoding transcript
                        query = f"""
                            SELECT DISTINCT tx_ac
                            FROM {self.uta_db.schema}.tx_exon_aln_v
                            WHERE hgnc = '{gene}'
                            AND alt_ac = '{genomic_accession}'
                            """  # noqa: S608
                        result = await self.uta_db.execute_query(query)

                        if result:
                            transcript = result[0]["tx_ac"]
                        else:
                            return _GenomicTxSeg(
                                errors=[
                                    f"Could not find a transcript for {gene} on {genomic_accession}"
                                ]
                            )

            tx_exons = await self._get_exons_coords(
                tx_ac=transcript, genomic_accession=genomic_accession
            )
            if not tx_exons:
                return _GenomicTxSeg(errors=[f"Unable to get exons for {transcript}"])

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

                use_start = (
                    strand == Strand.POSITIVE if is_start else strand != Strand.POSITIVE
                )

                offset = self._get_exon_offset(
                    start=tx_exons[exon_num].alt_start_i,  # Start exon coordinate
                    end=tx_exons[exon_num].alt_end_i,  # End exon coordinate
                    pos=genomic_pos,
                    is_start=use_start,
                )

                genomic_location_data = self._get_seq_loc(
                    genomic_accession, genomic_pos, is_start, strand
                )
                if genomic_location_data.errors:
                    return _GenomicTxSeg(errors=genomic_location_data.errors)

                return _GenomicTxSeg(
                    gene=gene,
                    genomic_accession=genomic_accession,
                    tx_ac=transcript,
                    seg=TxSegment(
                        exon_ord=exon_num,
                        offset=offset,
                        genomic_location=genomic_location_data.sequence_location,
                    ),
                )

        if genomic_accession:
            # Check if valid accession is given
            if not await self.uta_db.validate_genomic_ac(genomic_accession):
                return _GenomicTxSeg(
                    errors=[f"Invalid genomic accession: {genomic_accession}"]
                )

            genes_alt_acs, warning = await self.uta_db.get_genes_and_alt_acs(
                genomic_pos, alt_ac=genomic_accession, gene=gene
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

        gene_alt_ac, warning = self._get_gene_and_alt_ac(genes_alt_acs, gene)
        if not gene_alt_ac:
            return _GenomicTxSeg(errors=[warning])
        gene, genomic_accession = gene_alt_ac

        if transcript is None:
            return await self._get_mane_genomic_tx_seg(
                gene, genomic_accession, genomic_pos, is_start
            )

        return await self._get_genomic_tx_seg(
            genomic_accession, genomic_pos, transcript, is_start, gene
        )

    def _get_seq_loc(
        self, genomic_accession: str, genomic_pos: int, is_start: bool, strand: Strand
    ) -> _SeqLocAndError:
        ga4gh_seq_id, err_msg = self.seqrepo_access.translate_identifier(
            genomic_accession, "ga4gh"
        )
        if err_msg:
            return _SeqLocAndError(errors=[err_msg])

        use_start = strand == Strand.POSITIVE if is_start else strand != Strand.POSITIVE

        return _SeqLocAndError(
            sequence_location=models.SequenceLocation(
                sequenceReference=models.SequenceReference(
                    refgetAccession=ga4gh_seq_id[0].split("ga4gh:")[-1]
                ),
                start=genomic_pos if use_start else None,
                end=genomic_pos if not use_start else None,
            )
        )

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
        genomic_accession = next(iter(alt_acs))

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
        return (gene, genomic_accession), None

    @staticmethod
    def _get_exon_number(tx_exons: list[ExonCoord], tx_pos: int) -> int:
        """Find related exon number for a position

        :param tx_exons: List of exon coordinates for a transcript
        :param tx_pos: Transcript position change
        :return: Exon number associated to transcript position change. Will be 1-based
        """
        for coords in tx_exons:
            if coords.tx_start_i <= tx_pos < coords.tx_end_i:
                return coords.ord

        return -1

    async def _get_mane_genomic_tx_seg(
        self,
        gene: str,
        genomic_accession: str,
        pos: int,
        is_start: bool,
    ) -> _GenomicTxSeg:
        """Set genomic data in `params` found from MANE.

        :param gene: Gene symbol
        :param genomic_accession: Genomic accession
        :param pos: Genomic position
        :param is_start: `True` if `pos` is start position. `False` if `pos` is end
            position.
        :return: Warnings if found
        """
        mane_data: (
            CdnaRepresentation | None
        ) = await self.mane_transcript.get_mane_transcript(
            genomic_accession,
            pos,
            pos,
            AnnotationLayer.GENOMIC,
            gene=gene,
            try_longest_compatible=True,
            coordinate_type=CoordinateType.INTER_RESIDUE,
        )
        if not mane_data:
            err_msg = (
                f"Unable to find mane data for {genomic_accession} with position {pos}"
            )
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
        tx_exons = await self._get_exons_coords(tx_ac, genomic_accession)
        if not tx_exons:
            return _GenomicTxSeg(errors=[f"Unable to get exons for {tx_ac}"])

        tx_pos = mane_data.pos[0] + mane_data.coding_start_site
        exon_num = self._get_exon_number(tx_exons, tx_pos)
        exon = tx_exons[exon_num]

        if exon_num == -1:
            err_msg = (
                f"{tx_ac} with position {tx_pos} does not exist on exons: {tx_exons}"
            )
            _logger.warning(err_msg)
            return _GenomicTxSeg(errors=[err_msg])

        strand = Strand(mane_data.strand)
        if strand == Strand.POSITIVE:
            offset = tx_pos - exon.tx_start_i if is_start else tx_pos - exon.tx_end_i
        else:
            offset = exon.tx_end_i - tx_pos if is_start else exon.tx_start_i - tx_pos

        # Need to check if we need to change pos for liftover
        genomic_data, err_msg = await self.uta_db.get_alt_ac_start_or_end(
            tx_ac, tx_pos, tx_pos, gene
        )
        if genomic_data is None:
            return _GenomicTxSeg(errors=[err_msg])

        genomic_pos = (
            genomic_data.alt_start_i + offset
            if is_start
            else genomic_data.alt_end_i + offset
        )
        genomic_accession = genomic_data.alt_ac
        genomic_location_data = self._get_seq_loc(
            genomic_accession, genomic_pos, is_start, strand
        )
        if genomic_location_data.errors:
            return _GenomicTxSeg(errors=genomic_location_data.errors)

        return _GenomicTxSeg(
            gene=mane_gene,
            genomic_accession=genomic_accession,
            tx_ac=tx_ac,
            seg=TxSegment(
                exon_ord=exon_num,
                offset=offset,
                genomic_location=genomic_location_data.sequence_location,
            ),
        )

    async def _get_genomic_tx_seg(
        self,
        genomic_accession: str,
        genomic_pos: int,
        tx_ac: str,
        is_start: bool,
        gene: str,
    ) -> _GenomicTxSeg:
        """Set genomic data in ``params``

        :param params: Parameters for response
        :param is_start: ``True`` if ``pos`` is start position. ``False`` if ``pos`` is
            end position.
        :return: Warnings if found
        """
        # We should always try to liftover
        grch38_ac = await self.uta_db.get_newest_assembly_ac(genomic_accession)
        if not grch38_ac:
            return _GenomicTxSeg(
                errors=[f"Invalid genomic accession: {genomic_accession}"]
            )

        grch38_ac = grch38_ac[0]
        if grch38_ac != genomic_accession:  # alt_ac is genomic accession
            # Liftover to 38
            descr = await self.uta_db.get_chr_assembly(genomic_accession)
            if descr is None:
                return _GenomicTxSeg(
                    errors=[
                        f"Unable to get chromosome and assembly for "
                        f"{genomic_accession}"
                    ]
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
            genomic_accession = grch38_ac

        tx_exons = await self._get_exons_coords(tx_ac, genomic_accession=grch38_ac)
        if not tx_exons:
            return _GenomicTxSeg(errors=[f"Unable to get exons for {tx_ac}"])

        tx_exon_aln_data = await self.uta_db.get_tx_exon_aln_v_data(
            tx_ac,
            genomic_pos,
            genomic_pos,
            alt_ac=genomic_accession,
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
            tx_exon_aln_data.alt_start_i,
            tx_exon_aln_data.alt_end_i,
            genomic_pos,
            is_start=is_start,
        )

        genomic_location_data = self._get_seq_loc(
            genomic_accession, genomic_pos, is_start, tx_exon_aln_data.alt_strand
        )
        if genomic_location_data.errors:
            return _GenomicTxSeg(errors=genomic_location_data.errors)

        return _GenomicTxSeg(
            gene=tx_exon_aln_data.hgnc,
            genomic_accession=genomic_accession,
            tx_ac=tx_exon_aln_data.tx_ac,
            seg=TxSegment(
                exon_ord=tx_exon_aln_data.ord,
                offset=offset,
                genomic_location=genomic_location_data.sequence_location,
            ),
        )

    @staticmethod
    def _get_exon_offset(start: int, end: int, pos: int, is_start: bool) -> None:
        """Set value for ``exon_offset`` in ``params``.

        :param params: Parameters for response
        :param start: Start exon coord (can be transcript or aligned genomic)
        :param end: End exon coord (can be transcript or aligned genomic)
        :param pos: Position change (can be transcript or genomic)
        :param is_start: ``True`` if ``pos`` is start position. ``False`` if ``pos`` is
            end position
        """
        return pos - start if is_start else pos - end

    async def _get_exons_coords(
        self, tx_ac: str, genomic_accession: str | None = None
    ) -> list[ExonCoord]:
        """Get list of transcript exons start/end coordinates.

        :param tx_ac: Transcript accession
        :param genomic_accession: Genomic accession
        :return: List of exon coordinate data
        """
        if genomic_accession:
            # We know what assembly we're looking for since we have the
            # genomic accession
            query = f"""
                SELECT DISTINCT ord, tx_start_i, tx_end_i, alt_start_i, alt_end_i, alt_strand
                FROM {self.uta_db.schema}.tx_exon_aln_v
                WHERE tx_ac = '{tx_ac}'
                AND alt_aln_method = 'splign'
                AND alt_ac = '{genomic_accession}'
                """  # noqa: S608
        else:
            # Use GRCh38 by default if no genomic accession is provided
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
        result = await self.uta_db.execute_query(query)

        if not result:
            msg = f"Unable to get exons for {tx_ac}"
            _logger.warning(msg)
            return []
        return [ExonCoord(**r) for r in result]

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

        :param: tx_exons_genomic_coords: List of tuples describing exons and genomic
            coordinates for a transcript. Each tuple contains the transcript number
            (0-indexed), the transcript coordinates for the exon, and the genomic
            coordinates for the exon. Pos 0 in the tuple corresponds to the exon
            number, pos 1 and pos 2 refer to the start and end transcript coordinates,
            respectively, and pos 3 and 4 refer to the start and end genomic
            coordinates, respectively.
        :param strand: Strand
        :param: start: Genomic coordinate of breakpoint
        :param: end: Genomic coordinate of breakpoint
        :return: Exon number corresponding to adjacent exon. Will be 0-based.
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
        return exon.ord if end else exon.ord + 1

    @staticmethod
    def _is_exonic_breakpoint(pos: int, tx_genomic_coords: list[ExonCoord]) -> bool:
        """Check if a breakpoint occurs on an exon

        :param pos: Genomic breakpoint
        :param tx_genomic_coords: A list of genomic coordinates for a transcript
        :return: True if the breakpoint occurs on an exon
        """
        return any(
            pos >= exon.alt_start_i and pos <= exon.alt_end_i
            for exon in tx_genomic_coords
        )
