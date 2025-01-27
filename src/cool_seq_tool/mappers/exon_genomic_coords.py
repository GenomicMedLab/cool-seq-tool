"""Provide mapping capabilities between transcript exon and genomic coordinates."""

import logging

from ga4gh.vrs.models import SequenceLocation, SequenceReference
from pydantic import ConfigDict, Field, StrictInt, StrictStr, model_validator

from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.mappers.liftover import LiftOver
from cool_seq_tool.schemas import (
    Assembly,
    BaseModelForbidExtra,
    CoordinateType,
    ServiceMeta,
    Strand,
)
from cool_seq_tool.sources.mane_transcript_mappings import ManeTranscriptMappings
from cool_seq_tool.sources.uta_database import GenomicAlnData, UtaDatabase
from cool_seq_tool.utils import service_meta

_logger = logging.getLogger(__name__)


class _ExonCoord(BaseModelForbidExtra):
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

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "ord": 0,
                "tx_start_i": 0,
                "tx_end_i": 234,
                "alt_start_i": 154191901,
                "alt_end_i": 154192135,
                "alt_strand": Strand.NEGATIVE,
            }
        }
    )


class TxSegment(BaseModelForbidExtra):
    """Model for representing transcript segment data."""

    exon_ord: StrictInt = Field(..., description="Exon number. 0-based.")
    offset: StrictInt = Field(
        0,
        description="The value added to or subtracted from the `genomic_location` to find the start or end of an exon.",
    )
    genomic_location: SequenceLocation = Field(
        ..., description="The genomic position of a transcript segment."
    )

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "exon_ord": 0,
                "offset": 0,
                "genomic_location": {
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "type": "SequenceReference",
                        "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                    },
                    "end": 154192135,
                },
            }
        }
    )


class GenomicTxSeg(BaseModelForbidExtra):
    """Model for representing a boundary for a transcript segment."""

    seg: TxSegment | None = Field(None, description="Transcript segment.")
    gene: StrictStr | None = Field(
        None, description="Valid, case-sensitive HGNC gene symbol."
    )
    genomic_ac: StrictStr | None = Field(None, description="RefSeq genomic accession.")
    tx_ac: StrictStr | None = Field(None, description="RefSeq transcript accession.")
    errors: list[StrictStr] = Field([], description="Error messages.")

    @model_validator(mode="before")
    def check_errors(cls, values: dict) -> dict:  # noqa: N805
        """Ensure that fields are (un)set depending on errors

        :param values: Values in model
        :raises ValueError: If `seg`, `genomic_ac` and `tx_ac` are not
        provided when there are no errors
        :return: Values in model
        """
        if not values.get("errors") and not all(
            (
                values.get("seg"),
                values.get("genomic_ac"),
                values.get("tx_ac"),
            )
        ):
            err_msg = "`seg`, `genomic_ac` and `tx_ac` must be provided"
            raise ValueError(err_msg)
        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "gene": "TPM3",
                "genomic_ac": "NC_000001.11",
                "tx_ac": "NM_152263.3",
                "seg": {
                    "exon_ord": 0,
                    "offset": 0,
                    "genomic_location": {
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        },
                        "end": 154192135,
                    },
                },
                "errors": [],
            }
        }
    )


class GenomicTxSegService(BaseModelForbidExtra):
    """Service model for genomic and transcript data."""

    gene: StrictStr | None = Field(
        None, description="Valid, case-sensitive HGNC gene symbol."
    )
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
        :raises ValueError: If `genomic_ac`, `tx_ac` and `seg_start` or `seg_end`
            not provided when there are no errors
        :return: Values in model, including service metadata
        """
        values["service_meta"] = service_meta()
        if not values.get("errors") and not all(
            (
                values.get("genomic_ac"),
                values.get("tx_ac"),
                values.get("seg_start") or values.get("seg_end"),
            )
        ):
            err_msg = (
                "`genomic_ac`, `tx_ac` and `seg_start` or `seg_end` must be provided"
            )
            raise ValueError(err_msg)

        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "gene": "TPM3",
                "genomic_ac": "NC_000001.11",
                "tx_ac": "NM_152263.3",
                "seg_start": {
                    "exon_ord": 0,
                    "offset": 0,
                    "genomic_location": {
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        },
                        "end": 154192135,
                    },
                },
                "seg_end": {
                    "exon_ord": 7,
                    "offset": 0,
                    "genomic_location": {
                        "type": "SequenceLocation",
                        "sequenceReference": {
                            "type": "SequenceReference",
                            "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                        },
                        "start": 154170399,
                    },
                },
            }
        }
    )


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
        :param mane_transcript_mappings: Instance to provide access to ManeTranscriptMappings class
        :param liftover: Instance to provide mapping between human genome assemblies
        """
        self.seqrepo_access = seqrepo_access
        self.uta_db = uta_db
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
        """Get aligned genomic data given transcript segment data.

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
        >>> (
        ...     tpm3.genomic_ac,
        ...     tpm3.seg_start.genomic_location.end,
        ...     tpm3.seg_end.genomic_location.start,
        ... )
        ('NC_000001.11', 154192135, 154170399)

        :param transcript: RefSeq transcript accession
        :param gene: Valid, case-sensitive HGNC gene symbol
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

        # Get aligned genomic data (hgnc gene, alt_ac, alt_start_i, alt_end_i, strand)
        # for exon(s)
        (
            genomic_aln_start,
            genomic_aln_end,
            err_msg,
        ) = await self._get_genomic_aln_coords(
            transcript, tx_exon_start_coords, tx_exon_end_coords, gene=gene
        )
        if err_msg:
            return _return_service_errors([err_msg])

        # Get gene and chromosome data, check that at least one was retrieved
        gene = genomic_aln_start.hgnc if genomic_aln_start else genomic_aln_end.hgnc
        genomic_ac = (
            genomic_aln_start.alt_ac if genomic_aln_start else genomic_aln_end.alt_ac
        )
        if gene is None or genomic_ac is None:
            return _return_service_errors(
                [
                    "Unable to retrieve `gene` or `genomic_ac` from genomic start and genomic end data"
                ],
            )

        strand = (
            Strand(genomic_aln_start.alt_strand)
            if genomic_aln_start
            else Strand(genomic_aln_end.alt_strand)
        )

        if exon_start_exists:
            seg_start, err_msg = self._get_tx_segment(
                genomic_ac,
                strand,
                exon_start_offset,
                genomic_aln_start,
                is_seg_start=True,
            )
            if err_msg:
                return _return_service_errors([err_msg])
        else:
            seg_start = None

        if exon_end_exists:
            seg_end, err_msg = self._get_tx_segment(
                genomic_ac,
                strand,
                exon_end_offset,
                genomic_aln_end,
                is_seg_start=False,
            )
            if err_msg:
                return _return_service_errors([err_msg])
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
        seg_start_genomic: int | None = None,
        seg_end_genomic: int | None = None,
        transcript: str | None = None,
        gene: str | None = None,
        coordinate_type: CoordinateType = CoordinateType.INTER_RESIDUE,
        starting_assembly: Assembly = Assembly.GRCH38,
    ) -> GenomicTxSegService:
        """Get transcript segment data for genomic data, lifted over to GRCh38.

        If liftover to GRCh38 is unsuccessful, will return errors.

        MANE Transcript data will be returned if and only if ``transcript`` is not
        supplied. ``gene`` must be given in order to retrieve MANE Transcript data.

        >>> import asyncio
        >>> from cool_seq_tool import CoolSeqTool
        >>> from cool_seq_tool.schemas import Strand
        >>> egc = CoolSeqTool().ex_g_coords_mapper
        >>> result = asyncio.run(
        ...     egc.genomic_to_tx_segment(
        ...         genomic_ac="NC_000001.11",
        ...         seg_start_genomic=154192135,
        ...         seg_end_genomic=154170399,
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
            ``genomic_ac`` will be used. If the genomic accession is from GRCh37, it
            will be lifted over to GRCh38 and the original accession version will be
            ignored
        :param seg_start_genomic: Genomic position where the transcript segment starts
        :param seg_end_genomic: Genomic position where the transcript segment ends
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript. See the :ref:`Transcript Selection policy <transcript_selection_policy>`
            page.
        :param gene: A valid, case-sensitive HGNC symbol. Must be given if no ``transcript``
            value is provided.
        :param coordinate_type: Coordinate type for ``seg_start_genomic`` and
            ``seg_end_genomic``. Expects inter-residue coordinates by default
        :param starting_assembly: The assembly that the supplied coordinate comes from. Set to
            GRCh38 by default. Will attempt to liftover if starting assembly is GRCh37
        :return: Genomic data (inter-residue coordinates)
        """
        errors = []
        if seg_start_genomic is None and seg_end_genomic is None:
            errors.append(
                "Must provide either `seg_start_genomic` or `seg_end_genomic`"
            )
        if chromosome is None and genomic_ac is None:
            errors.append("Must provide either `chromosome` or `alt_ac`")
        if transcript is None and gene is None:
            errors.append("Must provide either `gene` or `transcript`")
        if errors:
            return _return_service_errors(errors)

        params = {}

        if seg_start_genomic:
            start_tx_seg_data = await self._genomic_to_tx_segment(
                seg_start_genomic,
                chromosome=chromosome,
                genomic_ac=genomic_ac,
                transcript=transcript,
                gene=gene,
                is_seg_start=True,
                coordinate_type=coordinate_type,
                starting_assembly=starting_assembly,
            )
            if start_tx_seg_data.errors:
                return _return_service_errors(start_tx_seg_data.errors)

            params["gene"] = start_tx_seg_data.gene
            params["genomic_ac"] = start_tx_seg_data.genomic_ac
            params["tx_ac"] = start_tx_seg_data.tx_ac
            params["seg_start"] = start_tx_seg_data.seg
        else:
            start_tx_seg_data = None

        if seg_end_genomic:
            end_tx_seg_data = await self._genomic_to_tx_segment(
                seg_end_genomic,
                chromosome=chromosome,
                genomic_ac=genomic_ac,
                transcript=transcript,
                gene=gene,
                is_seg_start=False,
                coordinate_type=coordinate_type,
                starting_assembly=starting_assembly,
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

    async def _get_start_end_exon_coords(
        self,
        tx_ac: str,
        exon_start: int | None = None,
        exon_end: int | None = None,
        genomic_ac: str | None = None,
    ) -> tuple[_ExonCoord | None, _ExonCoord | None, list[str]]:
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
            return None, None, [f"Transcript does not exist in UTA: {tx_ac}"]

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

    async def _get_all_exon_coords(
        self, tx_ac: str, genomic_ac: str | None = None
    ) -> list[_ExonCoord]:
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
        return [_ExonCoord(**r) for r in results]

    async def _get_genomic_aln_coords(
        self,
        tx_ac: str,
        tx_exon_start: _ExonCoord | None = None,
        tx_exon_end: _ExonCoord | None = None,
        gene: str | None = None,
    ) -> tuple[GenomicAlnData | None, GenomicAlnData | None, str | None]:
        """Get aligned genomic coordinates for transcript exon start and end.

        ``tx_exon_start`` and ``tx_exon_end`` is expected to reference the same
        transcript and genomic accession.

        :param tx_ac: Transcript accession
        :param tx_exon_start: Transcript's exon start coordinates. If not provided,
            must provide ``tx_exon_end``
        :param tx_exon_end: Transcript's exon end coordinates. If not provided, must
            provide ``tx_exon_start``
        :param gene: A valid, case-sensitive HGNC gene symbol
        :return: Tuple containing aligned genomic data for start and end exon and
            warnings if found
        """
        if tx_exon_start is None and tx_exon_end is None:
            msg = "Must provide either `tx_exon_start` or `tx_exon_end` or both"
            _logger.warning(msg)
            return None, None, msg

        aligned_coords = {"start": None, "end": None}
        for exon, key in [(tx_exon_start, "start"), (tx_exon_end, "end")]:
            if exon:
                aligned_coord, warning = await self.uta_db.get_alt_ac_start_or_end(
                    tx_ac, exon.tx_start_i, exon.tx_end_i, gene=gene
                )
                if aligned_coord:
                    aligned_coords[key] = aligned_coord
                else:
                    return None, None, warning

        return *aligned_coords.values(), None

    def _get_tx_segment(
        self,
        genomic_ac: str,
        strand: Strand,
        offset: int,
        genomic_ac_data: _ExonCoord,
        is_seg_start: bool = False,
    ) -> tuple[TxSegment | None, str | None]:
        """Get transcript segment data given ``genomic_ac`` and offset data

        :param genomic_ac: Genomic RefSeq accession
        :param strand: Strand
        :param offset: Exon offset
        :param genomic_ac_data: Exon coordinate data for ``genomic_ac``
        :param is_seg_start: ``True`` if retrieving genomic data where the transcript
            segment starts, defaults to ``False``
        :return: Transcript segment data
        """
        if is_seg_start:
            if strand == Strand.POSITIVE:
                seg_genomic_pos = offset + genomic_ac_data.alt_start_i
            else:
                seg_genomic_pos = genomic_ac_data.alt_end_i - offset
        else:
            if strand == Strand.POSITIVE:
                seg_genomic_pos = offset + genomic_ac_data.alt_end_i
            else:
                seg_genomic_pos = genomic_ac_data.alt_start_i - offset

        genomic_loc, err_msg = self._get_vrs_seq_loc(
            genomic_ac,
            seg_genomic_pos,
            is_seg_start=is_seg_start,
            strand=strand,
        )
        if err_msg:
            return None, err_msg

        return TxSegment(
            exon_ord=genomic_ac_data.ord,
            genomic_location=genomic_loc,
            offset=offset,
        ), None

    def _get_vrs_seq_loc(
        self, genomic_ac: str, genomic_pos: int, is_seg_start: bool, strand: Strand
    ) -> tuple[SequenceLocation | None, str | None]:
        """Create VRS Sequence Location for genomic position where transcript segment
        occurs

        :param genomic_ac: RefSeq genomic accession
        :param genomic_pos: Genomic position where the transcript segment occurs
        :param is_seg_start: ``True`` if ``genomic_pos`` is where the transcript segment
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

        use_start = (
            strand == Strand.POSITIVE if is_seg_start else strand != Strand.POSITIVE
        )

        return SequenceLocation(
            sequenceReference=SequenceReference(
                refgetAccession=ga4gh_seq_id[0].split("ga4gh:")[-1]
            ),
            start=genomic_pos if use_start else None,
            end=genomic_pos if not use_start else None,
        ), None

    async def _genomic_to_tx_segment(
        self,
        genomic_pos: int,
        chromosome: str | None = None,
        genomic_ac: str | None = None,
        transcript: str | None = None,
        gene: str | None = None,
        is_seg_start: bool = True,
        coordinate_type: CoordinateType = CoordinateType.INTER_RESIDUE,
        starting_assembly: Assembly = Assembly.GRCH38,
    ) -> GenomicTxSeg:
        """Given genomic data, generate a boundary for a transcript segment.

        Will liftover to GRCh38 assembly. If liftover is unsuccessful, will return
        errors.

        Either an HGNC gene symbol or transcript accession must be provided to this
        method

        :param genomic_pos: Genomic position where the transcript segment starts or ends
        :param chromosome: Chromosome. Must give chromosome without a prefix
            (i.e. ``1`` or ``X``). If not provided, must provide ``genomic_ac``. If
            position maps to both GRCh37 and GRCh38, GRCh38 assembly will be used.
            If ``genomic_ac`` is also provided, ``genomic_ac`` will be used.
        :param genomic_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``genomic_ac``
            will be used. If the genomic accession is from GRCh37, it will be lifted
            over to GRCh38 and the original accession version will be ignored
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript
        :param gene: Valid, case-sensitive HGNC gene symbol
        :param is_seg_start: ``True`` if ``genomic_pos`` is where the transcript segment starts.
            ``False`` if ``genomic_pos`` is where the transcript segment ends.
        :param coordinate_type: Coordinate type for ``seg_start_genomic`` and
            ``seg_end_genomic``. Expects inter-residue coordinates by default
        :param starting_assembly: The assembly that the supplied coordinate comes from. Set to
            GRCh38 by default. Will attempt to liftover if starting assembly is GRCh37
        :return: Data for a transcript segment boundary (inter-residue coordinates)
        """
        params = {key: None for key in GenomicTxSeg.model_fields}

        # Validate inputs exist in UTA
        if gene:
            gene_validation = await self.uta_db.gene_exists(gene)
            if not gene_validation:
                return GenomicTxSeg(errors=[f"Gene does not exist in UTA: {gene}"])

        if transcript:
            transcript_validation = await self.uta_db.transcript_exists(transcript)
            if not transcript_validation:
                return GenomicTxSeg(
                    errors=[f"Transcript does not exist in UTA: {transcript}"]
                )

        if genomic_ac:
            grch38_ac = await self.uta_db.get_newest_assembly_ac(genomic_ac)
            if grch38_ac:
                genomic_ac = grch38_ac[0]
            else:
                return GenomicTxSeg(
                    errors=[f"Genomic accession does not exist in UTA: {genomic_ac}"]
                )
        else:
            genomic_acs, err_msg = self.seqrepo_access.chromosome_to_acs(chromosome)

            if not genomic_acs:
                return GenomicTxSeg(
                    errors=[err_msg],
                )
            genomic_ac = genomic_acs[0]

        # Liftover to GRCh38 if the provided assembly is GRCh37
        if starting_assembly == Assembly.GRCH37:
            genomic_pos = await self._get_grch38_pos(
                genomic_ac, genomic_pos, chromosome=chromosome if chromosome else None
            )
            if not genomic_pos:
                return GenomicTxSeg(
                    errors=[
                        f"Lifting over {genomic_pos} on {genomic_ac} from {Assembly.GRCH37.value} to {Assembly.GRCH38.value} was unsuccessful."
                    ]
                )

        # Select a transcript if not provided
        if not transcript:
            mane_transcripts = self.mane_transcript_mappings.get_gene_mane_data(gene)

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
                        return GenomicTxSeg(
                            errors=[
                                f"Could not find a transcript for {gene} on {genomic_ac}"
                            ]
                        )
        # gene is not required to liftover coordinates if tx_ac and genomic_ac are given, but we should set the associated gene
        if not gene:
            _gene, err_msg = await self._get_tx_ac_gene(transcript)
            if err_msg:
                return GenomicTxSeg(errors=[err_msg])
            gene = _gene

        tx_exons = await self._get_all_exon_coords(
            tx_ac=transcript, genomic_ac=genomic_ac
        )
        if not tx_exons:
            return GenomicTxSeg(
                errors=[f"No exons found given transcript: {transcript}"]
            )

        strand = Strand(tx_exons[0].alt_strand)
        params["strand"] = strand
        use_alt_start_i = self._use_alt_start_i(
            is_seg_start=is_seg_start, strand=strand
        )
        if use_alt_start_i and coordinate_type == CoordinateType.RESIDUE:
            genomic_pos = genomic_pos - 1  # Convert residue coordinate to inter-residue

        # Validate that the breakpoint between the first and last exon for the selected transcript
        coordinate_check = await self._validate_genomic_breakpoint(
            pos=genomic_pos, genomic_ac=genomic_ac, tx_ac=transcript
        )
        if not coordinate_check:
            return GenomicTxSeg(
                errors=[
                    f"{genomic_pos} on {genomic_ac} does not occur within the exons for {transcript}"
                ]
            )

        # Check if breakpoint occurs on an exon.
        # If not, determine the adjacent exon given the selected transcript
        if not self._is_exonic_breakpoint(genomic_pos, tx_exons):
            exon_num = self._get_adjacent_exon(
                tx_exons_genomic_coords=tx_exons,
                strand=strand,
                start=genomic_pos if is_seg_start else None,
                end=genomic_pos if not is_seg_start else None,
            )
        else:
            exon_data = await self.uta_db.get_tx_exon_aln_v_data(
                transcript,
                genomic_pos,
                genomic_pos,
                alt_ac=genomic_ac,
                use_tx_pos=False,
            )
            exon_num = exon_data[0].ord

        offset = self._get_exon_offset(
            genomic_pos=genomic_pos,
            exon_boundary=tx_exons[exon_num].alt_start_i
            if use_alt_start_i
            else tx_exons[exon_num].alt_end_i,
            strand=strand,
        )

        genomic_location, err_msg = self._get_vrs_seq_loc(
            genomic_ac, genomic_pos, is_seg_start, strand
        )
        if err_msg:
            return GenomicTxSeg(errors=[err_msg])

        return GenomicTxSeg(
            gene=gene,
            genomic_ac=genomic_ac,
            tx_ac=transcript,
            seg=TxSegment(
                exon_ord=exon_num,
                offset=offset,
                genomic_location=genomic_location,
            ),
        )

    async def _get_grch38_pos(
        self,
        genomic_ac: str,
        genomic_pos: int,
        chromosome: str | None = None,
    ) -> int | None:
        """Get GRCh38 genomic representation for accession and position

        :param genomic_pos: A genomic coordinate in GRCh37
        :param genomic_ac: The genomic accession in GRCh38
        :param chromosome: The chromosome that genomic_pos occurs on
        :return The genomic coordinate in GRCh38
        """
        if not chromosome:
            chromosome, _ = self.seqrepo_access.translate_identifier(
                genomic_ac, target_namespaces=Assembly.GRCH38.value
            )
            chromosome = chromosome[-1].split(":")[-1]
        liftover_data = self.liftover.get_liftover(
            chromosome, genomic_pos, Assembly.GRCH38
        )
        return liftover_data[1] if liftover_data else None

    async def _validate_genomic_breakpoint(
        self,
        pos: int,
        genomic_ac: str,
        tx_ac: str,
    ) -> bool:
        """Validate that a genomic coordinate falls within the first and last exon
            for a transcript on a given accession

        :param pos: Genomic position on ``genomic_ac``
        :param genomic_ac: RefSeq genomic accession, e.g. ``"NC_000007.14"``
        :param transcript: A transcript accession
        :return: ``True`` if the coordinate falls within the first and last exon
            for the transcript, ``False`` if not
        """
        query = f"""
            WITH tx_boundaries AS (
                SELECT
                MIN(alt_start_i) AS min_start,
                MAX(alt_end_i) AS max_end
                FROM {self.uta_db.schema}.tx_exon_aln_v
                WHERE tx_ac = '{tx_ac}'
                AND alt_ac = '{genomic_ac}'
            )
            SELECT * FROM tx_boundaries
            WHERE {pos} between tx_boundaries.min_start and tx_boundaries.max_end
            """  # noqa: S608
        results = await self.uta_db.execute_query(query)
        return bool(results)

    async def _get_tx_ac_gene(
        self,
        tx_ac: str,
    ) -> tuple[str | None, str | None]:
        """Get gene given a transcript.

        If multiple genes are found for a given ``tx_ac``, only one
        gene will be returned.

        :param tx_ac: RefSeq transcript, e.g. ``"NM_004333.6"``
        :return: HGNC gene symbol associated to transcript and
            warning
        """
        query = f"""
            SELECT DISTINCT hgnc
            FROM {self.uta_db.schema}.tx_exon_aln_v
            WHERE tx_ac = '{tx_ac}'
            ORDER BY hgnc
            LIMIT 1;
            """  # noqa: S608
        results = await self.uta_db.execute_query(query)
        if not results:
            return None, f"No gene(s) found given {tx_ac}"

        return results[0]["hgnc"], None

    @staticmethod
    def _is_exonic_breakpoint(pos: int, tx_genomic_coords: list[_ExonCoord]) -> bool:
        """Check if a breakpoint occurs on an exon

        :param pos: Genomic breakpoint
        :param tx_genomic_coords: A list of transcript exon coordinate data
        :return: ``True`` if the breakpoint occurs on an exon
        """
        return any(
            exon.alt_start_i <= pos <= exon.alt_end_i for exon in tx_genomic_coords
        )

    @staticmethod
    def _use_alt_start_i(is_seg_start: bool, strand: Strand) -> bool:
        """Determine whether to use alt_start_i or alt_end_i from UTA when computing
        exon offset

        :param is_seg_start: ``True`` if ``genomic_pos`` is where the transcript segment starts.
            ``False`` if ``genomic_pos`` is where the transcript segment ends.
        :param strand: The transcribed strand
        :return ``True`` if alt_start_i should be used, ``False`` if alt_end_i should
        be used
        """
        return (
            is_seg_start
            and strand == Strand.POSITIVE
            or not is_seg_start
            and strand == Strand.NEGATIVE
        )

    @staticmethod
    def _get_adjacent_exon(
        tx_exons_genomic_coords: list[_ExonCoord],
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
        # If a transcript has only one exon, return 0
        if len(tx_exons_genomic_coords) == 1:
            return 0

        # Check if a breakpoint occurs before/after the transcript boundaries
        bp = start if start else end
        exon_list_len = len(tx_exons_genomic_coords) - 1

        if strand == Strand.POSITIVE:
            if bp < tx_exons_genomic_coords[0].alt_start_i:
                return 0
            if bp > tx_exons_genomic_coords[exon_list_len].alt_end_i:
                return exon_list_len
        if strand == Strand.NEGATIVE:
            if bp > tx_exons_genomic_coords[0].alt_end_i:
                return 0
            if bp < tx_exons_genomic_coords[exon_list_len].alt_start_i:
                return exon_list_len

        for i in range(exon_list_len):
            exon = tx_exons_genomic_coords[i]
            if start == exon.alt_start_i:
                break
            if end == exon.alt_end_i:
                break
            next_exon = tx_exons_genomic_coords[i + 1]
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
    def _get_exon_offset(
        genomic_pos: int,
        exon_boundary: int,
        strand: Strand,
    ) -> int:
        """Compute offset from exon start or end index

        :param genomic_pos: The supplied genomic position. This can represent, for
            example, a fusion junction breakpoint. This position is represented using
            inter-residue coordinates
        :param exon_boundary: The genomic position for the exon boundary that the offset
            is being computed against
        :paran strand: The transcribed strand
        :return: Offset from exon start or end index
        """
        return (
            genomic_pos - exon_boundary
            if strand == Strand.POSITIVE
            else (genomic_pos - exon_boundary) * -1
        )
