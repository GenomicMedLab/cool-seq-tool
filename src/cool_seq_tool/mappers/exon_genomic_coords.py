"""Provide mapping capabilities between transcript exon and genomic coordinates."""
import logging
from typing import Dict, List, Optional, Tuple, TypeVar, Union

from cool_seq_tool.mappers.mane_transcript import CdnaRepresentation, ManeTranscript
from cool_seq_tool.schemas import (
    AnnotationLayer,
    Assembly,
    GenomicData,
    GenomicDataResponse,
    ResidueMode,
    Strand,
    TranscriptExonData,
    TranscriptExonDataResponse,
)
from cool_seq_tool.sources.uta_database import UtaDatabase
from cool_seq_tool.utils import get_inter_residue_pos, service_meta

CoordinatesResponseType = TypeVar(
    "CoordinatesResponseType", GenomicDataResponse, TranscriptExonDataResponse
)

logger = logging.getLogger(__name__)


class ExonGenomicCoordsMapper:
    """Provide capabilities for mapping transcript exon representation to/from genomic
    coordinate representation.
    """

    def __init__(self, uta_db: UtaDatabase, mane_transcript: ManeTranscript) -> None:
        """Initialize ExonGenomicCoordsMapper class.

        A lot of resources are required for initialization, so when defaults are enough,
        it's easiest to let the core CoolSeqTool class handle it for you:

        >>> from cool_seq_tool.app import CoolSeqTool
        >>> egc = CoolSeqTool().ex_g_coords_mapper

        Note that this class's public methods are all defined as ``async``, so they will
        need to be called with ``await`` when called from a function, or run from an
        event loop. See the :ref:`Usage section <async_note>` for more information.

        >>> import asyncio
        >>> result = asyncio.run(egc.transcript_to_genomic_coordinates(
        ...     "NM_002529.3",
        ...     exon_start=2,
        ...     exon_end=17
        ... ))
        >>> result.genomic_data.start, result.genomic_data.end
        (156864428, 156881456)

        :param uta_db: UtaDatabase instance to give access to query UTA database
        :param mane_transcript: Instance to align to MANE or compatible representation
        """
        self.uta_db = uta_db
        self.mane_transcript = mane_transcript

    @staticmethod
    def _return_warnings(
        resp: CoordinatesResponseType, warning_msg: str
    ) -> CoordinatesResponseType:
        """Add warnings to response object

        :param resp: Response object
        :param warning_msg: Warning message on why ``transcript_exon_data`` or
            ``genomic_data`` field is ``None``
        :return: Response object with warning message
        """
        logger.warning(warning_msg)
        resp.warnings.append(warning_msg)
        return resp

    async def transcript_to_genomic_coordinates(
        self,
        transcript: str,
        gene: Optional[str] = None,
        exon_start: Optional[int] = None,
        exon_start_offset: int = 0,
        exon_end: Optional[int] = None,
        exon_end_offset: int = 0,
    ) -> GenomicDataResponse:
        """Get genomic data given transcript data.

        By default, transcript data is aligned to the GRCh38 assembly.

        >>> import asyncio
        >>> from cool_seq_tool.app import CoolSeqTool
        >>> egc = CoolSeqTool().ex_g_coords_mapper
        >>> tpm3 = asyncio.run(egc.transcript_to_genomic_coordinates(
        ...     "NM_152263.3"
        ...     gene="TPM3", chr="NC_000001.11",
        ...     exon_start=1, exon_end=8,
        ... ))
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
        if not transcript:
            return self._return_warnings(resp, "Must provide `transcript`")
        else:
            transcript = transcript.strip()

        exon_start_exists, exon_end_exists = False, False
        if exon_start is not None:
            if exon_start < 1:
                return self._return_warnings(resp, "`exon_start` cannot be less than 1")
            exon_start_exists = True

        if exon_end is not None:
            if exon_end < 1:
                return self._return_warnings(resp, "`exon_end` cannot be less than 1")
            exon_end_exists = True

        if not exon_start_exists and not exon_end_exists:
            return self._return_warnings(
                resp, "Must provide either `exon_start` or `exon_end`"
            )
        elif exon_start_exists and exon_end_exists:
            if exon_start > exon_end:
                return self._return_warnings(
                    resp,
                    f"Start exon {exon_start} is greater than end exon {exon_end}",
                )

        # Get all exons and associated start/end coordinates for transcript
        tx_exons, warning = await self.uta_db.get_tx_exons(transcript)
        if not tx_exons:
            return self._return_warnings(resp, warning or "")

        # Get exon start and exon end coordinates
        tx_exon_coords, warning = self.get_tx_exon_coords(
            transcript, tx_exons, exon_start, exon_end
        )
        if not tx_exon_coords:
            return self._return_warnings(resp, warning or "")
        tx_exon_start_coords, tx_exon_end_coords = tx_exon_coords

        if gene:
            gene = gene.upper().strip()

        # Get aligned genomic data (hgnc gene, alt_ac, alt_start_i, alt_end_i, strand)
        # for exon(s)
        alt_ac_start_end, warning = await self._get_alt_ac_start_and_end(
            transcript, tx_exon_start_coords, tx_exon_end_coords, gene=gene
        )
        if not alt_ac_start_end:
            return self._return_warnings(resp, warning or "")
        alt_ac_start_data, alt_ac_end_data = alt_ac_start_end

        # Get gene and chromosome data, check that at least one was retrieved
        gene = alt_ac_start_data[0] if alt_ac_start_data else alt_ac_end_data[0]
        chromosome = alt_ac_start_data[1] if alt_ac_start_data else alt_ac_end_data[1]
        if gene is None or chromosome is None:
            return self._return_warnings(
                resp,
                "Unable to retrieve `gene` or `chromosome` from genomic start and "
                "genomic end data",
            )

        g_start = alt_ac_start_data[3] - 1 if alt_ac_start_data else None
        g_end = alt_ac_end_data[2] + 1 if alt_ac_end_data else None
        strand = (
            Strand(alt_ac_start_data[4])
            if alt_ac_start_data
            else Strand(alt_ac_end_data[4])
        )

        # Using none since could set to 0
        start_exits = g_start is not None
        end_exists = g_end is not None

        # Calculate offsets
        if strand == Strand.NEGATIVE:
            start_offset = exon_start_offset * -1 if start_exits else None
            end_offset = exon_end_offset * -1 if end_exists else None
        else:
            start_offset = exon_start_offset if start_exits else 0
            end_offset = exon_end_offset if end_exists else None

        # Get genomic coordinates with offsets included
        g_start = g_start + start_offset if start_exits else None
        g_end = g_end + end_offset if end_exists else None

        resp.genomic_data = GenomicData(
            gene=gene,
            chr=chromosome,
            start=g_start,
            end=g_end,
            exon_start=exon_start if start_exits else None,
            exon_start_offset=exon_start_offset if start_exits else 0,
            exon_end=exon_end if end_exists else None,
            exon_end_offset=exon_end_offset if end_exists else 0,
            transcript=transcript,
            strand=strand,
        )

        return resp

    async def genomic_to_transcript_exon_coordinates(
        self,
        chromosome: Optional[str] = None,
        alt_ac: Optional[str] = None,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[Strand] = None,
        transcript: Optional[str] = None,
        gene: Optional[str] = None,
        residue_mode: Union[
            ResidueMode.INTER_RESIDUE, ResidueMode.RESIDUE
        ] = ResidueMode.RESIDUE,
    ) -> GenomicDataResponse:
        """Get transcript data for genomic data, lifted over to GRCh38.

        MANE Transcript data will be returned if and only if ``transcript`` is not
        supplied. ``gene`` must be given in order to retrieve MANE Transcript data.

        >>> import asyncio
        >>> from cool_seq_tool.app import CoolSeqTool
        >>> from cool_seq_tool.schemas import Strand
        >>> egc = CoolSeqTool().ex_g_coords_mapper
        >>> result = asyncio.run(egc.genomic_to_transcript_exon_coordinates(
        ...     chromosome="NC_000001.11",
        ...     start=154192136,
        ...     end=154170400,
        ...     strand=Strand.NEGATIVE,
        ...     transcript="NM_152263.3"
        ... ))
        >>> result.genomic_data.exon_start, result.genomic_data.exon_end
        (1, 8)

        :param chromosome: Chromosome. Must give chromosome without a prefix
            (i.e. ``1`` or ``X``). If not provided, must provide ``alt_ac``.
            If ``alt_ac`` is also provided, ``alt_ac`` will be used.
        :param alt_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``alt_ac``
            will be used.
        :param start: Start genomic position
        :param end: End genomic position
        :param strand: Strand
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript. See the :ref:`Transcript Selection policy <transcript_selection_policy>`
            page.
        :param gene: HGNC gene symbol
        :param residue_mode: Residue mode for ``start`` and ``end``
        :return: Genomic data (inter-residue coordinates)
        """
        resp = GenomicDataResponse(
            genomic_data=None, warnings=[], service_meta=service_meta()
        )
        if start is None and end is None:
            return self._return_warnings(resp, "Must provide either `start` or `end`")

        params = {key: None for key in GenomicData.model_fields.keys()}
        if gene is not None:
            gene = gene.upper().strip()

        if start:
            if residue_mode == ResidueMode.RESIDUE:
                # zero-based for UTA
                start -= 1
            residue_mode = ResidueMode.ZERO
            start_data = await self._genomic_to_transcript_exon_coordinate(
                start,
                chromosome=chromosome,
                alt_ac=alt_ac,
                strand=strand,
                transcript=transcript,
                gene=gene,
                is_start=True,
            )
            if start_data.transcript_exon_data:
                start_data = start_data.transcript_exon_data.model_dump()
            else:
                return self._return_warnings(resp, start_data.warnings[0])
        else:
            start_data = None

        if end:
            end -= 1
            residue_mode = ResidueMode.ZERO
            end_data = await self._genomic_to_transcript_exon_coordinate(
                end,
                chromosome=chromosome,
                alt_ac=alt_ac,
                strand=strand,
                transcript=transcript,
                gene=gene,
                is_start=False,
            )
            if end_data.transcript_exon_data:
                end_data = end_data.transcript_exon_data.model_dump()
            else:
                return self._return_warnings(resp, end_data.warnings[0])
        else:
            end_data = None

        for field in ["transcript", "gene", "chr", "strand"]:
            if start_data:
                if end_data:
                    if start_data[field] != end_data[field]:
                        msg = (
                            f"Start `{field}`, {start_data[field]}, does "
                            f"not match End `{field}`, {end_data[field]}"
                        )
                        return self._return_warnings(resp, msg)
                params[field] = start_data[field]
            else:
                params[field] = end_data[field]

        if gene and gene != params["gene"]:
            msg = (
                f"Input gene, {gene}, does not match expected output"
                f"gene, {params['gene']}"
            )
            return self._return_warnings(resp, msg)

        for label, data in [("start", start_data), ("end", end_data)]:
            if data:
                params[label] = data["pos"]
                params[f"exon_{label}"] = data["exon"]
                params[f"exon_{label}_offset"] = data["exon_offset"]
        resp.genomic_data = GenomicData(**params)
        return resp

    @staticmethod
    def _validate_exon(
        transcript: str, tx_exons: List[Tuple[int, int]], exon_number: int
    ) -> Tuple[Optional[Tuple[int, int]], Optional[str]]:
        """Validate that exon number exists on a given transcript

        :param transcript: Transcript accession
        :param tx_exons: List of transcript's exons and associated coordinates
        :param exon_number: Exon number to validate
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

    def get_tx_exon_coords(
        self,
        transcript: str,
        tx_exons: List[Tuple[int, int]],
        exon_start: Optional[int] = None,
        exon_end: Optional[int] = None,
    ) -> Tuple[
        Optional[Tuple[Optional[Tuple[int, int]], Optional[Tuple[int, int]]]],
        Optional[str],
    ]:
        """Get exon coordinates for ``exon_start`` and ``exon_end``

        :param transcript: Transcript accession
        :param tx_exons: List of all transcript exons and coordinates
        :param exon_start: Start exon number
        :param exon_end: End exon number
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
        tx_exon_start: Optional[Tuple[int, int]] = None,
        tx_exon_end: Optional[Tuple[int, int]] = None,
        gene: Optional[str] = None,
    ) -> Tuple[Optional[Tuple[Tuple, Tuple]], Optional[str]]:
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
            logger.warning(msg)
            return None, msg

        alt_ac_data = {"start": None, "end": None}
        for exon, key in [(tx_exon_start, "start"), (tx_exon_end, "end")]:
            if exon:
                alt_ac_val, warning = await self.uta_db.get_alt_ac_start_or_end(
                    tx_ac, exon[0], exon[1], gene=gene
                )
                if alt_ac_val:
                    alt_ac_data[key] = alt_ac_val
                else:
                    return None, warning

        alt_ac_data_values = alt_ac_data.values()
        # Validate that start and end alignments have matching gene, genomic accession,
        # and strand
        if all(alt_ac_data_values):
            for i in (0, 1, 4):
                if alt_ac_data["start"][i] != alt_ac_data["end"][i]:
                    if i == 0:
                        error = "HGNC gene symbol does not match"
                    elif i == 1:
                        error = "Genomic accession does not match"
                    else:
                        error = "Strand does not match"
                    logger.warning(
                        f"{error}: {alt_ac_data['start'][i]} != {alt_ac_data['end'][i]}"
                    )
                    return None, error
        return tuple(alt_ac_data_values), None

    async def _genomic_to_transcript_exon_coordinate(
        self,
        pos: int,
        chromosome: Optional[str] = None,
        alt_ac: Optional[str] = None,
        strand: Optional[Strand] = None,
        transcript: Optional[str] = None,
        gene: Optional[str] = None,
        is_start: bool = True,
    ) -> TranscriptExonDataResponse:
        """Convert individual genomic data to transcript data

        :param pos: Genomic position (zero-based)
        :param chromosome: Chromosome. Must give chromosome without a prefix
            (i.e. ``1`` or ``X``). If not provided, must provide ``alt_ac``.
            If ``alt_ac`` is also provided, ``alt_ac`` will be used.
        :param alt_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``alt_ac``
            will be used.
        :param strand: Strand
        :param transcript: The transcript to use. If this is not given, we will try the
            following transcripts: MANE Select, MANE Clinical Plus, Longest Remaining
            Compatible Transcript
        :param gene: HGNC gene symbol
        :param is_start: ``True`` if ``pos`` is start position. ``False`` if ``pos`` is
            end position.
        :return: Transcript data (inter-residue coordinates)
        """
        resp = TranscriptExonDataResponse(
            transcript_exon_data=None, warnings=[], service_meta=service_meta()
        )

        if transcript is None and gene is None:
            return self._return_warnings(
                resp, "Must provide either `gene` or `transcript`"
            )

        params = {key: None for key in TranscriptExonData.model_fields.keys()}

        if alt_ac:
            # Check if valid accession is given
            if not await self.uta_db.validate_genomic_ac(alt_ac):
                return self._return_warnings(
                    resp, f"Invalid genomic accession: {alt_ac}"
                )

            genes_alt_acs, warning = await self.uta_db.get_genes_and_alt_acs(
                pos, strand=strand, alt_ac=alt_ac, gene=gene
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
                pos, strand=strand, chromosome=chromosome, gene=gene
            )
        else:
            genes_alt_acs = None

        if not genes_alt_acs:
            return self._return_warnings(resp, warning)

        gene_alt_ac, warning = self._get_gene_and_alt_ac(genes_alt_acs, gene)
        if not gene_alt_ac:
            return self._return_warnings(resp, warning)
        gene, alt_ac = gene_alt_ac

        if transcript is None:
            warnings = await self._set_mane_genomic_data(
                params, gene, alt_ac, pos, strand, is_start
            )
            if warnings:
                return self._return_warnings(resp, warnings)
        else:
            params["transcript"] = transcript
            params["gene"] = gene
            params["pos"] = pos
            params["chr"] = alt_ac
            warning = await self._set_genomic_data(params, strand, is_start)
            if warning:
                return self._return_warnings(resp, warning)

        resp.transcript_exon_data = TranscriptExonData(**params)
        return resp

    @staticmethod
    def _get_gene_and_alt_ac(
        genes_alt_acs: Dict, gene: Optional[str]
    ) -> Tuple[Optional[Tuple[str, str]], Optional[str]]:
        """Return gene genomic accession

        :param genes_alt_acs: Dictionary containing genes and genomic accessions
        :param gene: Gene symbol
        :return: (Gene, Genomic accession) if both exist
        """
        alt_acs = genes_alt_acs["alt_acs"]
        len_alt_acs = len(alt_acs)
        if len_alt_acs > 1:
            return None, f"Found more than one accessions: {alt_acs}"
        elif len_alt_acs == 0:
            return None, "No genomic accessions found"
        alt_ac = next(iter(alt_acs))

        genes = genes_alt_acs["genes"]
        len_genes = len(genes)
        input_gene = gene
        output_gene = None
        if len_genes == 1:
            output_gene = next(iter(genes))
        elif len_genes > 1:
            return None, f"Found more than one gene: {genes}"
        elif len_genes == 0:
            return None, "No genes found"

        if input_gene is not None:
            if output_gene != input_gene.upper():
                return (
                    None,
                    f"Input gene, {input_gene}, does not match "
                    f"expected output gene, {output_gene}",
                )

        gene = output_gene if output_gene else input_gene
        return (gene, alt_ac), None

    async def _set_mane_genomic_data(
        self,
        params: Dict,
        gene: str,
        alt_ac: str,
        pos: int,
        strand: Strand,
        is_start: bool,
    ) -> Optional[str]:
        """Set genomic data in `params` found from MANE.

        :param params: Parameters for response
        :param gene: Gene symbol
        :param alt_ac: Genomic accession
        :param pos: Genomic position
        :param strand: Strand
        :param is_start: `True` if `pos` is start position. `False` if `pos` is end
            position.
        :return: Warnings if found
        """
        start, end = get_inter_residue_pos(pos, pos, residue_mode=ResidueMode.ZERO)
        mane_data: Optional[
            CdnaRepresentation
        ] = await self.mane_transcript.get_mane_transcript(
            alt_ac,
            start,
            end,
            AnnotationLayer.GENOMIC,
            gene=gene,
            try_longest_compatible=True,
            residue_mode=ResidueMode.INTER_RESIDUE,
        )
        if not mane_data:
            msg = f"Unable to find mane data for {alt_ac} with position {pos}"
            if gene:
                msg += f" on gene {gene}"
            logger.warning(msg)
            return msg

        params["gene"] = mane_data.gene
        params["transcript"] = (
            mane_data.refseq
            if mane_data.refseq
            else mane_data.ensembl
            if mane_data.ensembl
            else None
        )
        tx_exons = await self._structure_exons(params["transcript"], alt_ac=alt_ac)
        if not tx_exons:
            return f"Unable to get exons for {params['transcript']}"
        tx_pos = mane_data.pos[0] + mane_data.coding_start_site
        params["exon"] = self._get_exon_number(tx_exons, tx_pos)

        try:
            tx_exon = tx_exons[params["exon"] - 1]
        except IndexError:
            msg = (
                f"{params['transcript']} with position {tx_pos} "
                f"does not exist on exons: {tx_exons}"
            )
            logger.warning(msg)
            return msg

        strand_to_use = strand if strand is not None else mane_data.strand
        params["strand"] = strand_to_use
        self._set_exon_offset(
            params,
            tx_exon[0],
            tx_exon[1],
            tx_pos,
            is_start=is_start,
            strand=strand_to_use,
        )

        # Need to check if we need to change pos for liftover
        genomic_data, warnings = await self.uta_db.get_alt_ac_start_or_end(
            params["transcript"], tx_pos, tx_pos, gene
        )
        if genomic_data is None:
            return warnings

        params["chr"] = genomic_data[1]
        genomic_coords = genomic_data[2], genomic_data[3]
        genomic_pos = genomic_coords[1] - 1 if is_start else genomic_coords[0] + 1
        params["pos"] = (
            genomic_pos - params["exon_offset"]
            if strand_to_use == -1
            else genomic_pos + params["exon_offset"]
        )
        return None

    async def _set_genomic_data(
        self, params: Dict, strand: Strand, is_start: bool
    ) -> Optional[str]:
        """Set genomic data in ``params``

        :param params: Parameters for response
        :param strand: Strand
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
            liftover_data = self.uta_db.get_liftover(
                chromosome_number, params["pos"], Assembly.GRCH38
            )
            if liftover_data is None:
                return (
                    f"Position {params['pos']} does not exist on "
                    f"chromosome {chromosome_number}"
                )

            params["pos"] = liftover_data[1]
            params["chr"] = grch38_ac

        tx_exons = await self._structure_exons(params["transcript"], alt_ac=grch38_ac)
        if not tx_exons:
            return f"Unable to get exons for {params['transcript']}"

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
        data_exons = data[2], data[3]
        i = 1
        found_tx_exon = False
        for exon in tx_exons:
            if data_exons == exon:
                found_tx_exon = True
                break
            i += 1
        if not found_tx_exon:
            # Either first or last
            i = 1 if data_exons == (0, tx_exons[0][1]) else i - 1
        params["exon"] = i

        strand_to_use = strand if strand is not None else Strand(data[7])
        params["strand"] = strand_to_use
        if not is_start:
            # convert back to inter-residue for end position
            params["pos"] += 1
        self._set_exon_offset(
            params,
            data[5] if is_start else data[5] + 1,  # need to convert to inter-residue
            data[6] - 1 if is_start else data[6],  # need to convert to inter-residue
            params["pos"],
            is_start=is_start,
            strand=strand_to_use,
        )
        return None

    @staticmethod
    def _set_exon_offset(
        params: Dict, start: int, end: int, pos: int, is_start: bool, strand: Strand
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

    async def _structure_exons(
        self, transcript: str, alt_ac: Optional[str] = None
    ) -> List[Tuple[int, int]]:
        """Structure exons as list of tuples.

        :param transcript: Transcript accession
        :param alt_ac: Genomic accession
        :return: List of tuples containing transcript exon coordinates
        """
        result = []
        tx_exons, _ = await self.uta_db.get_tx_exons(transcript, alt_ac=alt_ac)

        if not tx_exons:
            return result

        for coords in tx_exons:
            result.append((coords[0], coords[1]))
        return result

    @staticmethod
    def _get_exon_number(tx_exons: List, tx_pos: int) -> int:
        """Find related exon number for a position

        :param tx_exons: List of exon coordinates for a transcript
        :param tx_pos: Transcript position change
        :return: Exon number associated to transcript position change. Will be 1-based
        """
        i = 1
        for coords in tx_exons:
            if coords[0] <= tx_pos <= coords[1]:
                break
            i += 1
        return i
