"""Module for initializing data sources."""
from datetime import datetime
from typing import Optional, Union, List, Tuple, Dict
from pathlib import Path

from cool_seq_tool import logger
from cool_seq_tool.data_sources.alignment_mapper import AlignmentMapper
from cool_seq_tool.schemas import Assembly, GenomicData, TranscriptExonData, \
    ResidueMode, GenomicDataResponse, ServiceMeta, TranscriptExonDataResponse
from cool_seq_tool.data_sources import MANETranscript, MANETranscriptMappings,\
    SeqRepoAccess, TranscriptMappings, UTADatabase, GeneNormalizer
from cool_seq_tool import SEQREPO_DATA_PATH, \
    TRANSCRIPT_MAPPINGS_PATH, LRG_REFSEQGENE_PATH, MANE_SUMMARY_PATH, \
    UTA_DB_URL
from cool_seq_tool.version import __version__


class CoolSeqTool:
    """Class to initialize data sources."""

    def __init__(self, seqrepo_data_path: str = SEQREPO_DATA_PATH,
                 transcript_file_path: str = TRANSCRIPT_MAPPINGS_PATH,
                 lrg_refseqgene_path: str = LRG_REFSEQGENE_PATH,
                 mane_data_path: str = MANE_SUMMARY_PATH,
                 db_url: str = UTA_DB_URL, db_pwd: str = "",
                 gene_db_url: str = "", gene_db_region: str = "us-east-2"
                 ) -> None:
        """Initialize CoolSeqTool class

        :param str seqrepo_data_path: The path to the seqrepo directory.
        :param str transcript_file_path: The path to transcript_mappings.tsv
        :param str lrg_refseqgene_path: The path to LRG_RefSeqGene
        :param str mane_data_path: Path to RefSeq MANE summary data
        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        :param str gene_db_url: URL to gene normalizer dynamodb
        :param str gene_db_region: AWS region for gene normalizer db
        """
        self.seqrepo_access = SeqRepoAccess(
            seqrepo_data_path=seqrepo_data_path)
        self.transcript_mappings = TranscriptMappings(
            transcript_file_path=transcript_file_path,
            lrg_refseqgene_path=lrg_refseqgene_path)
        self.mane_transcript_mappings = MANETranscriptMappings(
            mane_data_path=mane_data_path)
        self.uta_db = UTADatabase(db_url=db_url, db_pwd=db_pwd)
        gene_normalizer = GeneNormalizer(gene_db_url, gene_db_region)
        self.alignment_mapper = AlignmentMapper(
            self.seqrepo_access, self.transcript_mappings, self.uta_db)
        self.mane_transcript = MANETranscript(
            self.seqrepo_access, self.transcript_mappings,
            self.mane_transcript_mappings, self.uta_db, gene_normalizer)

    @staticmethod
    def service_meta() -> ServiceMeta:
        """Return ServiceMeta for cool_seq_tool

        :return: ServiceMeta object
        """
        return ServiceMeta(
            version=__version__,
            response_datetime=datetime.now()
        )

    @staticmethod
    def _return_warnings(
            resp: Union[GenomicDataResponse, TranscriptExonDataResponse],
            warning_msg: str) -> Union[GenomicDataResponse, TranscriptExonDataResponse]:  # noqa: E501
        """Add warnings to response object

        :param Union[GenomicDataResponse, TranscriptExonDataResponse] resp:
            Response object
        :param str warning_msg: Warning message on why `transcript_exon_data`
            or `genomic_data` field is None
        :return: Response object with warning message
        """
        logger.warning(warning_msg)
        resp.warnings.append(warning_msg)
        return resp

    async def transcript_to_genomic_coordinates(
            self, gene: Optional[str] = None, transcript: Optional[str] = None,
            exon_start: Optional[int] = None, exon_start_offset: Optional[int] = 0,  # noqa: E501
            exon_end: Optional[int] = None, exon_end_offset: Optional[int] = 0,
            **kwargs) -> GenomicDataResponse:
        """Get genomic data given transcript data.
        Will use GRCh38 coordinates if possible

        :param Optional[str] gene: Gene symbol
        :param Optional[str] transcript: Transcript accession
        :param Optional[int] exon_start: Starting transcript exon number
        :param Optional[int] exon_end: Ending transcript exon number
        :param Optional[int] exon_start_offset: Starting exon offset
        :param Optional[int] exon_end_offset: Ending exon offset
        :return: GRCh38 genomic data (inter-residue coordinates)
        """
        resp = GenomicDataResponse(
            genomic_data=None,
            warnings=[],
            service_meta=self.service_meta()
        )

        if not transcript:
            return self._return_warnings(resp, "Must provide `transcript`")
        else:
            transcript = transcript.strip()

        if exon_start is None and exon_end is None:
            return self._return_warnings(
                resp, "Must provide either `exon_start` or `exon_end`")

        if gene:
            gene = gene.upper().strip()

        if exon_start and exon_end:
            if exon_start > exon_end:
                return self._return_warnings(
                    resp,
                    f"Start exon {exon_start} is greater than end exon {exon_end}"  # noqa: E501
                )

        tx_exons, warning = await self.uta_db.get_tx_exons(transcript)
        if not tx_exons:
            return self._return_warnings(resp, warning)

        tx_exon_coords, warning = self.uta_db.get_tx_exon_coords(
            transcript, tx_exons, exon_start, exon_end)
        if not tx_exon_coords:
            return self._return_warnings(resp, warning)
        tx_exon_start, tx_exon_end = tx_exon_coords

        alt_ac_start_end, warning = await self.uta_db.get_alt_ac_start_and_end(
            transcript, tx_exon_start, tx_exon_end, gene=gene)
        if not alt_ac_start_end:
            return self._return_warnings(resp, warning)
        alt_ac_start, alt_ac_end = alt_ac_start_end

        gene = alt_ac_start[0] if alt_ac_start else alt_ac_end[0]
        chromosome = alt_ac_start[1] if alt_ac_start else alt_ac_end[1]
        if gene is None or chromosome is None:
            return self._return_warnings(
                resp, "Unable to retrieve `gene` or `chromosome` from "
                      "genomic start or end data")

        start = alt_ac_start[3] if alt_ac_start else None
        end = alt_ac_end[2] if alt_ac_end else None
        strand = alt_ac_start[4] if alt_ac_start else alt_ac_end[4]

        # Using none since could set to 0
        start_exits = start is not None
        end_exists = end is not None

        if strand == -1:
            start_offset = exon_start_offset * -1 if start_exits else None
            end_offset = exon_end_offset * -1 if end_exists else None
        else:
            start_offset = exon_start_offset if start_exits else None
            end_offset = exon_end_offset if end_exists else None

        start = start + start_offset if start_exits else None
        end = end + end_offset if end_exists else None

        resp.genomic_data = GenomicData(
            gene=gene,
            chr=chromosome,
            start=start,
            end=end,
            exon_start=exon_start if start_exits else None,
            exon_start_offset=exon_start_offset if start_exits else None,
            exon_end=exon_end if end_exists else None,
            exon_end_offset=exon_end_offset if end_exists else None,
            transcript=transcript,
            strand=strand
        )

        return resp

    async def genomic_to_transcript_exon_coordinates(
            self, chromosome: Union[str, int], start: Optional[int] = None,
            end: Optional[int] = None, strand: Optional[int] = None,
            transcript: Optional[str] = None, gene: Optional[str] = None,
            residue_mode: ResidueMode = ResidueMode.RESIDUE,
            **kwargs) -> GenomicDataResponse:
        """Get transcript data for genomic data.
        MANE Transcript data will be returned iff `transcript` is not supplied.
            `gene` must be supplied in order to retrieve MANE Transcript data.
        Liftovers genomic coordinates to GRCh38

        :param str chromosome: Chromosome. Must either give chromosome number
            (i.e. `1`) or accession (i.e. `NC_000001.11`).
        :param int start: Start genomic position
        :param int end: End genomic position
        :param str strand: Strand. Must be either `-1` or `1`.
        :param str transcript: The transcript to use. If this is not given,
            we will try the following transcripts: MANE Select, MANE Clinical
            Plus, Longest Remaining Compatible Transcript
        :param str gene: Gene symbol
        :param str residue_mode: Default is `resiude` (1-based).
            Must be either `residue` or `inter-residue` (0-based).
        :return: Genomic data (inter-residue coordinates)
        """
        resp = GenomicDataResponse(
            genomic_data=None,
            warnings=[],
            service_meta=self.service_meta()
        )
        if start is None and end is None:
            return self._return_warnings(
                resp, "Must provide either `start` or `end`")

        params = {key: None for key in GenomicData.__fields__.keys()}
        if gene is not None:
            gene = gene.upper().strip()

        if start:
            if residue_mode == ResidueMode.RESIDUE:
                start -= 1
            start_data = await self._genomic_to_transcript_exon_coordinate(
                chromosome, start, strand=strand, transcript=transcript,
                gene=gene, is_start=True,
                residue_mode=ResidueMode.INTER_RESIDUE
            )
            if start_data.transcript_exon_data:
                start_data = start_data.transcript_exon_data.dict()
            else:
                return self._return_warnings(resp, start_data.warnings[0])
        else:
            start_data = None

        if end:
            if residue_mode == ResidueMode.RESIDUE:
                end -= 1
            end_data = await self._genomic_to_transcript_exon_coordinate(
                chromosome, end, strand=strand, transcript=transcript,
                gene=gene, is_start=False,
                residue_mode=ResidueMode.INTER_RESIDUE
            )
            if end_data.transcript_exon_data:
                end_data = end_data.transcript_exon_data.dict()
            else:
                return self._return_warnings(resp, end_data.warnings[0])
        else:
            end_data = None

        for field in ["transcript", "gene", "chr", "strand"]:
            if start_data:
                if end_data:
                    if start_data[field] != end_data[field]:
                        msg = f"Start `{field}`, {start_data[field]}, does " \
                              f"not match End `{field}`, {end_data[field]}"
                        return self._return_warnings(resp, msg)
                params[field] = start_data[field]
            else:
                params[field] = end_data[field]

        if gene and gene != params["gene"]:
            msg = f"Input gene, {gene}, does not match expected output" \
                  f"gene, {params['gene']}"
            return self._return_warnings(resp, msg)

        for label, data in [("start", start_data), ("end", end_data)]:
            if data:
                params[label] = data["pos"]
                params[f"exon_{label}"] = data["exon"]
                params[f"exon_{label}_offset"] = data["exon_offset"]
        resp.genomic_data = GenomicData(**params)
        return resp

    async def _genomic_to_transcript_exon_coordinate(
            self, chromosome: Union[str, int], pos: int, strand: int = None,
            transcript: str = None, gene: str = None, is_start: bool = True,
            residue_mode: ResidueMode = ResidueMode.RESIDUE) -> TranscriptExonDataResponse:  # noqa: E501
        """Convert individual genomic data to transcript data

        :param str chromosome: Chromosome. Must either give chromosome number
            (i.e. `1`) or accession (i.e. `NC_000001.11`).
        :param int pos: Genomic position
        :param str strand: Strand. Must be either `-1` or `1`.
        :param str transcript: The transcript to use. If this is not given,
            we will try the following transcripts: MANE Select, MANE Clinical
            Plus, Longest Remaining Compatible Transcript
        :param str gene: Gene symbol
        :param bool is_start: `True` if `pos` is start position. `False` if
            `pos` is end position.
        :param str residue_mode: Default is `resiude` (1-based).
            Must be either `residue` or `inter-residue` (0-based).
        :return: Transcript data (inter-residue coordinates)
        """
        resp = TranscriptExonDataResponse(
            transcript_exon_data=None,
            warnings=[],
            service_meta=self.service_meta()
        )

        if transcript is None and gene is None:
            return self._return_warnings(
                resp, "Must provide either `gene` or `transcript`"
            )

        params = {key: None for key in TranscriptExonData.__fields__.keys()}

        try:
            # Check if just chromosome is given. If it is, we should
            # convert this to the correct accession version
            if chromosome == "X":
                chromosome = 23
            elif chromosome == "Y":
                chromosome = 24
            else:
                chromosome = int(chromosome)
        except ValueError:
            # Check if valid accession is given
            if not await self.uta_db.validate_genomic_ac(chromosome):
                return self._return_warnings(
                    resp, f"Invalid chromosome: {chromosome}")

        if isinstance(chromosome, str):
            # Accession given
            genes_alt_acs, warning = \
                await self.uta_db.chr_to_gene_and_accessions(
                    chromosome, pos, strand=strand, alt_ac=chromosome, gene=gene)
        else:
            # Number given
            genes_alt_acs, warning = \
                await self.uta_db.chr_to_gene_and_accessions(
                    chromosome, pos, strand=strand, alt_ac=None, gene=gene)
        if not genes_alt_acs:
            return self._return_warnings(resp, warning)

        gene_alt_ac, warning = self._get_gene_and_alt_ac(genes_alt_acs, gene)
        if not gene_alt_ac:
            return self._return_warnings(resp, warning)
        gene, alt_ac = gene_alt_ac

        if transcript is None:
            warnings = await self._set_mane_genomic_data(
                params, gene, alt_ac, pos, strand, is_start, residue_mode)
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

        :param Dict genes_alt_acs: Dictionary containing genes and
            genomic accessions
        :param Optional[str] gene: Gene symbol
        :return: [Gene, Genomic accession] if both exist
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
                return None, f"Input gene, {input_gene}, does not match " \
                             f"expected output gene, {output_gene}"

        gene = output_gene if output_gene else input_gene
        return (gene, alt_ac), None

    async def _set_mane_genomic_data(
            self, params: Dict, gene: str, alt_ac: str, pos: int, strand: int,
            is_start: bool, residue_mode: str
    ) -> Optional[str]:
        """Set genomic data in `params` found from MANE.

        :param Dict params: Parameters for response
        :param str gene: Gene symbol
        :param str alt_ac: Genomic accession
        :param int pos: Genomic position
        :param int strand: Strand
        :param bool is_start: `True` if `pos` is start position. `False` if
            `pos` is end position.
        :param str residue_mode: Residue mode for start/end positions
            Must be either `inter-residue` or `residue`
        :return: Warnings if found
        """
        mane_data = await self.mane_transcript.get_mane_transcript(
            alt_ac, pos, "g", gene=gene,
            try_longest_compatible=True, residue_mode=residue_mode
        )
        if not mane_data:
            msg = f"Unable to find mane data for {alt_ac} with position {pos}"
            if gene:
                msg += f" on gene {gene}"
            logger.warning(msg)
            return msg

        if mane_data["strand"] == "-":
            mane_data["strand"] = -1
        elif mane_data["strand"] == "+":
            mane_data["strand"] = 1

        params["gene"] = mane_data["gene"]
        params["transcript"] = mane_data["refseq"] if mane_data["refseq"] \
            else mane_data["ensembl"] if mane_data["ensembl"] else None
        tx_exons = await self._structure_exons(params["transcript"], alt_ac=alt_ac)
        if not tx_exons:
            return f"Unable to get exons for {params['transcript']}"
        tx_pos = mane_data["pos"][0] + mane_data["coding_start_site"]
        params["exon"] = self._get_exon_number(tx_exons, tx_pos)

        try:
            tx_exon = tx_exons[params["exon"] - 1]
        except IndexError:
            msg = f"{params['transcript']} with position {tx_pos} "\
                  f"does not exist on exons: {tx_exons}"
            logger.warning(msg)
            return msg

        strand_to_use = strand if strand is not None else mane_data["strand"]
        params["strand"] = strand_to_use
        self._set_exon_offset(params, tx_exon[0], tx_exon[1], tx_pos,
                              is_start=is_start, strand=strand_to_use)

        # Need to check if we need to change pos for liftover
        genomic_data, warnings = await self.uta_db.get_alt_ac_start_or_end(
            params["transcript"], tx_pos, tx_pos, gene)
        if genomic_data is None:
            return warnings

        params["chr"] = genomic_data[1]
        genomic_coords = genomic_data[2], genomic_data[3]
        genomic_pos = genomic_coords[1] if is_start else genomic_coords[0]
        params["pos"] = genomic_pos - params["exon_offset"] if \
            strand_to_use == -1 else genomic_pos + params["exon_offset"]
        return None

    async def _set_genomic_data(self, params: Dict, strand: int,
                                is_start: bool) -> Optional[str]:
        """Set genomic data in `params`.

        :param Dict params: Parameters for response
        :param int strand: Strand
        :param bool is_start: `True` if `pos` is start position. `False` if
            `pos` is end position.
        :return: Warnings if found
        """
        # We should always try to liftover
        grch38_ac = await self.uta_db.get_newest_assembly_ac(params["chr"])
        if not grch38_ac:
            return f"Invalid genomic accession: {params['chr']}"

        grch38_ac = grch38_ac[0][0]
        if grch38_ac != params["chr"]:  # params["chr"] is genomic accession
            # Liftover to 38
            descr = await self.uta_db.get_chr_assembly(params["chr"])
            if descr is None:
                return f"Unable to get chromosome and assembly for " \
                       f"{params['chr']}"

            chromosome_number, assembly = descr
            liftover_data = self.uta_db.get_liftover(
                chromosome_number, params["pos"], Assembly.GRCH38)
            if liftover_data is None:
                return f"Position {params['pos']} does not exist on " \
                       f"chromosome {chromosome_number}"

            params["pos"] = liftover_data[1]
            params["chr"] = grch38_ac

        tx_exons = await self._structure_exons(params["transcript"], alt_ac=grch38_ac)
        if not tx_exons:
            return f"Unable to get exons for {params['transcript']}"
        data = await self.uta_db.get_tx_exon_aln_v_data(
            params["transcript"], params["pos"], params["pos"],
            alt_ac=params["chr"], use_tx_pos=False)
        if len(data) != 1:
            return f"Must find exactly one row for genomic data, " \
                   f"but found: {len(data)}"

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

        strand_to_use = strand if strand is not None else data[7]
        params["strand"] = strand_to_use
        self._set_exon_offset(params, data[5], data[6], params["pos"],
                              is_start=is_start, strand=strand_to_use)
        return None

    @staticmethod
    def _set_exon_offset(params: Dict, start: int, end: int, pos: int,
                         is_start: bool, strand: int) -> None:
        """Set `exon_offset` in params.

        :param Dict params: Parameters for response
        :param int start: Start exon coord (can be transcript or genomic)
        :param int end: End exon coord (can be transcript or genomic)
        :param int pos: Position change (can be transcript or genomic)
        :param bool is_start: `True` if `pos` is start position.
            `False` if `pos` is end position
        :param int strand: Strand
        """
        if is_start:
            if strand == -1:
                params["exon_offset"] = end - pos
            else:
                params["exon_offset"] = pos - end
        else:
            if strand == -1:
                params["exon_offset"] = start - pos
            else:
                params["exon_offset"] = pos - start

    async def _structure_exons(
        self, transcript: str, alt_ac: Optional[str] = None
    ) -> List[Tuple[int, int]]:
        """Structure exons as list of tuples.

        :param str transcript: Transcript accession
        :param Optional[str] alt_ac: Genomic accession
        :return: List of tuples containing transcript exon coordinates
        """
        result = list()
        tx_exons, _ = await self.uta_db.get_tx_exons(transcript, alt_ac=alt_ac)
        if not tx_exons:
            return result
        for coords in tx_exons:
            result.append((coords[0], coords[1]))
        return result

    @staticmethod
    def _get_exon_number(tx_exons: List, tx_pos: int) -> int:
        """Find exon number.

        :param List tx_exons: List of exon coordinates
        :param int tx_pos: Transcript position change
        :return: Exon number associated to transcript position change
        """
        i = 1
        for coords in tx_exons:
            if coords[0] <= tx_pos <= coords[1]:
                break
            i += 1
        return i

    def get_fasta_file(
        self, sequence_id: str, outfile_path: Path
    ) -> None:
        """Retrieve FASTA file containing sequence for requested sequence ID.
        :param sequence_id: accession ID, sans namespace, eg `NM_152263.3`
        :param outfile_path: path to save file to
        :return: None, but saves sequence data to `outfile_path` if successful
        :raise: KeyError if SeqRepo doesn't have sequence data for the given ID
        """
        sequence = self.seqrepo_access.get_reference_sequence(sequence_id)[0]
        if not sequence:
            raise KeyError

        REFSEQ_PREFIXES = [
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
            "WP_"
        ]
        ENSEMBL_PREFIXES = [
            "ENSE",
            "ENSFM",
            "ENSG",
            "ENSGT",
            "ENSP",
            "ENSR",
            "ENST"
        ]

        if sequence_id[:3] in REFSEQ_PREFIXES:
            aliases = self.seqrepo_access.translate_identifier(
                sequence_id, ["ensembl", "ga4gh"]
            )
            header = f">ref|refseq:{sequence_id}|{'|'.join(aliases[0])}"
        elif sequence_id[:4] in ENSEMBL_PREFIXES:
            aliases = self.seqrepo_access.translate_identifier(
                sequence_id, ["refseq", "ga4gh"]
            )
            header = f">emb|ensembl:{sequence_id}|{'|'.join(aliases[0])}"
        else:
            aliases = self.seqrepo_access.translate_identifier(
                sequence_id, ["ensembl", "refseq", "ga4gh"]
            )
            header = f">gnl|ID|{sequence_id}|{'|'.join(aliases[0])}"

        LINE_LENGTH = 60
        file_data = [header] + [sequence[i: i + LINE_LENGTH]
                                for i in range(0, len(sequence), LINE_LENGTH)]
        text = "\n".join(file_data)
        outfile_path.write_text(text)
