"""Module for initializing data sources."""
from typing import Optional, Union, List, Tuple
from uta_tools import logger
from uta_tools.schemas import GenomicData, TranscriptExonData, ResidueMode
from uta_tools.data_sources import MANETranscript, MANETranscriptMappings,\
    SeqRepoAccess, TranscriptMappings, UTADatabase
from uta_tools import SEQREPO_DATA_PATH, \
    TRANSCRIPT_MAPPINGS_PATH, LRG_REFSEQGENE_PATH, MANE_SUMMARY_PATH, \
    UTA_DB_URL


class UTATools:
    """Class to initialize data sources."""

    def __init__(self, seqrepo_data_path: str = SEQREPO_DATA_PATH,
                 transcript_file_path: str = TRANSCRIPT_MAPPINGS_PATH,
                 lrg_refseqgene_path: str = LRG_REFSEQGENE_PATH,
                 mane_data_path: str = MANE_SUMMARY_PATH,
                 db_url: str = UTA_DB_URL, db_pwd: str = '',
                 liftover_from: str = 'hg19',
                 liftover_to: str = 'hg38'
                 ):
        """Initialize UTATools class

        :param str seqrepo_data_path: The path to the seqrepo directory.
        :param str transcript_file_path: The path to transcript_mappings.tsv
        :param str lrg_refseqgene_path: The path to LRG_RefSeqGene
        :param str mane_data_path: Path to RefSeq MANE summary data
        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        :param str liftover_from: Assembly to liftover from
        :param str liftover_to: Assembly to liftover to
        """
        self.seqrepo_access = SeqRepoAccess(
            seqrepo_data_path=seqrepo_data_path)
        self.transcript_mappings = TranscriptMappings(
            transcript_file_path=transcript_file_path,
            lrg_refseqgene_path=lrg_refseqgene_path)
        self.mane_transcript_mappings = MANETranscriptMappings(
            mane_data_path=mane_data_path)
        self.uta_db = UTADatabase(
            db_url=db_url, db_pwd=db_pwd, liftover_from=liftover_from,
            liftover_to=liftover_to)
        self.mane_transcript = MANETranscript(
            self.seqrepo_access, self.transcript_mappings,
            self.mane_transcript_mappings, self.uta_db)

    async def transcript_to_genomic(
            self, exon_start: int, exon_end: int,
            exon_start_offset: int = 0, exon_end_offset: int = 0,
            gene: Optional[str] = None, transcript: str = None,
            *args, **kwargs) -> Optional[GenomicData]:
        """Get genomic data given transcript data.
        Will liftover to GRCh38 coordinates if possible.

        :param int exon_start: Starting transcript exon number
        :param int exon_end: Ending transcript exon number
        :param int exon_start_offset: Starting exon offset
        :param int exon_end_offset: Ending exon offset
        :param str gene: Gene symbol
        :param str transcript: Transcript accession
        :return: Genomic data (inter-residue coordinates)
        """
        if gene:
            gene = gene.upper().strip()

        if transcript:
            transcript = transcript.strip()

        tx_exon_start_end = await self.uta_db.get_tx_exon_start_end(
            transcript, exon_start, exon_end)
        if not tx_exon_start_end:
            return None
        (tx_exons, exon_start, exon_end) = tx_exon_start_end

        tx_exon_coords = self.uta_db.get_tx_exon_coords(
            tx_exons, exon_start, exon_end)
        if not tx_exon_coords:
            return None
        tx_exon_start, tx_exon_end = tx_exon_coords

        alt_ac_start_end = await self.uta_db.get_alt_ac_start_and_end(
            transcript, tx_exon_start, tx_exon_end, gene=gene)
        if not alt_ac_start_end:
            return None
        alt_ac_start, alt_ac_end = alt_ac_start_end

        gene = alt_ac_start[0]
        chromosome = alt_ac_start[1]
        if gene is None or chromosome is None:
            return None

        start = alt_ac_start[3]
        end = alt_ac_end[2]
        strand = alt_ac_start[4]
        if strand == -1:
            start_offset = exon_start_offset * -1
            end_offset = exon_end_offset * -1
        else:
            start_offset = exon_start_offset
            end_offset = exon_end_offset
        start += start_offset
        end += end_offset

        return GenomicData(
            gene=gene,
            chr=chromosome,
            start=start,
            end=end,
            exon_start=exon_start,
            exon_start_offset=exon_start_offset,
            exon_end=exon_end,
            exon_end_offset=exon_end_offset,
            transcript=transcript
        )

    async def genomic_to_transcript(self, chromosome: Union[str, int],
                                    start: int, end: int,
                                    strand: Optional[int] = None,
                                    transcript: Optional[str] = None,
                                    gene: Optional[str] = None,
                                    residue_mode: ResidueMode = ResidueMode.RESIDUE,  # noqa: E501
                                    *args, **kwargs) -> Optional[GenomicData]:
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
        params = {key: None for key in GenomicData.__fields__.keys()}
        if gene is not None:
            gene = gene.upper().strip()

        start_data = await self._individual_genomic_to_transcript(
            chromosome, start, strand=strand, transcript=transcript,
            gene=gene, is_start=True, residue_mode=residue_mode
        )
        if not start_data:
            return None

        end_data = await self._individual_genomic_to_transcript(
            chromosome, end, strand=strand, transcript=transcript,
            gene=gene, is_start=False, residue_mode=residue_mode
        )
        if not end_data:
            return None

        start_data = start_data.dict()
        end_data = end_data.dict()

        for field in ["transcript", "gene", "chr"]:
            if start_data[field] != end_data[field]:
                logger.warning(f"Start `{field}`, {start_data[field]}, "
                               f"does not match End `{field}`, "
                               f"{end_data[field]}")
                return None
            else:
                params[field] = start_data[field]

        if gene and gene != params["gene"]:
            logger.warning(f"Input gene, {gene}, does not match output "
                           f"gene, {params['gene']}")
            return None

        for label, data in [("start", start_data), ("end", end_data)]:
            params[label] = data["pos"]
            params[f"exon_{label}"] = data["exon"]
            params[f"exon_{label}_offset"] = data["exon_offset"]
        return GenomicData(**params)

    async def _individual_genomic_to_transcript(
            self, chromosome: Union[str, int], pos: int, strand: int = None,
            transcript: str = None, gene: str = None, is_start: bool = True,
            residue_mode: ResidueMode = ResidueMode.RESIDUE) -> Optional[TranscriptExonData]:  # noqa: E501
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
        params = {key: None for key in TranscriptExonData.__fields__.keys()}

        try:
            # Check if just chromosome number is given. If it is, we should
            # convert this to the correct accession version
            chromosome = int(chromosome)
        except ValueError:
            # Check if valid accession is given
            if not await self.uta_db.validate_genomic_ac(chromosome):
                logger.warning(f"Chromosome is not valid: {chromosome}")
                return None

        if isinstance(chromosome, str):
            # Accession given
            genes_alt_acs = await self.uta_db.chr_to_gene_and_accessions(
                chromosome, pos, strand=strand,
                alt_ac=chromosome
            )
        else:
            # Number given
            genes_alt_acs = await self.uta_db.chr_to_gene_and_accessions(
                chromosome, pos, strand=strand, alt_ac=None
            )
        gene_alt_ac = self._get_gene_and_alt_ac(genes_alt_acs, gene)
        if not gene_alt_ac:
            return None
        gene, alt_ac = gene_alt_ac

        if transcript is None:
            await self._set_mane_genomic_data(params, gene, alt_ac, pos,
                                              strand, is_start)
        else:
            params["transcript"] = transcript
            params["gene"] = gene
            params["pos"] = pos
            params["chr"] = alt_ac
            await self._set_genomic_data(params, strand, is_start)

        if residue_mode.lower().strip() == ResidueMode.RESIDUE and is_start:
            params["pos"] -= 1
        return TranscriptExonData(**params)

    def _get_gene_and_alt_ac(self,
                             genes_alt_acs: dict,
                             gene: Optional[str]) -> Optional[Tuple[str, str]]:
        """Return gene genomic accession

        :param dict genes_alt_acs: Dictionary containing genes and
            genomic accessions
        :param dict gene: Gene symbol
        :return: [Gene, Genomic accession] if both exist
        """
        alt_acs = genes_alt_acs["alt_acs"]
        len_alt_acs = len(alt_acs)
        if len_alt_acs > 1:
            logger.warning(f"Found more than one accessions: {alt_acs}")
            return None
        elif len_alt_acs == 0:
            logger.warning("No accessions found")
            return None
        alt_ac = next(iter(alt_acs))

        genes = genes_alt_acs["genes"]
        len_genes = len(genes)
        input_gene = gene
        output_gene = None
        if len_genes == 1:
            output_gene = next(iter(genes))
        elif len_genes > 1:
            logger.warning(f"Found more than one gene: {genes}")
            return None
        elif len_genes == 0:
            logger.warning("No genes found")
            return None

        if input_gene is not None:
            input_gene = input_gene.upper()
            if output_gene != input_gene:
                logger.warning(f"Input gene, {input_gene}, does not match "
                               f"output gene, {output_gene}")
                return None

        gene = output_gene if output_gene else input_gene
        return gene, alt_ac

    async def _set_mane_genomic_data(self, params: dict, gene: str,
                                     alt_ac: str, pos: int, strand: int,
                                     is_start: bool) -> None:
        """Set genomic data in `params` found from MANE.

        :param dict params: Parameters for response
        :param str gene: Gene symbol
        :param str alt_ac: Genomic accession
        :param int pos: Genomic position
        :param int strand: Strand
        :param bool is_start: `True` if `pos` is start position. `False` if
            `pos` is end position.
        """
        mane_data = await self.mane_transcript.get_mane_transcript(
            alt_ac, pos, None, 'g', gene=gene,
            normalize_endpoint=True
        )
        if not mane_data:
            logger.warning("Could not find mane data")
            return None

        params["gene"] = mane_data["gene"]
        params["transcript"] = mane_data["refseq"] if mane_data["refseq"] \
            else mane_data["ensembl"] if mane_data["ensembl"] else None
        tx_exons = await self._structure_exons(params["transcript"])
        tx_pos = mane_data["pos"][0] + mane_data["coding_start_site"]
        params["exon"] = self._get_exon_number(tx_exons, tx_pos)
        tx_exon = tx_exons[params["exon"] - 1]
        strand_to_use = strand if strand is not None else mane_data["strand"]  # noqa: E501
        self._set_exon_offset(params, tx_exon[0], tx_exon[1], tx_pos,
                              is_start=is_start, strand=strand_to_use)

        # Need to check if we need to change pos for liftover
        genomic_data = await self.uta_db.get_alt_ac_start_or_end(
            params["transcript"], tx_pos, tx_pos, gene)
        if genomic_data is None:
            return None
        params["chr"] = genomic_data[1]
        genomic_coords = genomic_data[2], genomic_data[3]
        genomic_pos = genomic_coords[1] if is_start else genomic_coords[0]
        params["pos"] = genomic_pos - params["exon_offset"] if \
            strand_to_use == -1 else genomic_pos + params["exon_offset"]

    async def _set_genomic_data(self, params: dict, strand: int,
                                is_start: bool) -> None:
        """Set genomic data in `params`.

        :param dict params: Parameters for response
        :param int strand: Strand
        :param bool is_start: `True` if `pos` is start position. `False` if
            `pos` is end position.
        """
        # We should always try to liftover
        grch38_ac = await self.uta_db.get_newest_assembly_ac(params["chr"])
        if not grch38_ac:
            logger.warning(f"Invalid genomic accession: {params['chr']}")
            return None

        grch38_ac = grch38_ac[0][0]
        if grch38_ac != params["chr"]:
            # Liftover
            descr = await self.uta_db.get_chr_assembly(params["chr"])
            if descr is None:
                return None

            chromosome_number, assembly = descr
            liftover_data = self.uta_db.get_liftover(
                chromosome_number, params["pos"])
            if liftover_data is None:
                return None

            params["pos"] = liftover_data[1]
            params["chr"] = grch38_ac

        tx_exons = await self._structure_exons(params["transcript"])
        data = await self.uta_db.get_tx_exon_aln_v_data(
            params["transcript"], params["pos"], params["pos"],
            alt_ac=params["chr"], use_tx_pos=False)
        if len(data) != 1:
            logger.warning(f"Must find exactly one row for genomic data, "
                           f"but found: {len(data)}")
            return None

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
        self._set_exon_offset(params, data[5], data[6], params["pos"],
                              is_start=is_start, strand=strand_to_use)

    def _set_exon_offset(self, params: dict, start: int, end: int, pos: int,
                         is_start: bool, strand: int) -> None:
        """Set `exon_offset` in params.

        :param dict params: Parameters for response
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

    async def _structure_exons(self, transcript: str) -> List[Tuple[int, int]]:
        """Structure exons as list of tuples.

        :param str transcript: Transcript accession
        :return: List of tuples containing transcript exon coordinates
        """
        result = list()
        tx_exons = await self.uta_db.get_tx_exons(transcript)
        for tx_exon in tx_exons:
            coords = tx_exon.split(",")
            result.append((int(coords[0]), int(coords[1])))
        return result

    def _get_exon_number(self, tx_exons: list, tx_pos: int) -> int:
        """Find exon number.

        :param list tx_exons: List of exon coordinates
        :param int tx_pos: Transcript position change
        :return: Exon number associated to transcript position change
        """
        i = 1
        for coords in tx_exons:
            if coords[0] <= tx_pos <= coords[1]:
                break
            i += 1
        return i
