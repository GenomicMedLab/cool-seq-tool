"""Module for initializing data sources."""
from typing import Optional, Dict, Union, List, Tuple
from uta_tools import logger
from uta_tools.schemas import GenomicData
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
            self, transcript: str, exon_start: int, exon_end: int,
            exon_start_offset: int = 0, exon_end_offset: int = 0,
            gene: str = None, *args, **kwargs) -> Optional[GenomicData]:
        """Get genomic data given transcript data.
        Will liftover to GRCh38 coordinates if possible.

        :param str transcript: Transcript accession
        :param int exon_start: Starting exon number
        :param int exon_end: Ending exon number
        :param int exon_start_offset: Starting exon offset
        :param int exon_end_offset: Ending exon offset
        :param str gene: Gene symbol
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

        gene = alt_ac_start[0]
        chromosome = alt_ac_start[1]
        if gene is None or chromosome is None:
            return None
        return GenomicData(
            gene=gene,
            chr=chromosome,
            start=start - 1,
            end=end,
            exon_start=exon_start,
            exon_start_offset=exon_start_offset,
            exon_end=exon_end,
            exon_end_offset=exon_end_offset,
            transcript=transcript
        )

    async def genomic_to_transcript(self, chromosome: Union[str, int],
                                    start: int, end: int,
                                    strand: int = None, transcript: str = None,
                                    gene: str = None,
                                    residue_mode: str = "residue",
                                    *args, **kwargs) -> Optional[GenomicData]:
        """Get transcript data for genomic range data.
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
        :param str gene: Gene
        :param str residue_mode: Default is `resiude` (1-based).
            Must be either `residue` or `inter-residue` (0-based).
        :return: Genomic data (inter-residue coordinates)
        """
        params = {
            "start": None,
            "exon_start": None,
            "exon_start_offset": 0,
            "end": None,
            "exon_end": None,
            "exon_end_offset": 0,
            "gene": None,
            "transcript": None,
            "chr": None
        }
        if gene is not None:
            gene = gene.upper()
        start_data = await self._genomic_to_transcript(
            chromosome, start, strand=strand, transcript=transcript,
            gene=gene, is_start=True
        )
        if not start_data:
            return None

        end_data = await self._genomic_to_transcript(
            chromosome, end, strand=strand, transcript=transcript,
            gene=gene, is_start=False
        )
        if not end_data:
            return None

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

        params["start"] = start_data["pos"] - 1 if residue_mode.lower() == "residue" else start_data["pos"]  # noqa: E501
        params["exon_start"] = start_data["exon"]
        params["exon_start_offset"] = start_data["exon_offset"] if start_data["exon_offset"] is not None else 0  # noqa: E501

        params["end"] = end_data["pos"]
        params["exon_end"] = end_data["exon"]
        params["exon_end_offset"] = end_data["exon_offset"] if end_data["exon_offset"] is not None else 0  # noqa: E501
        return GenomicData(**params)

    async def _genomic_to_transcript(self, chromosome: Union[str, int],
                                     pos: int, strand: int = None,
                                     transcript: str = None,
                                     gene: str = None,
                                     is_start: bool = True) -> Optional[Dict]:
        """Get transcript data given genomic data.

        :param str chromosome: Chromosome. Must either give chromosome number
            (i.e. `1`) or accession (i.e. `NC_000001.11`).
        :param int pos: Position
        :param str strand: Strand. Must be either `-1` or `1`.
        :param str transcript: The transcript to use. If this is not given,
            we will try the following transcripts: MANE Select, MANE Clinical
            Plus, Longest Remaining Compatible Transcript
        :param str gene: Gene
        :param bool is_start: `True` if `pos` is start position. `False` if
            `pos` is end position.
        :return: Transcript data (inter-residue coordinates)
        """
        result = {
            "transcript": transcript,
            "pos": pos,
            "exon": None,
            "exon_offset": 0,
            "gene": None,
            "chr": None
        }
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
            genes_alt_acs = await self.uta_db.chr_to_accession(
                chromosome, pos, strand=strand,
                alt_ac=chromosome
            )
        else:
            # Number given
            genes_alt_acs = await self.uta_db.chr_to_accession(
                chromosome, pos, strand=strand, alt_ac=None
            )
        alt_acs = genes_alt_acs["alt_acs"]
        len_alt_acs = len(alt_acs)
        if len_alt_acs > 1:
            logger.warning(f"Found more than one accessions: {alt_acs}")
            return None
        elif len_alt_acs == 0:
            logger.warning("No accessions found")
            return None
        alt_ac = next(iter(alt_acs))
        result["chr"] = alt_ac

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
            if gene and gene != input_gene:
                logger.warning(f"Input gene, {input_gene}, does not match "
                               f"output gene, {output_gene}")
                return None

        result["gene"] = output_gene if output_gene else input_gene

        if transcript is None:
            if gene is None:
                logger.warning("Must provide gene for MANE transcript")
                return None
            mane_data = await self.mane_transcript.get_mane_transcript( # noqa
                alt_ac, pos, None, 'g', gene=gene,
                normalize_endpoint=True
            )

            if mane_data["refseq"]:
                transcript = mane_data["refseq"]
            elif mane_data["ensembl"]:
                transcript = mane_data["ensembl"]
            else:
                transcript = None
            result["transcript"] = transcript

            tx_exons = await self._structure_exons(transcript)
            tx_pos = mane_data["pos"][0] + mane_data["coding_start_site"]
            exon_pos = self._get_exon_number(tx_exons, tx_pos)
            tx_exon = tx_exons[exon_pos - 1]

            strand_to_use = strand if strand is not None else mane_data["strand"]  # noqa: E501
            if is_start:
                if strand_to_use == -1:
                    result["exon_offset"] = tx_exon[1] - tx_pos
                else:
                    result["exon_offset"] = tx_pos - tx_exon[1]
            else:
                if strand_to_use == -1:
                    result["exon_offset"] = tx_exon[0] - tx_pos
                else:
                    result["exon_offset"] = tx_pos - tx_exon[0]

            # Need to check if we need to change pos
            genomic_data = await self.uta_db.get_alt_ac_start_or_end(
                transcript, tx_pos, tx_pos, gene)
            if genomic_data is None:
                return None
            result["chr"] = genomic_data[1]

            genomic_coords = genomic_data[2], genomic_data[3]
            if is_start:
                genomic_pos = genomic_coords[1]
            else:
                genomic_pos = genomic_coords[0]

            if strand_to_use == -1:
                result["pos"] = genomic_pos - result["exon_offset"]
            else:
                result["pos"] = genomic_pos + result["exon_offset"]
            result["exon"] = exon_pos
        else:
            # We should always try to liftover
            grch38_ac = await self.uta_db.get_newest_assembly_ac(alt_ac)
            if not grch38_ac:
                logger.warning(f"Invalid genomic accession: {alt_ac}")
                return None
            grch38_ac = grch38_ac[0][0]
            if grch38_ac != alt_ac:
                # Liftover
                descr = await self.uta_db.get_chr_assembly(alt_ac)
                if descr is None:
                    return None
                chromosome_number, assembly = descr
                liftover_data = self.uta_db.get_liftover(
                    chromosome_number, pos)
                if liftover_data is None:
                    return None
                pos = liftover_data[1]
                result["pos"] = pos

                alt_ac = grch38_ac
                result["chr"] = alt_ac

            tx_exons = await self._structure_exons(transcript)
            data = await self.uta_db.get_tx_exon_aln_v_data(
                transcript, pos, pos, alt_ac=alt_ac, use_tx_pos=False
            )
            if len(data) == 1:
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
                    first_exon = 0, tx_exons[0][1]
                    if data_exons == first_exon:
                        i = 1
                    else:
                        i -= 1
                result["exon"] = i

                # Find exon offset
                alt_start_i = data[5]
                alt_end_i = data[6]
                uta_strand = data[7]
                strand_to_use = strand if strand is not None else uta_strand
                if is_start:
                    if strand_to_use == - 1:
                        result["exon_offset"] = alt_end_i - pos
                    else:
                        result["exon_offset"] = pos - alt_end_i
                else:
                    if strand_to_use == -1:
                        result["exon_offset"] = alt_start_i - pos
                    else:
                        result["exon_offset"] = pos - alt_start_i

        return result

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
