"""Module for initializing data sources."""
from typing import Optional, Dict, Union
from uta_tools import logger
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

    async def genomic_to_transcript(self, chromosome: Union[str, int],
                                    start_pos: int, end_pos: int,
                                    strand: int = None, transcript: str = None,
                                    gene: str = None,
                                    residue_mode: str = "residue")\
            -> Optional[Dict]:
        """Get transcript data for genomic range data.

        :param str chromosome: Chromosome. Must either give chromosome number
            (i.e. `1`) or accession (i.e. `NC_000001.11`).
        :param int start_pos: Start position
        :param int end_pos: End position
        :param str strand: Strand. Must be either `-1` or `1`.
        :param str transcript: The transcript to use. If this is not given,
            we will try the following transcripts: MANE Select, MANE Clinical
            Plus, Longest Remaining Compatible Transcript
        :param str gene: Gene
        :param str residue_mode: Default is `resiude` (1-based).
            Must be either `residue` or `inter-residue` (0-based).
        :return: Transcript data (inter-residue coordinates)
        """
        params = {
            "start": None,
            "exon_start": None,
            "exon_start_offset": 0,
            "end": None,
            "exon_end": None,
            "exon_end_offset": 0,
            "gene": None,
            "transcript": None
        }
        start = await self._genomic_to_transcript(
            chromosome, start_pos, strand=strand, transcript=transcript,
            gene=gene, residue_mode=residue_mode
        )
        if not start:
            return None

        end = await self._genomic_to_transcript(
            chromosome, end_pos, strand=strand, transcript=transcript,
            gene=gene, residue_mode=residue_mode
        )
        if not end:
            return None

        for field in ["transcript", "gene"]:
            if start[field] != end[field]:
                logger.warning(f"Start `{field}`, {start[field]}, does not "
                               f"match End `{field}`, {end[field]}")
                return None
            else:
                params[field] = start[field]

        params["start"] = start["pos"] - 1 if residue_mode.lower() == "residue" else start["pos"]  # noqa: E501
        params["exon_start"] = start["exon"]
        params["exon_start_offset"] = start["exon_offset"] if start["exon_offset"] is not None else 0  # noqa: E501

        params["end"] = end["pos"]
        params["exon_end"] = end["exon"]
        params["exon_end_offset"] = end["exon_offset"] if end["exon_offset"] is not None else 0  # noqa: E501
        return params

    async def _genomic_to_transcript(self, chromosome: Union[str, int],
                                     pos: int, strand: int = None,
                                     transcript: str = None, gene: str = None,
                                     residue_mode: str = 'residue')\
            -> Optional[Dict]:
        """Get transcript data given genomic data.

        :param str chromosome: Chromosome. Must either give chromosome number
            (i.e. `1`) or accession (i.e. `NC_000001.11`).
        :param int pos: Position
        :param str strand: Strand. Must be either `-1` or `1`.
        :param str transcript: The transcript to use. If this is not given,
            we will try the following transcripts: MANE Select, MANE Clinical
            Plus, Longest Remaining Compatible Transcript
        :param str gene: Gene
        :param str residue_mode: Default is `resiude` (1-based).
            Must be either `residue` or `inter-residue` (0-based).
        :return: Transcript data (inter-residue coordinates)
        """
        result = {
            "transcript": transcript,
            "pos": pos,
            "exon": None,
            "exon_offset": 0,
            "gene": None
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
        genes = genes_alt_acs["genes"]
        len_alt_acs = len(alt_acs)

        if len_alt_acs > 1:
            logger.warning(f"Found more than one accessions: {alt_acs}")
            return None
        elif len_alt_acs == 0:
            logger.warning("No accessions found")
            return None

        alt_ac = next(iter(alt_acs))
        if gene is None:
            if len(genes) == 1:
                gene = next(iter(genes))
                result["gene"] = gene

        if transcript is None:
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

            tx_exons = await self._get_exons_tuple(transcript)

            tx_pos = mane_data["pos"][0] + mane_data["coding_start_site"]
            exon_pos = self._find_exon(tx_exons, tx_pos)
            result["exon"] = exon_pos
        else:
            tx_exons = await self._get_exons_tuple(transcript)
            data = await self.uta_db.get_tx_exon_aln_v_data(
                transcript, pos, pos, alt_ac=alt_ac, use_tx_pos=False
            )
            if len(data) == 1:
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
                    i = 1
                result["exon"] = i
        return result

    async def _get_exons_tuple(self, transcript: str):
        """Reformat exons as tuples"""
        result = list()
        tx_exons = await self.uta_db.get_tx_exons(transcript)
        for tx_exon in tx_exons:
            coords = tx_exon.split(",")
            result.append((int(coords[0]), int(coords[1])))
        return result

    def _find_exon(self, tx_exons: list, tx_pos: int):
        """Find exon number"""
        i = 0
        for coords in tx_exons:
            if coords[0] <= tx_pos <= coords[1]:
                break
            i += 1
        return i
