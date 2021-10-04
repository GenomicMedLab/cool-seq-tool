"""Module for initializing data sources."""
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
