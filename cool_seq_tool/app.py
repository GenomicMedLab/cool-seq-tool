"""Module for initializing data sources."""
from typing import Optional
from pathlib import Path
import logging

from biocommons.seqrepo import SeqRepo

from cool_seq_tool.mappers import (
    MANETranscript, AlignmentMapper, ExonGenomicCoordsMapper
)
from cool_seq_tool.sources.uta_database import UTA_DB_URL, UTADatabase
from cool_seq_tool.sources.mane_transcript_mappings import MANETranscriptMappings
from cool_seq_tool.sources.transcript_mappings import TranscriptMappings
from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.paths import LRG_REFSEQGENE_PATH, MANE_SUMMARY_PATH, \
    SEQREPO_ROOT_DIR, TRANSCRIPT_MAPPINGS_PATH


logger = logging.getLogger(__name__)


class CoolSeqTool:
    """Class to initialize data sources."""

    def __init__(
        self,
        transcript_file_path: Path = TRANSCRIPT_MAPPINGS_PATH,
        lrg_refseqgene_path: Path = LRG_REFSEQGENE_PATH,
        mane_data_path: Path = MANE_SUMMARY_PATH,
        db_url: str = UTA_DB_URL,
        sr: Optional[SeqRepo] = None
    ) -> None:
        """Initialize CoolSeqTool class

        :param transcript_file_path: The path to transcript_mapping.tsv
        :param lrg_refseqgene_path: The path to LRG_RefSeqGene
        :param mane_data_path: Path to RefSeq MANE summary data
        :param db_url: PostgreSQL connection URL
            Format: `driver://user:password@host/database/schema`
        :param sr: SeqRepo instance. If this is not provided, will create a new instance
        """
        if not sr:
            sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
        self.seqrepo_access = SeqRepoAccess(sr)
        self.transcript_mappings = TranscriptMappings(
            transcript_file_path=transcript_file_path,
            lrg_refseqgene_path=lrg_refseqgene_path)
        self.mane_transcript_mappings = MANETranscriptMappings(
            mane_data_path=mane_data_path)
        self.uta_db = UTADatabase(db_url=db_url)
        self.alignment_mapper = AlignmentMapper(
            self.seqrepo_access, self.transcript_mappings, self.uta_db)
        self.mane_transcript = MANETranscript(
            self.seqrepo_access, self.transcript_mappings,
            self.mane_transcript_mappings, self.uta_db)
        self.exon_genomic_coords_mapper = ExonGenomicCoordsMapper(self.uta_db,
                                                                  self.mane_transcript)
