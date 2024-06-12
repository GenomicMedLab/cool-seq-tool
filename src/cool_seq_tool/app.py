"""Provides core CoolSeqTool class, which non-redundantly initializes all Cool-Seq-Tool
data handler and mapping resources for straightforward access.
"""
import logging
from pathlib import Path
from typing import Optional

from biocommons.seqrepo import SeqRepo

from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.mappers import (
    AlignmentMapper,
    ExonGenomicCoordsMapper,
    ManeTranscript,
)
from cool_seq_tool.resources import (
    SEQREPO_ROOT_DIR,
    TRANSCRIPT_MAPPINGS_PATH,
)
from cool_seq_tool.sources.mane_transcript_mappings import ManeTranscriptMappings
from cool_seq_tool.sources.transcript_mappings import TranscriptMappings
from cool_seq_tool.sources.uta_database import UTA_DB_URL, UtaDatabase

logger = logging.getLogger(__name__)


class CoolSeqTool:
    """Non-redundantly initialize all Cool-Seq-Tool data resources, available under the
    following attribute names:

    * ``self.seqrepo_access``: :py:class:`SeqRepoAccess <cool_seq_tool.handlers.seqrepo_access.SeqRepoAccess>`
    * ``self.transcript_mappings``: :py:class:`TranscriptMappings <cool_seq_tool.sources.transcript_mappings.TranscriptMappings>`
    * ``self.mane_transcript_mappings``: :py:class:`ManeTranscriptMappings <cool_seq_tool.sources.mane_transcript_mappings.ManeTranscriptMappings>`
    * ``self.uta_db``: :py:class:`UtaDatabase <cool_seq_tool.sources.uta_database.UtaDatabase>`
    * ``self.alignment_mapper``: :py:class:`AlignmentMapper <cool_seq_tool.mappers.alignment.AlignmentMapper>`
    * ``self.mane_transcript``: :py:class:`ManeTranscript <cool_seq_tool.mappers.mane_transcript.ManeTranscript>`
    * ``self.ex_g_coords_mapper``: :py:class:`ExonGenomicCoordsMapper <cool_seq_tool.mappers.exon_genomic_coords.ExonGenomicCoordsMapper>`

    Initialization with default resource locations is straightforward:

    .. code-block:: pycon

       >>> from cool_seq_tool.app import CoolSeqTool
       >>> cst = CoolSeqTool()

    See the :ref:`configuration <configuration>` section for more information.
    """

    def __init__(
        self,
        transcript_file_path: Path = TRANSCRIPT_MAPPINGS_PATH,
        lrg_refseqgene_path: Path | None = None,
        mane_data_path: Path | None = None,
        db_url: str = UTA_DB_URL,
        sr: Optional[SeqRepo] = None,
    ) -> None:
        """Initialize CoolSeqTool class

        :param transcript_file_path: The path to ``transcript_mapping.tsv``
        :param lrg_refseqgene_path: The path to the LRG_RefSeqGene file
        :param mane_data_path: Path to RefSeq MANE summary data
        :param db_url: PostgreSQL connection URL
            Format: ``driver://user:password@host/database/schema``
        :param sr: SeqRepo instance. If this is not provided, will create a new instance
        """
        if not sr:
            sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
        self.seqrepo_access = SeqRepoAccess(sr)
        self.transcript_mappings = TranscriptMappings(
            transcript_file_path=transcript_file_path,
            lrg_refseqgene_path=lrg_refseqgene_path,
        )
        self.mane_transcript_mappings = ManeTranscriptMappings(
            mane_data_path=mane_data_path
        )
        self.uta_db = UtaDatabase(db_url=db_url)
        self.alignment_mapper = AlignmentMapper(
            self.seqrepo_access, self.transcript_mappings, self.uta_db
        )
        self.mane_transcript = ManeTranscript(
            self.seqrepo_access,
            self.transcript_mappings,
            self.mane_transcript_mappings,
            self.uta_db,
        )
        self.ex_g_coords_mapper = ExonGenomicCoordsMapper(
            self.seqrepo_access,
            self.uta_db,
            self.mane_transcript,
            self.mane_transcript_mappings,
        )
