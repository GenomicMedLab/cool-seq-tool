"""Provides core CoolSeqTool class, which non-redundantly initializes all Cool-Seq-Tool
data handler and mapping resources for straightforward access.
"""

import logging
from pathlib import Path

from biocommons.seqrepo import SeqRepo

from cool_seq_tool.handlers.seqrepo_access import SEQREPO_ROOT_DIR, SeqRepoAccess
from cool_seq_tool.mappers import (
    AlignmentMapper,
    ExonGenomicCoordsMapper,
    LiftOver,
    ManeTranscript,
)
from cool_seq_tool.sources.mane_transcript_mappings import ManeTranscriptMappings
from cool_seq_tool.sources.transcript_mappings import TranscriptMappings
from cool_seq_tool.sources.uta_database import UTA_DB_URL, UtaDatabase

_logger = logging.getLogger(__name__)


class CoolSeqTool:
    """Non-redundantly initialize all Cool-Seq-Tool data resources, available under the
    following attribute names:

    * ``self.seqrepo_access``: :py:class:`SeqRepoAccess <cool_seq_tool.handlers.seqrepo_access.SeqRepoAccess>`
    * ``self.transcript_mappings``: :py:class:`TranscriptMappings <cool_seq_tool.sources.transcript_mappings.TranscriptMappings>`
    * ``self.mane_transcript_mappings``: :py:class:`ManeTranscriptMappings <cool_seq_tool.sources.mane_transcript_mappings.ManeTranscriptMappings>`
    * ``self.uta_db``: :py:class:`UtaDatabase <cool_seq_tool.sources.uta_database.UtaDatabase>`
    * ``self.alignment_mapper``: :py:class:`AlignmentMapper <cool_seq_tool.mappers.alignment.AlignmentMapper>`
    * ``self.liftover``: :py:class:`LiftOver <cool_seq_tool.mappers.liftover.LiftOver>`
    * ``self.mane_transcript``: :py:class:`ManeTranscript <cool_seq_tool.mappers.mane_transcript.ManeTranscript>`
    * ``self.ex_g_coords_mapper``: :py:class:`ExonGenomicCoordsMapper <cool_seq_tool.mappers.exon_genomic_coords.ExonGenomicCoordsMapper>`
    """

    def __init__(
        self,
        transcript_file_path: Path | None = None,
        lrg_refseqgene_path: Path | None = None,
        mane_data_path: Path | None = None,
        db_url: str = UTA_DB_URL,
        sr: SeqRepo | None = None,
        force_local_files: bool = False,
    ) -> None:
        """Initialize CoolSeqTool class.

        Initialization with default resource locations is straightforward:

        >>> from cool_seq_tool import CoolSeqTool
        >>> cst = CoolSeqTool()

        By default, this will attempt to fetch the latest versions of static resources,
        which means brief FTP and HTTPS requests to NCBI servers upon initialization.
        To suppress this check and simply rely on the most recent locally-available
        data:

        >>> cst = CoolSeqTool(force_local_files=True)

        Note that this will raise a FileNotFoundError if no locally-available data exists.

        Paths to those files can also be explicitly passed to avoid checks as well:

        >>> from pathlib import Path
        >>> cst = CoolSeqTool(
        ...     lrg_refseqgene_path=Path("lrg_refseqgene_20240625.tsv"),
        ...     mane_data_path=Path("ncbi_mane_summary_1.3.txt"),
        ... )

        If not passed explicit arguments, these locations can also be set via
        environment variables. See the :ref:`configuration <configuration>` section of
        the docs for more information.

        :param transcript_file_path: The path to ``transcript_mapping.tsv``
        :param lrg_refseqgene_path: The path to the LRG_RefSeqGene file
        :param mane_data_path: Path to RefSeq MANE summary data
        :param db_url: PostgreSQL connection URL
            Format: ``driver://user:password@host/database/schema``
        :param sr: SeqRepo instance. If this is not provided, will create a new instance
        :param force_local_files: if ``True``, don't check for or try to acquire latest
            versions of static data files -- just use most recently available, if any
        """
        if not sr:
            sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
        self.seqrepo_access = SeqRepoAccess(sr)
        self.transcript_mappings = TranscriptMappings(
            transcript_file_path=transcript_file_path,
            lrg_refseqgene_path=lrg_refseqgene_path,
            from_local=force_local_files,
        )
        self.mane_transcript_mappings = ManeTranscriptMappings(
            mane_data_path=mane_data_path, from_local=force_local_files
        )
        self.uta_db = UtaDatabase(db_url=db_url)
        self.alignment_mapper = AlignmentMapper(
            self.seqrepo_access, self.transcript_mappings, self.uta_db
        )
        self.liftover = LiftOver()
        self.mane_transcript = ManeTranscript(
            self.seqrepo_access,
            self.transcript_mappings,
            self.mane_transcript_mappings,
            self.uta_db,
            self.liftover,
        )
        self.ex_g_coords_mapper = ExonGenomicCoordsMapper(
            self.seqrepo_access,
            self.uta_db,
            self.mane_transcript_mappings,
            self.liftover,
        )
