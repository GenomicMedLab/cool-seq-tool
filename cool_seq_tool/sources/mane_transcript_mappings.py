"""The module for loading MANE Transcript mappings to genes."""
import logging
from pathlib import Path
from typing import Dict, List

import polars as pl

from cool_seq_tool.paths import MANE_SUMMARY_PATH

logger = logging.getLogger(__name__)


class MANETranscriptMappings:
    """The MANE Transcript mappings class."""

    def __init__(self, mane_data_path: Path = MANE_SUMMARY_PATH) -> None:
        """Initialize the MANE Transcript mappings class.
        :param Path mane_data_path: Path to RefSeq MANE summary data
        """
        self.mane_data_path = mane_data_path
        self.df = self._load_mane_transcript_data()

    def _load_mane_transcript_data(self) -> pl.DataFrame:
        """Load RefSeq MANE data file into DataFrame.
        :return: DataFrame containing RefSeq MANE Transcript data
        """
        return pl.read_csv(self.mane_data_path, separator="\t")

    def get_gene_mane_data(self, gene_symbol: str) -> List[Dict]:
        """Return MANE Transcript data for a gene.
        :param str gene_symbol: HGNC Gene Symbol
        :return: List of MANE Transcript data (Transcript accessions,
            gene, and location information). Sorted list: MANE Select and then MANE Plus
            Clinical
        """
        data = self.df.filter(pl.col("symbol") == gene_symbol.upper())

        if len(data) == 0:
            logger.warning(
                f"Unable to get MANE Transcript data for gene: " f"{gene_symbol}"
            )
            return []

        data = data.sort(by="MANE_status", descending=True)
        return data.to_dicts()

    def get_mane_from_transcripts(self, transcripts: List[str]) -> List[Dict]:
        """Get mane transcripts from a list of transcripts

        :param List[str] transcripts: RefSeq transcripts on c. coordinate
        :return: MANE data
        """
        mane_rows = self.df.filter(pl.col("RefSeq_nuc").is_in(transcripts))
        if len(mane_rows) == 0:
            return []
        return mane_rows.to_dicts()

    def get_mane_data_from_chr_pos(
        self, alt_ac: str, start: int, end: int
    ) -> List[Dict]:
        """Get MANE data given chromosome, start pos, end end pos. Assumes GRCh38.
        :param str alt_ac: NC Accession
        :param int start: Start genomic position. Assumes residue coordinates.
        :param int end: End genomic position. Assumes residue coordinates.
        :return: List of MANE data. Will return sorted list:
            MANE Select then MANE Plus Clinical.
        """
        mane_rows = self.df.filter(
            (start >= pl.col("chr_start"))
            & (end <= pl.col("chr_end"))
            & (pl.col("GRCh38_chr") == alt_ac)
        )
        if len(mane_rows) == 0:
            return []

        mane_rows = mane_rows.sort(by="MANE_status", descending=True)
        return mane_rows.to_dicts()
