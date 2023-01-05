"""The module for loading MANE Transcript mappings to genes."""
from typing import Dict, Optional, List

import pandas as pd

from cool_seq_tool import MANE_SUMMARY_PATH, logger


class MANETranscriptMappings:
    """The MANE Transcript mappings class."""

    def __init__(self, mane_data_path: str = MANE_SUMMARY_PATH) -> None:
        """Initialize the MANE Transcript mappings class.
        :param str mane_data_path: Path to RefSeq MANE summary data
        """
        self.mane_data_path = mane_data_path
        self.df = self._load_mane_transcript_data()

    def _load_mane_transcript_data(self) -> pd.core.frame.DataFrame:
        """Load RefSeq MANE data file into DataFrame.
        :return: DataFrame containing RefSeq MANE Transcript data
        """
        return pd.read_csv(self.mane_data_path, delimiter="\t")

    def get_gene_mane_data(self, gene_symbol: str) -> Optional[List[Dict]]:
        """Return MANE Transcript data for a gene.
        :param str gene_symbol: HGNC Gene Symbol
        :return: MANE Transcript data (Transcript accessions,
            gene, and location information)
        """
        data = self.df.loc[self.df["symbol"] == gene_symbol.upper()]

        if len(data) == 0:
            logger.warning(f"Unable to get MANE Transcript data for gene: "
                           f"{gene_symbol}")
            return None

        # Ordering: MANE Plus Clinical (If it exists), MANE Select
        data = data.sort_values("MANE_status")
        return data.to_dict("records")

    def get_mane_from_transcripts(self, transcripts: List[str]) -> List[Dict]:
        """Get mane transcripts from a list of transcripts

        :param List[str] transcripts: RefSeq transcripts on c. coordinate
        :return: MANE data
        """
        mane_rows = self.df["RefSeq_nuc"].isin(transcripts)
        result = self.df[mane_rows]
        if len(result) == 0:
            return []
        return result.to_dict("records")

    def get_mane_data_from_chr_pos(self, chromosome: str, start: int,
                                   end: int) -> List[Dict]:
        """Get MANE data given chromosome, start pos, end end pos. Assumes GRCh38.
        :param str chromosome: Chromosome. Case sensitive for X, Y chromosomes.
            Do not prefix. Valid examples: "1" ... "22", "X", "Y"
        :param int start: Start genomic position. Assumes residue coordinates.
        :param int end: End genomic position. Assumes residue coordinates.
        :return: List of MANE data. Will return sorted list:
            MANE Select then MANE Plus Clinical.
        """
        mane_rows = self.df[(start >= self.df["chr_start"].astype(int)) & (end <= self.df["chr_end"].astype(int)) & (self.df["GRCh38_chr"] == chromosome)]  # noqa: E501
        if len(mane_rows) == 0:
            return []
        mane_rows = mane_rows.sort_values("MANE_status", ascending=False)
        return mane_rows.to_dict("records")
