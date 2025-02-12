"""Provide fast tabular access to MANE summary file. Enables retrieval of associated
MANE transcripts for gene symbols, genomic positions, or transcript accessions.
"""

import logging
from pathlib import Path

import polars as pl

from cool_seq_tool.resources.data_files import DataFile, get_data_file
from cool_seq_tool.schemas import ManeGeneData

_logger = logging.getLogger(__name__)


class ManeTranscriptMappings:
    """Provide fast tabular access to MANE summary file.

    By default, acquires data from `NCBI FTP server <ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/>`_
    if unavailable locally. The local data location can be passed as an argument or
    given under the environment variable ``MANE_SUMMARY_PATH``.

    See the `NCBI MANE page <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ for more information.
    """

    def __init__(
        self, mane_data_path: Path | None = None, from_local: bool = False
    ) -> None:
        """Initialize the MANE Transcript mappings class.

        :param mane_data_path: Path to RefSeq MANE summary data
        :param from_local: if ``True``, don't check for or acquire latest version --
            just provide most recent locally available file, if possible, and raise
            error otherwise
        """
        if not mane_data_path:
            mane_data_path = get_data_file(DataFile.MANE_SUMMARY, from_local)
        self.mane_data_path = mane_data_path
        self.df = self._load_mane_transcript_data()

    def _load_mane_transcript_data(self) -> pl.DataFrame:
        """Load RefSeq MANE data file into DataFrame.

        :return: DataFrame containing RefSeq MANE Transcript data
        """
        return pl.read_csv(self.mane_data_path, separator="\t")

    def get_gene_mane_data(self, gene_symbol: str) -> list[dict]:
        """Return MANE Transcript data for a gene.

        >>> from cool_seq_tool.sources import ManeTranscriptMappings
        >>> m = ManeTranscriptMappings()
        >>> braf_mane = m.get_gene_mane_data("BRAF")
        >>> braf_mane[0]["RefSeq_nuc"], braf_mane[0]["MANE_status"]
        ('NM_004333.6', 'MANE Select')
        >>> braf_mane[1]["RefSeq_nuc"], braf_mane[1]["MANE_status"]
        ('NM_001374258.1', 'MANE Plus Clinical')

        :param str gene_symbol: HGNC Gene Symbol
        :return: List of MANE Transcript data (Transcript accessions, gene, and
            location information). The list is sorted so that a MANE Select entry comes
            first, followed by a MANE Plus Clinical entry, if available.
        """
        data = self.df.filter(
            pl.col("symbol").str.to_uppercase() == gene_symbol.upper()
        )

        if len(data) == 0:
            _logger.warning(
                "Unable to get MANE Transcript data for gene: %s", gene_symbol
            )
            return []

        data = data.sort(by="MANE_status", descending=True)
        return data.to_dicts()

    def get_mane_from_transcripts(self, transcripts: list[str]) -> list[dict]:
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
    ) -> list[dict]:
        """Get MANE data given a GRCh38 genomic position.

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

    def get_genomic_mane_genes(
        self, ac: str, start: int, end: int
    ) -> list[ManeGeneData]:
        """Get MANE gene(s) for genomic location

        :param ac: RefSeq genomic accession
        :param start: Genomic start position. Assumes residue coordinates.
        :param end: Genomic end position. Assumes residue coordinates.
        :return: Unique MANE gene(s) found for a genomic location
        """
        # Only interested in rows where genomic location lives
        mane_rows = self.df.filter(
            (start >= pl.col("chr_start"))
            & (end <= pl.col("chr_end"))
            & (pl.col("GRCh38_chr") == ac)
        )

        if mane_rows.is_empty():
            return []

        # Group rows by NCBI ID, transform values to representation we want, MANE status
        # will be converted to list with DESC order
        mane_rows = mane_rows.group_by("#NCBI_GeneID").agg(
            [
                pl.col("#NCBI_GeneID")
                .first()
                .str.split_exact(":", 1)
                .struct.field("field_1")
                .cast(pl.Int32)
                .alias("ncbi_gene_id"),
                pl.col("HGNC_ID")
                .first()
                .str.split_exact(":", 1)
                .struct.field("field_1")
                .cast(pl.Int32)
                .alias("hgnc_id"),
                pl.col("MANE_status")
                .unique()
                .str.to_lowercase()
                .str.replace_all(" ", "_")
                .alias("status")
                .sort(descending=True),
                pl.col("symbol").first(),
            ]
        )

        # Sort final rows based on MANE status
        # First by length (which means gene has both select and plus clinical)
        # Then by DESC order
        # Then by NCBI ID ASC order
        mane_rows = (
            mane_rows.with_columns(
                [
                    pl.col("status").list.len().alias("status_count"),
                    pl.col("status").list.join("_").alias("status_str"),
                    pl.col("ncbi_gene_id"),
                ]
            )
            .sort(
                ["status_count", "status_str", "ncbi_gene_id"],
                descending=[True, True, False],
            )
            .drop(["status_count", "status_str", "#NCBI_GeneID"])
        )
        return [ManeGeneData(**mane_gene) for mane_gene in mane_rows.to_dicts()]
