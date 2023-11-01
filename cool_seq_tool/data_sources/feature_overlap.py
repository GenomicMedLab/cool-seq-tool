"""Module for getting feature (gene/exon) overlap"""
from pathlib import Path
from typing import Optional, Dict
import re

import pandas as pd

from cool_seq_tool.data_sources import SeqRepoAccess
from cool_seq_tool.paths import MANE_REFSEQ_GFF_PATH


class FeatureOverlapError(Exception):
    """Custom exception for the Feature Overlap class"""


class FeatureOverlap:
    """The class for getting feature overlap"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        mane_refseq_gff_path: Path = MANE_REFSEQ_GFF_PATH
    ) -> None:
        """Initialize the FeatureOverlap class

        :param mane_refseq_gff_path: Path to the MANE RefSeq GFF file
        """
        self.seqrepo_access = seqrepo_access
        self.mane_refseq_gff_path = mane_refseq_gff_path
        self.df = self._load_mane_refseq_gff_data()

    def _load_mane_refseq_gff_data(self) -> pd.core.frame.DataFrame:
        """Load MANE RefSeq GFF data file into DataFrame.
        Does transformations on the data.

        :return: DataFrame containing MANE RefSeq GFF data
        """
        df = pd.read_csv(
            self.mane_refseq_gff_path,
            sep="\t",
            header=None,
            skiprows=9,
            usecols=[0, 1, 2, 3, 4, 8]
        )
        df.columns = ["chromosome", "source", "type", "exon_start", "exon_stop", "info"]

        # Restrict to only feature of interest: coding exons (which has gene info)
        df = df[df["type"].isin(["CDS"])].copy()

        # Get name from the info field
        df["info_name"] = df["info"].apply(
            lambda info: re.findall("Name=([^;]+)", info)[0]
        )

        # Get gene from the info field
        df["gene"] = df["info"].apply(
            lambda info: re.findall("gene=([^;]+)", info)[0]
        )

        # Get chromosome names without prefix and without suffix for alternate
        # transcripts
        df["chrom_normalized"] = df["chromosome"].apply(
            lambda chromosome: chromosome.strip("chr").split("_")[0]
        )

        # Convert start and stop to ints
        df["exon_start"] = df["exon_start"].astype(int)
        df["exon_stop"] = df["exon_stop"].astype(int)

        # Only retain certain columns
        df = df[
            [
                "type",
                "source",
                "chromosome",
                "chrom_normalized",
                "exon_start",
                "exon_stop",
                "info_name",
                "gene"
            ]
        ]

        return df

    def _get_chr_from_alt_ac(self, alt_ac: str) -> str:
        """Get chromosome given RefSeq genomic accession

        :param alt_ac: RefSeq genomic accession on GRCh38 assembly
        :return: Chromosome. 1..22, X, Y. No 'chr' prefix.
        """
        aliases, error_msg = self.seqrepo_access.translate_identifier(alt_ac, "GRCh38")

        if error_msg:
            raise FeatureOverlapError(str(error_msg))

        chr_pattern = r"^GRCh38:(?P<chromosome>X|Y|([1-9]|1[0-9]|2[0-2]))$"
        for a in aliases:
            chr_match = re.match(chr_pattern, a)
            if chr_match:
                break

        if not chr_match:
            raise FeatureOverlapError(f"Unable to find GRCh38 chromosome for: {alt_ac}")

        chr_groupdict = chr_match.groupdict()
        return chr_groupdict["chromosome"]

    def get_grch38_overlap(
        self,
        start: int,
        end: int,
        chromosome: Optional[str] = None,
        alt_ac: Optional[str] = None
    ) -> Optional[Dict]:
        """Get feature overlap for GRCh38 genomic data

        :param start: GRCh38 start position
        :param end: GRCh38 end position
        :param chromosome: Chromosome. 1..22, X, or Y. If not provided, must provide
            `alt_ac`
        :param alt_ac: RefSeq genomic accession on GRCh38 assembly. If not provided,
            must provide `chromosome`
        :raise FeatureOverlapError: If missing required fields
        :return: Feature overlap dictionary where the key is the gene name and the value
            is the list of CDS overlap (exon_start, exon_stop, overlap_start,
            overlap_stop)
        """
        if chromosome:
            if not re.match(r"^X|Y|([1-9]|1[0-9]|2[0-2])$", chromosome):
                raise FeatureOverlapError("`chromosome` must be 1, ..., 22, X, or Y")
        else:
            if alt_ac:
                chromosome = self._get_chr_from_alt_ac(alt_ac)
            else:
                raise FeatureOverlapError(
                    "Must provide either `chromosome` or `alt_ac`"
                )

        # Get feature dataframe
        feature_df = self.df[
            (self.df["chrom_normalized"] == chromosome) & (self.df["exon_start"] <= end) & (self.df["exon_stop"] >= start)  # noqa: E501
        ].copy()

        if feature_df.empty:
            return None

        # Add overlap columns
        feature_df["overlap_start"] = feature_df["exon_start"].apply(
            lambda x: x if start <= x else x + (start - x)
        )
        feature_df["overlap_stop"] = feature_df["exon_stop"].apply(
            lambda x: x - (x - end) if end <= x else x
        )

        feature_df = feature_df.sort_values(["exon_start"], ascending=True)

        return (
            feature_df.groupby(["gene"])[["info_name", "exon_start", "exon_stop", "overlap_start", "overlap_stop"]]  # noqa: E501
            .apply(lambda x: x.set_index("info_name").to_dict(orient="records"))
            .to_dict()
        )
