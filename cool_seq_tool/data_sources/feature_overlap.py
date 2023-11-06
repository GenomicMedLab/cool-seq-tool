"""Module for getting feature (gene/exon) overlap"""
import re
from pathlib import Path
from typing import Dict, Optional

import pandas as pd
from ga4gh.core import ga4gh_identify
from ga4gh.vrs import models

from cool_seq_tool.data_sources import SeqRepoAccess
from cool_seq_tool.paths import MANE_REFSEQ_GFF_PATH
from cool_seq_tool.schemas import Assembly, CdsOverlap, ResidueMode


# Pattern for chromosome
CHR_PATTERN = r"X|Y|([1-9]|1[0-9]|2[0-2])"


class FeatureOverlapError(Exception):
    """Custom exception for the Feature Overlap class"""


class FeatureOverlap:
    """The class for getting feature overlap"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        mane_refseq_gff_path: Path = MANE_REFSEQ_GFF_PATH,
    ) -> None:
        """Initialize the FeatureOverlap class. Will load RefSeq data and store as df.

        :param seqrepo_access: Client for accessing SeqRepo data
        :param mane_refseq_gff_path: Path to the MANE RefSeq GFF file
        """
        self.seqrepo_access = seqrepo_access
        self.mane_refseq_gff_path = mane_refseq_gff_path
        self.df = self._load_mane_refseq_gff_data()

    def _load_mane_refseq_gff_data(self) -> pd.core.frame.DataFrame:
        """Load MANE RefSeq GFF data file into DataFrame.

        :return: DataFrame containing MANE RefSeq GFF data for CDS. Columns include
            `type`, `chromosome` (chromosome without 'chr' prefix), `cds_start`,
            `cds_stop`, `info_name` (name of record), and `gene`. `cds_start` and
            `cds_stop` use inter-residue coordinates.
        """
        df = pd.read_csv(
            self.mane_refseq_gff_path,
            sep="\t",
            header=None,
            skiprows=9,
            usecols=[0, 2, 3, 4, 8],
        )
        df.columns = ["chromosome", "type", "cds_start", "cds_stop", "info"]

        # Restrict to only feature of interest: CDS (which has gene info)
        df = df[df["type"] == "CDS"].copy()

        # Get name from the info field
        df["info_name"] = df["info"].apply(
            lambda info: re.findall("Name=([^;]+)", info)[0]
        )

        # Get gene from the info field
        df["gene"] = df["info"].apply(lambda info: re.findall("gene=([^;]+)", info)[0])

        # Get chromosome names without prefix and without suffix for alternate
        # transcripts
        df["chromosome"] = df["chromosome"].apply(
            lambda chromosome: chromosome.strip("chr").split("_")[0]
        )
        df["chromosome"] = df["chromosome"].astype(str)

        # Convert start and stop to ints
        df["cds_start"] = df["cds_start"].astype(int)
        df["cds_stop"] = df["cds_stop"].astype(int)

        # Convert to inter-residue coordinates
        df["cds_start"] -= 1

        # Only retain certain columns
        df = df[["type", "chromosome", "cds_start", "cds_stop", "info_name", "gene"]]

        return df

    def _get_chr_from_alt_ac(self, identifier: str) -> str:
        """Get chromosome given genomic identifier

        :param identifier: Genomic identifier on GRCh38 assembly
        :raises FeatureOverlapError: If unable to find associated GRCh38 chromosome
        :return: Chromosome. 1..22, X, Y. No 'chr' prefix.
        """
        aliases, error_msg = self.seqrepo_access.translate_identifier(
            identifier, Assembly.GRCH38.value
        )

        if error_msg:
            raise FeatureOverlapError(str(error_msg))

        if not aliases:
            raise FeatureOverlapError(
                f"Unable to find {Assembly.GRCH38.value} aliases for: {identifier}"
            )

        assembly_chr_pattern = rf"^{Assembly.GRCH38.value}:(?P<chromosome>{CHR_PATTERN})$"  # noqa: E501
        for a in aliases:
            chr_match = re.match(assembly_chr_pattern, a)
            if chr_match:
                break

        if not chr_match:
            raise FeatureOverlapError(
                f"Unable to find {Assembly.GRCH38.value} chromosome for: {identifier}"
            )

        chr_groupdict = chr_match.groupdict()
        return chr_groupdict["chromosome"]

    def get_grch38_mane_gene_cds_overlap(
        self,
        start: int,
        end: int,
        chromosome: Optional[str] = None,
        identifier: Optional[str] = None,
        residue_mode: ResidueMode = ResidueMode.RESIDUE,
    ) -> Optional[Dict[str, CdsOverlap]]:
        """Given GRCh38 genomic data, find the overlapping MANE features (gene and cds)

        :param start: GRCh38 start position
        :param end: GRCh38 end position
        :param chromosome: Chromosome. 1..22, X, or Y. If not provided, must provide
            `identifier`. If both `chromosome` and `identifier` are provided,
            `chromosome` will be used.
        :param identifier: Genomic identifier on GRCh38 assembly. If not provided, must
            provide `chromosome`. If both `chromosome` and `identifier` are provided,
            `chromosome` will be used.
        :param residue_mode: Residue mode for `start` and `end`
        :raise FeatureOverlapError: If missing required fields or unable to find
            associated ga4gh identifier
        :return: MANE feature (gene/cds) overlap data represented as a dict
            {
                gene: {
                    'cds': VRS Sequence Location
                    'overlap': VRS Sequence Location
                }
            }
        """
        ga4gh_seq_id = None
        if chromosome:
            if not re.match(f"^{CHR_PATTERN}$", chromosome):
                raise FeatureOverlapError("`chromosome` must be 1, ..., 22, X, or Y")
        else:
            if identifier:
                chromosome = self._get_chr_from_alt_ac(identifier)
                if identifier.startswith("ga4gh:SQ."):
                    ga4gh_seq_id = identifier
            else:
                raise FeatureOverlapError(
                    "Must provide either `chromosome` or `identifier`"
                )

        # Convert residue to inter-residue
        if residue_mode == ResidueMode.RESIDUE:
            start -= 1

        # Get feature dataframe (df uses inter-residue)
        feature_df = self.df[
            (self.df["chromosome"] == chromosome)
            & (self.df["cds_start"] <= end)  # noqa: W503
            & (self.df["cds_stop"] >= start)  # noqa: W503
        ].copy()

        if feature_df.empty:
            return None

        # Add overlap columns
        feature_df["overlap_start"] = feature_df["cds_start"].apply(
            lambda x: x if start <= x else start
        )
        feature_df["overlap_stop"] = feature_df["cds_stop"].apply(
            lambda x: end if end <= x else x
        )

        # Get ga4gh identifier for chromosome
        if not ga4gh_seq_id:
            grch38_chr = f"{Assembly.GRCH38.value}:{chromosome}"
            ga4gh_aliases, error_msg = self.seqrepo_access.translate_identifier(
                grch38_chr, "ga4gh"
            )

            # Errors should never happen but catching just in case
            if error_msg:
                raise FeatureOverlapError(str(error_msg))
            elif not ga4gh_aliases:
                raise FeatureOverlapError(
                    f"Unable to find ga4gh identifier for: {grch38_chr}"
                )

            ga4gh_seq_id = ga4gh_aliases[0]

        def _get_seq_loc(start_pos: int, stop_pos: int, ga4gh_sequence_id: str) -> Dict:
            """Get VRS Sequence Location represented as a dict

            :param start_pos: Start position
            :param stop_pos: Stop position
            :param ga4gh_sequence_id: ga4gh sequence identifier
            :return: VRS Sequence Location represented as dictionary with the ga4gh ID
                included
            """
            _sl = models.SequenceLocation(
                sequence_id=ga4gh_sequence_id,
                interval=models.SequenceInterval(
                    start=models.Number(value=start_pos),
                    end=models.Number(value=stop_pos),
                ),
            )
            _sl._id = ga4gh_identify(_sl)
            return _sl.as_dict()

        resp = {}
        for gene, group in feature_df.groupby("gene"):
            _gene_overlap_data = []

            for cds_row in group.itertuples():
                _gene_overlap_data.append(
                    CdsOverlap(
                        cds=_get_seq_loc(
                            cds_row.cds_start, cds_row.cds_stop, ga4gh_seq_id
                        ),
                        overlap=_get_seq_loc(
                            cds_row.overlap_start, cds_row.overlap_stop, ga4gh_seq_id
                        ),
                    ).dict(by_alias=True)
                )
            resp[gene] = _gene_overlap_data

        return resp
