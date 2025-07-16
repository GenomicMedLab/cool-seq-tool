"""Module for getting feature (gene/exon) overlap"""

import re
from pathlib import Path

import polars as pl
from ga4gh.core import ga4gh_identify
from ga4gh.vrs.models import SequenceLocation, SequenceReference

from cool_seq_tool.handlers import SeqRepoAccess
from cool_seq_tool.resources.data_files import DataFile, get_data_file
from cool_seq_tool.schemas import Assembly, CdsOverlap, CoordinateType

# Pattern for chromosome
CHR_PATTERN = r"X|Y|([1-9]|1[0-9]|2[0-2])"


class FeatureOverlapError(Exception):
    """Custom exception for the Feature Overlap class"""


class FeatureOverlap:
    """The class for getting feature overlap"""

    def __init__(
        self,
        seqrepo_access: SeqRepoAccess,
        mane_refseq_genomic_path: Path | None = None,
        from_local: bool = False,
    ) -> None:
        """Initialize the FeatureOverlap class. Will load RefSeq data and store as df.

        :param seqrepo_access: Client for accessing SeqRepo data
        :param mane_refseq_genomic_path: Path to MANE RefSeq Genomic GFF data
        :param from_local: if ``True``, don't check for or acquire latest version --
            just provide most recent locally available file, if possible, and raise
            error otherwise
        """
        if not mane_refseq_genomic_path:
            mane_refseq_genomic_path = get_data_file(
                DataFile.MANE_REFSEQ_GENOMIC, from_local
            )
        self.seqrepo_access = seqrepo_access
        self.mane_refseq_genomic_path = mane_refseq_genomic_path
        self.df = self._load_mane_refseq_gff_data()

    def _load_mane_refseq_gff_data(self) -> pl.DataFrame:
        """Load MANE RefSeq GFF data file into DataFrame.

        :return: DataFrame containing MANE RefSeq Genomic GFF data for CDS. Columns
            include `type`, `chromosome` (chromosome without 'chr' prefix), `cds_start`,
            `cds_stop`, `info_name` (name of record), and `gene`. `cds_start` and
            `cds_stop` use inter-residue coordinates.
        """
        df = pl.read_csv(
            self.mane_refseq_genomic_path,
            separator="\t",
            has_header=False,
            skip_rows=9,
            columns=[0, 2, 3, 4, 8],
        )
        df.columns = ["chromosome", "type", "cds_start", "cds_stop", "info"]

        # Restrict to only feature of interest: CDS (which has gene info)
        df = df.filter(pl.col("type") == "CDS")

        # Get name from the info field
        # Get gene from the info field
        # Get chromosome names without prefix and without suffix for alternate transcripts
        # Convert start and stop to ints
        # Convert to inter-residue coordinates
        # Only return certain columns
        return df.with_columns(
            (pl.col("info").str.extract(r"Name=([^;]+)", 1).alias("info_name")),
            (pl.col("info").str.extract(r"gene=([^;]+)", 1).alias("gene")),
            (pl.col("chromosome").str.extract(r"^chr?([^_]+)", 1).alias("chromosome")),
            (pl.col("cds_start").cast(pl.Int64) - 1).alias("cds_start"),
            (pl.col("cds_stop").cast(pl.Int64).alias("cds_stop")),
        ).select(
            [
                pl.col("type"),
                pl.col("chromosome"),
                pl.col("cds_start"),
                pl.col("cds_stop"),
                pl.col("info_name"),
                pl.col("gene"),
            ]
        )

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
            error_msg = (
                f"Unable to find {Assembly.GRCH38.value} aliases for: {identifier}"
            )
            raise FeatureOverlapError(error_msg)

        assembly_chr_pattern = (
            rf"^{Assembly.GRCH38.value}:(?P<chromosome>{CHR_PATTERN})$"
        )
        for a in aliases:
            chr_match = re.match(assembly_chr_pattern, a)
            if chr_match:
                break

        if not chr_match:
            error_msg = (
                f"Unable to find {Assembly.GRCH38.value} chromosome for: {identifier}"
            )
            raise FeatureOverlapError(error_msg)

        chr_groupdict = chr_match.groupdict()
        return chr_groupdict["chromosome"]

    def get_grch38_mane_gene_cds_overlap(
        self,
        start: int,
        end: int,
        chromosome: str | None = None,
        identifier: str | None = None,
        coordinate_type: CoordinateType = CoordinateType.RESIDUE,
    ) -> dict[str, list[CdsOverlap]] | None:
        """Given GRCh38 genomic data, find the overlapping MANE features (gene and cds).
        The genomic data is specified as a sequence location by `chromosome`, `start`,
        `end`. All CDS regions with which the input sequence location has nonzero base
        pair overlap will be returned.

        :param start: GRCh38 start position
        :param end: GRCh38 end position
        :param chromosome: Chromosome. 1..22, X, or Y. If not provided, must provide
            `identifier`. If both `chromosome` and `identifier` are provided,
            `chromosome` will be used.
        :param identifier: Genomic identifier on GRCh38 assembly. If not provided, must
            provide `chromosome`. If both `chromosome` and `identifier` are provided,
            `chromosome` will be used.
        :param coordinate_type: Coordinate type for ``start`` and ``end``
        :raise FeatureOverlapError: If missing required fields or unable to find
            associated ga4gh identifier
        :return: MANE feature (gene/cds) overlap data represented as a dict. The
            dictionary will be keyed by genes which overlap the input sequence location.
            Each gene contains a list of the overlapping CDS regions with the beginning
            and end of the input sequence location's overlap with each
        """
        ga4gh_seq_id = None
        if chromosome:
            if not re.match(f"^{CHR_PATTERN}$", chromosome):
                error_msg = "`chromosome` must be 1, ..., 22, X, or Y"
                raise FeatureOverlapError(error_msg)
        else:
            if identifier:
                chromosome = self._get_chr_from_alt_ac(identifier)
                if identifier.startswith("ga4gh:SQ."):
                    ga4gh_seq_id = identifier
            else:
                error_msg = "Must provide either `chromosome` or `identifier`"
                raise FeatureOverlapError(error_msg)

        # Convert residue to inter-residue
        if coordinate_type == CoordinateType.RESIDUE:
            start -= 1

        # Get feature dataframe (df uses inter-residue)
        feature_df = self.df.filter(
            (pl.col("chromosome") == chromosome)
            & (pl.col("cds_start") <= end)
            & (pl.col("cds_stop") >= start)
        )

        if feature_df.is_empty():
            return None

        # Add overlap columns
        feature_df = feature_df.with_columns(
            [
                pl.when(pl.col("cds_start") < start)
                .then(start)
                .otherwise(pl.col("cds_start"))
                .alias("overlap_start"),
                pl.when(pl.col("cds_stop") > end)
                .then(end)
                .otherwise(pl.col("cds_stop"))
                .alias("overlap_stop"),
            ]
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

            if not ga4gh_aliases:
                error_msg = f"Unable to find ga4gh identifier for: {grch38_chr}"
                raise FeatureOverlapError(error_msg)

            ga4gh_seq_id = ga4gh_aliases[0]

        def _get_seq_loc(start_pos: int, stop_pos: int, refget_ac: str) -> dict:
            """Get VRS Sequence Location represented as a dict

            :param start_pos: Start position
            :param stop_pos: Stop position
            :param refget_ac: Refget Accession (SQ.)
            :return: VRS Sequence Location represented as dictionary with the ga4gh ID
                included
            """
            _sl = SequenceLocation(
                sequenceReference=SequenceReference(
                    refgetAccession=refget_ac,
                ),
                start=start_pos,
                end=stop_pos,
            )
            ga4gh_identify(_sl)
            return _sl.model_dump(exclude_none=True)

        resp = {}
        refget_ac = ga4gh_seq_id.split("ga4gh:")[-1]
        for gene, group in feature_df.group_by("gene"):
            gene = gene[0]
            _gene_overlap_data = [
                CdsOverlap(
                    cds=_get_seq_loc(
                        cds_row["cds_start"], cds_row["cds_stop"], refget_ac
                    ),
                    overlap=_get_seq_loc(
                        cds_row["overlap_start"], cds_row["overlap_stop"], refget_ac
                    ),
                ).model_dump(by_alias=True, exclude_none=True)
                for cds_row in group.iter_rows(named=True)
            ]
            resp[gene] = _gene_overlap_data

        return resp
