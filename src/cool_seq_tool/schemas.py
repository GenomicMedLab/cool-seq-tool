"""Defines attribute constants, useful object structures, and API response schemas."""

import datetime
from enum import Enum, IntEnum
from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    StrictInt,
    StrictStr,
    model_validator,
)

from cool_seq_tool import __version__

_now = str(datetime.datetime.now(tz=datetime.timezone.utc))


class AnnotationLayer(str, Enum):
    """Create enum for supported annotation layers"""

    PROTEIN = "p"
    CDNA = "c"
    GENOMIC = "g"


class Strand(IntEnum):
    """Create enum for positive and negative strand"""

    POSITIVE = 1
    NEGATIVE = -1


class Assembly(str, Enum):
    """Define supported genomic assemblies. Must be defined in ascending order"""

    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"

    @classmethod
    def values(cls) -> list[str]:
        """Return list of values in enum (ascending assembly order)"""
        return [item.value for item in cls]


class TranscriptPriority(str, Enum):
    """Create Enum for Transcript Priority labels"""

    MANE_SELECT = "mane_select"
    MANE_PLUS_CLINICAL = "mane_plus_clinical"
    LONGEST_COMPATIBLE_REMAINING = "longest_compatible_remaining"
    GRCH38 = "grch38"


class ResidueMode(str, Enum):
    """Create Enum for residue modes.

    We typically prefer to operate in inter-residue coordinates, but users should be
    careful to define the coordinate mode of their data when calling ``cool-seq-tool``
    functions.

                      |   | C |   | T |   | G |   |
    ZERO              |   | 0 |   | 1 |   | 2 |   |
    RESIDUE           |   | 1 |   | 2 |   | 3 |   |
    INTER_RESIDUE     | 0 |   | 1 |   | 2 |   | 3 |

    .. tabularcolumns:: |L|C|C|C|C|C|C|C|
    .. list-table::
       :header-rows: 1

       * -
         -
         - C
         -
         - T
         -
         - G
         -
       * - ``ZERO``
         -
         - 0
         -
         - 1
         -
         - 2
         -
       * - ``RESIDUE``
         -
         - 1
         -
         - 2
         -
         - 3
         -
       * - ``INTER_RESIDUE``
         - 0
         -
         - 1
         -
         - 2
         -
         - 3


    See "Conventions that promote reliable data sharing" and figure 3 within the
    `Variation Representation Schema (VRS) paper <https://www.ncbi.nlm.nih.gov/pmc/articles/pmid/35311178/>`_ for further discussion.
    """

    ZERO = "zero"
    RESIDUE = "residue"
    INTER_RESIDUE = "inter-residue"


class BaseModelForbidExtra(BaseModel, extra="forbid"):
    """Base Pydantic model class with extra values forbidden."""


class GenesGenomicAcs(BaseModelForbidExtra):
    """Represent HGNC gene symbols and genomic accessions"""

    genes: set[str]
    alt_acs: set[str]


class GenomicTxData(BaseModelForbidExtra):
    """Represent aligned genomic/transcript exon data"""

    gene: str
    strand: Strand
    tx_pos_range: tuple[int, int]
    alt_pos_range: tuple[int, int]
    alt_aln_method: str
    tx_exon_id: int
    alt_exon_id: int


class GenomicTxMetadata(GenomicTxData):
    """Store relevant metadata for genomic and transcript accessions"""

    tx_ac: str
    alt_ac: str
    coding_start_site: int = 0
    coding_end_site: int = 0
    alt_pos_change_range: tuple[int, int]
    pos_change: tuple[int, int] | None


class ManeGeneData(BaseModel, extra="forbid"):
    """Define minimal object model for representing a MANE gene"""

    ncbi_gene_id: StrictInt
    hgnc_id: StrictInt | None
    symbol: StrictStr


class TranscriptExonData(BaseModelForbidExtra):
    """Model containing transcript exon data."""

    transcript: StrictStr
    pos: StrictInt
    exon: StrictInt
    exon_offset: StrictInt = 0
    gene: StrictStr
    chr: StrictStr
    strand: Strand

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "chr": "NC_000001.11",
                "gene": "TPM3",
                "pos": 154192135,
                "exon": 1,
                "exon_offset": 0,
                "transcript": "NM_152263.3",
                "strand": Strand.NEGATIVE,
            }
        }
    )


class GenomicData(BaseModelForbidExtra):
    """Model containing genomic and transcript exon data."""

    gene: StrictStr
    chr: StrictStr
    start: StrictInt | None = None  # Genomic start position
    end: StrictInt | None = None  # Genomic end position
    exon_start: StrictInt | None = None
    exon_start_offset: StrictInt | None = 0
    exon_end: StrictInt | None = None
    exon_end_offset: StrictInt | None = 0
    transcript: StrictStr
    strand: Strand

    @model_validator(mode="after")
    def check_start_end(cls, values):
        """Check that at least one of {``start``, ``end``} is set.
        Check that at least one of {``exon_start``, ``exon_end``} is set.
        If not set, set corresponding offset to ``None``
        """
        start = values.start
        end = values.end
        if not start and not end:
            msg = "Missing values for `start` or `end`"
            raise ValueError(msg)

        if start:
            if not values.exon_start:
                msg = "Missing value `exon_start`"
                raise ValueError(msg)
        else:
            values.exon_start_offset = None

        if end:
            if not values.exon_end:
                msg = "Missing value `exon_end`"
                raise ValueError(msg)
        else:
            values.exon_end_offset = None
        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "gene": "TPM3",
                "chr": "NC_000001.11",
                "start": 154192135,
                "end": None,
                "exon_start": 1,
                "exon_end": None,
                "exon_start_offset": 0,
                "exon_end_offset": None,
                "transcript": "NM_152263.3",
                "strand": Strand.NEGATIVE,
            }
        }
    )


class ServiceMeta(BaseModelForbidExtra):
    """Metadata for cool_seq_tool service"""

    name: Literal["cool_seq_tool"] = "cool_seq_tool"
    version: StrictStr
    response_datetime: datetime.datetime
    url: Literal["https://github.com/GenomicMedLab/cool-seq-tool"] = (
        "https://github.com/GenomicMedLab/cool-seq-tool"
    )

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "name": "cool_seq_tool",
                "version": __version__,
                "response_datetime": _now,
                "url": "https://github.com/GenomicMedLab/cool-seq-tool",
            }
        }
    )


class TranscriptExonDataResponse(BaseModelForbidExtra):
    """Response model for Transcript Exon Data"""

    transcript_exon_data: TranscriptExonData | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "transcript_exon_data": {
                    "chr": "NC_000001.11",
                    "gene": "TPM3",
                    "pos": 154192135,
                    "exon": 1,
                    "exon_offset": 0,
                    "transcript": "NM_152263.3",
                    "strand": Strand.NEGATIVE,
                },
                "warnings": [],
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": _now,
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool",
                },
            }
        }
    )


class GenomicDataResponse(BaseModelForbidExtra):
    """Response model for Genomic Data"""

    genomic_data: GenomicData | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "genomic_data": {
                    "gene": "TPM3",
                    "chr": "NC_000001.11",
                    "start": 154192135,
                    "end": None,
                    "exon_start": 1,
                    "exon_end": None,
                    "exon_start_offset": 0,
                    "exon_end_offset": None,
                    "transcript": "NM_152263.3",
                    "strand": Strand.NEGATIVE,
                },
                "warnings": [],
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": _now,
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool",
                },
            }
        }
    )
