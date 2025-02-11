"""Defines attribute constants, useful object structures, and API response schemas."""

import datetime
from enum import Enum, IntEnum
from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    StrictInt,
    StrictStr,
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


class ManeStatus(str, Enum):
    """Define constraints for mane status"""

    SELECT = "mane_select"
    PLUS_CLINICAL = "mane_plus_clinical"


class TranscriptPriority(str, Enum):
    """Create Enum for Transcript Priority labels"""

    MANE_SELECT = ManeStatus.SELECT.value
    MANE_PLUS_CLINICAL = ManeStatus.PLUS_CLINICAL.value
    LONGEST_COMPATIBLE_REMAINING = "longest_compatible_remaining"
    GRCH38 = "grch38"


class CoordinateType(str, Enum):
    """Create Enum for coordinate types.

    It is preferred to operate in inter-residue coordinates, but users should be
    careful to define the coordinate mode of their data when calling ``cool-seq-tool``
    functions.

    ``RESIDUE`` means 1-indexed, residue coordinates and ``INTER_RESIDUE`` means
    0-indexed, inter-residue coordinates.

                      |   | C |   | T |   | G |   |
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

    RESIDUE = "residue"
    INTER_RESIDUE = "inter-residue"


class BaseModelForbidExtra(BaseModel, extra="forbid"):
    """Base Pydantic model class with extra values forbidden."""


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
    status: list[ManeStatus]


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
