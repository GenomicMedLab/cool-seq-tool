"""Defines attribute constants, useful object structures, and API response schemas."""

import datetime
import re
from enum import Enum, IntEnum
from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    StrictInt,
    StrictStr,
    field_validator,
    model_validator,
)

from cool_seq_tool.version import __version__

_now = str(datetime.datetime.now(tz=datetime.timezone.utc))


class AnnotationLayer(str, Enum):
    """Create enum for supported annotation layers"""

    PROTEIN: Literal["p"] = "p"
    CDNA: Literal["c"] = "c"
    GENOMIC: Literal["g"] = "g"


class Strand(IntEnum):
    """Create enum for positive and negative strand"""

    POSITIVE = 1
    NEGATIVE = -1


class Assembly(str, Enum):
    """Create Enum for supported genomic assemblies"""

    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


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


class GenomicRequestBody(BaseModelForbidExtra):
    """Define constraints for genomic to transcript exon coordinates request body"""

    chromosome: StrictStr | StrictInt
    start: StrictInt | None = None
    end: StrictInt | None = None
    strand: Strand | None = None
    transcript: StrictStr | None = None
    gene: StrictStr | None = None
    residue_mode: ResidueMode = ResidueMode.RESIDUE

    @model_validator(mode="after")
    def check_start_and_end(cls, values):
        """Check that at least one of {``start``, ``end``} is set"""
        start, end = values.start, values.end
        if not start or end:
            msg = "Must provide either `start` or `end`"
            raise ValueError(msg)
        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "chromosome": "NC_000001.11",
                "start": 154192135,
                "end": None,
                "strand": Strand.NEGATIVE,
                "transcript": "NM_152263.3",
                "gene": "TPM3",
                "residue_mode": "residue",
            }
        }
    )


class TranscriptRequestBody(BaseModelForbidExtra):
    """Define constraints for transcript exon to genomic coordinates request body"""

    transcript: StrictStr
    gene: StrictStr | None = None
    exon_start: StrictInt | None = None
    exon_start_offset: StrictInt | None = 0
    exon_end: StrictInt | None = None
    exon_end_offset: StrictInt | None = 0

    @model_validator(mode="after")
    def check_exon_start_and_exon_end(cls, values):
        """Check that at least one of {``exon_start``, ``exon_end``} is set"""
        exon_start, exon_end = values.exon_start, values.exon_end
        if not exon_start or exon_end:
            msg = "Must provide either `exon_start` or `exon_end`"
            raise ValueError(msg)
        return values

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "gene": "TPM3",
                "transcript": "NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 1,
                "exon_end": None,
                "exon_end_offset": None,
            }
        }
    )


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

    @field_validator("version")
    def validate_version(cls, v):
        """Check version matches semantic versioning regex pattern.
        https://semver.org/#is-there-a-suggested-regular-expression-regex-to-check-a-semver-string
        """
        version_regex = r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"
        if not re.match(version_regex, v):
            msg = f"Invalid version {v}"
            raise ValueError(msg)
        return v

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


class MappedManeData(BaseModel):
    """Define mapped mane data fields"""

    gene: StrictStr
    refseq: StrictStr
    ensembl: StrictStr | None = None
    strand: Strand
    status: TranscriptPriority
    alt_ac: StrictStr
    assembly: Assembly

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "gene": "BRAF",
                "refseq": "NM_001374258.1",
                "ensembl": "ENST00000644969.2",
                "strand": Strand.NEGATIVE,
                "status": TranscriptPriority.MANE_PLUS_CLINICAL,
                "alt_ac": "NC_000007.13",
                "assembly": "GRCh37",
            }
        }
    )


class MappedManeDataService(BaseModelForbidExtra):
    """Service model response for mapped mane data"""

    mapped_mane_data: MappedManeData | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "mapped_mane_data": {
                    "gene": "BRAF",
                    "refseq": "NM_001374258.1",
                    "ensembl": "ENST00000644969.2",
                    "strand": Strand.NEGATIVE,
                    "status": TranscriptPriority.MANE_PLUS_CLINICAL,
                    "alt_ac": "NC_000007.13",
                    "assembly": "GRCh37",
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


class ManeData(BaseModel):
    """Define mane data fields"""

    gene: StrictStr | None = None
    refseq: StrictStr | None = None
    ensembl: StrictStr | None = None
    pos: tuple[int, int]
    strand: Strand
    status: TranscriptPriority

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "gene": "BRAF",
                "refseq": "NP_004324.2",
                "ensembl": "ENSP00000493543.1",
                "pos": (598, 598),
                "strand": Strand.NEGATIVE,
                "status": TranscriptPriority.MANE_SELECT,
            }
        }
    )


class ManeDataService(BaseModelForbidExtra):
    """Service model response for getting mane data"""

    mane_data: ManeData | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "mane_data": {
                    "gene": "BRAF",
                    "refseq": "NP_004324.2",
                    "ensembl": "ENSP00000493543.1",
                    "pos": (598, 598),
                    "strand": Strand.NEGATIVE,
                    "status": TranscriptPriority.MANE_SELECT,
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


# ALIGNMENT MAPPER SERVICE SCHEMAS


class CdnaRepresentation(BaseModelForbidExtra):
    """Model response for cDNA representation"""

    c_ac: StrictStr
    c_start_pos: StrictInt
    c_end_pos: StrictInt
    cds_start: StrictInt
    residue_mode: Literal[ResidueMode.INTER_RESIDUE] = ResidueMode.INTER_RESIDUE.value

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "c_ac": "NM_004333.6",
                "c_start_pos": 1797,
                "c_end_pos": 1800,
                "cds_start": 226,
                "residue_mode": ResidueMode.INTER_RESIDUE,
            }
        }
    )


class ToCdnaService(BaseModelForbidExtra):
    """Service model response for protein -> cDNA"""

    c_data: CdnaRepresentation | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "c_data": {
                    "c_ac": "NM_004333.6",
                    "c_start_pos": 1797,
                    "c_end_pos": 1800,
                    "cds_start": 226,
                    "residue_mode": ResidueMode.INTER_RESIDUE,
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


class GenomicRepresentation(BaseModelForbidExtra):
    """Model response for genomic representation"""

    g_ac: str
    g_start_pos: int
    g_end_pos: int
    residue_mode: Literal[ResidueMode.INTER_RESIDUE] = ResidueMode.INTER_RESIDUE.value

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "g_ac": "NC_000007.13",
                "g_start_pos": 140453134,
                "g_end_pos": 140453137,
                "residue_mode": ResidueMode.INTER_RESIDUE,
            }
        }
    )


class ToGenomicService(BaseModelForbidExtra):
    """Service model response for cDNA -> genomic"""

    g_data: GenomicRepresentation | None = None
    warnings: list[StrictStr] = []
    service_meta: ServiceMeta

    model_config = ConfigDict(
        json_schema_extra={
            "example": {
                "g_data": {
                    "g_ac": "NC_000007.13",
                    "g_start_pos": 140453134,
                    "g_end_pos": 140453137,
                    "residue_mode": ResidueMode.INTER_RESIDUE,
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
