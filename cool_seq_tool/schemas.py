"""Module for data models."""
from datetime import datetime
from enum import Enum
import re
from typing import Literal, Optional, List, Tuple, Union, Dict, Any, Type

from pydantic import BaseModel, root_validator, validator
from pydantic.main import Extra
from pydantic.types import StrictStr, StrictInt

from cool_seq_tool.version import __version__


class AnnotationLayer(str, Enum):
    """Create enum for supported annotation layers"""

    PROTEIN = "p"
    CDNA = "c"
    GENOMIC = "g"


class Strand(str, Enum):
    """Create enum for positive and negative strand"""

    POSITIVE = "+"
    NEGATIVE = "-"


class Assembly(str, Enum):
    """Create Enum for supported genomic assemblies"""

    GRCH37 = "GRCh37"
    GRCH38 = "GRCh38"


class TranscriptPriorityLabel(str, Enum):
    """Create Enum for Transcript Priority labels"""

    MANESelect = "mane_select"
    MANEPlusClinical = "mane_plus_clinical"
    LongestCompatibleRemaining = "longest_compatible_remaining"


class ResidueMode(str, Enum):
    """Create Enum for residue modes."""

    RESIDUE = "residue"
    INTER_RESIDUE = "inter-residue"


class BaseModelForbidExtra(BaseModel):
    """Base Pydantic model class with extra values forbidden."""

    class Config:
        """Class configs."""

        extra = Extra.forbid


class GenomicRequestBody(BaseModelForbidExtra):
    """Define constraints for genomic to transcript exon coordinates request body"""

    chromosome: Union[StrictStr, StrictInt]
    start: Optional[StrictInt] = None
    end: Optional[StrictInt] = None
    strand: Optional[StrictInt] = None
    transcript: Optional[StrictStr] = None
    gene: Optional[StrictStr] = None
    residue_mode: ResidueMode = ResidueMode.RESIDUE

    @root_validator(pre=False)
    def check_start_and_end(cls, values):
        """Check that at least one of {`start`, `end`} is set"""
        msg = "Must provide either `start` or `end`"
        start, end = values.get("start"), values.get("end")
        assert start or end, msg
        return values

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["GenomicRequestBody"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "chromosome": "NC_000001.11",
                "start": 154192135,
                "end": None,
                "strand": -1,
                "transcript": "NM_152263.3",
                "gene": "TPM3",
                "residue_mode": "residue"
            }


class TranscriptRequestBody(BaseModelForbidExtra):
    """Define constraints for transcript exon to genomic coordinates request body"""

    gene: Optional[StrictStr] = None
    transcript: Optional[StrictStr] = None
    exon_start: Optional[StrictInt] = None
    exon_start_offset: Optional[StrictInt] = 0
    exon_end: Optional[StrictInt] = None
    exon_end_offset: Optional[StrictInt] = 0

    @root_validator(pre=False)
    def check_exon_start_and_exon_end(cls, values):
        """Check that at least one of {`exon_start`, `exon_end`} is set"""
        msg = "Must provide either `exon_start` or `exon_end`"
        exon_start, exon_end = values.get("exon_start"), values.get("exon_end")
        assert exon_start or exon_end, msg
        return values

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["TranscriptRequestBody"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "gene": "TPM3",
                "transcript": "NM_152263.3",
                "exon_start": 1,
                "exon_start_offset": 1,
                "exon_end": None,
                "exon_end_offset": None,
            }


class TranscriptExonData(BaseModelForbidExtra):
    """Model containing transcript exon data."""

    transcript: StrictStr
    pos: StrictInt
    exon: StrictInt
    exon_offset: StrictInt = 0
    gene: StrictStr
    chr: StrictStr
    strand: StrictInt

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["TranscriptExonData"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "chr": "NC_000001.11",
                "gene": "TPM3",
                "pos": 154192135,
                "exon": 1,
                "exon_offset": 0,
                "transcript": "NM_152263.3",
                "strand": -1
            }


class GenomicData(BaseModelForbidExtra):
    """Model containing genomic and transcript exon data."""

    gene: StrictStr
    chr: StrictStr
    start: Optional[StrictInt] = None  # Genomic start position
    end: Optional[StrictInt] = None  # Genomic end position
    exon_start: Optional[StrictInt] = None
    exon_start_offset: Optional[StrictInt] = 0
    exon_end: Optional[StrictInt] = None
    exon_end_offset: Optional[StrictInt] = 0
    transcript: StrictStr
    strand: StrictInt

    @root_validator(pre=True)
    def check_start_end(cls, values):
        """
        Check that at least one of {`start`, `end`} is set.
        Check that at least one of {`exon_start`, `exon_end`} is set.
        If not set, set corresponding offset to `None`
        """
        msg = "Missing values for `start` or `end`"
        start = values.get("start")
        end = values.get("end")
        assert start or end, msg

        if start:
            msg = "Missing value `exon_start`"
            assert values.get("exon_start"), msg
        else:
            values["exon_start_offset"] = None

        if end:
            msg = "Missing value `exon_end`"
            assert values.get("exon_end"), msg
        else:
            values["exon_end_offset"] = None
        return values

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["GenomicData"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "gene": "TPM3",
                "chr": "NC_000001.11",
                "start": 154192135,
                "end": None,
                "exon_start": 1,
                "exon_end": None,
                "exon_start_offset": 0,
                "exon_end_offset": None,
                "transcript": "NM_152263.3",
                "strand": -1
            }


class ServiceMeta(BaseModelForbidExtra):
    """Metadata for cool_seq_tool service"""

    name: Literal["cool_seq_tool"] = "cool_seq_tool"
    version: StrictStr
    response_datetime: datetime
    url: Literal["https://github.com/GenomicMedLab/cool-seq-tool"] = "https://github.com/GenomicMedLab/cool-seq-tool"  # noqa: E501

    @validator("version")
    def validate_version(cls, v):
        """Check version matches semantic versioning regex pattern.
        https://semver.org/#is-there-a-suggested-regular-expression-regex-to-check-a-semver-string
        """
        version_regex = r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"  # noqa: E501
        assert bool(re.match(version_regex, v))
        return v

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ServiceMeta"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "name": "cool_seq_tool",
                "version": __version__,
                "response_datetime": datetime.now(),
                "url": "https://github.com/GenomicMedLab/cool-seq-tool"
            }


class TranscriptExonDataResponse(BaseModelForbidExtra):
    """Response model for Transcript Exon Data"""

    transcript_exon_data: Optional[TranscriptExonData] = None
    warnings: Optional[List[StrictStr]] = []
    service_meta: ServiceMeta

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["TranscriptExonDataResponse"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "transcript_exon_data": {
                    "chr": "NC_000001.11",
                    "gene": "TPM3",
                    "pos": 154192135,
                    "exon": 1,
                    "exon_offset": 0,
                    "transcript": "NM_152263.3",
                    "strand": -1
                },
                "warnings": list(),
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": datetime.now(),
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool"
                }
            }


class GenomicDataResponse(BaseModelForbidExtra):
    """Response model for Genomic Data"""

    genomic_data: Optional[GenomicData] = None
    warnings: Optional[List[StrictStr]] = []
    service_meta: ServiceMeta

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["GenomicDataResponse"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
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
                    "strand": -1
                },
                "warnings": list(),
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": datetime.now(),
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool"
                }
            }


class MappedManeData(BaseModel):
    """Define mapped mane data fields"""

    gene: StrictStr
    refseq: StrictStr
    ensembl: Optional[StrictStr] = None
    strand: Strand
    status: TranscriptPriorityLabel
    alt_ac: StrictStr
    assembly: Assembly

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["MappedManeData"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "gene": "BRAF",
                "refseq": "NM_001374258.1",
                "ensembl": "ENST00000644969.2",
                "strand": "-",
                "status": "mane_plus_clinical",
                "alt_ac": "NC_000007.13",
                "assembly": "GRCh37"
            }


class MappedManeDataService(BaseModelForbidExtra):
    """Service model response for mapped mane data"""

    mapped_mane_data: Optional[MappedManeData] = None
    warnings: List[StrictStr] = []
    service_meta: ServiceMeta

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["MappedManeDataService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "mapped_mane_data": {
                    "gene": "BRAF",
                    "refseq": "NM_001374258.1",
                    "ensembl": "ENST00000644969.2",
                    "strand": "-",
                    "status": "mane_plus_clinical",
                    "alt_ac": "NC_000007.13",
                    "assembly": "GRCh37"
                },
                "warnings": list(),
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": datetime.now(),
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool"
                }
            }


class ManeData(BaseModel):
    """Define mane data fields"""

    gene: Optional[StrictStr] = None
    refseq: Optional[StrictStr] = None
    ensembl: Optional[StrictStr] = None
    pos: Tuple[int, int]
    strand: Strand
    status: TranscriptPriorityLabel

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ManeData"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "gene": "BRAF",
                "refseq": "NP_004324.2",
                "ensembl": "ENSP00000493543.1",
                "pos": (598, 598),
                "strand": "-",
                "status": "mane_select"
            }


class ManeDataService(BaseModelForbidExtra):
    """Service model response for getting mane data"""

    mane_data: Optional[ManeData] = None
    warnings: List[StrictStr] = []
    service_meta: ServiceMeta

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ManeDataService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "mane_data": {
                    "gene": "BRAF",
                    "refseq": "NP_004324.2",
                    "ensembl": "ENSP00000493543.1",
                    "pos": (598, 598),
                    "strand": "-",
                    "status": "mane_select"
                },
                "warnings": list(),
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": datetime.now(),
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool"
                }
            }


# ALIGNMENT MAPPER SERVICE SCHEMAS


class CdnaRepresentation(BaseModelForbidExtra):
    """Model response for cDNA representation"""

    c_ac: str
    c_start_pos: str
    c_end_pos: str
    cds_start: int
    residue_mode = ResidueMode.INTER_RESIDUE.value

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["CdnaRepresentation"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "c_ac": "NM_004333.6",
                "c_start_pos": 1797,
                "c_end_pos": 1800,
                "cds_start": 226,
                "residue_mode": "inter-residue"
            }


class ToCdnaService(BaseModelForbidExtra):
    """Service model response for protein -> cDNA"""

    c_data: Optional[CdnaRepresentation] = None
    warnings: List[StrictStr] = []
    service_meta: ServiceMeta

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ToCdnaService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "c_data": {
                    "c_ac": "NM_004333.6",
                    "c_start_pos": 1797,
                    "c_end_pos": 1800,
                    "cds_start": 226,
                    "residue_mode": "inter-residue"
                },
                "warnings": list(),
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": datetime.now(),
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool"
                }
            }


class GenomicRepresentation(BaseModelForbidExtra):
    """Model response for genomic representation"""

    g_ac: str
    g_start_pos: int
    g_end_pos: int
    residue_mode = ResidueMode.INTER_RESIDUE.value

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["CdnaRepresentation"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "g_ac": "NC_000007.13",
                "g_start_pos": 140453134,
                "g_end_pos": 140453137,
                "residue_mode": "inter-residue"
            }


class ToGenomicService(BaseModelForbidExtra):
    """Service model response for cDNA -> genomic"""

    g_data: Optional[GenomicRepresentation] = None
    warnings: List[StrictStr] = []
    service_meta: ServiceMeta

    class Config(BaseModelForbidExtra.Config):
        """Configure model."""

        @staticmethod
        def schema_extra(schema: Dict[str, Any],
                         model: Type["ToGenomicService"]) -> None:
            """Configure OpenAPI schema."""
            if "title" in schema.keys():
                schema.pop("title", None)
            for prop in schema.get("properties", {}).values():
                prop.pop("title", None)
            schema["example"] = {
                "g_data": {
                    "g_ac": "NC_000007.13",
                    "g_start_pos": 140453134,
                    "g_end_pos": 140453137,
                    "residue_mode": "inter-residue"
                },
                "warnings": list(),
                "service_meta": {
                    "name": "cool_seq_tool",
                    "version": __version__,
                    "response_datetime": datetime.now(),
                    "url": "https://github.com/GenomicMedLab/cool-seq-tool"
                }
            }
