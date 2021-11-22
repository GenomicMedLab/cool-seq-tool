"""Module for data models."""
from datetime import datetime
from enum import Enum
import re
from typing import Literal, Optional, List

from pydantic import BaseModel, root_validator, validator
from pydantic.main import Extra
from pydantic.types import StrictStr, StrictInt


class ResidueMode(str, Enum):
    """Create Enum for residue modes."""

    RESIDUE: Literal["residue"] = "residue"
    INTER_RESIDUE: Literal["inter-residue"] = "inter-residue"


class TranscriptExonData(BaseModel):
    """Model containing transcript exon data."""

    class Config:
        """Class configs."""

        extra = Extra.forbid

    transcript: StrictStr
    pos: StrictInt
    exon: StrictInt
    exon_offset: StrictInt = 0
    gene: StrictStr
    chr: StrictStr
    strand: StrictInt


class GenomicData(BaseModel):
    """Model containing genomic and transcript exon data."""

    class Config:
        """Class configs."""

        extra = Extra.forbid

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


class ServiceMeta(BaseModel):
    """Metadata for uta_tools service"""

    name: Literal["uta_tools"] = "uta_tools"
    version: StrictStr
    response_datetime: datetime
    url: Literal["https://github.com/cancervariants/uta_tools"] = "https://github.com/cancervariants/uta_tools"  # noqa: E501

    @validator("version")
    def validate_version(cls, v):
        """Check version matches semantic versioning regex pattern.
        https://semver.org/#is-there-a-suggested-regular-expression-regex-to-check-a-semver-string
        """
        version_regex = r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"  # noqa: E501
        assert bool(re.match(version_regex, v))
        return v


class TranscriptExonDataResponse(BaseModel):
    """Response model for Transcript Exon Data"""

    transcript_exon_data: Optional[TranscriptExonData] = None
    warnings: Optional[List[StrictStr]] = []
    service_meta: ServiceMeta

    class Config:
        """Class configs."""

        extra = Extra.forbid


class GenomicDataResponse(BaseModel):
    """Response model for Genomic Data"""

    genomic_data: Optional[GenomicData] = None
    warnings: Optional[List[StrictStr]] = []
    service_meta: ServiceMeta

    class Config:
        """Class configs."""

        extra = Extra.forbid
