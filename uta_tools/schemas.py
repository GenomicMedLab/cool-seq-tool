"""Module for data models."""
from enum import Enum
from pydantic import BaseModel, root_validator
from pydantic.main import Extra
from typing import Literal, Optional


class ResidueMode(str, Enum):
    """Create Enum for residue modes."""

    RESIDUE: Literal["residue"] = "residue"
    INTER_RESIDUE: Literal["inter-residue"] = "inter-residue"


class GenomicData(BaseModel):
    """Model containing genomic and transcript exon data."""

    class Config:
        """Class configs."""

        extra = Extra.forbid

    gene: str
    chr: str
    start: Optional[int] = None  # Genomic start position
    end: Optional[int] = None  # Genomic end position
    exon_start: Optional[int] = None
    exon_start_offset: Optional[int] = 0
    exon_end: Optional[int] = None
    exon_end_offset: Optional[int] = 0
    transcript: str

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


class TranscriptExonData(BaseModel):
    """Model containing transcript exon data."""

    class Config:
        """Class configs."""

        extra = Extra.forbid

    transcript: str
    pos: int
    exon: int
    exon_offset: int = 0
    gene: str
    chr: str
