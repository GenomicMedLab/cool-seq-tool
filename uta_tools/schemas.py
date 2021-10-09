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
    start: int  # Genomic start position
    end: int  # Genomic end position
    exon_start: Optional[int] = None
    exon_start_offset: Optional[int] = 0
    exon_end: Optional[int] = None
    exon_end_offset: Optional[int] = 0
    transcript: str

    @root_validator(pre=True)
    def check_exons(cls, values):
        """Check that at least one of {`exon_start`, `exon_end`} is set.
        If not set, set corresponding offset to `None`
        """
        msg = "Must give values for either `exon_start`, `exon_end`, or both"
        exon_start = values.get("exon_start")
        exon_end = values.get("exon_end")
        assert exon_start or exon_end, msg

        for exon, exon_offset in [(exon_start, "exon_start_offset"),
                                  (exon_end, "exon_end_offset")]:
            if exon is None:
                values[exon_offset] = None
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
