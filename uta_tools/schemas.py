"""Module for data models."""
from enum import Enum
from pydantic import BaseModel
from pydantic.main import Extra
from typing import Literal


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
    exon_start: int
    exon_start_offset: int = 0
    exon_end: int
    exon_end_offset: int = 0
    transcript: str


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
