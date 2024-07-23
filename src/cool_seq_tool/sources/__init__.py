"""Module for providing basic acquisition/setup for the various resources"""

from .mane_transcript_mappings import ManeGeneData, ManeTranscriptMappings
from .transcript_mappings import TranscriptMappings
from .uta_database import UtaDatabase

__all__ = [
    "ManeGeneData",
    "ManeTranscriptMappings",
    "TranscriptMappings",
    "UtaDatabase",
]
