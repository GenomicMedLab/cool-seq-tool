"""Module for providing basic acquisition/setup for the various resources"""

from .mane_transcript_mappings import ManeTranscriptMappings
from .transcript_mappings import TranscriptMappings
from .uta_database import UtaDatabase

__all__ = ["ManeTranscriptMappings", "TranscriptMappings", "UtaDatabase"]
