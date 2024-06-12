"""Module for data"""
from .resources import (
    SEQREPO_ROOT_DIR,
    TRANSCRIPT_MAPPINGS_PATH,
    get_lrg_refseqgene,
    get_mane_summary,
)

__all__ = [
    "TRANSCRIPT_MAPPINGS_PATH",
    "SEQREPO_ROOT_DIR",
    "get_mane_summary",
    "get_lrg_refseqgene",
]
