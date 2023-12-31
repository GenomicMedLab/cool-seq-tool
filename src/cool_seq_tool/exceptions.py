"""Module for custom exceptions"""


class CoolSeqToolError(Exception):
    """Custom exception for Cool-Seq-Tool"""


class SeqRepoAccessError(CoolSeqToolError):
    """Custom exception for SeqRepoAccess class"""
