"""Module for routers"""

from enum import Enum

from cool_seq_tool.app import CoolSeqTool

cool_seq_tool = CoolSeqTool()
SERVICE_NAME = "cool_seq_tool"
RESP_DESCR = "A response to a validly-formed query."
UNHANDLED_EXCEPTION_MSG = "Unhandled exception occurred. Check logs for more details."


class Tags(str, Enum):
    """Define tags for endpoints"""

    MANE_TRANSCRIPT = "MANE Transcript"
    ALIGNMENT_MAPPER = "Alignment Mapper"
