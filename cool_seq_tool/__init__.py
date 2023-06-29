"""The cool_seq_tool package"""
from os import environ
from pathlib import Path
import logging

APP_ROOT = Path(__file__).resolve().parents[0]

logging.basicConfig(
    filename="cool_seq_tool.log",
    format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s"
)
logger = logging.getLogger("cool_seq_tool")
logger.setLevel(logging.DEBUG)

LOG_FN = "cool_seq_tool.log"

from .app import CoolSeqTool  # noqa: E402, F401, I202
