"""The cool_seq_tool package"""
import logging
from pathlib import Path

APP_ROOT = Path(__file__).resolve().parents[0]

logging.basicConfig(
    filename="cool_seq_tool.log",
    format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s",
)
logger = logging.getLogger("cool_seq_tool")
logger.setLevel(logging.DEBUG)
