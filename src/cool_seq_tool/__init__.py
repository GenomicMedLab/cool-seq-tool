"""The cool_seq_tool package"""
import logging

logging.basicConfig(
    filename="cool_seq_tool.log",
    format="[%(asctime)s] - %(name)s - %(levelname)s : %(message)s",
)
logger = logging.getLogger("cool_seq_tool")
logger.setLevel(logging.DEBUG)

LOG_FN = "cool_seq_tool.log"
