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

UTA_DB_URL = environ.get("UTA_DB_URL",
                         "postgresql://uta_admin@localhost:5433/uta/uta_20210129")
SEQREPO_DATA_PATH = Path(environ.get("SEQREPO_DATA_PATH",
                                     "/usr/local/share/seqrepo/latest"))
TRANSCRIPT_MAPPINGS_PATH = Path(environ.get("TRANSCRIPT_MAPPINGS_PATH",
                                            f"{APP_ROOT}/data/transcript_mapping.tsv"))


MANE_SUMMARY_PATH = environ.get("MANE_SUMMARY_PATH")
LRG_REFSEQGENE_PATH = environ.get("LRG_REFSEQGENE_PATH")
if not all((MANE_SUMMARY_PATH, LRG_REFSEQGENE_PATH)):
    from cool_seq_tool.data import DataDownload  # noqa: E402, I202
    d = DataDownload()

    if not MANE_SUMMARY_PATH:
        MANE_SUMMARY_PATH = d._mane_summary_path

    if not LRG_REFSEQGENE_PATH:
        LRG_REFSEQGENE_PATH = d._lrg_refseqgene_path
MANE_SUMMARY_PATH = Path(MANE_SUMMARY_PATH)
LRG_REFSEQGENE_PATH = Path(LRG_REFSEQGENE_PATH)


from cool_seq_tool.cool_seq_tool import CoolSeqTool  # noqa: E402, F401, I202
