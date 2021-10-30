"""The uta_tools package"""
from os import environ
from pathlib import Path
import logging

APP_ROOT = Path(__file__).resolve().parents[0]

logging.basicConfig(
    filename='uta_tools.log',
    format='[%(asctime)s] - %(name)s - %(levelname)s : %(message)s'
)
logger = logging.getLogger('uta_tools')
logger.setLevel(logging.DEBUG)

if "UTA_DB_URL" in environ:
    UTA_DB_URL = environ["UTA_DB_URL"]
else:
    UTA_DB_URL = "postgresql://uta_admin@localhost:5433/uta/uta_20210129"

if "SEQREPO_DATA_PATH" in environ:
    SEQREPO_DATA_PATH = environ["SEQREPO_DATA_PATH"]
else:
    SEQREPO_DATA_PATH = "/usr/local/share/seqrepo/latest"

TRANSCRIPT_MAPPINGS_PATH = f"{APP_ROOT}/data/transcript_mapping.tsv"

from uta_tools.data import DataDownload  # noqa: E402
d = DataDownload()
MANE_SUMMARY_PATH = d._mane_summary_path
LRG_REFSEQGENE_PATH = d._lrg_refseqgene_path

from uta_tools.uta_tools import UTATools  # noqa: E402, F401
