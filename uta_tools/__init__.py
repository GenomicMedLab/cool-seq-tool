"""The uta_tools package"""
from os import environ
import logging

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
