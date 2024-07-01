"""Fetch data files regarding transcript mapping and annotation."""

import logging
from enum import Enum
from importlib import resources
from os import environ
from pathlib import Path

from wags_tails import NcbiLrgRefSeqGeneData, NcbiManeSummaryData

_logger = logging.getLogger(__name__)


class DataFile(str, Enum):
    """Constrain legal values for file resource fetching in :py:meth:`get_data_file() <cool_seq_tool.resources.data_files.get_data_file>`."""

    TRANSCRIPT_MAPPINGS = "transcript_mappings"
    MANE_SUMMARY = "mane_summary"
    LRG_REFSEQGENE = "lrg_refseqgene"

    def lower(self) -> str:
        """Return lower-cased value

        :return: lower case string
        """
        return self.value.lower()


_resource_acquisition_params = {
    DataFile.TRANSCRIPT_MAPPINGS: (
        "TRANSCRIPT_MAPPINGS_PATH",
        lambda _: resources.files(__package__) / "transcript_mapping.tsv",
    ),
    DataFile.MANE_SUMMARY: (
        "MANE_SUMMARY_PATH",
        lambda from_local: NcbiManeSummaryData(silent=True).get_latest(
            from_local=from_local
        )[0],
    ),
    DataFile.LRG_REFSEQGENE: (
        "LRG_REFSEQGENE_PATH",
        lambda from_local: NcbiLrgRefSeqGeneData(silent=True).get_latest(
            from_local=from_local
        )[0],
    ),
}


def get_data_file(resource: DataFile, from_local: bool = False) -> Path:
    """Acquire Cool-Seq-Tool file dependency.

    Each resource can be defined using an environment variable:

    * ``Resource.TRANSCRIPT_MAPPINGS`` -> ``TRANSCRIPT_MAPPINGS_PATH``
    * ``Resource.MANE_SUMMARY`` -> ``MANE_SUMMARY_PATH``
    * ``Resource.LRG_REFSEQGENE`` -> ``LRG_REFSEQGENE_PATH``

    Otherwise, this function falls back on default expected locations:

    * ``transcript_mappings.tsv`` is bundled with this library.
    * LRG RefseqGene and MANE summary files are acquired from NCBI using the `wags-tails <https://wags-tails.readthedocs.io/stable/>`_ if unavailable locally, or out of date.

    :param resource: resource to fetch
    :param from_local: if ``True``, don't check for or acquire latest version -- just
        provide most recent locally available file and raise FileNotFoundError otherwise
    :return: path to file. Consuming functions can assume that it exists and is a file.
    :raise FileNotFoundError: if file location configured by env var doesn't exist
    :raise ValueError: if file location configured by env var isn't a file
    """
    params = _resource_acquisition_params[resource]
    configured_path = environ.get(params[0])
    if configured_path:
        _logger.debug(
            "Acquiring %s via env var %s:%s", resource, params[0], configured_path
        )
        path = Path(configured_path)
        loc_descr = (
            "the default file bundled with Cool-Seq-Tool"
            if resource == DataFile.TRANSCRIPT_MAPPINGS
            else "the the default file pattern and possibly acquire from source via the `wags-tails` package"
        )
        msg = f'No {params[0].replace("_", " ").title()} file exists at path {configured_path} defined under env var {params[0]}. Either unset to use {loc_descr}, or ensure that it is available at this location. See the "Environment configuration" section under the Usage page within the documentation for more: https://coolseqtool.readthedocs.io/stable/usage.html#environment-configuration'
        if not path.exists():
            raise FileNotFoundError(msg)
        if not path.is_file():
            raise ValueError(msg)
    else:
        _logger.debug("Acquiring %s from default location/method.", resource)
        # param[1] is the resource fetcher function -- use `from_local` param to
        # optionally avoid unnecessary fetches
        path = params[1](from_local)
    _logger.debug("Acquired %s at %s", resource, path)
    return path
