"""Provide paths to shared files, and trigger data acquisition if unavailable."""
from importlib import resources
from os import environ
from pathlib import Path

from wags_tails.ncbi_lrg_refseqgene import NcbiLrgRefSeqGeneData
from wags_tails.ncbi_mane_summary import NcbiManeSummaryData

_configured_transcript_mappings_path = environ.get("TRANSCRIPT_MAPPINGS_PATH")
if _configured_transcript_mappings_path:
    TRANSCRIPT_MAPPINGS_PATH = Path(_configured_transcript_mappings_path)
    if not TRANSCRIPT_MAPPINGS_PATH.exists():
        msg = f'No transcript mappings file exists at path {_configured_transcript_mappings_path} defined under env var TRANSCRIPT_MAPPINGS_PATH. Either unset to use default file bundled with cool-seq-tool, or ensure that is available at this location. See the "Environment configuration" section under the Usage page within the documentation for more.'
        raise FileNotFoundError(msg)
else:
    TRANSCRIPT_MAPPINGS_PATH: Path = (
        resources.files(__package__) / "transcript_mapping.tsv"
    )  # type: ignore[reportAssignmentType]


def get_mane_summary() -> Path:
    """Acquire NCBI MANE summary file.

    Exact location can be user-configured with the environment variable MANE_SUMMARY_PATH.
    Otherwise, uses `wags-tails <https://wags-tails.readthedocs.io/stable/>`_ to acquire
    the latest version (either locally, if already available, or from source if out of
    date).

    :return: path to existing MANE summary file.
    :raise FileNotFoundError: if MANE_SUMMARY_PATH location doesn't point to a file that
        exists
    """
    configured_mane_summary_path = environ.get("MANE_SUMMARY_PATH")
    if configured_mane_summary_path:
        mane_summary_path = Path(configured_mane_summary_path)
        if not mane_summary_path.exists() or not mane_summary_path.is_file():
            msg = f'No MANE summary file exists at path {configured_mane_summary_path} defined under env var MANE_SUMMARY_PATH. Either unset to use the default file pattern and possibly acquire from source via `wags-tails` package, or ensure that it is available at this location. See the "Environment configuration" section under the Usage page within the documentation for more.'
            raise FileNotFoundError(msg)
    else:
        _provider = NcbiManeSummaryData(silent=True)
        mane_summary_path, _ = _provider.get_latest()
    return mane_summary_path


def get_lrg_refseqgene() -> Path:
    """Acquire NCBI LRG Refseq Gene summary file.

    Exact location can be user-configured with the environment variable LRG_REFSEQGENE_PATH.
    Otherwise, uses `wags-tails <https://wags-tails.readthedocs.io/stable/>`_ to acquire
    the latest version (either locally, if already available, or from source if out of
    date).

    :return: path to existing MANE summary file.
    :raise FileNotFoundError: if MANE_SUMMARY_PATH location doesn't point to a file that
        exists
    """
    configured_lrg_refseqgene_path = environ.get("LRG_REFSEQGENE_PATH")
    if configured_lrg_refseqgene_path:
        lrg_refseqgene_path = Path(configured_lrg_refseqgene_path)
        if not lrg_refseqgene_path.exists() or not lrg_refseqgene_path.is_file():
            msg = f'No LRG Refseq Gene exists at path {configured_lrg_refseqgene_path} defined under env var LRG_REFSEQGENE_PATH. Either unset to use the default file pattern and possibly acquire from source via `wags-tails` package, or ensure that it is available at this location. See the "Environment configuration" section under the Usage page within the documentation for more.'
            raise FileNotFoundError(msg)
    else:
        _provider = NcbiLrgRefSeqGeneData(silent=True)
        lrg_refseqgene_path, _ = _provider.get_latest()
    return lrg_refseqgene_path


SEQREPO_ROOT_DIR = environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest")
