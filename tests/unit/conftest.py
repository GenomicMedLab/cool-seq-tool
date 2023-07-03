"""Provide utilities for test cases."""
from biocommons.seqrepo import SeqRepo
import pytest

from cool_seq_tool.data_sources import SeqRepoAccess
from cool_seq_tool.paths import SEQREPO_ROOT_DIR


@pytest.fixture(scope="session")
def test_seqrepo_access():
    """Create SeqRepoAccess test fixture"""
    return SeqRepoAccess(SeqRepo(root_dir=SEQREPO_ROOT_DIR))
