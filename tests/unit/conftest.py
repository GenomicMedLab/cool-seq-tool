"""Provide utilities for test cases."""
from biocommons.seqrepo import SeqRepo
import pytest

from cool_seq_tool import SEQREPO_ROOT_DIR
from cool_seq_tool.data_sources import SeqRepoAccess


@pytest.fixture(scope="session")
def test_seqrepo_access():
    """Create SeqRepoAccess test fixture"""
    return SeqRepoAccess(SeqRepo(root_dir=SEQREPO_ROOT_DIR))
