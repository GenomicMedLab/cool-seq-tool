"""Provide utilities for test cases."""
from biocommons.seqrepo import SeqRepo
import pytest

from cool_seq_tool import SEQREPO_ROOT_DIR
from cool_seq_tool.data_sources import MANETranscriptMappings,\
    SeqRepoAccess, TranscriptMappings, UTADatabase


@pytest.fixture(scope="session")
def test_seqrepo_access():
    """Create SeqRepoAccess test fixture"""
    return SeqRepoAccess(SeqRepo(root_dir=SEQREPO_ROOT_DIR))


@pytest.fixture(scope="session")
def test_transcript_mappings():
    """Create TranscriptMappings test fixture"""
    return TranscriptMappings()


@pytest.fixture(scope="session")
async def test_uta_db():
    """Create UTADatabase test fixture."""
    test_uta_db = UTADatabase()
    await test_uta_db._create_genomic_table()
    return test_uta_db


@pytest.fixture(scope="session")
def test_mane_transcript_mappings():
    """Create MANETranscriptMappings test fixture"""
    return MANETranscriptMappings()
