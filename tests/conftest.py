"""Provide utilities for test cases."""
import asyncio

import pytest

from cool_seq_tool.app import CoolSeqTool


@pytest.fixture(scope="session")
def event_loop(request):
    """Create an instance of the default event loop for each test case."""
    loop = asyncio.get_event_loop_policy().new_event_loop()
    yield loop
    loop.close()


@pytest.fixture(scope="session")
def test_cool_seq_tool():
    """Create CoolSeqTool test fixture"""
    return CoolSeqTool()


@pytest.fixture(scope="session")
def test_seqrepo_access(test_cool_seq_tool):
    """Create SeqRepoAccess test fixture"""
    return test_cool_seq_tool.seqrepo_access


@pytest.fixture(scope="session")
def test_db(test_cool_seq_tool):
    """Create UTA Database test fixture"""
    return test_cool_seq_tool.uta_db


@pytest.fixture(scope="session")
def test_transcript_mappings(test_cool_seq_tool):
    """Create Transcript Mappings test fixture"""
    return test_cool_seq_tool.transcript_mappings


@pytest.fixture(scope="session")
def test_mane_transcript_mappings(test_cool_seq_tool):
    """Create MANE Transcript Mappings test fixture"""
    return test_cool_seq_tool.mane_transcript_mappings
