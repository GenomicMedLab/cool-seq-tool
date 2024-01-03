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


@pytest.fixture(scope="session")
def nm_152263_exons():
    """Create test fixture for NM_152263.3 exons."""
    return [
        (0, 234),
        (234, 360),
        (360, 494),
        (494, 612),
        (612, 683),
        (683, 759),
        (759, 822),
        (822, 892),
        (892, 971),
        (971, 7099),
    ]


@pytest.fixture(scope="session")
def tpm3_1_8_start_genomic():
    """Create test fixture for genomic data for exon 1, 8"""
    return "TPM3", "NC_000001.11", 154191901, 154192135, -1


@pytest.fixture(scope="session")
def tpm3_1_8_end_genomic():
    """Create test fixture for genomic data for exon 1, 8"""
    return "TPM3", "NC_000001.11", 154170399, 154170469, -1
