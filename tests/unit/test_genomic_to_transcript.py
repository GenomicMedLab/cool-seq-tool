"""Test that genomic_to_transcript method works correctly."""
import pytest
from uta_tools import UTATools


@pytest.fixture(scope="session")
async def test_uta_tools():
    """Create a UTATools test fixture"""
    test_uta_tools = UTATools()
    await test_uta_tools.uta_db._create_genomic_table()
    return test_uta_tools


@pytest.fixture(scope="module")
def tpm3_exon1():
    """Create test fixture for TPM3 exon 1."""
    return {
        "gene": "TPM3",
        "pos": 154192135,
        "exon": 1,
        "exon_offset": 0,
        "transcript": "NM_152263.3"

    }


@pytest.fixture(scope="module")
def tpm3_exon8():
    """Create test fixture for TPM3 exon 8."""
    return {
        "gene": "TPM3",
        "pos": 154170399,
        "exon": 8,
        "exon_offset": 0,
        "transcript": "NM_152263.3"

    }


@pytest.fixture(scope="module")
def tpm3_exon1_exon8():
    """Create test fixture for TPM3."""
    return {
        "gene": "TPM3",
        "start": 154192134,
        "exon_start": 1,
        "exon_start_offset": 0,
        "end": 154170399,
        "exon_end": 8,
        "exon_end_offset": 0,
        "transcript": "NM_152263.3"

    }


@pytest.fixture(scope="module")
def mane_tpm3_exon1_exon8():
    """Create test fixture for TPM3."""
    return {
        "gene": "BRAF",
        "start": 140501359,
        "exon_start": 6,
        "exon_start_offset": 0,
        "end": 140453136,
        "exon_end": 16,
        "exon_end_offset": 0,
        "transcript": "NM_001374258.1"

    }


@pytest.mark.asyncio
async def test__genomic_to_transcript(test_uta_tools, tpm3_exon1, tpm3_exon8):
    """Test that _genomic_to_transcript method works correctly."""
    resp = await test_uta_tools._genomic_to_transcript(
        "NC_000001.11", 154192135, strand=-1, transcript="NM_152263.3",
        gene="TPM3"
    )
    assert resp == tpm3_exon1

    resp = await test_uta_tools._genomic_to_transcript(
        1, 154192135, strand=-1, transcript="NM_152263.3"
    )
    assert resp == tpm3_exon1

    resp = await test_uta_tools._genomic_to_transcript(
        1, 154192135, transcript="NM_152263.3"
    )
    assert resp == tpm3_exon1

    resp = await test_uta_tools._genomic_to_transcript(
        "NC_000001.11", 154170399, strand=-1, transcript="NM_152263.3"
    )
    assert resp == tpm3_exon8

    resp = await test_uta_tools._genomic_to_transcript(
        1, 154170399, strand=-1, transcript="NM_152263.3"
    )
    assert resp == tpm3_exon8

    resp = await test_uta_tools._genomic_to_transcript(
        1, 154170399, transcript="NM_152263.3"
    )
    assert resp == tpm3_exon8


@pytest.mark.asyncio
async def test_genomic_to_transcript(test_uta_tools, tpm3_exon1_exon8,
                                     mane_tpm3_exon1_exon8):
    """Test that genomic_to_transcript method works correctly."""
    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000001.11", 154192135, 154170399, strand=-1,
        transcript="NM_152263.3")
    assert resp == tpm3_exon1_exon8

    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000007.13", 140501360, 140453136, strand=-1, gene="BRAF")
    assert resp == mane_tpm3_exon1_exon8


@pytest.mark.asyncio
async def test_invalid(test_uta_tools):
    """Test that invalid queries return `None`."""
    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000001.11", 154192135, 154170399, strand=-1,
        transcript="NM_152263.3", gene="dummy gene")
    assert resp is None

    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000001.200", 154192135, 154170399, strand=-1,
        transcript="NM_152263.3")
    assert resp is None

    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000001.11", 9999999999999, 9999999999999, strand=-1,
        transcript="NM_152263.3")
    assert resp is None

    resp = await test_uta_tools._genomic_to_transcript(
        "NC_000001.11", 154192135, strand=1, transcript="NM_152263.3",
        gene="TPM3"
    )
    assert resp is None
