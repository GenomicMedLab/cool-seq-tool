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
def tpm3_exon1_exon8_offset():
    """Create test fixture for TPM3."""
    return {
        "gene": "TPM3",
        "start": 154192131,
        "exon_start": 1,
        "exon_start_offset": 3,
        "end": 154170404,
        "exon_end": 8,
        "exon_end_offset": -5,
        "transcript": "NM_152263.3"

    }


@pytest.fixture(scope="module")
def mane_BRAF():
    """Create test fixture for BRAF."""
    return {
        "gene": "BRAF",
        "start": 140801411,
        "exon_start": 6,
        "exon_start_offset": 148,
        "end": 140753332,
        "exon_end": 16,
        "exon_end_offset": -58,
        "transcript": "NM_001374258.1"

    }


@pytest.fixture(scope="module")
def wee1_exon2_exon11():
    """Create test fixture for WEE1."""
    return {
        "gene": "WEE1",
        "start": 9597639,
        "exon_start": 2,
        "exon_start_offset": 0,
        "end": 9609995,
        "exon_end": 11,
        "exon_end_offset": 0,
        "transcript": "NM_003390.3"

    }


@pytest.fixture(scope="module")
def mane_wee1_exon2_exon11():
    """Create test fixture for WEE1."""
    return {
        "gene": "WEE1",
        "start": 9576092,
        "exon_start": 2,
        "exon_start_offset": 0,
        "end": 9586856,
        "exon_end": 10,
        "exon_end_offset": 146,
        "transcript": "NM_003390.4"

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
        "NC_000001.11", 154170399, strand=-1, transcript="NM_152263.3",
        is_start=False
    )
    assert resp == tpm3_exon8

    resp = await test_uta_tools._genomic_to_transcript(
        1, 154170399, strand=-1, transcript="NM_152263.3", is_start=False
    )
    assert resp == tpm3_exon8

    resp = await test_uta_tools._genomic_to_transcript(
        1, 154170399, transcript="NM_152263.3", is_start=False
    )
    assert resp == tpm3_exon8


@pytest.mark.asyncio
async def test_genomic_to_transcript(test_uta_tools, tpm3_exon1_exon8,
                                     mane_BRAF,
                                     tpm3_exon1_exon8_offset,
                                     wee1_exon2_exon11,
                                     mane_wee1_exon2_exon11):
    """Test that genomic_to_transcript method works correctly."""
    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000001.11", 154192135, 154170399, strand=-1,
        transcript="NM_152263.3")
    assert resp == tpm3_exon1_exon8

    # Offset, no strand provided
    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000001.11", 154192132, 154170404, transcript="NM_152263.3")
    assert resp == tpm3_exon1_exon8_offset

    # Offset, strand provided
    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000001.11", 154192132, 154170404, strand=-1,
        transcript="NM_152263.3")
    assert resp == tpm3_exon1_exon8_offset

    # MANE Transcript Test
    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000007.13", 140501360, 140453136, strand=-1, gene="BRAF")
    assert resp == mane_BRAF

    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000011.9", 9597640, 9609995, strand=1, transcript="NM_003390.3"
    )
    assert resp == wee1_exon2_exon11

    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000011.9", 9597640, 9609995, transcript="NM_003390.3", gene="wee1"
    )
    assert resp == wee1_exon2_exon11

    # MANE Transcript
    resp = await test_uta_tools.genomic_to_transcript(
        "NC_000011.9", 9597640, 9609995, gene="wee1"
    )
    assert resp == mane_wee1_exon2_exon11


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
