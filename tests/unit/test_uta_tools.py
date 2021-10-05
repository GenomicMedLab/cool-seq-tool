"""Module for testing that UTATools works correctly."""
import pytest
from uta_tools import UTATools
import copy
from uta_tools.schemas import GenomicData, TranscriptExonData


@pytest.fixture(scope="session")
async def test_uta_tools():
    """Create a UTATools test fixture"""
    test_uta_tools = UTATools()
    await test_uta_tools.uta_db._create_genomic_table()
    return test_uta_tools


@pytest.fixture(scope="module")
def tpm3_exon1():
    """Create test fixture for TPM3 exon 1."""
    params = {
        "chr": "NC_000001.11",
        "gene": "TPM3",
        "pos": 154192135,
        "exon": 1,
        "exon_offset": 0,
        "transcript": "NM_152263.3"

    }
    return TranscriptExonData(**params)


@pytest.fixture(scope="module")
def tpm3_exon8():
    """Create test fixture for TPM3 exon 8."""
    params = {
        "chr": "NC_000001.11",
        "gene": "TPM3",
        "pos": 154170399,
        "exon": 8,
        "exon_offset": 0,
        "transcript": "NM_152263.3"

    }
    return TranscriptExonData(**params)


@pytest.fixture(scope='module')
def tpm3_exon1_exon8():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192134,
        "end": 154170399,
        "exon_start": 1,
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start_offset": 0,
        "transcript": "NM_152263.3"
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_exon8_offset():
    """Create test fixture for TPM3."""
    params = {
        "chr": "NC_000001.11",
        "gene": "TPM3",
        "start": 154192131,
        "exon_start": 1,
        "exon_start_offset": 3,
        "end": 154170404,
        "exon_end": 8,
        "exon_end_offset": -5,
        "transcript": "NM_152263.3"

    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def mane_BRAF():
    """Create test fixture for BRAF."""
    params = {
        "chr": "NC_000007.14",
        "gene": "BRAF",
        "start": 140801411,
        "exon_start": 6,
        "exon_start_offset": 148,
        "end": 140753332,
        "exon_end": 16,
        "exon_end_offset": -58,
        "transcript": "NM_001374258.1"

    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def wee1_exon2_exon11():
    """Create test fixture for WEE1."""
    params = {
        "chr": "NC_000011.10",
        "gene": "WEE1",
        "start": 9576092,
        "exon_start": 2,
        "exon_start_offset": 0,
        "end": 9588448,
        "exon_end": 11,
        "exon_end_offset": 0,
        "transcript": "NM_003390.3"

    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def mane_wee1_exon2_exon11():
    """Create test fixture for WEE1."""
    params = {
        "chr": "NC_000011.10",
        "gene": "WEE1",
        "start": 9576092,
        "exon_start": 2,
        "exon_start_offset": 0,
        "end": 9586856,
        "exon_end": 10,
        "exon_end_offset": 146,
        "transcript": "NM_003390.4"

    }
    return GenomicData(**params)


@pytest.fixture(scope='module')
def ntrk1_exon10_exon17():
    """Create test fixture for NTRK1."""
    params = {
        "gene": "NTRK1",
        "chr": "NC_000001.11",
        "start": 156874625,
        "end": 156881456,
        "exon_start": 10,
        "exon_end": 17,
        "exon_end_offset": 0,
        "exon_start_offset": 0,
        "transcript": "NM_002529.3"
    }
    return GenomicData(**params)


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
async def test_tpm3(test_uta_tools, tpm3_exon1_exon8,
                    tpm3_exon1_exon8_offset):
    """Test TPM3 genomic_to_transcript and transcript_to_genomic."""
    inputs = {
        "chromosome": "NC_000001.11",
        "start": 154192135,
        "end": 154170399,
        "strand": -1,
        "transcript": "NM_152263.3"
    }
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == tpm3_exon1_exon8
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == tpm3_exon1_exon8

    # No strand
    del inputs["strand"]
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == tpm3_exon1_exon8
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == tpm3_exon1_exon8

    # Offset, no strand
    inputs["start"] = 154192132
    inputs["end"] = 154170404
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == tpm3_exon1_exon8_offset
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == tpm3_exon1_exon8_offset

    # Offset, strand
    inputs["strand"] = -1
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == tpm3_exon1_exon8_offset
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == tpm3_exon1_exon8_offset


@pytest.mark.asyncio
async def test_braf(test_uta_tools, mane_BRAF):
    """Test BRAF genomic_to_transcript and transcript_to_genomic."""
    inputs = {
        "chromosome": "NC_000007.13",
        "start": 140501360,
        "end": 140453136,
        "strand": -1,
        "gene": "BRAF"
    }
    # MANE
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == mane_BRAF
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == mane_BRAF


@pytest.mark.asyncio
async def test_wee1(test_uta_tools, wee1_exon2_exon11, mane_wee1_exon2_exon11):
    """Test WEE1 genomic_to_transcript and transcript_to_genomic."""
    inputs = {
        "chromosome": "NC_000011.9",
        "start": 9597640,
        "end": 9609995,
        "strand": 1,
        "transcript": "NM_003390.3"
    }
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == wee1_exon2_exon11
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == wee1_exon2_exon11

    inputs["gene"] = "wee1"
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == wee1_exon2_exon11
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == wee1_exon2_exon11

    # MANE
    del inputs["transcript"]
    g_to_t_resp = await test_uta_tools.genomic_to_transcript(**inputs)
    assert g_to_t_resp == mane_wee1_exon2_exon11
    t_to_g_resp = await test_uta_tools.transcript_to_genomic(**g_to_t_resp.dict())  # noqa: E501
    assert t_to_g_resp == mane_wee1_exon2_exon11


@pytest.mark.asyncio
async def test_transcript_to_genomic(test_uta_tools, tpm3_exon1_exon8,
                                     ntrk1_exon10_exon17):
    """Test that transcript_to_genomic works correctly."""
    # TPM3
    resp = await test_uta_tools.transcript_to_genomic('NM_152263.3', 0, 8)
    assert resp == tpm3_exon1_exon8

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3       ', 0, 8)
    assert resp == tpm3_exon1_exon8

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 0, 8, gene="TPM3")
    assert resp == tpm3_exon1_exon8

    resp = await test_uta_tools.transcript_to_genomic(
        ' NM_152263.3 ', 0, 8, gene=" TPM3 ")
    assert resp == tpm3_exon1_exon8

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 0, 8, gene="tpm3")
    assert resp == tpm3_exon1_exon8

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 0, 0, gene="tpm3")
    expected = copy.deepcopy(tpm3_exon1_exon8.dict())
    expected["exon_end"] = 10
    expected["end"] = 154161812
    assert resp == expected

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 0, 8, exon_end_offset=-5)
    expected["exon_end"] = 8
    expected["exon_end_offset"] = -5
    expected["end"] = 154170404
    assert resp == expected

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 0, 8, exon_end_offset=5)
    expected["exon_end_offset"] = 5
    expected["end"] = 154170394
    assert resp == expected

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 3, 8, exon_start_offset=3, exon_end_offset=5)
    expected["exon_start"] = 3
    expected["exon_start_offset"] = 3
    expected["start"] = 154176244
    assert resp == expected

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 3, 8, exon_start_offset=-3, exon_end_offset=5)
    expected["exon_start_offset"] = -3
    expected["start"] = 154176250
    assert resp == expected

    # NTRK1
    resp = await test_uta_tools.transcript_to_genomic('NM_002529.3', 10, 0)
    assert resp == ntrk1_exon10_exon17

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_002529.3', 10, 0, gene="NTRK1")
    assert resp == ntrk1_exon10_exon17

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_002529.3', 10, 0, gene="NTRK1")
    assert resp == ntrk1_exon10_exon17

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_002529.3', 10, 0, exon_start_offset=3)
    expected = copy.deepcopy(ntrk1_exon10_exon17.dict())
    expected["exon_start_offset"] = 3
    expected["start"] = 156874628
    assert resp == expected

    resp = await test_uta_tools.transcript_to_genomic(
        'NM_002529.3', 10, 0, exon_start_offset=-3)
    expected["exon_start_offset"] = -3
    expected["start"] = 156874622
    assert resp == expected


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

    # Exon 22 does not exist
    resp = await test_uta_tools.transcript_to_genomic('NM_152263.3', 0, 22)
    assert resp is None

    # Start > End
    resp = await test_uta_tools.transcript_to_genomic('NM_152263.3', 8, 1)
    assert resp is None

    # End < Start
    resp = await test_uta_tools.transcript_to_genomic('NM_152263.3', 7, 6)
    assert resp is None

    # Transcript DNE
    resp = await test_uta_tools.transcript_to_genomic('NM_12345.6', 7, 0)
    assert resp is None

    # Index error for invalid exon
    resp = await test_uta_tools.transcript_to_genomic('NM_12345.6', -1, 0)
    assert resp is None

    # Gene that does not match transcript
    resp = await test_uta_tools.transcript_to_genomic(
        'NM_152263.3', 8, 1, gene='NTKR1')
    assert resp is None

    # No transcript given
    resp = await test_uta_tools.transcript_to_genomic(None, 8, 1, gene='NTKR1')
    assert resp is None

    resp = await test_uta_tools.transcript_to_genomic('', 8, 1, gene='NTKR1')
    assert resp is None
