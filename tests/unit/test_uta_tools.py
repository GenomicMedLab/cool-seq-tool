"""Module for testing that UTATools works correctly."""
from datetime import datetime
import re
import copy

import pytest

from uta_tools import UTATools
from uta_tools.schemas import GenomicData, TranscriptExonData


@pytest.fixture(scope="module")
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
        "transcript": "NM_152263.3",
        "strand": -1
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
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return TranscriptExonData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_g():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192135,
        "end": None,
        "exon_start": 1,
        "exon_end": None,
        "exon_start_offset": 0,
        "exon_end_offset": None,
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon8_g():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": None,
        "end": 154170399,
        "exon_start": None,
        "exon_end": 8,
        "exon_start_offset": None,
        "exon_end_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_exon8():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192135,
        "end": 154170399,
        "exon_start": 1,
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_t_to_g():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192135,
        "end": None,
        "exon_start": 1,
        "exon_end": None,
        "exon_end_offset": None,
        "exon_start_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon8_t_to_g():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": None,
        "end": 154170399,
        "exon_start": None,
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start_offset": None,
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_exon8_t_to_g():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192135,
        "end": 154170399,
        "exon_start": 1,
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_exon8_offset():
    """Create test fixture for TPM3."""
    params = {
        "chr": "NC_000001.11",
        "gene": "TPM3",
        "start": 154192132,
        "exon_start": 1,
        "exon_start_offset": 3,
        "end": 154170404,
        "exon_end": 8,
        "exon_end_offset": -5,
        "transcript": "NM_152263.3",
        "strand": -1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def mane_braf():
    """Create test fixture for BRAF."""
    params = {
        "chr": "NC_000007.14",
        "gene": "BRAF",
        "start": 140808062,
        "exon_start": 5,
        "exon_start_offset": 0,
        "end": 140753331,
        "exon_end": 15,
        "exon_end_offset": -57,
        "transcript": "NM_004333.6",
        "strand": -1
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
        "exon_start_offset": -1,
        "end": 9588448,
        "exon_end": 11,
        "exon_end_offset": 0,
        "transcript": "NM_003390.3",
        "strand": 1
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
        "exon_start_offset": -1,
        "end": 9586856,
        "exon_end": 10,
        "exon_end_offset": 146,
        "transcript": "NM_003390.4",
        "strand": 1
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def ntrk1_exon10_exon17():
    """Create test fixture for NTRK1."""
    params = {
        "gene": "NTRK1",
        "chr": "NC_000001.11",
        "start": 156874626,
        "end": 156881456,
        "exon_start": 10,
        "exon_end": 17,
        "exon_end_offset": 0,
        "exon_start_offset": 0,
        "transcript": "NM_002529.3",
        "strand": 1
    }
    return GenomicData(**params)


def check_service_meta(actual):
    """Check that service metadata matches expected

    :param ServiceMeta actual: Actual service metadata
    """
    assert actual.name == "uta_tools"
    version_regex = r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"  # noqa: E501
    assert bool(re.match(version_regex, actual.version))
    assert isinstance(actual.response_datetime, datetime)
    assert actual.url == "https://github.com/cancervariants/uta_tools"


def genomic_data_assertion_checks(actual, expected=None, is_valid=True):
    """Check that actual matches expected for both valid and invalid
    genomic data responses

    :param GenomicDataResponse actual: Actual data
    :param GenomicData expected: Expected GenomicData
    :param bool is_valid: `True` if expected is valid response.
        `False` otherwise.
    """
    if is_valid:
        assert actual.genomic_data == expected
        assert actual.warnings == []
    else:
        assert actual.genomic_data is None
        assert len(actual.warnings) > 0
    check_service_meta(actual.service_meta)


def transcript_exon_data_assertion_checks(actual, expected=None,
                                          is_valid=True):
    """Check that actual matches expected for both valid and invalid
    transcript exon data responses

    :param TranscriptExonDataResponse actual: Actual data
    :param TranscriptExonData expected: Expected TranscriptExonData
    :param bool is_valid: `True` if expected is valid response.
        `False` otherwise.
    """
    if is_valid:
        assert actual.transcript_exon_data == expected
        assert actual.warnings == []
    else:
        assert actual.transcript_exon_data is None
        assert len(actual.warnings) > 0
    check_service_meta(actual.service_meta)


@pytest.mark.asyncio
async def test__genomic_to_transcript(test_uta_tools, tpm3_exon1, tpm3_exon8):
    """Test that _genomic_to_transcript_exon_coordinate
    method works correctly.
    """
    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154192135, strand=-1, transcript="NM_152263.3",
        gene="TPM3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        1, 154192135, strand=-1, transcript="NM_152263.3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        1, 154192135, transcript="NM_152263.3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154170399, strand=-1, transcript="NM_152263.3",
        is_start=False
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)

    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        1, 154170399, strand=-1, transcript="NM_152263.3", is_start=False
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)

    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        1, 154170399, transcript="NM_152263.3", is_start=False
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)


@pytest.mark.asyncio
async def test_tpm3(test_uta_tools, tpm3_exon1_exon8,
                    tpm3_exon1_exon8_offset, tpm3_exon1_g, tpm3_exon8_g,
                    tpm3_exon1_exon8_t_to_g):
    """Test TPM3 genomic_to_transcript_exon_coordinates and
    transcript_to_genomic_coordinates.
    """
    inputs = {
        "chromosome": "NC_000001.11",
        "start": 154192136,
        "end": 154170400,
        "strand": -1,
        "transcript": "NM_152263.3"
    }
    tpm3_exon1_exon8_t_to_g = copy.deepcopy(tpm3_exon1_exon8)
    tpm3_exon1_exon8_t_to_g.start = 154192135

    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_t_to_g)

    inputs["residue_mode"] = "INTER-RESIDUE"
    inputs["start"] = 154192135
    inputs["end"] = 154170399
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8_t_to_g)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_t_to_g)

    # No strand
    del inputs["strand"]
    del inputs["residue_mode"]
    inputs["start"] = 154192136
    inputs["end"] = 154170400
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_t_to_g)

    # Offset, no strand
    inputs["start"] = 154192132
    inputs["end"] = 154170404
    inputs["residue_mode"] = "INTER-RESIDUE"
    tpm3_exon1_exon8_offset_t_to_g = copy.deepcopy(tpm3_exon1_exon8_offset)
    tpm3_exon1_exon8_offset_t_to_g.start = 154192132
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8_offset)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_offset_t_to_g)

    # Offset, strand
    inputs["strand"] = -1
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8_offset)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_offset_t_to_g)

    # Test only setting start
    inputs = {
        "chromosome": "NC_000001.11",
        "start": 154192135,
        "strand": -1,
        "transcript": "NM_152263.3",
        "residue_mode": "inter-residue"
    }
    tpm3_exon1_exon8_t_to_g = copy.deepcopy(tpm3_exon1_g)
    tpm3_exon1_exon8_t_to_g.start = 154192135

    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_g)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_t_to_g)

    # Test only setting end
    inputs = {
        "chromosome": "NC_000001.11",
        "end": 154170399,
        "strand": -1,
        "transcript": "NM_152263.3",
        "residue_mode": "inter-residue"
    }
    tpm3_exon1_exon8_t_to_g = copy.deepcopy(tpm3_exon8_g)

    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon8_g)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_t_to_g)


@pytest.mark.asyncio
async def test_braf(test_uta_tools, mane_braf):
    """Test BRAF genomic_to_transcript_exon_coordinates and
    transcript_to_genomic_coordinates.
    """
    inputs = {
        "chromosome": "NC_000007.13",
        "start": 140501360,
        "end": 140453136,
        "strand": -1,
        "gene": "BRAF"
    }
    # MANE
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_braf)

    del inputs["strand"]
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_braf)

    mane_braf_t_to_g = copy.deepcopy(mane_braf)
    t_to_g_resp = \
        await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    mane_braf_t_to_g.start = 140808062
    genomic_data_assertion_checks(t_to_g_resp, mane_braf_t_to_g)


@pytest.mark.asyncio
async def test_wee1(test_uta_tools, wee1_exon2_exon11, mane_wee1_exon2_exon11):
    """Test WEE1 genomic_to_transcript_exon_coordinates and
    transcript_to_genomic_coordinates.
    """
    inputs = {
        "chromosome": "NC_000011.9",
        "start": 9597640,
        "end": 9609996,
        "strand": 1,
        "transcript": "NM_003390.3"
    }
    wee1_exon2_exon11_t_to_g = copy.deepcopy(wee1_exon2_exon11)
    wee1_exon2_exon11_t_to_g.start = 9576092
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, wee1_exon2_exon11)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, wee1_exon2_exon11_t_to_g)

    inputs["gene"] = "wee1"
    del inputs["strand"]
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, wee1_exon2_exon11)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, wee1_exon2_exon11_t_to_g)

    # MANE
    del inputs["transcript"]
    mane_wee1_exon2_exon11_t_to_g = copy.deepcopy(mane_wee1_exon2_exon11)
    mane_wee1_exon2_exon11_t_to_g.start = 9576092
    g_to_t_resp = \
        await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_wee1_exon2_exon11)
    t_to_g_resp = await test_uta_tools.transcript_to_genomic_coordinates(**g_to_t_resp.genomic_data.dict())  # noqa: E501
    genomic_data_assertion_checks(t_to_g_resp, mane_wee1_exon2_exon11_t_to_g)


@pytest.mark.asyncio
async def test_transcript_to_genomic(test_uta_tools, tpm3_exon1_exon8_t_to_g,
                                     tpm3_exon1_t_to_g, tpm3_exon8_t_to_g,
                                     ntrk1_exon10_exon17):
    """Test that transcript_to_genomic_coordinates works correctly."""
    # TPM3
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, tpm3_exon8_t_to_g)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=None, transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, tpm3_exon1_t_to_g)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, transcript="NM_152263.3       ")
    genomic_data_assertion_checks(resp, tpm3_exon8_t_to_g)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, gene="TPM3", transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, tpm3_exon8_t_to_g)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, gene=" TPM3 ", transcript=" NM_152263.3 ")
    genomic_data_assertion_checks(resp, tpm3_exon8_t_to_g)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, gene="tpm3", transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, tpm3_exon8_t_to_g)

    expected = copy.deepcopy(tpm3_exon1_exon8_t_to_g)
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, exon_end_offset=-5, transcript="NM_152263.3")
    expected.exon_end = 8
    expected.exon_end_offset = -5
    expected.end = 154170404
    genomic_data_assertion_checks(resp, expected)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, exon_end_offset=5, transcript="NM_152263.3")
    expected.exon_end_offset = 5
    expected.end = 154170394
    genomic_data_assertion_checks(resp, expected)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=3, exon_end=8, exon_start_offset=3, exon_end_offset=5,
        transcript="NM_152263.3")
    expected.exon_start = 3
    expected.exon_start_offset = 3
    expected.start = 154176245
    genomic_data_assertion_checks(resp, expected)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=3, exon_end=8, exon_start_offset=-3, exon_end_offset=5,
        transcript="NM_152263.3")
    expected.exon_start_offset = -3
    expected.start = 154176251
    genomic_data_assertion_checks(resp, expected)

    # NTRK1
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, transcript="NM_002529.3")
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, gene="NTRK1", transcript="NM_002529.3")
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, gene="NTRK1", transcript="NM_002529.3")
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, exon_start_offset=3,
        transcript="NM_002529.3")
    expected = copy.deepcopy(ntrk1_exon10_exon17)
    expected.exon_start_offset = 3
    expected.start = 156874629
    genomic_data_assertion_checks(resp, expected)

    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, exon_start_offset=-3,
        transcript="NM_002529.3")
    expected.exon_start_offset = -3
    expected.start = 156874623
    genomic_data_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_valid_inputs(test_uta_tools):
    """Test that valid inputs don"t return any errors"""
    inputs = {
        "gene": "TPM3",
        "chromosome": "NC_000001.11",
        "start": 154171413
    }
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)  # noqa: E501
    assert resp.genomic_data

    inputs = {
        "gene": "WEE1",
        "chromosome": "NC_000011.9",
        "end": 9609996
    }
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)  # noqa: E501
    assert resp.genomic_data

    inputs["chromosome"] = "11"
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)  # noqa: E501
    assert resp.genomic_data

    inputs = {
        "transcript": "NM_003390.3",
        "exon_start": 2
    }
    resp = await test_uta_tools.transcript_to_genomic_coordinates(**inputs)
    assert resp.genomic_data

    inputs["gene"] = "WEE1"
    resp = await test_uta_tools.transcript_to_genomic_coordinates(**inputs)
    assert resp.genomic_data

    # Test X/Y chromosome bug
    inputs = {
        "chromosome": "X",
        "strand": 1,
        "start": 154437254,
        "end": 154437299,
        "gene": "GDI1",
        "residue_mode": "inter-residue"
    }
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(**inputs)
    assert resp.genomic_data


@pytest.mark.asyncio
async def test_invalid(test_uta_tools):
    """Test that invalid queries return `None`."""
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(
        transcript="NM_152263 3", start=154192135, end=154170399,
        chromosome="NC_000001.11"
    )
    assert resp.warnings == ["Unable to get exons for NM_152263 3"]

    # start and end not given
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(
        "NC_000001.11", start=None, end=None, strand=-1,
        transcript="NM_152263.3", gene="TPM3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `start` or `end`"]

    # Invalid gene
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(
        "NC_000001.11", start=154192135, end=154170399, strand=-1,
        transcript="NM_152263.3", gene="dummy gene")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Unable to find a result for chromosome NC_000001.11 "
                             "where genomic coordinate 154192134 is mapped between an "
                             "exon's start and end coordinates on the negative strand "
                             "and on gene DUMMY GENE"]

    # Invalid chromosome
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(
        "NC_000001.200", start=154192135, end=154170399, strand=-1,
        transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Invalid chromosome: NC_000001.200"]

    # Invalid coordinates
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(
        "NC_000001.11", start=9999999999999, end=9999999999999, strand=-1,
        transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result for chromosome NC_000001.11 where genomic "
        "coordinate 9999999999998 is mapped between an exon's start and end "
        "coordinates on the negative strand"]

    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(
        chromosome="1", start=154170400, strand=-1, transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Must find exactly one row for genomic data, but found: 0"
    ]

    # Strand does not match
    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154192135, strand=1, transcript="NM_152263.3",
        gene="TPM3"
    )
    transcript_exon_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result for chromosome NC_000001.11 where genomic "
        "coordinate 154192135 is mapped between an exon's start and end "
        "coordinates on the positive strand and on gene TPM3"
    ]

    # Must supply either gene or transcript
    resp = await test_uta_tools._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154192135, strand=1
    )
    transcript_exon_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `gene` or `transcript`"]

    # Exon 22 does not exist
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=22, transcript="NM_152263.3", )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Exon 22 does not exist on NM_152263.3"]

    # Start > End
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=8, exon_end=1, transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Start exon 8 is greater than end exon 1"]

    # Transcript DNE
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=7, exon_end=None, transcript="NM_12345.6")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Unable to get exons for NM_12345.6"]

    # Index error for invalid exon
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=-1, exon_end=0, transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Exon -1 does not exist on NM_152263.3"]

    # Cant supply 0 based exons
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=0, exon_end=1, transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Exon 0 does not exist on NM_152263.3"]

    # Gene that does not match transcript
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, gene="NTKR1", transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result where NM_152263.3 has transcript coordinates"
        " 117 and 234 between an exon's start and end coordinates on gene "
        "NTKR1"]

    # No transcript given
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, gene="NTKR1", transcript=None)
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide `transcript`"]

    # No transcript given
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, gene="NTKR1", transcript="")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide `transcript`"]

    # No exons given
    resp = await test_uta_tools.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=None, transcript="NM_152263.3")
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `exon_start` or `exon_end`"]

    # Does not live on exon
    params = {"chromosome": 1, "strand": 1, "start": 3629864, "gene": "TPRG1L",
              "residue_mode": "inter-residue"}
    resp = await test_uta_tools.genomic_to_transcript_exon_coordinates(
        **params)
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["NM_182752.4 with position 2138 does not exist"
                             " on exons: [(58, 259), (259, 351), (351, 528),"
                             " (528, 682), (682, 877)]"]


def test_get_fasta_file(test_uta_tools, tmp_path):
    """Test get_fasta_file method"""
    tpm3 = tmp_path / "NM_002529.3.fasta"
    test_uta_tools.get_fasta_file("NM_002529.3", tpm3)
    tpm3_expected = """>ref|refseq:NM_002529.3|ga4gh:SQ.RSkww1aYmsMiWbNdNnOTnVDAM3ZWp1uA
TGCAGCTGGGAGCGCACAGACGGCTGCCCCGCCTGAGCGAGGCGGGCGCCGCCGCGATGC
TGCGAGGCGGACGGCGCGGGCAGCTTGGCTGGCACAGCTGGGCTGCGGGGCCGGGCAGCC
TGCTGGCTTGGCTGATACTGGCATCTGCGGGCGCCGCACCCTGCCCCGATGCCTGCTGCC
CCCACGGCTCCTCGGGACTGCGATGCACCCGGGATGGGGCCCTGGATAGCCTCCACCACC
TGCCCGGCGCAGAGAACCTGACTGAGCTCTACATCGAGAACCAGCAGCATCTGCAGCATC
TGGAGCTCCGTGATCTGAGGGGCCTGGGGGAGCTGAGAAACCTCACCATCGTGAAGAGTG
GTCTCCGTTTCGTGGCGCCAGATGCCTTCCATTTCACTCCTCGGCTCAGTCGCCTGAATC
TCTCCTTCAACGCTCTGGAGTCTCTCTCCTGGAAAACTGTGCAGGGCCTCTCCTTACAGG
AACTGGTCCTGTCGGGGAACCCTCTGCACTGTTCTTGTGCCCTGCGCTGGCTACAGCGCT
GGGAGGAGGAGGGACTGGGCGGAGTGCCTGAACAGAAGCTGCAGTGTCATGGGCAAGGGC
CCCTGGCCCACATGCCCAATGCCAGCTGTGGTGTGCCCACGCTGAAGGTCCAGGTGCCCA
ATGCCTCGGTGGATGTGGGGGACGACGTGCTGCTGCGGTGCCAGGTGGAGGGGCGGGGCC
TGGAGCAGGCCGGCTGGATCCTCACAGAGCTGGAGCAGTCAGCCACGGTGATGAAATCTG
GGGGTCTGCCATCCCTGGGGCTGACCCTGGCCAATGTCACCAGTGACCTCAACAGGAAGA
ACGTGACGTGCTGGGCAGAGAACGATGTGGGCCGGGCAGAGGTCTCTGTTCAGGTCAACG
TCTCCTTCCCGGCCAGTGTGCAGCTGCACACGGCGGTGGAGATGCACCACTGGTGCATCC
CCTTCTCTGTGGATGGGCAGCCGGCACCGTCTCTGCGCTGGCTCTTCAATGGCTCCGTGC
TCAATGAGACCAGCTTCATCTTCACTGAGTTCCTGGAGCCGGCAGCCAATGAGACCGTGC
GGCACGGGTGTCTGCGCCTCAACCAGCCCACCCACGTCAACAACGGCAACTACACGCTGC
TGGCTGCCAACCCCTTCGGCCAGGCCTCCGCCTCCATCATGGCTGCCTTCATGGACAACC
CTTTCGAGTTCAACCCCGAGGACCCCATCCCTGTCTCCTTCTCGCCGGTGGACACTAACA
GCACATCTGGAGACCCGGTGGAGAAGAAGGACGAAACACCTTTTGGGGTCTCGGTGGCTG
TGGGCCTGGCCGTCTTTGCCTGCCTCTTCCTTTCTACGCTGCTCCTTGTGCTCAACAAAT
GTGGACGGAGAAACAAGTTTGGGATCAACCGCCCGGCTGTGCTGGCTCCAGAGGATGGGC
TGGCCATGTCCCTGCATTTCATGACATTGGGTGGCAGCTCCCTGTCCCCCACCGAGGGCA
AAGGCTCTGGGCTCCAAGGCCACATCATCGAGAACCCACAATACTTCAGTGATGCCTGTG
TTCACCACATCAAGCGCCGGGACATCGTGCTCAAGTGGGAGCTGGGGGAGGGCGCCTTTG
GGAAGGTCTTCCTTGCTGAGTGCCACAACCTCCTGCCTGAGCAGGACAAGATGCTGGTGG
CTGTCAAGGCACTGAAGGAGGCGTCCGAGAGTGCTCGGCAGGACTTCCAGCGTGAGGCTG
AGCTGCTCACCATGCTGCAGCACCAGCACATCGTGCGCTTCTTCGGCGTCTGCACCGAGG
GCCGCCCCCTGCTCATGGTCTTTGAGTATATGCGGCACGGGGACCTCAACCGCTTCCTCC
GATCCCATGGACCTGATGCCAAGCTGCTGGCTGGTGGGGAGGATGTGGCTCCAGGCCCCC
TGGGTCTGGGGCAGCTGCTGGCCGTGGCTAGCCAGGTCGCTGCGGGGATGGTGTACCTGG
CGGGTCTGCATTTTGTGCACCGGGACCTGGCCACACGCAACTGTCTAGTGGGCCAGGGAC
TGGTGGTCAAGATTGGTGATTTTGGCATGAGCAGGGATATCTACAGCACCGACTATTACC
GTGTGGGAGGCCGCACCATGCTGCCCATTCGCTGGATGCCGCCCGAGAGCATCCTGTACC
GTAAGTTCACCACCGAGAGCGACGTGTGGAGCTTCGGCGTGGTGCTCTGGGAGATCTTCA
CCTACGGCAAGCAGCCCTGGTACCAGCTCTCCAACACGGAGGCAATCGACTGCATCACGC
AGGGACGTGAGTTGGAGCGGCCACGTGCCTGCCCACCAGAGGTCTACGCCATCATGCGGG
GCTGCTGGCAGCGGGAGCCCCAGCAACGCCACAGCATCAAGGATGTGCACGCCCGGCTGC
AAGCCCTGGCCCAGGCACCTCCTGTCTACCTGGATGTCCTGGGCTAGGGGGCCGGCCCAG
GGGCTGGGAGTGGTTAGCCGGAATACTGGGGCCTGCCCTCAGCATCCCCCATAGCTCCCA
GCAGCCCCAGGGTGATCTCAAAGTATCTAATTCACCCTCAGCATGTGGGAAGGGACAGGT
GGGGGCTGGGAGTAGAGGATGTTCCTGCTTCTCTAGGCAAGGTCCCGTCATAGCAATTAT
ATTTATTATCCCTTGAAAAAAAA"""
    assert tpm3.read_text() == tpm3_expected

    limk2 = tmp_path / "ENST00000331728.9.fasta"
    test_uta_tools.get_fasta_file("ENST00000331728.9", limk2)
    limk2_expected = """>emb|ensembl:ENST00000331728.9|refseq:NM_005569.4|ga4gh:SQ.7_mlQyDN-uWH0RlxTQFvFEv6ykd2D-xF
GTCTTCCCGCGCCTGAGGCGGCGGCGGCAGGAGCTGAGGGGAGTTGTAGGGAACTGAGGG
GAGCTGCTGTGTCCCCCGCCTCCTCCTCCCCATTTCCGCGCTCCCGGGACCATGTCCGCG
CTGGCGGGTGAAGATGTCTGGAGGTGTCCAGGCTGTGGGGACCACATTGCTCCAAGCCAG
ATATGGTACAGGACTGTCAACGAAACCTGGCACGGCTCTTGCTTCCGGTGTTCAGAATGC
CAGGATTCCCTCACCAACTGGTACTATGAGAAGGATGGGAAGCTCTACTGCCCCAAGGAC
TACTGGGGGAAGTTTGGGGAGTTCTGTCATGGGTGCTCCCTGCTGATGACAGGGCCTTTT
ATGGTGGCTGGGGAGTTCAAGTACCACCCAGAGTGCTTTGCCTGTATGAGCTGCAAGGTG
ATCATTGAGGATGGGGATGCATATGCACTGGTGCAGCATGCCACCCTCTACTGTGGGAAG
TGCCACAATGAGGTGGTGCTGGCACCCATGTTTGAGAGACTCTCCACAGAGTCTGTTCAG
GAGCAGCTGCCCTACTCTGTCACGCTCATCTCCATGCCGGCCACCACTGAAGGCAGGCGG
GGCTTCTCCGTGTCCGTGGAGAGTGCCTGCTCCAACTACGCCACCACTGTGCAAGTGAAA
GAGGTCAACCGGATGCACATCAGTCCCAACAATCGAAACGCCATCCACCCTGGGGACCGC
ATCCTGGAGATCAATGGGACCCCCGTCCGCACACTTCGAGTGGAGGAGGTGGAGGATGCA
ATTAGCCAGACGAGCCAGACACTTCAGCTGTTGATTGAACATGACCCCGTCTCCCAACGC
CTGGACCAGCTGCGGCTGGAGGCCCGGCTCGCTCCTCACATGCAGAATGCCGGACACCCC
CACGCCCTCAGCACCCTGGACACCAAGGAGAATCTGGAGGGGACACTGAGGAGACGTTCC
CTAAGGCGCAGTAACAGTATCTCCAAGTCCCCTGGCCCCAGCTCCCCAAAGGAGCCCCTG
CTGTTCAGCCGTGACATCAGCCGCTCAGAATCCCTTCGTTGTTCCAGCAGCTATTCACAG
CAGATCTTCCGGCCCTGTGACCTAATCCATGGGGAGGTCCTGGGGAAGGGCTTCTTTGGG
CAGGCTATCAAGGTGACACACAAAGCCACGGGCAAAGTGATGGTCATGAAAGAGTTAATT
CGATGTGATGAGGAGACCCAGAAAACTTTTCTGACTGAGGTGAAAGTGATGCGCAGCCTG
GACCACCCCAATGTGCTCAAGTTCATTGGTGTGCTGTACAAGGATAAGAAGCTGAACCTC
CTGACAGAGTACATTGAGGGGGGCACACTGAAGGACTTTCTGCGCAGTATGGATCCGTTC
CCCTGGCAGCAGAAGGTCAGGTTTGCCAAAGGAATCGCCTCCGGAATGGCCTATTTGCAC
TCTATGTGCATCATCCACCGGGATCTGAACTCGCACAACTGCCTCATCAAGTTGGACAAG
ACTGTGGTGGTGGCAGACTTTGGGCTGTCACGGCTCATAGTGGAAGAGAGGAAAAGGGCC
CCCATGGAGAAGGCCACCACCAAGAAACGCACCTTGCGCAAGAACGACCGCAAGAAGCGC
TACACGGTGGTGGGAAACCCCTACTGGATGGCCCCTGAGATGCTGAACGGAAAGAGCTAT
GATGAGACGGTGGATATCTTCTCCTTTGGGATCGTTCTCTGTGAGATCATTGGGCAGGTG
TATGCAGATCCTGACTGCCTTCCCCGAACACTGGACTTTGGCCTCAACGTGAAGCTTTTC
TGGGAGAAGTTTGTTCCCACAGATTGTCCCCCGGCCTTCTTCCCGCTGGCCGCCATCTGC
TGCAGACTGGAGCCTGAGAGCAGACCAGCATTCTCGAAATTGGAGGACTCCTTTGAGGCC
CTCTCCCTGTACCTGGGGGAGCTGGGCATCCCGCTGCCTGCAGAGCTGGAGGAGTTGGAC
CACACTGTGAGCATGCAGTACGGCCTGACCCGGGACTCACCTCCCTAGCCCTGGCCCAGC
CCCCTGCAGGGGGGTGTTCTACAGCCAGCATTGCCCCTCTGTGCCCCATTCCTGCTGTGA
GCAGGGCCGTCCGGGCTTCCTGTGGATTGGCGGAATGTTTAGAAGCAGAACAAGCCATTC
CTATTACCTCCCCAGGAGGCAAGTGGGCGCAGCACCAGGGAAATGTATCTCCACAGGTTC
TGGGGCCTAGTTACTGTCTGTAAATCCAATACTTGCCTGAAAGCTGTGAAGAAGAAAAAA
ACCCCTGGCCTTTGGGCCAGGAGGAATCTGTTACTCGAATCCACCCAGGAACTCCCTGGC
AGTGGATTGTGGGAGGCTCTTGCTTACACTAATCAGCGTGACCTGGACCTGCTGGGCAGG
ATCCCAGGGTGAACCTGCCTGTGAACTCTGAAGTCACTAGTCCAGCTGGGTGCAGGAGGA
CTTCAAGTGTGTGGACGAAAGAAAGACTGATGGCTCAAAGGGTGTGAAAAAGTCAGTGAT
GCTCCCCCTTTCTACTCCAGATCCTGTCCTTCCTGGAGCAAGGTTGAGGGAGTAGGTTTT
GAAGAGTCCCTTAATATGTGGTGGAACAGGCCAGGAGTTAGAGAAAGGGCTGGCTTCTGT
TTACCTGCTCACTGGCTCTAGCCAGCCCAGGGACCACATCAATGTGAGAGGAAGCCTCCA
CCTCATGTTTTCAAACTTAATACTGGAGACTGGCTGAGAACTTACGGACAACATCCTTTC
TGTCTGAAACAAACAGTCACAAGCAAAGGAAGAGGCTGGGGGACTAGAAAGAGGCCCTGC
CCTCTAGAAAGCTCAGATCTTGGCTTCTGTTACTCATACTCGGGTGGGCTCCTTAGTCAG
ATGCCTAAAACATTTTGCCTAAAGCTCGATGGGTTCTGGAGGACAGTGTGGCTTGTCACA
GGCCTAGAGTCTGAGGGAGGGGAGTGGGAGTCTCAGCAATCTCTTGGTCTTGGCTTCATG
GCAACCACTGCTCACCCTTCAACATGCCTGGTTTAGGCAGCAGCTTGGGCTGGGAAGAGG
TGGTGGCAGAGTCTCAAAGCTGAGATGCTGAGAGAGATAGCTCCCTGAGCTGGGCCATCT
GACTTCTACCTCCCATGTTTGCTCTCCCAACTCATTAGCTCCTGGGCAGCATCCTCCTGA
GCCACATGTGCAGGTACTGGAAAACCTCCATCTTGGCTCCCAGAGCTCTAGGAACTCTTC
ATCACAACTAGATTTGCCTCTTCTAAGTGTCTATGAGCTTGCACCATATTTAATAAATTG
GGAATGGGTTTGGGGTATTAATGCAATGTGTGGTGGTTGTATTGGAGCAGGGGGAATTGA
TAAAGGAGAGTGGTTGCTGTTAATATTATCTTATCTATTGGGTGGTATGTGAAATATTGT
ACATAGACCTGATGAGTTGTGGGACCAGATGTCATCTCTGGTCAGAGTTTACTTGCTATA
TAGACTGTACTTATGTGTGAAGTTTGCAAGCTTGCTTTAGGGCTGAGCCCTGGACTCCCA
GCAGCAGCACAGTTCAGCATTGTGTGGCTGGTTGTTTCCTGGCTGTCCCCAGCAAGTGTA
GGAGTGGTGGGCCTGAACTGGGCCATTGATCAGACTAAATAAATTAAGCAGTTAACATAA
CTGGCAA"""  # noqa: E501
    assert limk2.read_text() == limk2_expected

    limk2_seguid = tmp_path / "SEGUID_LIMK2.fasta"
    test_uta_tools.get_fasta_file("ugqOFdlaed2cnxrGa7zngGMrLlY", limk2_seguid)
    limk2_seguid_expected = """>gnl|ID|ugqOFdlaed2cnxrGa7zngGMrLlY|ensembl:ENST00000331728.9|refseq:NM_005569.4|ga4gh:SQ.7_mlQyDN-uWH0RlxTQFvFEv6ykd2D-xF
GTCTTCCCGCGCCTGAGGCGGCGGCGGCAGGAGCTGAGGGGAGTTGTAGGGAACTGAGGG
GAGCTGCTGTGTCCCCCGCCTCCTCCTCCCCATTTCCGCGCTCCCGGGACCATGTCCGCG
CTGGCGGGTGAAGATGTCTGGAGGTGTCCAGGCTGTGGGGACCACATTGCTCCAAGCCAG
ATATGGTACAGGACTGTCAACGAAACCTGGCACGGCTCTTGCTTCCGGTGTTCAGAATGC
CAGGATTCCCTCACCAACTGGTACTATGAGAAGGATGGGAAGCTCTACTGCCCCAAGGAC
TACTGGGGGAAGTTTGGGGAGTTCTGTCATGGGTGCTCCCTGCTGATGACAGGGCCTTTT
ATGGTGGCTGGGGAGTTCAAGTACCACCCAGAGTGCTTTGCCTGTATGAGCTGCAAGGTG
ATCATTGAGGATGGGGATGCATATGCACTGGTGCAGCATGCCACCCTCTACTGTGGGAAG
TGCCACAATGAGGTGGTGCTGGCACCCATGTTTGAGAGACTCTCCACAGAGTCTGTTCAG
GAGCAGCTGCCCTACTCTGTCACGCTCATCTCCATGCCGGCCACCACTGAAGGCAGGCGG
GGCTTCTCCGTGTCCGTGGAGAGTGCCTGCTCCAACTACGCCACCACTGTGCAAGTGAAA
GAGGTCAACCGGATGCACATCAGTCCCAACAATCGAAACGCCATCCACCCTGGGGACCGC
ATCCTGGAGATCAATGGGACCCCCGTCCGCACACTTCGAGTGGAGGAGGTGGAGGATGCA
ATTAGCCAGACGAGCCAGACACTTCAGCTGTTGATTGAACATGACCCCGTCTCCCAACGC
CTGGACCAGCTGCGGCTGGAGGCCCGGCTCGCTCCTCACATGCAGAATGCCGGACACCCC
CACGCCCTCAGCACCCTGGACACCAAGGAGAATCTGGAGGGGACACTGAGGAGACGTTCC
CTAAGGCGCAGTAACAGTATCTCCAAGTCCCCTGGCCCCAGCTCCCCAAAGGAGCCCCTG
CTGTTCAGCCGTGACATCAGCCGCTCAGAATCCCTTCGTTGTTCCAGCAGCTATTCACAG
CAGATCTTCCGGCCCTGTGACCTAATCCATGGGGAGGTCCTGGGGAAGGGCTTCTTTGGG
CAGGCTATCAAGGTGACACACAAAGCCACGGGCAAAGTGATGGTCATGAAAGAGTTAATT
CGATGTGATGAGGAGACCCAGAAAACTTTTCTGACTGAGGTGAAAGTGATGCGCAGCCTG
GACCACCCCAATGTGCTCAAGTTCATTGGTGTGCTGTACAAGGATAAGAAGCTGAACCTC
CTGACAGAGTACATTGAGGGGGGCACACTGAAGGACTTTCTGCGCAGTATGGATCCGTTC
CCCTGGCAGCAGAAGGTCAGGTTTGCCAAAGGAATCGCCTCCGGAATGGCCTATTTGCAC
TCTATGTGCATCATCCACCGGGATCTGAACTCGCACAACTGCCTCATCAAGTTGGACAAG
ACTGTGGTGGTGGCAGACTTTGGGCTGTCACGGCTCATAGTGGAAGAGAGGAAAAGGGCC
CCCATGGAGAAGGCCACCACCAAGAAACGCACCTTGCGCAAGAACGACCGCAAGAAGCGC
TACACGGTGGTGGGAAACCCCTACTGGATGGCCCCTGAGATGCTGAACGGAAAGAGCTAT
GATGAGACGGTGGATATCTTCTCCTTTGGGATCGTTCTCTGTGAGATCATTGGGCAGGTG
TATGCAGATCCTGACTGCCTTCCCCGAACACTGGACTTTGGCCTCAACGTGAAGCTTTTC
TGGGAGAAGTTTGTTCCCACAGATTGTCCCCCGGCCTTCTTCCCGCTGGCCGCCATCTGC
TGCAGACTGGAGCCTGAGAGCAGACCAGCATTCTCGAAATTGGAGGACTCCTTTGAGGCC
CTCTCCCTGTACCTGGGGGAGCTGGGCATCCCGCTGCCTGCAGAGCTGGAGGAGTTGGAC
CACACTGTGAGCATGCAGTACGGCCTGACCCGGGACTCACCTCCCTAGCCCTGGCCCAGC
CCCCTGCAGGGGGGTGTTCTACAGCCAGCATTGCCCCTCTGTGCCCCATTCCTGCTGTGA
GCAGGGCCGTCCGGGCTTCCTGTGGATTGGCGGAATGTTTAGAAGCAGAACAAGCCATTC
CTATTACCTCCCCAGGAGGCAAGTGGGCGCAGCACCAGGGAAATGTATCTCCACAGGTTC
TGGGGCCTAGTTACTGTCTGTAAATCCAATACTTGCCTGAAAGCTGTGAAGAAGAAAAAA
ACCCCTGGCCTTTGGGCCAGGAGGAATCTGTTACTCGAATCCACCCAGGAACTCCCTGGC
AGTGGATTGTGGGAGGCTCTTGCTTACACTAATCAGCGTGACCTGGACCTGCTGGGCAGG
ATCCCAGGGTGAACCTGCCTGTGAACTCTGAAGTCACTAGTCCAGCTGGGTGCAGGAGGA
CTTCAAGTGTGTGGACGAAAGAAAGACTGATGGCTCAAAGGGTGTGAAAAAGTCAGTGAT
GCTCCCCCTTTCTACTCCAGATCCTGTCCTTCCTGGAGCAAGGTTGAGGGAGTAGGTTTT
GAAGAGTCCCTTAATATGTGGTGGAACAGGCCAGGAGTTAGAGAAAGGGCTGGCTTCTGT
TTACCTGCTCACTGGCTCTAGCCAGCCCAGGGACCACATCAATGTGAGAGGAAGCCTCCA
CCTCATGTTTTCAAACTTAATACTGGAGACTGGCTGAGAACTTACGGACAACATCCTTTC
TGTCTGAAACAAACAGTCACAAGCAAAGGAAGAGGCTGGGGGACTAGAAAGAGGCCCTGC
CCTCTAGAAAGCTCAGATCTTGGCTTCTGTTACTCATACTCGGGTGGGCTCCTTAGTCAG
ATGCCTAAAACATTTTGCCTAAAGCTCGATGGGTTCTGGAGGACAGTGTGGCTTGTCACA
GGCCTAGAGTCTGAGGGAGGGGAGTGGGAGTCTCAGCAATCTCTTGGTCTTGGCTTCATG
GCAACCACTGCTCACCCTTCAACATGCCTGGTTTAGGCAGCAGCTTGGGCTGGGAAGAGG
TGGTGGCAGAGTCTCAAAGCTGAGATGCTGAGAGAGATAGCTCCCTGAGCTGGGCCATCT
GACTTCTACCTCCCATGTTTGCTCTCCCAACTCATTAGCTCCTGGGCAGCATCCTCCTGA
GCCACATGTGCAGGTACTGGAAAACCTCCATCTTGGCTCCCAGAGCTCTAGGAACTCTTC
ATCACAACTAGATTTGCCTCTTCTAAGTGTCTATGAGCTTGCACCATATTTAATAAATTG
GGAATGGGTTTGGGGTATTAATGCAATGTGTGGTGGTTGTATTGGAGCAGGGGGAATTGA
TAAAGGAGAGTGGTTGCTGTTAATATTATCTTATCTATTGGGTGGTATGTGAAATATTGT
ACATAGACCTGATGAGTTGTGGGACCAGATGTCATCTCTGGTCAGAGTTTACTTGCTATA
TAGACTGTACTTATGTGTGAAGTTTGCAAGCTTGCTTTAGGGCTGAGCCCTGGACTCCCA
GCAGCAGCACAGTTCAGCATTGTGTGGCTGGTTGTTTCCTGGCTGTCCCCAGCAAGTGTA
GGAGTGGTGGGCCTGAACTGGGCCATTGATCAGACTAAATAAATTAAGCAGTTAACATAA
CTGGCAA"""  # noqa: E501
    assert limk2_seguid.read_text() == limk2_seguid_expected

    invalid = tmp_path / "invalid.fasta"
    with pytest.raises(KeyError):
        test_uta_tools.get_fasta_file("NM_2529.3", invalid)
