"""Module for testing that Cool Seq Tool works correctly."""
import copy
import re
from datetime import datetime

import pytest

from cool_seq_tool.schemas import GenomicData, ResidueMode, TranscriptExonData


@pytest.fixture(scope="module")
def test_egc_mapper(test_cool_seq_tool):
    """Build mane ExonGenomicCoordsMapper test fixture."""
    return test_cool_seq_tool.ex_g_coords_mapper


@pytest.fixture(scope="module")
def tpm3_exon1():
    """Create test fixture for TPM3 exon 1."""
    params = {
        "chr": "NC_000001.11",
        "gene": "TPM3",
        "pos": 154192134,
        "exon": 1,
        "exon_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1,
    }
    return TranscriptExonData(**params)


@pytest.fixture(scope="module")
def tpm3_exon8():
    """Create test fixture for TPM3 exon 8."""
    params = {
        "chr": "NC_000001.11",
        "gene": "TPM3",
        "pos": 154170400,
        "exon": 8,
        "exon_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1,
    }
    return TranscriptExonData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_g():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192134,
        "end": None,
        "exon_start": 1,
        "exon_end": None,
        "exon_start_offset": 0,
        "exon_end_offset": None,
        "transcript": "NM_152263.3",
        "strand": -1,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon8_g():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": None,
        "end": 154170400,
        "exon_start": None,
        "exon_end": 8,
        "exon_start_offset": None,
        "exon_end_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_exon8():
    """Create test fixture for TPM3."""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192134,
        "end": 154170400,
        "exon_start": 1,
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start_offset": 0,
        "transcript": "NM_152263.3",
        "strand": -1,
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
        "exon_start_offset": 2,
        "end": 154170404,
        "exon_end": 8,
        "exon_end_offset": -4,
        "transcript": "NM_152263.3",
        "strand": -1,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def mane_braf():
    """Create test fixture for BRAF."""
    params = {
        "chr": "NC_000007.14",
        "gene": "BRAF",
        "start": 140808061,
        "exon_start": 5,
        "exon_start_offset": 0,
        "end": 140753332,
        "exon_end": 15,
        "exon_end_offset": -57,
        "transcript": "NM_004333.6",
        "strand": -1,
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
        "end": 9588449,
        "exon_end": 11,
        "exon_end_offset": 0,
        "transcript": "NM_003390.3",
        "strand": 1,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def mane_wee1_exon2_exon11():
    """Create test fixture for WEE1."""
    params = {
        "chr": "NC_000011.10",
        "gene": "WEE1",
        "start": 9576091,
        "exon_start": 2,
        "exon_start_offset": -1,
        "end": 9586857,
        "exon_end": 10,
        "exon_end_offset": 146,
        "transcript": "NM_003390.4",
        "strand": 1,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def ntrk1_exon10_exon17():
    """Create test fixture for NTRK1."""
    params = {
        "gene": "NTRK1",
        "chr": "NC_000001.11",
        "start": 156874625,
        "end": 156881457,
        "exon_start": 10,
        "exon_end": 17,
        "exon_end_offset": 0,
        "exon_start_offset": 0,
        "transcript": "NM_002529.3",
        "strand": 1,
    }
    return GenomicData(**params)


def check_service_meta(actual):
    """Check that service metadata matches expected

    :param ServiceMeta actual: Actual service metadata
    """
    assert actual.name == "cool_seq_tool"
    version_regex = r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"  # noqa: E501
    assert bool(re.match(version_regex, actual.version))
    assert isinstance(actual.response_datetime, datetime)
    assert actual.url == "https://github.com/GenomicMedLab/cool-seq-tool"


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


def transcript_exon_data_assertion_checks(actual, expected=None, is_valid=True):
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
async def test_genomic_to_transcript(test_egc_mapper, tpm3_exon1, tpm3_exon8):
    """Test that _genomic_to_transcript_exon_coordinate
    method works correctly.
    """
    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154192134, strand=-1, transcript="NM_152263.3", gene="TPM3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        1, 154192134, strand=-1, transcript="NM_152263.3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        1, 154192134, transcript="NM_152263.3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154170399, strand=-1, transcript="NM_152263.3", is_start=False
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        1, 154170399, strand=-1, transcript="NM_152263.3", is_start=False
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        1, 154170399, transcript="NM_152263.3", is_start=False
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)


@pytest.mark.asyncio
async def test_tpm3(
    test_egc_mapper,
    tpm3_exon1_exon8,
    tpm3_exon1_exon8_offset,
    tpm3_exon1_g,
    tpm3_exon8_g,
):
    """Test TPM3 genomic_to_transcript_exon_coordinates and
    transcript_to_genomic_coordinates.
    """
    inputs = {
        "chromosome": "NC_000001.11",
        "start": 154192135,
        "end": 154170400,
        "strand": -1,
        "transcript": "NM_152263.3",
    }

    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8)

    inputs["residue_mode"] = ResidueMode.INTER_RESIDUE
    inputs["start"] = 154192134
    inputs["end"] = 154170400
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8)

    # No strand
    del inputs["strand"]
    del inputs["residue_mode"]
    inputs["start"] = 154192135
    inputs["end"] = 154170400
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8)

    # Offset, no strand
    inputs["start"] = 154192132
    inputs["end"] = 154170404
    inputs["residue_mode"] = ResidueMode.INTER_RESIDUE
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8_offset)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_offset)

    # Offset, strand
    inputs["strand"] = -1
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8_offset)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_offset)

    # Test only setting start
    inputs = {
        "chromosome": "NC_000001.11",
        "start": 154192134,
        "strand": -1,
        "transcript": "NM_152263.3",
        "residue_mode": ResidueMode.INTER_RESIDUE,
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_g)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_g)

    # Test only setting end
    inputs = {
        "chromosome": "NC_000001.11",
        "end": 154170400,
        "strand": -1,
        "transcript": "NM_152263.3",
        "residue_mode": ResidueMode.INTER_RESIDUE,
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon8_g)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon8_g)


@pytest.mark.asyncio
async def test_braf(test_egc_mapper, mane_braf):
    """Test BRAF genomic_to_transcript_exon_coordinates and
    transcript_to_genomic_coordinates.
    """
    inputs = {
        "chromosome": "NC_000007.13",
        "start": 140501360,
        "end": 140453136,
        "strand": -1,
        "gene": "BRAF",
    }
    # MANE
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_braf)

    del inputs["strand"]
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_braf)

    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, mane_braf)


@pytest.mark.asyncio
async def test_wee1(test_egc_mapper, wee1_exon2_exon11, mane_wee1_exon2_exon11):
    """Test WEE1 genomic_to_transcript_exon_coordinates and
    transcript_to_genomic_coordinates.
    """
    inputs = {
        "chromosome": "NC_000011.9",
        "start": 9597640,
        "end": 9609996,
        "strand": 1,
        "transcript": "NM_003390.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, wee1_exon2_exon11)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, wee1_exon2_exon11)

    inputs["gene"] = "wee1"
    del inputs["strand"]
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, wee1_exon2_exon11)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, wee1_exon2_exon11)

    # MANE
    del inputs["transcript"]
    g_to_t_resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_wee1_exon2_exon11)
    t_to_g_resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        **g_to_t_resp.genomic_data.model_dump()
    )
    genomic_data_assertion_checks(t_to_g_resp, mane_wee1_exon2_exon11)


@pytest.mark.asyncio
async def test_transcript_to_genomic(
    test_egc_mapper,
    tpm3_exon1_g,
    tpm3_exon8_g,
    tpm3_exon1_exon8,
    ntrk1_exon10_exon17,
):
    """Test that transcript_to_genomic_coordinates works correctly."""
    # TPM3
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=None, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon1_g)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, transcript="NM_152263.3       "
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, gene="TPM3", transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, gene=" TPM3 ", transcript=" NM_152263.3 "
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=8, gene="tpm3", transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    expected = copy.deepcopy(tpm3_exon1_exon8)
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, exon_end_offset=-5, transcript="NM_152263.3"
    )
    expected.exon_end = 8
    expected.exon_end_offset = -5
    expected.end = 154170405
    genomic_data_assertion_checks(resp, expected)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, exon_end_offset=5, transcript="NM_152263.3"
    )
    expected.exon_end_offset = 5
    expected.end = 154170395
    genomic_data_assertion_checks(resp, expected)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=3,
        exon_end=8,
        exon_start_offset=3,
        exon_end_offset=5,
        transcript="NM_152263.3",
    )
    expected.exon_start = 3
    expected.exon_start_offset = 3
    expected.start = 154176244
    genomic_data_assertion_checks(resp, expected)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=3,
        exon_end=8,
        exon_start_offset=-3,
        exon_end_offset=5,
        transcript="NM_152263.3",
    )
    expected.exon_start_offset = -3
    expected.start = 154176250
    genomic_data_assertion_checks(resp, expected)

    # NTRK1
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, gene="NTRK1", transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, gene="NTRK1", transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, exon_start_offset=3, transcript="NM_002529.3"
    )
    expected = copy.deepcopy(ntrk1_exon10_exon17)
    expected.exon_start_offset = 3
    expected.start = 156874628
    genomic_data_assertion_checks(resp, expected)

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=10, exon_end=17, exon_start_offset=-3, transcript="NM_002529.3"
    )
    expected.exon_start_offset = -3
    expected.start = 156874622
    genomic_data_assertion_checks(resp, expected)


@pytest.mark.asyncio
async def test_valid_inputs(test_egc_mapper):
    """Test that valid inputs don"t return any errors"""
    inputs = {"gene": "TPM3", "chromosome": "NC_000001.11", "start": 154171413}
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    assert resp.genomic_data

    inputs = {"gene": "WEE1", "chromosome": "NC_000011.9", "end": 9609996}
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    assert resp.genomic_data

    inputs["chromosome"] = "11"
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    assert resp.genomic_data

    inputs = {"transcript": "NM_003390.3", "exon_start": 2}
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(**inputs)
    assert resp.genomic_data

    inputs["gene"] = "WEE1"
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(**inputs)
    assert resp.genomic_data

    # Test X/Y chromosome bug
    inputs = {
        "chromosome": "X",
        "strand": 1,
        "start": 154437254,
        "end": 154437299,
        "gene": "GDI1",
        "residue_mode": ResidueMode.INTER_RESIDUE,
    }
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(**inputs)
    assert resp.genomic_data

    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        gene="PDGFRB", transcript="NM_002609.4", exon_start=11, exon_end=23
    )
    assert resp.genomic_data


@pytest.mark.asyncio
async def test_invalid(test_egc_mapper):
    """Test that invalid queries return `None`."""
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(
        transcript="NM_152263 3",
        start=154192135,
        end=154170399,
        chromosome="NC_000001.11",
    )
    assert resp.warnings == ["Unable to get exons for NM_152263 3"]

    # start and end not given
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(
        "NC_000001.11",
        start=None,
        end=None,
        strand=-1,
        transcript="NM_152263.3",
        gene="TPM3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `start` or `end`"]

    # Invalid gene
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(
        "NC_000001.11",
        start=154192135,
        end=154170399,
        strand=-1,
        transcript="NM_152263.3",
        gene="dummy gene",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result for chromosome NC_000001.11 "
        "where genomic coordinate 154192134 is mapped between an "
        "exon's start and end coordinates on the negative strand "
        "and on gene DUMMY GENE"
    ]

    # Invalid chromosome
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(
        "NC_000001.200",
        start=154192135,
        end=154170399,
        strand=-1,
        transcript="NM_152263.3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Invalid chromosome: NC_000001.200"]

    # Invalid coordinates
    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(
        "NC_000001.11",
        start=9999999999999,
        end=9999999999999,
        strand=-1,
        transcript="NM_152263.3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result for chromosome NC_000001.11 where genomic "
        "coordinate 9999999999998 is mapped between an exon's start and end "
        "coordinates on the negative strand"
    ]

    resp = await test_egc_mapper.genomic_to_transcript_exon_coordinates(
        chromosome="1", start=154170400, strand=-1, transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must find exactly one row for genomic data, but found: 0"]

    # Strand does not match
    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154192135, strand=1, transcript="NM_152263.3", gene="TPM3"
    )
    transcript_exon_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result for chromosome NC_000001.11 where genomic "
        "coordinate 154192135 is mapped between an exon's start and end "
        "coordinates on the positive strand and on gene TPM3"
    ]

    # Must supply either gene or transcript
    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        "NC_000001.11", 154192135, strand=1
    )
    transcript_exon_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `gene` or `transcript`"]

    # Exon 22 does not exist
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=None,
        exon_end=22,
        transcript="NM_152263.3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Exon 22 does not exist on NM_152263.3"]

    # Start > End
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=8, exon_end=1, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Start exon 8 is greater than end exon 1"]

    # Transcript DNE
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=7, exon_end=None, transcript="NM_12345.6"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Unable to get exons for NM_12345.6"]

    # Index error for invalid exon
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=-1, exon_end=0, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Exon -1 does not exist on NM_152263.3"]

    # Cant supply 0 based exons
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=0, exon_end=1, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Exon 0 does not exist on NM_152263.3"]

    # Gene that does not match transcript
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, gene="NTKR1", transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result where NM_152263.3 has transcript coordinates"
        " 0 and 234 between an exon's start and end coordinates on gene "
        "NTKR1"
    ]

    # No transcript given
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, gene="NTKR1", transcript=None
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide `transcript`"]

    # No transcript given
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=1, exon_end=8, gene="NTKR1", transcript=""
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide `transcript`"]

    # No exons given
    resp = await test_egc_mapper.transcript_to_genomic_coordinates(
        exon_start=None, exon_end=None, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `exon_start` or `exon_end`"]
