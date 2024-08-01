"""Module for testing that Cool Seq Tool works correctly."""

import copy
from datetime import datetime

import pytest

from cool_seq_tool.schemas import (
    CoordinateType,
    GenomicData,
    Strand,
    TranscriptExonData,
)


@pytest.fixture(scope="module")
def test_egc_mapper(test_cool_seq_tool):
    """Build mane ExonGenomicCoordsMapper test fixture."""
    return test_cool_seq_tool.ex_g_coords_mapper


@pytest.fixture(scope="module")
def nm_152263_exons_genomic_coords():
    """Create test fixture for NM_152263.4 exons and genomic coordinates."""
    return [
        (0, 0, 199, 154191901, 154192100, -1),
        (1, 199, 325, 154191185, 154191311, -1),
        (2, 325, 459, 154176114, 154176248, -1),
        (3, 459, 577, 154173083, 154173201, -1),
        (4, 577, 648, 154172907, 154172978, -1),
        (5, 648, 724, 154171412, 154171488, -1),
        (6, 724, 787, 154170648, 154170711, -1),
        (7, 787, 857, 154170399, 154170469, -1),
        (8, 857, 936, 154169304, 154169383, -1),
        (9, 936, 7064, 154161812, 154167940, -1),
    ]


@pytest.fixture(scope="module")
def nm_001105539_exons_genomic_coords():
    """Create test fixture for NM_001105539.3 exons and genomic coordinates."""
    return [
        (0, 0, 1557, 80486225, 80487782, -1),
        (1, 1557, 2446, 80499493, 80500382, -1),
        (2, 2446, 2545, 80513909, 80514008, -1),
        (3, 2545, 2722, 80518402, 80518579, -1),
        (4, 2722, 2895, 80518781, 80518954, -1),
        (5, 2895, 9938, 80519222, 80526265, -1),
    ]


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
        "strand": Strand.NEGATIVE,
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
        "strand": Strand.NEGATIVE,
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
        "strand": Strand.NEGATIVE,
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
        "strand": Strand.NEGATIVE,
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
        "strand": Strand.NEGATIVE,
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
        "strand": Strand.NEGATIVE,
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
        "strand": Strand.NEGATIVE,
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
        "strand": Strand.POSITIVE,
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
        "strand": Strand.POSITIVE,
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
        "strand": Strand.POSITIVE,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def zbtb10_exon3_end():
    """Create test fixture for ZBTB10, end of exon 3"""
    params = {
        "gene": "ZBTB10",
        "chr": "NC_000008.11",
        "start": None,
        "end": 80514009,
        "exon_start": None,
        "exon_end": 3,
        "exon_end_offset": 100,
        "exon_start_offset": None,
        "transcript": "NM_001105539.3",
        "strand": Strand.POSITIVE,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def zbtb10_exon5_start():
    """Create test fixture for ZBTB10, start of exon 5"""
    params = {
        "gene": "ZBTB10",
        "chr": "NC_000008.11",
        "start": 80518580,
        "end": None,
        "exon_start": 5,
        "exon_start_offset": -374,
        "exon_end": None,
        "exon_end_offset": None,
        "transcript": "NM_001105539.3",
        "strand": Strand.POSITIVE,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon6_end():
    """Create test fixture for TPM3, end of exon 6"""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": None,
        "end": 154171409,
        "exon_start": None,
        "exon_start_offset": None,
        "exon_end": 6,
        "exon_end_offset": 3,
        "transcript": "NM_152263.4",
        "strand": Strand.NEGATIVE,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def tpm3_exon5_start():
    """Create test fixture for TPM3, start of exon 5"""
    params = {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154173080,
        "end": None,
        "exon_start": 5,
        "exon_start_offset": -102,
        "exon_end": None,
        "exon_end_offset": None,
        "transcript": "NM_152263.4",
        "strand": Strand.NEGATIVE,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def gusbp3_exon2_end():
    """Create test fixture for GUSBP3, end of exon 2"""
    params = {
        "gene": "GUSBP3",
        "chr": "NC_000005.10",
        "start": None,
        "end": 69680763,
        "exon_start": None,
        "exon_start_offset": None,
        "exon_end": 2,
        "exon_end_offset": 2,
        "transcript": "NR_027386.2",
        "strand": Strand.NEGATIVE,
    }
    return GenomicData(**params)


@pytest.fixture(scope="module")
def gusbp3_exon5_start():
    """Create test fixture for GUSBP3, start of exon 5"""
    params = {
        "gene": "GUSBP3",
        "chr": "NC_000005.10",
        "start": 69645878,
        "end": None,
        "exon_start": 5,
        "exon_start_offset": -3589,
        "exon_end": None,
        "exon_end_offset": None,
        "transcript": "NR_027386.2",
        "strand": Strand.NEGATIVE,
    }
    return GenomicData(**params)


def check_service_meta(actual):
    """Check that service metadata matches expected

    :param ServiceMeta actual: Actual service metadata
    """
    assert actual.name == "cool_seq_tool"
    # pydantic checks that version is a string
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


@pytest.mark.asyncio()
async def test_get_tx_exon_coords(test_egc_mapper, nm_152263_exons):
    """Test that get_tx_exon_coords works correctly."""
    resp = test_egc_mapper.get_tx_exon_coords("NM_152263.3", nm_152263_exons, 1, 8)
    assert resp[0] == ((0, 234), (822, 892))
    assert resp[1] is None

    resp = test_egc_mapper.get_tx_exon_coords("NM_152263.3", nm_152263_exons, 1, 11)
    assert resp[0] is None
    assert resp[1] == "Exon 11 does not exist on NM_152263.3"


@pytest.mark.asyncio()
async def test_get_adjacent_exon(
    test_egc_mapper, nm_152263_exons_genomic_coords, nm_001105539_exons_genomic_coords
):
    """Test that get_adjacent_exon works properly"""
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_152263_exons_genomic_coords,
        end=154191901,
        strand=Strand.NEGATIVE,
    )
    assert resp == 1
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_152263_exons_genomic_coords,
        end=154191184,
        strand=Strand.NEGATIVE,
    )
    assert resp == 2
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_152263_exons_genomic_coords,
        start=154191184,
        strand=Strand.NEGATIVE,
    )
    assert resp == 3
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_001105539_exons_genomic_coords,
        end=80500385,
        strand=Strand.POSITIVE,
    )
    assert resp == 2
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_001105539_exons_genomic_coords,
        start=80518580,
        strand=Strand.POSITIVE,
    )
    assert resp == 5


def test_is_exonic_breakpoint(test_egc_mapper, nm_001105539_exons_genomic_coords):
    """Test is breakpoint occurs on exon"""
    resp = test_egc_mapper._is_exonic_breakpoint(
        80514010, nm_001105539_exons_genomic_coords
    )
    assert resp is False  # Breakpoint does not occur on an exon

    resp = test_egc_mapper._is_exonic_breakpoint(
        80499495, nm_001105539_exons_genomic_coords
    )
    assert resp is True  # Breakpoint does occur on an exon


@pytest.mark.asyncio()
async def test_genomic_to_transcript_fusion_context(
    test_egc_mapper,
    zbtb10_exon3_end,
    zbtb10_exon5_start,
    tpm3_exon6_end,
    tpm3_exon5_start,
    gusbp3_exon2_end,
    gusbp3_exon5_start,
):
    """Test that genomic to transcript works correctly for non-exonic breakpoints"""
    inputs = {
        "chromosome": "8",
        "seg_end": 80514010,
        "gene": "ZBTB10",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, zbtb10_exon3_end)

    inputs = {
        "chromosome": "chr8",
        "seg_end": 80514010,
        "gene": "ZBTB10",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, zbtb10_exon3_end)

    inputs = {
        "chromosome": "8",
        "seg_start": 80518581,
        "gene": "ZBTB10",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, zbtb10_exon5_start)

    inputs = {
        "chromosome": "1",
        "seg_end": 154171410,
        "gene": "TPM3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, tpm3_exon6_end)

    inputs = {
        "chromosome": "1",
        "seg_start": 154173081,
        "gene": "TPM3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, tpm3_exon5_start)

    inputs = {
        "chromosome": "5",
        "seg_end": 69680764,
        "gene": "GUSBP3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, gusbp3_exon2_end)

    inputs = {
        "chromosome": "5",
        "seg_start": 69645879,
        "gene": "GUSBP3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, gusbp3_exon5_start)

    inputs = {  # Test when gene and strand are not provided
        "chromosome": "5",
        "seg_start": 69645879,
        "transcript": "NR_027386.2",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert (
        resp.warnings[0]
        == "Gene must be provided to select the adjacent transcript junction"
    )

    inputs = {  # Test when transcript is provided
        "chromosome": "5",
        "seg_start": 69645879,
        "gene": "GUSBP3",
        "transcript": "NR_027386.2",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(resp, gusbp3_exon5_start)


@pytest.mark.asyncio()
async def test_get_alt_ac_start_and_end(
    test_egc_mapper, tpm3_1_8_start_genomic, tpm3_1_8_end_genomic
):
    """Test that _get_alt_ac_start_and_end works correctly."""
    resp = await test_egc_mapper._get_alt_ac_start_and_end(
        "NM_152263.3", ["117", "234"], ["822", "892"], "TPM3"
    )
    assert resp[0] == (tpm3_1_8_start_genomic, tpm3_1_8_end_genomic)
    assert resp[1] is None

    resp = await test_egc_mapper._get_alt_ac_start_and_end("NM_152263.3", gene="TPM3")
    assert resp[0] is None
    assert resp[1] == "Must provide either `tx_exon_start` or `tx_exon_end` or both"


@pytest.mark.asyncio()
async def test_get_tx_exons_genomic_coords(
    test_egc_mapper, nm_152263_exons_genomic_coords
):
    """Test that _get_tx_exons_genomic_coords works correctly."""
    resp = await test_egc_mapper._get_tx_exons_genomic_coords(
        "NM_152263.4", "NC_000001.11"
    )
    assert resp[0] == nm_152263_exons_genomic_coords
    assert resp[1] is None

    # Invalid transcript accession given chromosome accession
    resp = await test_egc_mapper._get_tx_exons_genomic_coords(
        "NM_001105539.3", "NC_000001.11"
    )
    assert resp[0] is None
    assert (
        resp[1]
        == "Unable to get exons and genomic coordinates for NM_001105539.3 on NC_000001.11"
    )


@pytest.mark.asyncio()
async def test_genomic_to_transcript(test_egc_mapper, tpm3_exon1, tpm3_exon8):
    """Test that _genomic_to_transcript_exon_coordinate
    method works correctly.
    """
    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        154192134,
        alt_ac="NC_000001.11",
        transcript="NM_152263.3",
        gene="TPM3",
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        154192134, chromosome="1", transcript="NM_152263.3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        154192134, chromosome="1", transcript="NM_152263.3"
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        154170399,
        alt_ac="NC_000001.11",
        transcript="NM_152263.3",
        is_start=False,
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        154170399,
        chromosome="1",
        transcript="NM_152263.3",
        is_start=False,
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)

    resp = await test_egc_mapper._genomic_to_transcript_exon_coordinate(
        154170399, chromosome="1", transcript="NM_152263.3", is_start=False
    )
    transcript_exon_data_assertion_checks(resp, tpm3_exon8)


@pytest.mark.asyncio()
async def test_tpm3(
    test_egc_mapper,
    tpm3_exon1_exon8,
    tpm3_exon1_exon8_offset,
    tpm3_exon1_g,
    tpm3_exon8_g,
):
    """Test TPM3 genomic_to_tx_segment and
    tx_segment_to_genomic.
    """
    inputs = {
        "alt_ac": "NC_000001.11",
        "seg_start": 154192135,
        "seg_end": 154170400,
        "transcript": "NM_152263.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8)

    inputs = {
        "alt_ac": "NC_000001.11",
        "seg_start": 154192134,
        "seg_end": 154170400,
        "transcript": "NM_152263.3",
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8)

    # Offset
    inputs = {
        "alt_ac": "NC_000001.11",
        "seg_start": 154192132,
        "seg_end": 154170404,
        "transcript": "NM_152263.3",
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_exon8_offset)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_exon8_offset)

    # Test only setting start
    inputs = {
        "alt_ac": "NC_000001.11",
        "seg_start": 154192134,
        "transcript": "NM_152263.3",
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon1_g)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon1_g)

    # Test only setting end
    inputs = {
        "alt_ac": "NC_000001.11",
        "seg_end": 154170400,
        "transcript": "NM_152263.3",
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, tpm3_exon8_g)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, tpm3_exon8_g)


@pytest.mark.asyncio()
async def test_braf(test_egc_mapper, mane_braf):
    """Test BRAF genomic_to_tx_segment and
    tx_segment_to_genomic.
    """
    inputs = {
        "alt_ac": "NC_000007.13",
        "seg_start": 140501360,
        "seg_end": 140453136,
        "gene": "BRAF",
    }
    # MANE
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_braf)

    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, mane_braf)


@pytest.mark.asyncio()
async def test_wee1(test_egc_mapper, wee1_exon2_exon11, mane_wee1_exon2_exon11):
    """Test WEE1 genomic_to_tx_segment and
    tx_segment_to_genomic.
    """
    inputs = {
        "alt_ac": "NC_000011.9",
        "seg_start": 9597640,
        "seg_end": 9609996,
        "transcript": "NM_003390.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, wee1_exon2_exon11)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, wee1_exon2_exon11)

    # add gene
    inputs = {
        "alt_ac": "NC_000011.9",
        "seg_start": 9597640,
        "seg_end": 9609996,
        "transcript": "NM_003390.3",
        "gene": "wee1",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, wee1_exon2_exon11)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, wee1_exon2_exon11)

    # MANE since no transcript provided
    inputs = {
        "alt_ac": "NC_000011.9",
        "seg_start": 9597640,
        "seg_end": 9609996,
        "gene": "wee1",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_data_assertion_checks(g_to_t_resp, mane_wee1_exon2_exon11)
    params = g_to_t_resp.genomic_data
    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        params.transcript,
        gene=params.gene,
        exon_start=params.exon_start,
        exon_start_offset=params.exon_start_offset,
        exon_end=params.exon_end,
        exon_end_offset=params.exon_end_offset,
    )
    genomic_data_assertion_checks(t_to_g_resp, mane_wee1_exon2_exon11)


@pytest.mark.asyncio()
async def test_transcript_to_genomic(
    test_egc_mapper,
    tpm3_exon1_g,
    tpm3_exon8_g,
    tpm3_exon1_exon8,
    ntrk1_exon10_exon17,
):
    """Test that tx_segment_to_genomic works correctly."""
    # TPM3
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None, exon_end=8, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=None, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon1_g)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None, exon_end=8, transcript="NM_152263.3       "
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None, exon_end=8, gene="TPM3", transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None, exon_end=8, gene=" TPM3 ", transcript=" NM_152263.3 "
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None, exon_end=8, gene="tpm3", transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, tpm3_exon8_g)

    expected = copy.deepcopy(tpm3_exon1_exon8)
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=8, exon_end_offset=-5, transcript="NM_152263.3"
    )
    expected.exon_end = 8
    expected.exon_end_offset = -5
    expected.end = 154170405
    genomic_data_assertion_checks(resp, expected)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=8, exon_end_offset=5, transcript="NM_152263.3"
    )
    expected.exon_end_offset = 5
    expected.end = 154170395
    genomic_data_assertion_checks(resp, expected)

    resp = await test_egc_mapper.tx_segment_to_genomic(
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

    resp = await test_egc_mapper.tx_segment_to_genomic(
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
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, gene="NTRK1", transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, gene="NTRK1", transcript="NM_002529.3"
    )
    genomic_data_assertion_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, exon_start_offset=3, transcript="NM_002529.3"
    )
    expected = copy.deepcopy(ntrk1_exon10_exon17)
    expected.exon_start_offset = 3
    expected.start = 156874628
    genomic_data_assertion_checks(resp, expected)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, exon_start_offset=-3, transcript="NM_002529.3"
    )
    expected.exon_start_offset = -3
    expected.start = 156874622
    genomic_data_assertion_checks(resp, expected)


@pytest.mark.asyncio()
async def test_valid_inputs(test_egc_mapper):
    """Test that valid inputs don"t return any errors"""
    inputs = {"gene": "TPM3", "alt_ac": "NC_000001.11", "seg_start": 154171413}
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert resp.genomic_data

    inputs = {"gene": "WEE1", "alt_ac": "NC_000011.9", "seg_end": 9609996}
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert resp.genomic_data

    inputs = {"gene": "WEE1", "chromosome": "11", "seg_end": 9609996}
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert resp.genomic_data

    inputs = {"transcript": "NM_003390.3", "exon_start": 2}
    resp = await test_egc_mapper.tx_segment_to_genomic(**inputs)
    assert resp.genomic_data

    # add gene
    inputs = {"transcript": "NM_003390.3", "exon_start": 2, "gene": "WEE1"}
    resp = await test_egc_mapper.tx_segment_to_genomic(**inputs)
    assert resp.genomic_data

    # Test X/Y chromosome bug
    inputs = {
        "chromosome": "X",
        "seg_start": 154437254,
        "seg_end": 154437299,
        "gene": "GDI1",
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert resp.genomic_data

    resp = await test_egc_mapper.tx_segment_to_genomic(
        gene="PDGFRB", transcript="NM_002609.4", exon_start=11, exon_end=23
    )
    assert resp.genomic_data


@pytest.mark.asyncio()
async def test_invalid(test_egc_mapper):
    """Test that invalid queries return `None`."""
    resp = await test_egc_mapper.genomic_to_tx_segment(
        transcript="NM_152263 3",
        seg_start=154192135,
        seg_end=154170399,
        alt_ac="NC_000001.11",
    )
    assert resp.warnings == ["Unable to get exons for NM_152263 3"]

    # start and end not given
    resp = await test_egc_mapper.genomic_to_tx_segment(
        alt_ac="NC_000001.11",
        seg_start=None,
        seg_end=None,
        transcript="NM_152263.3",
        gene="TPM3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `seg_start` or `seg_end`"]

    # Invalid gene
    resp = await test_egc_mapper.genomic_to_tx_segment(
        alt_ac="NC_000001.11",
        seg_start=154192135,
        seg_end=154170399,
        transcript="NM_152263.3",
        gene="dummy gene",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result for chromosome NC_000001.11 "
        "where genomic coordinate 154192134 is mapped between an "
        "exon's start and end coordinates and on gene DUMMY GENE"
    ]

    # Invalid accession
    resp = await test_egc_mapper.genomic_to_tx_segment(
        alt_ac="NC_000001.200",
        seg_start=154192135,
        seg_end=154170399,
        transcript="NM_152263.3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Invalid genomic accession: NC_000001.200"]

    # Invalid coordinates
    resp = await test_egc_mapper.genomic_to_tx_segment(
        alt_ac="NC_000001.11",
        seg_start=9999999999999,
        seg_end=9999999999999,
        transcript="NM_152263.3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result for chromosome NC_000001.11 where genomic "
        "coordinate 9999999999998 is mapped between an exon's start and end coordinates"
    ]

    resp = await test_egc_mapper.genomic_to_tx_segment(
        chromosome="1",
        seg_start=154170400,
        transcript="NM_002529.3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must find exactly one row for genomic data, but found: 0"]

    # Must supply either gene or transcript
    resp = await test_egc_mapper.genomic_to_tx_segment(
        seg_start=154192135, alt_ac="NC_000001.11"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `gene` or `transcript`"]

    # Exon 22 does not exist
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None,
        exon_end=22,
        transcript="NM_152263.3",
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Exon 22 does not exist on NM_152263.3"]

    # Start > End
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=8, exon_end=1, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Start exon 8 is greater than end exon 1"]

    # Transcript DNE
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=7, exon_end=None, transcript="NM_12345.6"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Unable to get exons for NM_12345.6"]

    # Index error for invalid exon
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=-1, exon_end=0, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "`exon_start` cannot be less than 1",
        "`exon_end` cannot be less than 1",
    ]

    # Cant supply 0 based exons
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_end=0, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["`exon_end` cannot be less than 1"]

    # Gene that does not match transcript
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=8, gene="NTKR1", transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == [
        "Unable to find a result where NM_152263.3 has transcript coordinates"
        " 0 and 234 between an exon's start and end coordinates on gene "
        "NTKR1"
    ]

    # No transcript given
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=8, gene="NTKR1", transcript=""
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide `transcript`"]

    # No exons given
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None, exon_end=None, transcript="NM_152263.3"
    )
    genomic_data_assertion_checks(resp, is_valid=False)
    assert resp.warnings == ["Must provide either `exon_start` or `exon_end`"]
