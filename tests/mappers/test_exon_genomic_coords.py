"""Module for testing that Cool Seq Tool works correctly."""

from datetime import datetime

import pytest

from cool_seq_tool.mappers.exon_genomic_coords import (
    ExonCoord,
    GenomicTxSegService,
    _GenomicTxSeg,
)
from cool_seq_tool.schemas import (
    Strand,
)


@pytest.fixture(scope="module")
def test_egc_mapper(test_cool_seq_tool):
    """Build mane ExonGenomicCoordsMapper test fixture."""
    return test_cool_seq_tool.ex_g_coords_mapper


@pytest.fixture(scope="module")
def nm_152263_exons_genomic_coords():
    """Create test fixture for NM_152263.4 exons and genomic coordinates."""
    return [
        ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=199,
            alt_start_i=154191901,
            alt_end_i=154192100,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=1,
            tx_start_i=199,
            tx_end_i=325,
            alt_start_i=154191185,
            alt_end_i=154191311,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=2,
            tx_start_i=325,
            tx_end_i=459,
            alt_start_i=154176114,
            alt_end_i=154176248,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=3,
            tx_start_i=459,
            tx_end_i=577,
            alt_start_i=154173083,
            alt_end_i=154173201,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=4,
            tx_start_i=577,
            tx_end_i=648,
            alt_start_i=154172907,
            alt_end_i=154172978,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=5,
            tx_start_i=648,
            tx_end_i=724,
            alt_start_i=154171412,
            alt_end_i=154171488,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=6,
            tx_start_i=724,
            tx_end_i=787,
            alt_start_i=154170648,
            alt_end_i=154170711,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=7,
            tx_start_i=787,
            tx_end_i=857,
            alt_start_i=154170399,
            alt_end_i=154170469,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=8,
            tx_start_i=857,
            tx_end_i=936,
            alt_start_i=154169304,
            alt_end_i=154169383,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=9,
            tx_start_i=936,
            tx_end_i=7064,
            alt_start_i=154161812,
            alt_end_i=154167940,
            alt_strand=Strand.NEGATIVE,
        ),
    ]


@pytest.fixture(scope="module")
def nm_001105539_exons_genomic_coords():
    """Create test fixture for NM_001105539.3 exons and genomic coordinates."""
    return [
        ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=1557,
            alt_start_i=80486225,
            alt_end_i=80487782,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=1,
            tx_start_i=1557,
            tx_end_i=2446,
            alt_start_i=80499493,
            alt_end_i=80500382,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=2,
            tx_start_i=2446,
            tx_end_i=2545,
            alt_start_i=80513909,
            alt_end_i=80514008,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=3,
            tx_start_i=2545,
            tx_end_i=2722,
            alt_start_i=80518402,
            alt_end_i=80518579,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=4,
            tx_start_i=2722,
            tx_end_i=2895,
            alt_start_i=80518781,
            alt_end_i=80518954,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=5,
            tx_start_i=2895,
            tx_end_i=9938,
            alt_start_i=80519222,
            alt_end_i=80526265,
            alt_strand=Strand.NEGATIVE,
        ),
    ]


@pytest.fixture(scope="module")
def tpm3_exon1():
    """Create test fixture for TPM3 exon 1 (negative strand)."""
    params = {
        "gene": "TPM3",
        "genomic_ac": "NC_000001.11",
        "tx_ac": "NM_152263.3",
        "seg": {
            "exon_ord": 0,
            "offset": 0,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "end": 154191901,
            },
        },
    }
    return _GenomicTxSeg(**params)


@pytest.fixture(scope="module")
def tpm3_exon8():
    """Create test fixture for TPM3 exon 8 (negative strand)."""
    params = {
        "gene": "TPM3",
        "genomic_ac": "NC_000001.11",
        "tx_ac": "NM_152263.3",
        "seg": {
            "exon_ord": 7,
            "offset": 0,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 154170469,
            },
        },
    }
    return _GenomicTxSeg(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_g(tpm3_exon1):
    """Create test fixture for TPM3."""
    params = {
        "gene": tpm3_exon1.gene,
        "genomic_ac": tpm3_exon1.genomic_ac,
        "tx_ac": tpm3_exon1.tx_ac,
        "seg_start": tpm3_exon1.seg,
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def tpm3_exon8_g(tpm3_exon8):
    """Create test fixture for TPM3."""
    params = {
        "gene": tpm3_exon8.gene,
        "genomic_ac": tpm3_exon8.genomic_ac,
        "tx_ac": tpm3_exon8.tx_ac,
        "seg_end": tpm3_exon8.seg,
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_exon8(tpm3_exon1, tpm3_exon8):
    """Create test fixture for TPM3."""
    params = {
        "gene": tpm3_exon8.gene,
        "genomic_ac": tpm3_exon8.genomic_ac,
        "tx_ac": tpm3_exon8.tx_ac,
        "seg_start": tpm3_exon1.seg,
        "seg_end": tpm3_exon8.seg,
    }

    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def tpm3_exon1_exon8_offset(tpm3_exon1, tpm3_exon8):
    """Create test fixture for TPM3."""
    tpm3_exon1_cpy = tpm3_exon1.model_copy(deep=True)
    tpm3_exon1_cpy.seg.genomic_location.end = 154192133
    tpm3_exon1_cpy.seg.offset = -232
    tpm3_exon8_cpy = tpm3_exon8.model_copy(deep=True)
    tpm3_exon8_cpy.seg.genomic_location.start = 154170403
    tpm3_exon8_cpy.seg.offset = 66
    params = {
        "gene": "TPM3",
        "genomic_ac": "NC_000001.11",
        "tx_ac": "NM_152263.3",
        "seg_start": tpm3_exon1_cpy.seg,
        "seg_end": tpm3_exon8_cpy.seg,
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def mane_braf():
    """Create test fixture for BRAF (negative strand)."""
    params = {
        "gene": "BRAF",
        "genomic_ac": "NC_000007.14",
        "tx_ac": "NM_004333.6",
        "seg_start": {
            "exon_ord": 5,
            "offset": -148,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                },
                "end": 140801559,
            },
        },
        "seg_end": {
            "exon_ord": 14,
            "offset": 57,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                },
                "start": 140753336,
            },
        },
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def wee1_exon2_exon11():
    """Create test fixture for WEE1 (positive strand)."""
    params = {
        "gene": "WEE1",
        "genomic_ac": "NC_000011.10",
        "tx_ac": "NM_003390.3",
        "seg_start": {
            "exon_ord": 1,
            "offset": 205,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1",
                },
                "start": 9576092,
            },
        },
        "seg_end": {
            "exon_ord": 10,
            "offset": -1318,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1",
                },
                "end": 9588449,
            },
        },
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def mane_wee1_exon2_exon11():
    """Create test fixture for WEE1 (positive strand)."""
    params = {
        "gene": "WEE1",
        "genomic_ac": "NC_000011.10",
        "tx_ac": "NM_003390.4",
        "seg_start": {
            "exon_ord": 1,
            "offset": 205,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1",
                },
                "start": 9576092,
            },
        },
        "seg_end": {
            "exon_ord": 10,
            "offset": -1536,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.2NkFm8HK88MqeNkCgj78KidCAXgnsfV1",
                },
                "end": 9588449,
            },
        },
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def ntrk1_exon10_exon17():
    """Create test fixture for NTRK1 (positive strand)."""
    params = {
        "gene": "NTRK1",
        "genomic_ac": "NC_000001.11",
        "tx_ac": "NM_002529.3",
        "seg_start": {
            "exon_ord": 9,
            "offset": 0,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 156874570,
            },
        },
        "seg_end": {
            "exon_ord": 16,
            "offset": 0,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "end": 156881850,
            },
        },
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def zbtb10_exon3_end():
    """Create test fixture for ZBTB10, end of exon 3 (positive strand)"""
    params = {
        "gene": "ZBTB10",
        "genomic_ac": "NC_000008.11",
        "tx_ac": "NM_001105539.3",
        "seg_start": None,
        "seg_end": {
            "exon_ord": 2,
            "offset": 2,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.209Z7zJ-mFypBEWLk4rNC6S_OxY5p7bs",
                },
                "end": 80514010,
            },
        },
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def zbtb10_exon5_start():
    """Create test fixture for ZBTB10, start of exon 5 (positive strand)"""
    params = {
        "gene": "ZBTB10",
        "genomic_ac": "NC_000008.11",
        "tx_ac": "NM_001105539.3",
        "seg_start": {
            "exon_ord": 4,
            "offset": -201,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.209Z7zJ-mFypBEWLk4rNC6S_OxY5p7bs",
                },
                "start": 80518580,
            },
        },
        "seg_end": None,
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def tpm3_exon6_end():
    """Create test fixture for TPM3, end of exon 6 (negative strand)"""
    params = {
        "gene": "TPM3",
        "genomic_ac": "NC_000001.11",
        "tx_ac": "NM_152263.4",
        "seg_start": None,
        "seg_end": {
            "exon_ord": 5,
            "offset": 2,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "start": 154171410,
            },
        },
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def tpm3_exon5_start():
    """Create test fixture for TPM3, start of exon 5 (negative strand)"""
    params = {
        "gene": "TPM3",
        "genomic_ac": "NC_000001.11",
        "tx_ac": "NM_152263.4",
        "seg_start": {
            "exon_ord": 4,
            "offset": -102,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO",
                },
                "end": 154173080,
            },
        },
        "seg_end": None,
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def gusbp3_exon2_end():
    """Create test fixture for GUSBP3, end of exon 2 (negative strand)"""
    params = {
        "gene": "GUSBP3",
        "genomic_ac": "NC_000005.10",
        "tx_ac": "NR_027386.2",
        "seg_start": None,
        "seg_end": {
            "exon_ord": 1,
            "offset": 1,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
                },
                "start": 69680764,
            },
        },
    }
    return GenomicTxSegService(**params)


@pytest.fixture(scope="module")
def gusbp3_exon5_start():
    """Create test fixture for GUSBP3, start of exon 5 (negative strand)"""
    params = {
        "gene": "GUSBP3",
        "genomic_ac": "NC_000005.10",
        "tx_ac": "NR_027386.2",
        "seg_start": {
            "exon_ord": 4,
            "offset": -3589,
            "genomic_location": {
                "type": "SequenceLocation",
                "sequenceReference": {
                    "type": "SequenceReference",
                    "refgetAccession": "SQ.aUiQCzCPZ2d0csHbMSbh2NzInhonSXwI",
                },
                "end": 69645878,
            },
        },
        "seg_end": None,
    }
    return GenomicTxSegService(**params)


def check_service_meta(actual):
    """Check that service metadata matches expected

    :param ServiceMeta actual: Actual service metadata
    """
    assert actual.name == "cool_seq_tool"
    # pydantic checks that version is a string
    assert isinstance(actual.response_datetime, datetime)
    assert actual.url == "https://github.com/GenomicMedLab/cool-seq-tool"


def genomic_tx_seg_service_checks(actual, expected=None, is_valid=True):
    """Check that actual matches expected for both valid and invalid
    genomic data responses

    :param actual: Actual data
    :param expected: Expected data
    :param is_valid: `True` if expected is valid response. `False` otherwise.
    """
    if is_valid:
        assert actual.gene == expected.gene
        assert actual.genomic_ac == expected.genomic_ac
        assert actual.tx_ac == expected.tx_ac

        for seg_attr in ["seg_start", "seg_end"]:
            expected_seg = getattr(expected, seg_attr)
            if expected_seg:
                actual_seg = getattr(actual, seg_attr)
                assert actual_seg

                assert actual_seg.exon_ord == expected_seg.exon_ord
                assert actual_seg.offset == expected_seg.offset
                assert (
                    actual_seg.genomic_location.sequenceReference.refgetAccession
                    == expected_seg.genomic_location.sequenceReference.refgetAccession
                )
                assert (
                    actual_seg.genomic_location.start
                    == expected_seg.genomic_location.start
                )
                assert (
                    actual_seg.genomic_location.end == expected_seg.genomic_location.end
                )

        assert actual.errors == expected.errors
    else:
        assert actual.gene is None
        assert actual.genomic_ac is None
        assert actual.tx_ac is None
        assert actual.seg_start is None
        assert actual.seg_end is None
        assert len(actual.errors) > 0
    check_service_meta(actual.service_meta)


def get_t_to_g_args(genomic_tx_seg_service: GenomicTxSegService) -> dict:
    """Get arguments for tx_segment_to_genomic given genomic_to_tx_segment response

    :param genomic_tx_seg_service: Response from genomic_to_tx_segment
    :return: Arguments for tx_segment_to_genomic method
    """
    return {
        "transcript": genomic_tx_seg_service.tx_ac,
        "gene": genomic_tx_seg_service.gene,
        "exon_start": genomic_tx_seg_service.seg_start.exon_ord + 1
        if genomic_tx_seg_service.seg_start
        else None,
        "exon_start_offset": genomic_tx_seg_service.seg_start.offset
        if genomic_tx_seg_service.seg_start
        else 0,
        "exon_end": genomic_tx_seg_service.seg_end.exon_ord + 1
        if genomic_tx_seg_service.seg_end
        else None,
        "exon_end_offset": genomic_tx_seg_service.seg_end.offset
        if genomic_tx_seg_service.seg_end
        else 0,
    }


def genomic_tx_seg_checks(actual, expected=None, is_valid=True):
    """Check that actual matches expected for both valid and invalid
    transcript exon data responses

    :param TranscriptExonDataResponse actual: Actual data
    :param TranscriptExonData expected: Expected TranscriptExonData
    :param bool is_valid: `True` if expected is valid response.
        `False` otherwise.
    """
    if is_valid:
        assert actual.gene == expected.gene
        assert actual.genomic_ac == expected.genomic_ac
        assert actual.tx_ac == expected.tx_ac

        expected_seg = expected.seg
        if expected_seg:
            actual_seg = actual.seg
            assert actual_seg

            assert actual_seg.exon_ord == expected_seg.exon_ord
            assert actual_seg.offset == expected_seg.offset
            assert (
                actual_seg.genomic_location.sequenceReference.refgetAccession
                == expected_seg.genomic_location.sequenceReference.refgetAccession
            )
            assert (
                actual_seg.genomic_location.start == expected_seg.genomic_location.start
            )
            assert actual_seg.genomic_location.end == expected_seg.genomic_location.end

        assert actual.errors == expected.errors
    else:
        assert actual.gene is None
        assert actual.genomic_ac is None
        assert actual.tx_ac is None
        assert actual.seg is None
        assert len(actual.errors) > 0


@pytest.mark.asyncio()
async def test_get_all_exon_coords(
    test_egc_mapper, nm_152263_exons, nm_152263_exons_genomic_coords
):
    """Test that _get_all_exon_coords works correctly."""
    resp = await test_egc_mapper._get_all_exon_coords("NM_152263.3")
    assert resp == nm_152263_exons

    # Invalid transcript accession
    resp = await test_egc_mapper._get_all_exon_coords("NM_152263.36")
    assert resp == []

    resp = await test_egc_mapper._get_all_exon_coords("NM_152263.4", "NC_000001.11")
    assert resp == nm_152263_exons_genomic_coords

    # Invalid transcript accession given chromosome accession
    resp = await test_egc_mapper._get_all_exon_coords("NM_001105539.3", "NC_000001.11")
    assert resp == []


@pytest.mark.asyncio()
async def test_get_start_end_exon_coords(test_egc_mapper):
    """Test that _get_start_end_exon_coords works correctly."""
    resp = await test_egc_mapper._get_start_end_exon_coords(
        "NM_152263.3", exon_start=1, exon_end=8
    )
    assert resp == (
        ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=234,
            alt_start_i=154191901,
            alt_end_i=154192135,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=7,
            tx_start_i=822,
            tx_end_i=892,
            alt_start_i=154170399,
            alt_end_i=154170469,
            alt_strand=Strand.NEGATIVE,
        ),
        [],
    )

    resp = await test_egc_mapper._get_start_end_exon_coords(
        "NM_152263.3", exon_start=1, exon_end=11
    )
    assert resp == (None, None, ["Exon 11 does not exist on NM_152263.3"])

    resp = await test_egc_mapper._get_start_end_exon_coords(
        "NM_1234.5", exon_start=1, exon_end=11
    )
    assert resp == (None, None, ["No exons found given NM_1234.5"])


@pytest.mark.asyncio()
async def test_get_adjacent_exon(
    test_egc_mapper, nm_152263_exons_genomic_coords, nm_001105539_exons_genomic_coords
):
    """Test that get_adjacent_exon works properly"""
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_152263_exons_genomic_coords,
        end=154192100,
        strand=Strand.NEGATIVE,
    )
    assert resp == 0
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_152263_exons_genomic_coords,
        end=154191184,
        strand=Strand.NEGATIVE,
    )
    assert resp == 1
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_152263_exons_genomic_coords,
        start=154191184,
        strand=Strand.NEGATIVE,
    )
    assert resp == 2
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_001105539_exons_genomic_coords,
        end=80500385,
        strand=Strand.POSITIVE,
    )
    assert resp == 1
    resp = test_egc_mapper._get_adjacent_exon(
        tx_exons_genomic_coords=nm_001105539_exons_genomic_coords,
        start=80518580,
        strand=Strand.POSITIVE,
    )
    assert resp == 4


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
        "genomic_end": 80514010,
        "gene": "ZBTB10",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, zbtb10_exon3_end)

    inputs = {
        "chromosome": "chr8",
        "genomic_end": 80514010,
        "gene": "ZBTB10",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, zbtb10_exon3_end)

    inputs = {
        "chromosome": "8",
        "genomic_start": 80518580,
        "gene": "ZBTB10",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, zbtb10_exon5_start)

    inputs = {
        "chromosome": "1",
        "genomic_end": 154171410,
        "gene": "TPM3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, tpm3_exon6_end)

    inputs = {
        "chromosome": "1",
        "genomic_start": 154173080,
        "gene": "TPM3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, tpm3_exon5_start)

    inputs = {
        "chromosome": "5",
        "genomic_end": 69680764,
        "gene": "GUSBP3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, gusbp3_exon2_end)

    inputs = {
        "chromosome": "5",
        "genomic_start": 69645878,
        "gene": "GUSBP3",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, gusbp3_exon5_start)

    inputs = {  # Test when gene and strand are not provided
        "chromosome": "5",
        "genomic_start": 69645878,
        "transcript": "NR_027386.2",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert (
        resp.errors[0]
        == "`gene` must be provided to select the adjacent transcript junction"
    )

    inputs = {  # Test when transcript is provided
        "chromosome": "5",
        "genomic_start": 69645878,
        "gene": "GUSBP3",
        "transcript": "NR_027386.2",
        "get_nearest_transcript_junction": True,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(resp, gusbp3_exon5_start)


@pytest.mark.asyncio()
async def test_get_alt_ac_start_and_end(
    test_egc_mapper, tpm3_1_8_start_genomic, tpm3_1_8_end_genomic
):
    """Test that _get_alt_ac_start_and_end works correctly."""
    resp = await test_egc_mapper._get_alt_ac_start_and_end(
        "NM_152263.3",
        ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=234,
            alt_start_i=154191901,
            alt_end_i=154192135,
            alt_strand=Strand.NEGATIVE,
        ),
        ExonCoord(
            ord=7,
            tx_start_i=822,
            tx_end_i=892,
            alt_start_i=154170399,
            alt_end_i=154170469,
            alt_strand=Strand.NEGATIVE,
        ),
        "TPM3",
    )
    assert resp[0] == (tpm3_1_8_start_genomic, tpm3_1_8_end_genomic)
    assert resp[1] is None

    resp = await test_egc_mapper._get_alt_ac_start_and_end("NM_152263.3", gene="TPM3")
    assert resp[0] is None
    assert resp[1] == "Must provide either `tx_exon_start` or `tx_exon_end` or both"


@pytest.mark.asyncio()
async def test_get_get_exons_coords(test_egc_mapper, nm_152263_exons_genomic_coords):
    """Test that _get_all_exon_coords works correctly."""
    resp = await test_egc_mapper._get_all_exon_coords("NM_152263.4", "NC_000001.11")
    assert resp == nm_152263_exons_genomic_coords

    # Invalid transcript accession given chromosome accession
    resp = await test_egc_mapper._get_all_exon_coords("NM_001105539.3", "NC_000001.11")
    assert resp == []


@pytest.mark.asyncio()
async def test_genomic_to_transcript(test_egc_mapper, tpm3_exon1, tpm3_exon8):
    """Test that _genomic_to_tx_segment
    method works correctly.
    """
    resp = await test_egc_mapper._genomic_to_tx_segment(
        154191901,
        genomic_ac="NC_000001.11",
        transcript="NM_152263.3",
        gene="TPM3",
    )
    genomic_tx_seg_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_tx_segment(
        154191901, chromosome="1", transcript="NM_152263.3"
    )
    genomic_tx_seg_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_tx_segment(
        154191901, chromosome="1", transcript="NM_152263.3"
    )
    genomic_tx_seg_checks(resp, tpm3_exon1)

    resp = await test_egc_mapper._genomic_to_tx_segment(
        154170469,
        genomic_ac="NC_000001.11",
        transcript="NM_152263.3",
        is_start=False,
    )
    genomic_tx_seg_checks(resp, tpm3_exon8)

    resp = await test_egc_mapper._genomic_to_tx_segment(
        154170469,
        chromosome="1",
        transcript="NM_152263.3",
        is_start=False,
    )
    genomic_tx_seg_checks(resp, tpm3_exon8)

    resp = await test_egc_mapper._genomic_to_tx_segment(
        154170469, chromosome="1", transcript="NM_152263.3", is_start=False
    )
    genomic_tx_seg_checks(resp, tpm3_exon8)


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
        "genomic_ac": "NC_000001.11",
        "genomic_start": 154191901,
        "genomic_end": 154170469,
        "transcript": "NM_152263.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, tpm3_exon1_exon8)

    # t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
    #     **get_t_to_g_args(g_to_t_resp)
    # )
    # genomic_tx_seg_service_checks(t_to_g_resp, tpm3_exon1_exon8)

    # Offset
    inputs = {
        "genomic_ac": "NC_000001.11",
        "genomic_start": 154192133,
        "genomic_end": 154170403,
        "transcript": "NM_152263.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, tpm3_exon1_exon8_offset)

    # t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
    #     **get_t_to_g_args(g_to_t_resp)
    # )
    # genomic_tx_seg_service_checks(t_to_g_resp, tpm3_exon1_exon8_offset)

    # Test only setting start
    inputs = {
        "genomic_ac": "NC_000001.11",
        "genomic_start": 154191901,
        "transcript": "NM_152263.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, tpm3_exon1_g)

    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        **get_t_to_g_args(g_to_t_resp)
    )
    genomic_tx_seg_service_checks(t_to_g_resp, tpm3_exon1_g)

    # Test only setting end
    inputs = {
        "genomic_ac": "NC_000001.11",
        "genomic_end": 154170469,
        "transcript": "NM_152263.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, tpm3_exon8_g)

    # t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
    #     **get_t_to_g_args(g_to_t_resp)
    # )
    # genomic_tx_seg_service_checks(t_to_g_resp, tpm3_exon8_g)


@pytest.mark.asyncio()
async def test_braf(test_egc_mapper, mane_braf):
    """Test BRAF genomic_to_tx_segment and
    tx_segment_to_genomic.
    """
    inputs = {
        "genomic_ac": "NC_000007.13",
        "genomic_start": 140501359,  # GRCh38 coords: 140801559
        "genomic_end": 140453136,  # GRCh38 coords: 140753336
        "gene": "BRAF",
    }
    # MANE
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, mane_braf)

    # t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
    #     **get_t_to_g_args(g_to_t_resp)
    # )
    # genomic_tx_seg_service_checks(t_to_g_resp, mane_braf)


@pytest.mark.asyncio()
async def test_wee1(test_egc_mapper, wee1_exon2_exon11, mane_wee1_exon2_exon11):
    """Test WEE1 genomic_to_tx_segment and
    tx_segment_to_genomic.
    """
    inputs = {
        "genomic_ac": "NC_000011.9",
        "genomic_start": 9597639,
        "genomic_end": 9609996,
        "transcript": "NM_003390.3",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, wee1_exon2_exon11)

    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        **get_t_to_g_args(g_to_t_resp)
    )
    genomic_tx_seg_service_checks(t_to_g_resp, wee1_exon2_exon11)

    # add gene
    inputs = {
        "genomic_ac": "NC_000011.9",
        "genomic_start": 9597639,
        "genomic_end": 9609996,
        "transcript": "NM_003390.3",
        "gene": "wee1",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, wee1_exon2_exon11)

    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        **get_t_to_g_args(g_to_t_resp)
    )
    genomic_tx_seg_service_checks(t_to_g_resp, wee1_exon2_exon11)

    # MANE since no transcript provided
    inputs = {
        "genomic_ac": "NC_000011.9",
        "genomic_start": 9597639,  # GRCh38 coords: 9576092
        "genomic_end": 9609996,  # GRCh38 coords: 9588449
        "gene": "wee1",
    }
    g_to_t_resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    genomic_tx_seg_service_checks(g_to_t_resp, mane_wee1_exon2_exon11)

    t_to_g_resp = await test_egc_mapper.tx_segment_to_genomic(
        **get_t_to_g_args(g_to_t_resp)
    )
    genomic_tx_seg_service_checks(t_to_g_resp, mane_wee1_exon2_exon11)


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
    expected = tpm3_exon8_g.model_copy(deep=True)
    # resp = await test_egc_mapper.tx_segment_to_genomic(
    #     exon_start=None, exon_end=8, transcript="NM_152263.3"
    # )
    # expected.seg_end.genomic_location.start = 154170399
    # genomic_tx_seg_service_checks(resp, expected)

    # resp = await test_egc_mapper.tx_segment_to_genomic(
    #     exon_start=1, exon_end=None, transcript="NM_152263.3"
    # )
    # genomic_tx_seg_service_checks(resp, tpm3_exon1_g)

    # resp = await test_egc_mapper.tx_segment_to_genomic(
    #     exon_start=None, exon_end=8, gene="TPM3", transcript="NM_152263.3"
    # )
    # genomic_tx_seg_service_checks(resp, expected)

    # resp = await test_egc_mapper.tx_segment_to_genomic(
    #     exon_start=None, exon_end=8, gene="tpm3", transcript="NM_152263.3"
    # )
    # genomic_tx_seg_service_checks(resp, expected)

    expected = tpm3_exon1_exon8.model_copy(deep=True)
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=8, exon_end_offset=5, transcript="NM_152263.3"
    )
    expected.seg_end.offset = 5
    expected.seg_end.genomic_location.start = 154170404
    genomic_tx_seg_service_checks(resp, expected)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=3,
        exon_end=8,
        exon_start_offset=3,
        exon_end_offset=5,
        transcript="NM_152263.3",
    )
    expected.seg_start.exon_ord = 2
    expected.seg_start.offset = 3
    expected.seg_start.genomic_location.end = 154176117
    genomic_tx_seg_service_checks(resp, expected)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=3,
        exon_end=8,
        exon_start_offset=-3,
        exon_end_offset=5,
        transcript="NM_152263.3",
    )
    expected.seg_start.offset = -3
    expected.seg_start.genomic_location.end = 154176111
    genomic_tx_seg_service_checks(resp, expected)

    # NTRK1
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, transcript="NM_002529.3"
    )
    genomic_tx_seg_service_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, gene="NTRK1", transcript="NM_002529.3"
    )
    genomic_tx_seg_service_checks(resp, ntrk1_exon10_exon17)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, exon_start_offset=3, transcript="NM_002529.3"
    )
    expected = ntrk1_exon10_exon17.model_copy(deep=True)
    expected.seg_start.offset = 3
    expected.seg_start.genomic_location.start = 156874573
    genomic_tx_seg_service_checks(resp, expected)

    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=10, exon_end=17, exon_start_offset=-3, transcript="NM_002529.3"
    )
    expected.seg_start.offset = -3
    expected.seg_start.genomic_location.start = 156874567
    genomic_tx_seg_service_checks(resp, expected)


@pytest.mark.asyncio()
async def test_valid_inputs(test_egc_mapper):
    """Test that valid inputs don"t return any errors"""
    inputs = {
        "gene": "TPM3",
        "genomic_ac": "NC_000001.11",
        "genomic_start": 154171412,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert all((resp.gene, resp.genomic_ac, resp.tx_ac, resp.seg_start))

    inputs = {
        "gene": "WEE1",
        "genomic_ac": "NC_000011.9",
        "genomic_end": 9609996,
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert all((resp.gene, resp.genomic_ac, resp.tx_ac, resp.seg_end))

    inputs = {"gene": "WEE1", "chromosome": "11", "genomic_end": 9609996}
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert all((resp.gene, resp.genomic_ac, resp.tx_ac, resp.seg_end))

    inputs = {"transcript": "NM_003390.3", "exon_start": 2}
    resp = await test_egc_mapper.tx_segment_to_genomic(**inputs)
    assert all((resp.gene, resp.genomic_ac, resp.tx_ac, resp.seg_start))

    # add gene
    inputs = {"transcript": "NM_003390.3", "exon_start": 2, "gene": "WEE1"}
    resp = await test_egc_mapper.tx_segment_to_genomic(**inputs)
    assert all((resp.gene, resp.genomic_ac, resp.tx_ac, resp.seg_start))

    # Test X/Y chromosome bug
    inputs = {
        "chromosome": "X",
        "genomic_start": 154437253,
        "genomic_end": 154437299,
        "gene": "GDI1",
    }
    resp = await test_egc_mapper.genomic_to_tx_segment(**inputs)
    assert all((resp.gene, resp.genomic_ac, resp.tx_ac, resp.seg_start, resp.seg_end))

    resp = await test_egc_mapper.tx_segment_to_genomic(
        gene="PDGFRB", transcript="NM_002609.4", exon_start=11, exon_end=23
    )
    assert all((resp.gene, resp.genomic_ac, resp.tx_ac, resp.seg_start, resp.seg_end))


@pytest.mark.asyncio()
async def test_invalid(test_egc_mapper):
    """Test that invalid queries return `None`."""
    resp = await test_egc_mapper.genomic_to_tx_segment(
        transcript="NM_152263 3",
        genomic_start=154191901,
        genomic_end=154170469,
        genomic_ac="NC_000001.11",
    )
    assert resp.errors == ["No exons found given NM_152263 3"]

    # start and end not given
    resp = await test_egc_mapper.genomic_to_tx_segment(
        genomic_ac="NC_000001.11",
        genomic_start=None,
        genomic_end=None,
        transcript="NM_152263.3",
        gene="TPM3",
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Must provide either `genomic_start` or `genomic_end`"]

    # Invalid gene
    resp = await test_egc_mapper.genomic_to_tx_segment(
        genomic_ac="NC_000001.11",
        genomic_start=154191901,
        genomic_end=154170469,
        transcript="NM_152263.3",
        gene="dummy gene",
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Expected gene, DUMMY GENE, but found TPM3"]

    # Invalid accession
    resp = await test_egc_mapper.genomic_to_tx_segment(
        genomic_ac="NC_000001.200",
        genomic_start=154191901,
        genomic_end=154170469,
        transcript="NM_152263.3",
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Invalid genomic accession: NC_000001.200"]

    # Invalid coordinates
    resp = await test_egc_mapper.genomic_to_tx_segment(
        genomic_ac="NC_000001.11",
        genomic_start=9999999999998,
        genomic_end=9999999999999,
        transcript="NM_152263.3",
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == [
        "No gene(s) found given NC_000001.11 on position 9999999999998"
    ]

    resp = await test_egc_mapper.genomic_to_tx_segment(
        chromosome="1",
        genomic_start=154170469,
        transcript="NM_002529.3",
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Must find exactly one row for genomic data, but found: 0"]

    # Must supply either gene or transcript
    resp = await test_egc_mapper.genomic_to_tx_segment(
        genomic_start=154191901, genomic_ac="NC_000001.11"
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Must provide either `gene` or `transcript`"]

    # Exon 22 does not exist
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None,
        exon_end=22,
        transcript="NM_152263.3",
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Exon 22 does not exist on NM_152263.3"]

    # Start > End
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=8, exon_end=1, transcript="NM_152263.3"
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Start exon 8 is greater than end exon 1"]

    # Transcript DNE
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=7, exon_end=None, transcript="NM_12345.6"
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["No exons found given NM_12345.6"]

    # Index error for invalid exon
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=-1, exon_end=0, transcript="NM_152263.3"
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == [
        "`exon_start` cannot be less than 1",
        "`exon_end` cannot be less than 1",
    ]

    # Cant supply 0 based exons
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_end=0, transcript="NM_152263.3"
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["`exon_end` cannot be less than 1"]

    # Gene that does not match transcript
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=8, gene="NTKR1", transcript="NM_152263.3"
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == [
        "Unable to find a result where NM_152263.3 has transcript coordinates"
        " 0 and 234 between an exon's start and end coordinates on gene "
        "NTKR1"
    ]

    # No transcript given
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=1, exon_end=8, gene="NTKR1", transcript=""
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Must provide `transcript`"]

    # No exons given
    resp = await test_egc_mapper.tx_segment_to_genomic(
        exon_start=None, exon_end=None, transcript="NM_152263.3"
    )
    genomic_tx_seg_service_checks(resp, is_valid=False)
    assert resp.errors == ["Must provide either `exon_start` or `exon_end`"]
