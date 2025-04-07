"""Provide utilities for test cases."""

import asyncio

import pytest

from cool_seq_tool import CoolSeqTool
from cool_seq_tool.mappers.exon_genomic_coords import _ExonCoord
from cool_seq_tool.schemas import ManeGeneData, Strand
from cool_seq_tool.sources.uta_database import GenomicAlnData, GenomicTxMetadata


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
        _ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=234,
            alt_start_i=154191901,
            alt_end_i=154192135,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=1,
            tx_start_i=234,
            tx_end_i=360,
            alt_start_i=154191185,
            alt_end_i=154191311,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=2,
            tx_start_i=360,
            tx_end_i=494,
            alt_start_i=154176114,
            alt_end_i=154176248,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=3,
            tx_start_i=494,
            tx_end_i=612,
            alt_start_i=154173083,
            alt_end_i=154173201,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=4,
            tx_start_i=612,
            tx_end_i=683,
            alt_start_i=154172907,
            alt_end_i=154172978,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=5,
            tx_start_i=683,
            tx_end_i=759,
            alt_start_i=154171412,
            alt_end_i=154171488,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=6,
            tx_start_i=759,
            tx_end_i=822,
            alt_start_i=154170648,
            alt_end_i=154170711,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=7,
            tx_start_i=822,
            tx_end_i=892,
            alt_start_i=154170399,
            alt_end_i=154170469,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=8,
            tx_start_i=892,
            tx_end_i=971,
            alt_start_i=154169304,
            alt_end_i=154169383,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=9,
            tx_start_i=971,
            tx_end_i=7099,
            alt_start_i=154161812,
            alt_end_i=154167940,
            alt_strand=Strand.NEGATIVE,
        ),
    ]


@pytest.fixture(scope="session")
def nm_152263_exons_genomic_coords():
    """Create test fixture for NM_152263.4 exons and genomic coordinates."""
    return [
        _ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=199,
            alt_start_i=154191901,
            alt_end_i=154192100,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=1,
            tx_start_i=199,
            tx_end_i=325,
            alt_start_i=154191185,
            alt_end_i=154191311,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=2,
            tx_start_i=325,
            tx_end_i=459,
            alt_start_i=154176114,
            alt_end_i=154176248,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=3,
            tx_start_i=459,
            tx_end_i=577,
            alt_start_i=154173083,
            alt_end_i=154173201,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=4,
            tx_start_i=577,
            tx_end_i=648,
            alt_start_i=154172907,
            alt_end_i=154172978,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=5,
            tx_start_i=648,
            tx_end_i=724,
            alt_start_i=154171412,
            alt_end_i=154171488,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=6,
            tx_start_i=724,
            tx_end_i=787,
            alt_start_i=154170648,
            alt_end_i=154170711,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=7,
            tx_start_i=787,
            tx_end_i=857,
            alt_start_i=154170399,
            alt_end_i=154170469,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=8,
            tx_start_i=857,
            tx_end_i=936,
            alt_start_i=154169304,
            alt_end_i=154169383,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=9,
            tx_start_i=936,
            tx_end_i=7064,
            alt_start_i=154161812,
            alt_end_i=154167940,
            alt_strand=Strand.NEGATIVE,
        ),
    ]


@pytest.fixture(scope="session")
def nm_001105539_exons_genomic_coords():
    """Create test fixture for NM_001105539.3 exons and genomic coordinates."""
    return [
        _ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=1557,
            alt_start_i=80486225,
            alt_end_i=80487782,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=1,
            tx_start_i=1557,
            tx_end_i=2446,
            alt_start_i=80499493,
            alt_end_i=80500382,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=2,
            tx_start_i=2446,
            tx_end_i=2545,
            alt_start_i=80513909,
            alt_end_i=80514008,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=3,
            tx_start_i=2545,
            tx_end_i=2722,
            alt_start_i=80518402,
            alt_end_i=80518579,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=4,
            tx_start_i=2722,
            tx_end_i=2895,
            alt_start_i=80518781,
            alt_end_i=80518954,
            alt_strand=Strand.NEGATIVE,
        ),
        _ExonCoord(
            ord=5,
            tx_start_i=2895,
            tx_end_i=9938,
            alt_start_i=80519222,
            alt_end_i=80526265,
            alt_strand=Strand.NEGATIVE,
        ),
    ]


@pytest.fixture(scope="session")
def mm_001005183_1_exons():
    """Create test fixture for NM_001005183.1 exons and genomic coordinates"""
    return [
        _ExonCoord(
            ord=0,
            tx_start_i=0,
            tx_end_i=939,
            alt_start_i=55426253,
            alt_end_i=55427192,
            alt_strand=Strand.POSITIVE,
        )
    ]


@pytest.fixture(scope="session")
def tpm3_1_8_start_genomic():
    """Create test fixture for genomic data for exon 1, 8"""
    return GenomicAlnData(
        hgnc="TPM3",
        ord=0,
        alt_ac="NC_000001.11",
        alt_start_i=154191901,
        alt_end_i=154192135,
        alt_strand=Strand.NEGATIVE,
    )


@pytest.fixture(scope="session")
def tpm3_1_8_end_genomic():
    """Create test fixture for genomic data for exon 1, 8"""
    return GenomicAlnData(
        hgnc="TPM3",
        ord=7,
        alt_ac="NC_000001.11",
        alt_start_i=154170399,
        alt_end_i=154170469,
        alt_strand=Strand.NEGATIVE,
    )


@pytest.fixture(scope="session")
def genomic_tx_data():
    """Create test fixture for genomic_tx_data"""
    params = {
        "gene": "BRAF",
        "strand": Strand.NEGATIVE,
        "tx_pos_range": (2053, 2188),
        "alt_pos_range": (140439611, 140439746),
        "alt_aln_method": "splign",
        "tx_exon_id": 780496,
        "alt_exon_id": 1927265,
        "pos_change": (92, 43),
        "alt_pos_change_range": (140439703, 140439703),
        "tx_ac": "NM_004333.4",
        "alt_ac": "NC_000007.13",
    }
    return GenomicTxMetadata(**params)


@pytest.fixture(scope="session")
def egfr_mane_gene():
    """Create test fixture for EGFR MANE gene"""
    return ManeGeneData(
        ncbi_gene_id=1956, hgnc_id=3236, symbol="EGFR", status=["mane_select"]
    )


@pytest.fixture(scope="session")
def braf_mane_genes():
    """Create test fixture for BRAF MANE gene"""
    return [
        ManeGeneData(
            ncbi_gene_id=673,
            hgnc_id=1097,
            symbol="BRAF",
            status=["mane_select", "mane_plus_clinical"],
        ),
    ]
