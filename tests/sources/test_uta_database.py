"""Test UTA data source."""

import pytest

from cool_seq_tool.schemas import Strand
from cool_seq_tool.sources.uta_database import (
    GenomicTxData,
    GenomicTxMetadata,
    TxExonAlnData,
)


@pytest.fixture(scope="module")
def tx_exon_aln_v_data():
    """Create test fixture for tx_aln_v_data test."""
    return TxExonAlnData(
        hgnc="BRAF",
        ord=14,
        tx_ac="NM_004333.4",
        tx_start_i=1802,
        tx_end_i=1921,
        alt_ac="NC_000007.13",
        alt_start_i=140453074,
        alt_end_i=140453193,
        alt_strand=Strand.NEGATIVE,
        alt_aln_method="splign",
        tx_exon_id=780494,
        alt_exon_id=1927263,
    )


@pytest.fixture(scope="module")
def data_from_result():
    """Create test fixture for data from result"""
    params = {
        "gene": "BRAF",
        "strand": Strand.NEGATIVE,
        "tx_pos_range": (1802, 1921),
        "alt_pos_range": (140453074, 140453193),
        "alt_aln_method": "splign",
        "tx_exon_id": 780494,
        "alt_exon_id": 1927263,
    }
    return GenomicTxData(**params)


@pytest.mark.asyncio()
async def test_get_cds_start_end(test_db):
    """Test that get_cds_start_end works correctly."""
    expected = (61, 2362)
    resp = await test_db.get_cds_start_end("NM_004333.4")
    assert resp == expected

    resp = await test_db.get_cds_start_end("ENST00000288602.6")
    assert resp == expected

    resp = await test_db.get_cds_start_end("NM_004333.999")
    assert resp is None


@pytest.mark.asyncio()
async def test_get_newest_assembly_ac(test_db):
    """Test that get_newest_assembly_ac works correctly."""
    resp = await test_db.get_newest_assembly_ac("NC_000007.13")
    assert resp == ["NC_000007.14"]

    resp = await test_db.get_newest_assembly_ac("NC_000011.9")
    assert resp == ["NC_000011.10"]

    resp = await test_db.get_newest_assembly_ac("NC_000011.10")
    assert resp == ["NC_000011.10"]

    resp = await test_db.get_newest_assembly_ac("ENST00000288602")
    assert resp == ["ENST00000288602"]

    resp = await test_db.get_newest_assembly_ac("NC_0000077.1")
    assert resp == []


@pytest.mark.asyncio()
async def test_validate_genomic_ac(test_db):
    """Test that validate_genomic_ac"""
    resp = await test_db.validate_genomic_ac("NC_000007.13")
    assert resp is True

    resp = await test_db.validate_genomic_ac("NC_000007.17")
    assert resp is False


@pytest.mark.asyncio()
async def test_validate_gene_exists(test_db):
    """Test validate_gene_symbol"""
    resp = await test_db.gene_exists("TPM3")
    assert resp is True

    resp = await test_db.gene_exists("dummy gene")
    assert resp is False


@pytest.mark.asyncio()
async def test_validate_transcript_exists(test_db):
    """Tests validate_transcript"""
    resp = await test_db.transcript_exists("NM_152263.3")
    assert resp is True

    resp = await test_db.transcript_exists("NM_152263 3")
    assert resp is False


@pytest.mark.asyncio()
async def test_get_ac_descr(test_db):
    """Test that get_ac_descr works correctly."""
    resp = await test_db.get_ac_descr("NC_000007.13")
    assert resp is not None

    resp = await test_db.get_ac_descr("NC_000007.14")
    assert resp is None


@pytest.mark.asyncio()
async def test_get_tx_exon_aln_v_data(test_db, tx_exon_aln_v_data):
    """Test that get_tx_exon_aln_v_data"""
    resp = await test_db.get_tx_exon_aln_v_data(
        "NM_004333.4", 140453136, 140453136, alt_ac="NC_000007.13", use_tx_pos=False
    )
    assert resp == [tx_exon_aln_v_data]

    resp = await test_db.get_tx_exon_aln_v_data(
        "NM_004333.4", 140453136, 140453136, alt_ac=None, use_tx_pos=False
    )
    assert resp == [tx_exon_aln_v_data]

    resp = await test_db.get_tx_exon_aln_v_data(
        "NM_004333.4", 1860, 1860, alt_ac=None, use_tx_pos=True
    )
    assert resp == [
        TxExonAlnData(
            hgnc="BRAF",
            ord=14,
            tx_ac="NM_004333.4",
            tx_start_i=1802,
            tx_end_i=1921,
            alt_ac="NC_000007.13",
            alt_start_i=140453074,
            alt_end_i=140453193,
            alt_strand=Strand.NEGATIVE,
            alt_aln_method="splign",
            tx_exon_id=780494,
            alt_exon_id=1927263,
        ),
        TxExonAlnData(
            hgnc="BRAF",
            ord=14,
            tx_ac="NM_004333.4",
            tx_start_i=1802,
            tx_end_i=1921,
            alt_ac="NC_000007.14",
            alt_start_i=140753274,
            alt_end_i=140753393,
            alt_strand=Strand.NEGATIVE,
            alt_aln_method="splign",
            tx_exon_id=780494,
            alt_exon_id=6619850,
        ),
    ]


@pytest.mark.asyncio()
async def test_data_from_result(test_db, tx_exon_aln_v_data, data_from_result):
    """Test that data_from_result works correctly."""
    resp = test_db.data_from_result(tx_exon_aln_v_data)
    assert resp == data_from_result


@pytest.mark.asyncio()
async def test_mane_c_genomic_data(test_db):
    """Test that get_mane_c_genomic_data works correctly."""
    resp = await test_db.get_mane_c_genomic_data(
        "NM_001374258.1", None, 140753335, 140753335
    )
    expected_params = {
        "gene": "BRAF",
        "strand": Strand.NEGATIVE,
        "tx_pos_range": (2087, 2206),
        "alt_pos_range": (140753274, 140753393),
        "alt_aln_method": "splign",
        "tx_exon_id": 8439617,
        "alt_exon_id": 9353476,
        "coding_start_site": 226,
        "coding_end_site": 2650,
        "pos_change": (58, 61),
        "alt_pos_change_range": (140753335, 140753335),
        "tx_ac": "NM_001374258.1",
        "alt_ac": "NC_000007.14",
    }
    assert resp == GenomicTxMetadata(**expected_params)


@pytest.mark.asyncio()
async def test_get_genomic_tx_data(test_db, genomic_tx_data):
    """Test that get_genomic_tx_data works correctly."""
    resp = await test_db.get_genomic_tx_data("NM_004333.4", (2145, 2145))
    expected_params = {
        "gene": "BRAF",
        "strand": Strand.NEGATIVE,
        "tx_pos_range": (2053, 2188),
        "alt_pos_range": (140739811, 140739946),
        "alt_aln_method": "splign",
        "tx_exon_id": 780496,
        "alt_exon_id": 6619852,
        "tx_ac": "NM_004333.4",
        "alt_ac": "NC_000007.14",
        "pos_change": (92, 43),
        "alt_pos_change_range": (140739854, 140739854),
    }
    assert resp == GenomicTxMetadata(**expected_params)


@pytest.mark.asyncio()
async def test_get_ac_from_gene(test_db):
    """Test that get_ac_from_gene works correctly."""
    resp = await test_db.get_ac_from_gene("BRAF")
    assert resp == ["NC_000007.14", "NC_000007.13"]

    resp = await test_db.get_ac_from_gene("HRAS")
    assert resp == ["NC_000011.10", "NC_000011.9"]

    resp = await test_db.get_ac_from_gene("dummy")
    assert resp == []


@pytest.mark.asyncio()
async def test_get_gene_from_ac(test_db):
    """Tet that get_gene_from_ac works correctly."""
    resp = await test_db.get_gene_from_ac("NC_000007.13", 140453136, None)
    assert resp == ["BRAF"]

    resp = await test_db.get_gene_from_ac("NC_000007.14", 140753336, None)
    assert resp == ["BRAF"]

    resp = await test_db.get_gene_from_ac("NC_000007.13", 55249071, None)
    assert resp == ["EGFR", "EGFR-AS1"]

    resp = await test_db.get_gene_from_ac("NC_0000078.1", 140453136, None)
    assert resp is None


@pytest.mark.asyncio()
async def test_get_transcripts_from_gene(test_db):
    """Test that get_transcripts works correctly."""
    resp = await test_db.get_transcripts(start_pos=2145, end_pos=2145, gene="BRAF")
    assert len(resp) == 32

    # using no start/end pos
    resp = await test_db.get_transcripts(gene="BRAF")
    assert len(resp) == 32

    # using 0 start/end pos
    resp = await test_db.get_transcripts(gene="BRAF", start_pos=0, end_pos=0)
    assert len(resp) == 32

    # using 0 genomic start/end pos
    resp = await test_db.get_transcripts(
        gene="BRAF", start_pos=0, end_pos=0, use_tx_pos=False
    )
    assert len(resp) == 0

    # using gene with genomic pos
    resp = await test_db.get_transcripts(
        gene="BRAF", start_pos=140753336, end_pos=140753336, use_tx_pos=False
    )
    assert len(resp) == 16

    resp = await test_db.get_transcripts(
        gene="BRAF", start_pos=140453136, end_pos=140453136
    )
    assert len(resp) == 0

    # No gene and no alt_ac
    resp = await test_db.get_transcripts(start_pos=140453136, end_pos=140453136)
    assert len(resp) == 0


@pytest.mark.asyncio()
async def test_get_chr_assembly(test_db):
    """Test that get_chr_assembly works correctly."""
    resp = await test_db.get_chr_assembly("NC_000007.13")
    assert resp == ("chr7", "GRCh37")

    resp = await test_db.get_chr_assembly("NC_000007.14")
    assert resp is None

    # Invalid ac
    resp = await test_db.get_chr_assembly("NC_00000714")
    assert resp is None


@pytest.mark.asyncio()
async def test_p_to_c_ac(test_db):
    """Test that p_to_c_ac works correctly."""
    resp = await test_db.p_to_c_ac("NP_004324.2")
    assert resp == ["NM_004333.4", "NM_004333.5", "NM_004333.6"]

    resp = await test_db.p_to_c_ac("NP_064502.9")
    assert resp == ["NM_020117.9", "NM_020117.10", "NM_020117.11"]

    resp = await test_db.p_to_c_ac("NP_004324.22")
    assert resp == []


@pytest.mark.asyncio()
async def test_get_alt_ac_start_or_end(
    test_db, tpm3_1_8_start_genomic, tpm3_1_8_end_genomic
):
    """Test that get_alt_ac_start_or_end works correctly."""
    resp = await test_db.get_alt_ac_start_or_end("NM_152263.3", 117, 234, None)
    assert resp[0] == tpm3_1_8_start_genomic
    assert resp[1] is None

    resp = await test_db.get_alt_ac_start_or_end("NM_152263.3", 822, 892, None)
    assert resp[0] == tpm3_1_8_end_genomic
    assert resp[1] is None

    resp = await test_db.get_alt_ac_start_or_end("NM_152263.63", 822, 892, None)
    assert resp[0] is None
    assert (
        resp[1] == "Unable to find a result where NM_152263.63 has "
        "transcript coordinates 822 and 892 between an exon's "
        "start and end coordinates"
    )


@pytest.mark.asyncio()
async def test_get_mane_transcripts_from_genomic_pos(test_db):
    """Test that get_mane_transcripts_from_genomic_pos works correctly"""
    resp = await test_db.get_transcripts_from_genomic_pos("NC_000007.14", 140753336)
    assert set(resp) == {
        "NM_001354609.1",
        "NM_001354609.2",
        "NM_001374244.1",
        "NM_001374258.1",
        "NM_001378467.1",
        "NM_001378468.1",
        "NM_001378469.1",
        "NM_001378470.1",
        "NM_001378471.1",
        "NM_001378472.1",
        "NM_001378473.1",
        "NM_001378474.1",
        "NM_001378475.1",
        "NM_004333.4",
        "NM_004333.5",
        "NM_004333.6",
    }

    # invalid pos
    resp = await test_db.get_transcripts_from_genomic_pos("NC_000007.14", 150753336)
    assert resp == []

    # invalid ac
    resp = await test_db.get_transcripts_from_genomic_pos("NC_000007.14232", 140753336)
    assert resp == []
