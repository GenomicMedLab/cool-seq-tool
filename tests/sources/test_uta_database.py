"""Test UTA data source."""

from urllib.parse import urlparse

import pytest
import pytest_asyncio

from cool_seq_tool.schemas import Strand
from cool_seq_tool.sources.uta_database import (
    GenomicTxData,
    GenomicTxMetadata,
    NoMatchingAlignmentError,
    ParseResult,
    TxExonAlnData,
    UtaRepository,
    create_uta_connection_pool,
)


@pytest_asyncio.fixture(scope="function")
async def uta_repo():
    pool = await create_uta_connection_pool()
    async with pool.connection() as conn:
        yield UtaRepository(conn)
    await pool.close()


@pytest.fixture(scope="module")
def tx_exon_aln_data():
    """Create test fixture for tx_exon_aln_data test."""
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


@pytest.mark.asyncio
async def test_get_cds_start_end(uta_repo: UtaRepository):
    """Test that get_cds_start_end works correctly."""
    expected = (61, 2362)
    resp = await uta_repo.get_cds_start_end("NM_004333.4")
    assert resp == expected

    resp = await uta_repo.get_cds_start_end("ENST00000288602.6")
    assert resp == expected

    resp = await uta_repo.get_cds_start_end("NM_004333.999")
    assert resp is None


@pytest.mark.asyncio
async def test_get_newest_assembly_ac(uta_repo: UtaRepository):
    """Test that get_newest_assembly_ac works correctly."""
    resp = await uta_repo.get_newest_assembly_ac("NC_000007.13")
    assert resp == ["NC_000007.14"]

    resp = await uta_repo.get_newest_assembly_ac("NC_000011.9")
    assert resp == ["NC_000011.10"]

    resp = await uta_repo.get_newest_assembly_ac("NC_000011.10")
    assert resp == ["NC_000011.10"]

    resp = await uta_repo.get_newest_assembly_ac("ENST00000288602")
    assert resp == ["ENST00000288602"]

    resp = await uta_repo.get_newest_assembly_ac("NC_0000077.1")
    assert resp == []


@pytest.mark.asyncio
async def test_validate_genomic_ac(uta_repo: UtaRepository):
    """Test that validate_genomic_ac"""
    resp = await uta_repo.validate_genomic_ac("NC_000007.13")
    assert resp is True

    resp = await uta_repo.validate_genomic_ac("NC_000007.17")
    assert resp is False


@pytest.mark.asyncio
async def test_validate_gene_exists(uta_repo: UtaRepository):
    """Test validate_gene_symbol"""
    resp = await uta_repo.gene_exists("TPM3")
    assert resp is True

    resp = await uta_repo.gene_exists("dummy gene")
    assert resp is False


@pytest.mark.asyncio
async def test_validate_transcript_exists(uta_repo: UtaRepository):
    """Tests validate_transcript"""
    resp = await uta_repo.transcript_exists("NM_152263.3")
    assert resp is True

    resp = await uta_repo.transcript_exists("NM_152263 3")
    assert resp is False


@pytest.mark.asyncio
async def test_get_ac_descr(uta_repo: UtaRepository):
    """Test that get_ac_descr works correctly."""
    resp = await uta_repo.get_ac_descr("NC_000007.13")
    assert resp is not None

    resp = await uta_repo.get_ac_descr("NC_000007.14")
    assert resp is None


@pytest.mark.asyncio
async def test_get_tx_exon_aln_data(uta_repo: UtaRepository, tx_exon_aln_data):
    """Test that get_tx_exon_aln_data"""
    resp = await uta_repo.get_tx_exon_aln_data(
        "NM_004333.4", 140453136, 140453136, alt_ac="NC_000007.13", use_tx_pos=False
    )
    assert resp == [tx_exon_aln_data]

    resp = await uta_repo.get_tx_exon_aln_data(
        "NM_004333.4", 140453136, 140453136, alt_ac=None, use_tx_pos=False
    )
    assert resp == [tx_exon_aln_data]

    resp = await uta_repo.get_tx_exon_aln_data(
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


@pytest.mark.asyncio
async def test_mane_c_genomic_data(uta_repo: UtaRepository):
    """Test that get_mane_c_genomic_data works correctly."""
    resp = await uta_repo.get_mane_c_genomic_data(
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

    # Test example where sorting of tx_exon_aln_mv is needed
    resp = await uta_repo.get_mane_c_genomic_data(
        "NM_000077.5", "NC_000009.12", 21971186, 21971187
    )
    expected_params = {
        "gene": "CDKN2A",
        "strand": Strand.NEGATIVE,
        "tx_pos_range": (180, 487),
        "alt_pos_range": (21970901, 21971208),
        "alt_aln_method": "splign",
        "tx_exon_id": 8314723,
        "alt_exon_id": 8960507,
        "coding_start_site": 30,
        "coding_end_site": 501,
        "pos_change": (21, 285),
        "alt_pos_change_range": (21971187, 21971186),
        "tx_ac": "NM_000077.5",
        "alt_ac": "NC_000009.12",
    }
    assert resp == GenomicTxMetadata(**expected_params)

    # Test case where chromosomal accession is not provided
    resp = await uta_repo.get_mane_c_genomic_data(
        "NM_000077.5", None, 21971186, 21971187
    )
    assert resp == GenomicTxMetadata(**expected_params)


@pytest.mark.asyncio
async def test_get_genomic_tx_data(uta_repo: UtaRepository):
    """Test that get_genomic_tx_data works correctly."""
    # Positive strand transcript
    resp = await uta_repo.get_genomic_tx_data("NM_004327.3", (3595, 3596))
    expected_params = {
        "gene": "BCR",
        "strand": Strand.POSITIVE,
        "tx_pos_range": (3476, 3608),
        "alt_pos_range": (23295023, 23295155),
        "alt_aln_method": "splign",
        "tx_exon_id": 956565,
        "alt_exon_id": 6619783,
        "tx_ac": "NM_004327.3",
        "alt_ac": "NC_000022.11",
        "pos_change": (119, 12),
        "alt_pos_change_range": (23295142, 23295143),
    }
    assert resp == GenomicTxMetadata(**expected_params)

    # Negative strand transcript
    resp = await uta_repo.get_genomic_tx_data("NM_004333.4", (2144, 2145))
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
        "pos_change": (91, 43),
        "alt_pos_change_range": (140739855, 140739854),
    }
    assert resp == GenomicTxMetadata(**expected_params)


@pytest.mark.asyncio
async def test_get_ac_from_gene(uta_repo: UtaRepository):
    """Test that get_ac_from_gene works correctly."""
    resp = await uta_repo.get_ac_from_gene("BRAF")
    assert resp == ["NC_000007.14", "NC_000007.13"]

    resp = await uta_repo.get_ac_from_gene("HRAS")
    assert resp == ["NC_000011.10", "NC_000011.9"]

    resp = await uta_repo.get_ac_from_gene("dummy")
    assert resp == []


@pytest.mark.asyncio
async def test_get_gene_from_ac(uta_repo: UtaRepository):
    """Tet that get_gene_from_ac works correctly."""
    resp = await uta_repo.get_gene_from_ac("NC_000007.13", 140453136, None)
    assert resp == ["BRAF"]

    resp = await uta_repo.get_gene_from_ac("NC_000007.14", 140753336, None)
    assert resp == ["BRAF"]

    resp = await uta_repo.get_gene_from_ac("NC_000007.13", 55249071, None)
    assert resp == ["EGFR", "EGFR-AS1"]

    resp = await uta_repo.get_gene_from_ac("NC_0000078.1", 140453136, None)
    assert resp is None


@pytest.mark.asyncio
async def test_get_transcripts_from_gene(uta_repo: UtaRepository):
    """Test that get_transcripts works correctly."""
    resp = await uta_repo.get_transcripts(start_pos=2145, end_pos=2145, gene="BRAF")
    assert len(resp) == 32

    # using no start/end pos
    resp = await uta_repo.get_transcripts(gene="BRAF")
    assert len(resp) == 32

    # using 0 start/end pos
    resp = await uta_repo.get_transcripts(gene="BRAF", start_pos=0, end_pos=0)
    assert len(resp) == 32

    # using 0 genomic start/end pos
    resp = await uta_repo.get_transcripts(
        gene="BRAF", start_pos=0, end_pos=0, use_tx_pos=False
    )
    assert len(resp) == 0

    # using gene with genomic pos
    resp = await uta_repo.get_transcripts(
        gene="BRAF", start_pos=140753336, end_pos=140753336, use_tx_pos=False
    )
    assert len(resp) == 16

    resp = await uta_repo.get_transcripts(
        gene="BRAF", start_pos=140453136, end_pos=140453136
    )
    assert len(resp) == 0

    # No gene and no alt_ac
    resp = await uta_repo.get_transcripts(start_pos=140453136, end_pos=140453136)
    assert len(resp) == 0


@pytest.mark.asyncio
async def test_get_chr_assembly(uta_repo: UtaRepository):
    """Test that get_chr_assembly works correctly."""
    resp = await uta_repo.get_chr_assembly("NC_000007.13")
    assert resp == ("chr7", "GRCh37")

    resp = await uta_repo.get_chr_assembly("NC_000007.14")
    assert resp is None

    # Invalid ac
    resp = await uta_repo.get_chr_assembly("NC_00000714")
    assert resp is None


@pytest.mark.asyncio
async def test_p_to_c_ac(uta_repo: UtaRepository):
    """Test that p_to_c_ac works correctly."""
    resp = await uta_repo.p_to_c_ac("NP_004324.2")
    assert resp == ["NM_004333.4", "NM_004333.5", "NM_004333.6"]

    resp = await uta_repo.p_to_c_ac("NP_064502.9")
    assert resp == ["NM_020117.9", "NM_020117.10", "NM_020117.11"]

    resp = await uta_repo.p_to_c_ac("NP_004324.22")
    assert resp == []


@pytest.mark.asyncio
async def test_get_alt_ac_start_or_end(
    uta_repo: UtaRepository, tpm3_1_8_start_genomic, tpm3_1_8_end_genomic
):
    """Test that get_alt_ac_start_or_end works correctly."""
    resp = await uta_repo.get_alt_ac_start_or_end("NM_152263.3", 117, 234, None)
    assert resp == tpm3_1_8_start_genomic

    resp = await uta_repo.get_alt_ac_start_or_end("NM_152263.3", 822, 892, None)
    assert resp == tpm3_1_8_end_genomic

    with pytest.raises(NoMatchingAlignmentError):
        await uta_repo.get_alt_ac_start_or_end("NM_152263.63", 822, 892, None)


@pytest.mark.asyncio
async def test_get_mane_transcripts_from_genomic_pos(uta_repo: UtaRepository):
    """Test that get_mane_transcripts_from_genomic_pos works correctly"""
    resp = await uta_repo.get_transcripts_from_genomic_pos("NC_000007.14", 140753336)
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
    resp = await uta_repo.get_transcripts_from_genomic_pos("NC_000007.14", 150753336)
    assert resp == []

    # invalid ac
    resp = await uta_repo.get_transcripts_from_genomic_pos("NC_000007.14232", 140753336)
    assert resp == []


@pytest.mark.parametrize(
    ("raw_url", "expected"),
    [
        # Username + password
        (
            "postgresql://user:pass@localhost:5432/dbname",
            {
                "scheme": "postgresql",
                "username": "user",
                "password": "pass",
                "hostname": "localhost",
                "port": 5432,
                "database": "dbname",
                "sanitized_url": "postgresql://user:***@localhost:5432/dbname",
            },
        ),
        # Username with null password
        (
            "postgresql://user@localhost/dbname",
            {
                "scheme": "postgresql",
                "username": "user",
                "password": None,
                "hostname": "localhost",
                "port": None,
                "database": "dbname",
                "sanitized_url": "postgresql://user@localhost/dbname",
            },
        ),
        # Password is "0"
        (
            "postgresql://user:0@localhost/dbname",
            {
                "scheme": "postgresql",
                "username": "user",
                "password": "0",
                "hostname": "localhost",
                "port": None,
                "database": "dbname",
                "sanitized_url": "postgresql://user:***@localhost/dbname",
            },
        ),
        # Empty password
        (
            "postgresql://user:@localhost/dbname",
            {
                "scheme": "postgresql",
                "username": "user",
                "password": "",
                "hostname": "localhost",
                "port": None,
                "database": "dbname",
                "sanitized_url": "postgresql://user@localhost/dbname",
            },
        ),
        # No username
        (
            "postgresql://localhost:5432/dbname",
            {
                "scheme": "postgresql",
                "username": None,
                "password": None,
                "hostname": "localhost",
                "port": 5432,
                "database": "dbname",
                "sanitized_url": "postgresql://localhost:5432/dbname",
            },
        ),
        # With query params
        (
            "postgresql://user:secret@localhost/dbname?query#fragment",
            {
                "scheme": "postgresql",
                "username": "user",
                "password": "secret",
                "hostname": "localhost",
                "port": None,
                "database": "dbname",
                "sanitized_url": "postgresql://user:***@localhost/dbname?query#fragment",
            },
        ),
    ],
)
async def test_parsed_url(raw_url, expected):
    parsed_result = ParseResult(urlparse(raw_url))

    assert parsed_result.scheme == expected["scheme"]
    assert parsed_result.username == expected["username"]
    assert parsed_result.password == expected["password"]
    assert parsed_result.hostname == expected["hostname"]
    assert parsed_result.port == expected["port"]
    assert parsed_result.database == expected["database"]
    assert parsed_result.sanitized_url == expected["sanitized_url"]
