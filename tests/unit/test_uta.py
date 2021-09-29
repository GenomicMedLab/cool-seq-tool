"""Test UTA data source."""
import pytest
from uta_tools.uta import UTADatabase
import copy


@pytest.fixture(scope='session')
async def test_db():
    """Create uta db test fixture."""
    class TestUTADatabase:
        def __init__(self):
            self.test_db = UTADatabase()
            self.test_db.create_pool()

        async def transcript_to_genomic(
                self, tx_ac: str, start_exon: int,
                end_exon: int, start_exon_offset: int = 0,
                end_exon_offset: int = 0, gene: str = None):
            return await self.test_db.transcript_to_genomic(
                tx_ac, start_exon, end_exon, start_exon_offset,
                end_exon_offset, gene)
    return TestUTADatabase()


@pytest.fixture(scope='module')
def tpm3_exon1_exon8():
    """Create test fixture for TPM3."""
    return {
        "gene": "TPM3",
        "chr": "NC_000001.11",
        "start": 154192135,
        "end": 154170399,
        "exon_start": 1,
        "exon_end": 8,
        "exon_end_offset": 0,
        "exon_start_offset": 0
    }


@pytest.fixture(scope='module')
def ntrk1_exon10_exon17():
    """Create test fixture for NTRK1."""
    return {
        "gene": "NTRK1",
        "chr": "NC_000001.11",
        "start": 156874626,
        "end": 156881456,
        "exon_start": 10,
        "exon_end": 17,
        "exon_end_offset": 0,
        "exon_start_offset": 0
    }


@pytest.mark.asyncio
async def test_transcript_to_genomic(test_db, tpm3_exon1_exon8,
                                     ntrk1_exon10_exon17):
    """Test that transcript_to_genomic works correctly."""
    # TPM3
    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 8)
    assert resp == tpm3_exon1_exon8

    resp = await test_db.transcript_to_genomic('NM_152263.3       ', 0, 8)
    assert resp == tpm3_exon1_exon8

    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 8,
                                               gene="TPM3")
    assert resp == tpm3_exon1_exon8

    resp = await test_db.transcript_to_genomic(' NM_152263.3 ', 0, 8,
                                               gene=" TPM3 ")
    assert resp == tpm3_exon1_exon8

    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 8,
                                               gene="tpm3")
    assert resp == tpm3_exon1_exon8

    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 0,
                                               gene="tpm3")
    expected = copy.deepcopy(tpm3_exon1_exon8)
    expected["exon_end"] = 10
    expected["end"] = 154161812
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 8,
                                               end_exon_offset=-5)
    expected["exon_end"] = 8
    expected["exon_end_offset"] = -5
    expected["end"] = 154170404
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 8,
                                               end_exon_offset=5)
    expected["exon_end_offset"] = 5
    expected["end"] = 154170394
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_152263.3', 3, 8,
                                               start_exon_offset=3,
                                               end_exon_offset=5)
    expected["exon_start"] = 3
    expected["exon_start_offset"] = 3
    expected["start"] = 154176245
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_152263.3', 3, 8,
                                               start_exon_offset=-3,
                                               end_exon_offset=5)
    expected["exon_start_offset"] = -3
    expected["start"] = 154176251
    assert resp == expected

    # NTRK1
    resp = await test_db.transcript_to_genomic('NM_002529.3', 10, 0)
    assert resp == ntrk1_exon10_exon17

    resp = await test_db.transcript_to_genomic('NM_002529.3', 10, 0,
                                               gene="NTRK1")
    assert resp == ntrk1_exon10_exon17

    resp = await test_db.transcript_to_genomic('NM_002529.3', 10, 0,
                                               gene="NTRK1")
    assert resp == ntrk1_exon10_exon17

    resp = await test_db.transcript_to_genomic('NM_002529.3', 10, 0,
                                               start_exon_offset=3)
    expected = copy.deepcopy(ntrk1_exon10_exon17)
    expected["exon_start_offset"] = 3
    expected["start"] = 156874629
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_002529.3', 10, 0,
                                               start_exon_offset=-3)
    expected["exon_start_offset"] = -3
    expected["start"] = 156874623
    assert resp == expected


@pytest.mark.asyncio
async def test_no_matches(test_db):
    """Test that invalid queries return None."""
    # Exon 22 does not exist
    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 22)
    assert resp is None

    # Start > End
    resp = await test_db.transcript_to_genomic('NM_152263.3', 8, 1)
    assert resp is None

    # End < Start
    resp = await test_db.transcript_to_genomic('NM_152263.3', 7, 6)
    assert resp is None

    # Transcript DNE
    resp = await test_db.transcript_to_genomic('NM_12345.6', 7, 0)
    assert resp is None

    # Index error for invalid exon
    resp = await test_db.transcript_to_genomic('NM_12345.6', -1, 0)
    assert resp is None

    # Gene that does not match transcript
    resp = await test_db.transcript_to_genomic('NM_152263.3', 8, 1,
                                               gene='NTKR1')
    assert resp is None

    # No transcript given
    resp = await test_db.transcript_to_genomic(None, 8, 1, gene='NTKR1')
    assert resp is None

    resp = await test_db.transcript_to_genomic('', 8, 1, gene='NTKR1')
    assert resp is None
