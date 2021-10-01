"""Test UTA data source."""
import pytest
from uta_tools.data_sources import UTADatabase
import copy


@pytest.fixture(scope='session')
async def test_db():
    """Create uta db test fixture."""
    test_uta_db = UTADatabase()
    await test_uta_db._create_genomic_table()
    return test_uta_db


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


@pytest.fixture(scope='module')
def nm_152263_exons():
    """Create test fixture for NM_152263.3 exons."""
    return ["117,234", "234,360", "360,494", "494,612", "612,683", "683,759",
            "759,822", "822,892", "892,971", "971,975"]


@pytest.fixture(scope='module')
def tpm3_1_8_start_genomic():
    """Create test fixture for genomic data for exon 1, 8"""
    return 'TPM3', 'NC_000001.11', 154191901, 154192135, -1


@pytest.fixture(scope='module')
def tpm3_1_8_end_genomic():
    """Create test fixture for genomic data for exon 1, 8"""
    return 'TPM3', 'NC_000001.11', 154170399, 154170469, -1


@pytest.fixture(scope='module')
def tx_exon_aln_v_data():
    """Create test fixture for tx_aln_v_data test."""
    return ["BRAF", "NM_004333.4", 1802, 1921, "NC_000007.13", 140453074,
            140453193, -1, "splign", 780494, 1927263]


@pytest.fixture(scope='module')
def data_from_result():
    """Create test fixture for data from result"""
    return dict(
        gene="BRAF",
        strand="-",
        tx_pos_range=(1802, 1921),
        alt_pos_range=(140453074, 140453193),
        alt_aln_method="splign",
        tx_exon_id=780494,
        alt_exon_id=1927263
    )


@pytest.fixture(scope='module')
def genomic_tx_data():
    """Create test fixture for genomic_tx_data"""
    return dict(
        gene="BRAF",
        strand="-",
        tx_pos_range=(2053, 2188),
        alt_pos_range=(140439611, 140439746),
        alt_aln_method="splign",
        tx_exon_id=780496,
        alt_exon_id=1927265,
        pos_change=(92, 43),
        alt_pos_change_range=(140439703, 140439703),
        tx_ac="NM_004333.4",
        alt_ac="NC_000007.13"
    )


@pytest.mark.asyncio
async def test_get_tx_exons(test_db, nm_152263_exons):
    """Test that get_tx_exons works correctly."""
    resp = await test_db.get_tx_exons('NM_152263.3')
    assert resp == nm_152263_exons

    # Invalid transcript accession
    resp = await test_db.get_tx_exons('NM_152263.36')
    assert resp is None


@pytest.mark.asyncio
async def test_get_cds_start_end(test_db):
    """Test that get_cds_start_end works correctly."""
    expected = (61, 2362)
    resp = await test_db.get_cds_start_end('NM_004333.4')
    assert resp == expected

    resp = await test_db.get_cds_start_end('ENST00000288602.6')
    assert resp == expected

    resp = await test_db.get_cds_start_end('NM_004333.999')
    assert resp is None


@pytest.mark.asyncio
async def test_get_newest_assembly(test_db):
    """Test that get_newest_assembly works correctly."""
    resp = await test_db.get_newest_assembly_ac("NC_000007.13")
    assert resp == [["NC_000007.14"]]

    resp = await test_db.get_newest_assembly_ac("NC_0000077.1")
    assert resp == []


@pytest.mark.asyncio
async def test_validate_genomic_ac(test_db):
    """Test that validate_genomic_ac"""
    resp = await test_db.validate_genomic_ac("NC_000007.13")
    assert resp is True

    resp = await test_db.validate_genomic_ac("NC_000007.17")
    assert resp is False


@pytest.mark.asyncio
async def test_get_ac_descr(test_db):
    """Test that get_ac_descr works correctly."""
    resp = await test_db.get_ac_descr("NC_000007.13")
    assert resp is not None

    resp = await test_db.get_ac_descr("NC_000007.14")
    assert resp is None


@pytest.mark.asyncio
async def test_get_tx_exon_aln_v_data(test_db, tx_exon_aln_v_data):
    """Test that get_tx_exon_aln_v_data"""
    resp = await test_db.get_tx_exon_aln_v_data(
        "NM_004333.4", 140453136, 140453136, alt_ac="NC_000007.13",
        use_tx_pos=False)
    assert resp == [tx_exon_aln_v_data]

    resp = await test_db.get_tx_exon_aln_v_data(
        "NM_004333.4", 140453136, 140453136, alt_ac=None, use_tx_pos=False)
    assert resp == [tx_exon_aln_v_data]

    resp = await test_db.get_tx_exon_aln_v_data(
        "NM_004333.4", 140453136, None, alt_ac=None, use_tx_pos=False)
    assert resp == [tx_exon_aln_v_data]

    resp = await test_db.get_tx_exon_aln_v_data(
        "NM_004333.4", 1860, None, alt_ac=None, use_tx_pos=True)
    assert resp == [
        ['BRAF', 'NM_004333.4', 1802, 1921, 'NC_000007.13', 140453074,
         140453193, -1, 'splign', 780494, 1927263],
        ['BRAF', 'NM_004333.4', 1802, 1921, 'NC_000007.14', 140753274,
         140753393, -1, 'splign', 780494, 6619850]]


@pytest.mark.asyncio
async def test_data_from_result(test_db, tx_exon_aln_v_data, data_from_result):
    """Test that data_from_result works correctly."""
    resp = test_db.data_from_result(tx_exon_aln_v_data)
    assert resp == data_from_result


@pytest.mark.asyncio
async def test_mane_c_genomic_data(test_db):
    """Test that get_mane_c_genomic_data works correctly."""
    resp = await test_db.get_mane_c_genomic_data(
        "NM_001374258.1", None, 140753335, 140753335
    )
    expected = dict(
        gene="BRAF",
        strand="-",
        tx_pos_range=(2087, 2206),
        alt_pos_range=(140753274, 140753393),
        alt_aln_method="splign",
        tx_exon_id=8439617,
        alt_exon_id=9353476,
        coding_start_site=226,
        coding_end_site=2650,
        alt_pos_change=(61, 58),
        alt_pos_change_range=(140753335, 140753335),
        tx_ac="NM_001374258.1",
        alt_ac="NC_000007.14"
    )
    assert resp == expected


@pytest.mark.asyncio
async def test_get_genomic_tx_data(test_db, genomic_tx_data):
    """Test that get_genomic_tx_data works correctly."""
    resp = await test_db.get_genomic_tx_data("NM_004333.4", (2145, 2145))
    assert resp == {
        'gene': 'BRAF',
        'strand': '-',
        'tx_pos_range': (2053, 2188),
        'alt_pos_range': (140739811, 140739946),
        'alt_aln_method': 'splign',
        'tx_exon_id': 780496,
        'alt_exon_id': 6619852,
        'tx_ac': 'NM_004333.4',
        'alt_ac': 'NC_000007.14',
        'pos_change': (92, 43),
        'alt_pos_change_range': (140739903, 140739903)
    }


@pytest.mark.asyncio
async def test_get_ac_from_gene(test_db):
    """Test that get_ac_from_gene works correctly."""
    resp = await test_db.get_ac_from_gene('BRAF')
    assert resp == ['NC_000007.14', 'NC_000007.13']


@pytest.mark.asyncio
async def test_get_gene_from_ac(test_db):
    """Tet that get_gene_from_ac works correctly."""
    resp = await test_db.get_gene_from_ac('NC_000007.13', 140453136, None)
    assert resp == ["BRAF"]

    resp = await test_db.get_gene_from_ac('NC_000007.14', 140753336, None)
    assert resp == ["BRAF"]

    resp = await test_db.get_gene_from_ac('NC_000007.13', 55249071, None)
    assert resp == ["EGFR", "EGFR-AS1"]

    resp = await test_db.get_gene_from_ac('NC_0000078.1', 140453136, None)
    assert resp is None


@pytest.mark.asyncio
async def test_get_transcripts_from_gene(test_db):
    """Test that get_trasncripts_from_gene works correctly."""
    resp = await test_db.get_transcripts_from_gene("BRAF", 2145, 2145)
    assert len(resp) == 32

    resp = await test_db.get_transcripts_from_gene("BRAF", 140453136,
                                                   140453136)
    assert len(resp) == 0


@pytest.mark.asyncio
async def test_get_chr_assembly(test_db):
    """Test that get_chr_assembly works correctly."""
    resp = await test_db.get_chr_assembly("NC_000007.13")
    assert resp == ("chr7", "GRCh37")

    resp = await test_db.get_chr_assembly("NC_000007.14")
    assert resp is None


@pytest.mark.asyncio
async def test_liftover_to_38(test_db, genomic_tx_data):
    """Test that liftover_to_38 works correctly."""
    cpy = copy.deepcopy(genomic_tx_data)
    expected = copy.deepcopy(genomic_tx_data)
    await test_db.liftover_to_38(cpy)
    expected["alt_ac"] = "NC_000007.14"
    expected["alt_pos_change_range"] = (140739903, 140739903)
    expected["alt_pos_range"] = (140739811, 140739946)
    assert cpy == expected


def test_get_liftover(test_db):
    """Test that get_liftover works correctly."""
    resp = test_db.get_liftover("chr7", 140453136)
    assert resp == ('chr7', 140753336, '+', 14633688187)

    resp = test_db.get_liftover("chr17", 140453136)
    assert resp is None


def test_set_liftover(test_db, genomic_tx_data):
    """Test that _set_liftover works correctly."""
    cpy = copy.deepcopy(genomic_tx_data)
    expected = copy.deepcopy(genomic_tx_data)
    test_db._set_liftover(cpy, "alt_pos_range", "chr7")
    expected["alt_pos_range"] = (140739811, 140739946)
    assert cpy == expected
    test_db._set_liftover(cpy, "alt_pos_change_range", "chr7")
    expected["alt_pos_change_range"] = (140739903, 140739903)
    assert cpy == expected


@pytest.mark.asyncio
async def test_p_to_c_ac(test_db):
    """Test that p_to_c_ac works correctly."""
    resp = await test_db.p_to_c_ac("NP_004324.2")
    assert resp == ['NM_004333.4', 'NM_004333.5', 'NM_004333.6']

    resp = await test_db.p_to_c_ac("NP_004324.22")
    assert resp == []


@pytest.mark.asyncio
async def test_get_tx_exon_start_end(test_db, nm_152263_exons):
    """Test that get_tx_exon_start_end works correctly."""
    resp = await test_db.get_tx_exon_start_end('NM_152263.3', 1, 8)
    assert resp == (nm_152263_exons, 1, 8)

    resp = await test_db.get_tx_exon_start_end('NM_152263.3', 0, 8)
    assert resp == (nm_152263_exons, 1, 8)

    resp = await test_db.get_tx_exon_start_end('NM_152263.3', 0, 0)
    assert resp == (nm_152263_exons, 1, 10)

    resp = await test_db.get_tx_exon_start_end('NM_152263.3', 8, 1)
    assert resp is None


@pytest.mark.asyncio
async def test_get_tx_exon_coords(test_db, nm_152263_exons):
    """Test that get_tx_exon_coords works correctly."""
    resp = test_db.get_tx_exon_coords(nm_152263_exons, 1, 8)
    assert resp == (["117", "234"], ["822", "892"])

    resp = test_db.get_tx_exon_coords(nm_152263_exons, 1, 11)
    assert resp is None


@pytest.mark.asyncio
async def test_get_alt_ac_start_and_end(test_db, tpm3_1_8_start_genomic,
                                        tpm3_1_8_end_genomic):
    """Test that get_alt_ac_start_and_end works correctly."""
    resp = await test_db.get_alt_ac_start_and_end(
        'NM_152263.3', ["117", "234"], ["822", "892"], "TPM3")
    assert resp == (tpm3_1_8_start_genomic, tpm3_1_8_end_genomic)


@pytest.mark.asyncio
async def test_get_alt_ac_start_or_end(test_db, tpm3_1_8_start_genomic,
                                       tpm3_1_8_end_genomic):
    """Test that get_alt_ac_start_or_end works correctly."""
    resp = await test_db.get_alt_ac_start_or_end('NM_152263.3', 117, 234, None)
    assert resp == tpm3_1_8_start_genomic

    resp = await test_db.get_alt_ac_start_or_end('NM_152263.3', 822, 892, None)
    assert resp == tpm3_1_8_end_genomic


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
                                               exon_end_offset=-5)
    expected["exon_end"] = 8
    expected["exon_end_offset"] = -5
    expected["end"] = 154170404
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_152263.3', 0, 8,
                                               exon_end_offset=5)
    expected["exon_end_offset"] = 5
    expected["end"] = 154170394
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_152263.3', 3, 8,
                                               exon_start_offset=3,
                                               exon_end_offset=5)
    expected["exon_start"] = 3
    expected["exon_start_offset"] = 3
    expected["start"] = 154176245
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_152263.3', 3, 8,
                                               exon_start_offset=-3,
                                               exon_end_offset=5)
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
                                               exon_start_offset=3)
    expected = copy.deepcopy(ntrk1_exon10_exon17)
    expected["exon_start_offset"] = 3
    expected["start"] = 156874629
    assert resp == expected

    resp = await test_db.transcript_to_genomic('NM_002529.3', 10, 0,
                                               exon_start_offset=-3)
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
