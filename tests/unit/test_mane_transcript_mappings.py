"""Module for testing MANE Transcript Mapping class."""
import pytest


@pytest.fixture(scope="module")
def braf_select():
    """Create test fixture for BRAF MANE Select Transcript data."""
    return {
        "#NCBI_GeneID": "GeneID:673",
        "Ensembl_Gene": "ENSG00000157764.14",
        "HGNC_ID": "HGNC:1097",
        "symbol": "BRAF",
        "name": "B-Raf proto-oncogene, serine/threonine kinase",
        "RefSeq_nuc": "NM_004333.6",
        "RefSeq_prot": "NP_004324.2",
        "Ensembl_nuc": "ENST00000646891.2",
        "Ensembl_prot": "ENSP00000493543.1",
        "MANE_status": "MANE Select",
        "GRCh38_chr": "NC_000007.14",
        "chr_start": 140730665,
        "chr_end": 140924929,
        "chr_strand": "-"
    }


@pytest.fixture(scope="module")
def braf_plus_clinical():
    """Create test fixture for BRAF MANE Plus Clinical data."""
    return {
        "#NCBI_GeneID": "GeneID:673",
        "Ensembl_Gene": "ENSG00000157764.14",
        "HGNC_ID": "HGNC:1097",
        "symbol": "BRAF",
        "name": "B-Raf proto-oncogene, serine/threonine kinase",
        "RefSeq_nuc": "NM_001374258.1",
        "RefSeq_prot": "NP_001361187.1",
        "Ensembl_nuc": "ENST00000644969.2",
        "Ensembl_prot": "ENSP00000496776.1",
        "MANE_status": "MANE Plus Clinical",
        "GRCh38_chr": "NC_000007.14",
        "chr_start": 140719337,
        "chr_end": 140924929,
        "chr_strand": "-"
    }


@pytest.fixture(scope="module")
def ercc6_plus_clinical():
    """Create test fixture for ERCC6 MANE Plus Clinical Transcript data."""
    return {
        "#NCBI_GeneID": "GeneID:2074",
        "Ensembl_Gene": "ENSG00000225830.16",
        "HGNC_ID": "HGNC:3438",
        "symbol": "ERCC6",
        "name": "ERCC excision repair 6, chromatin remodeling factor",
        "RefSeq_nuc": "NM_001277058.2",
        "RefSeq_prot": "NP_001263987.1",
        "Ensembl_nuc": "ENST00000447839.7",
        "Ensembl_prot": "ENSP00000387966.2",
        "MANE_status": "MANE Plus Clinical",
        "GRCh38_chr": "NC_000010.11",
        "chr_start": 49515198,
        "chr_end": 49539121,
        "chr_strand": "-"
    }


@pytest.fixture(scope="module")
def ercc6_select():
    """Create test fixture for ERCC6 MANE Select Transcript data."""
    return {
        "#NCBI_GeneID": "GeneID:2074",
        "Ensembl_Gene": "ENSG00000225830.16",
        "HGNC_ID": "HGNC:3438",
        "symbol": "ERCC6",
        "name": "ERCC excision repair 6, chromatin remodeling factor",
        "RefSeq_nuc": "NM_000124.4",
        "RefSeq_prot": "NP_000115.1",
        "Ensembl_nuc": "ENST00000355832.10",
        "Ensembl_prot": "ENSP00000348089.5",
        "MANE_status": "MANE Select",
        "GRCh38_chr": "NC_000010.11",
        "chr_start": 49454470,
        "chr_end": 49539121,
        "chr_strand": "-"
    }


def test_get_gene_mane_data(test_mane_transcript_mappings, braf_select,
                            braf_plus_clinical, ercc6_select,
                            ercc6_plus_clinical):
    """Test that get_gene_mane_data method works correctly."""
    # MANE Select
    actual = test_mane_transcript_mappings.get_gene_mane_data("BRAF")
    assert len(actual) == 2
    assert actual[0] == braf_select
    assert actual[1] == braf_plus_clinical

    actual = test_mane_transcript_mappings.get_gene_mane_data("braf")
    assert len(actual) == 2
    assert actual[0] == braf_select
    assert actual[1] == braf_plus_clinical

    # MANE Select and MANE Plus Clinical
    actual = test_mane_transcript_mappings.get_gene_mane_data("ERCC6")
    assert len(actual) == 2
    assert actual[0] == ercc6_select
    assert actual[1] == ercc6_plus_clinical

    actual = test_mane_transcript_mappings.get_gene_mane_data("ercc6")
    assert actual[0] == ercc6_select
    assert actual[1] == ercc6_plus_clinical

    # No Matches
    actual = test_mane_transcript_mappings.get_gene_mane_data("BRAFF")
    assert actual == []

    actual = test_mane_transcript_mappings.get_gene_mane_data("")
    assert actual == []


def test_get_mane_from_transcripts(test_mane_transcript_mappings,
                                   braf_plus_clinical,
                                   braf_select, ercc6_plus_clinical):
    """Test that get get_mane_from_transcripts method works correctly"""
    transcripts = [
        "NM_001354609.1", "NM_001354609.2", "NM_001374244.1", "NM_001374258.1",
        "NM_001378467.1", "NM_001378468.1", "NM_001378469.1", "NM_001378470.1",
        "NM_001378471.1", "NM_001378472.1", "NM_001378473.1", "NM_001378474.1",
        "NM_001378475.1", "NM_004333.4", "NM_004333.5", "NM_004333.6"
    ]
    resp = test_mane_transcript_mappings.get_mane_from_transcripts(transcripts)
    assert len(resp) == 2
    assert braf_select in resp
    assert braf_plus_clinical in resp

    transcripts.append("NM_001277058.2")
    resp = test_mane_transcript_mappings.get_mane_from_transcripts(transcripts)
    assert len(resp) == 3
    assert braf_select in resp
    assert braf_plus_clinical in resp
    assert ercc6_plus_clinical in resp

    # Invalid transcripts
    resp = test_mane_transcript_mappings.get_mane_from_transcripts(
        ["NM_012334.34"])
    assert resp == []


def test_get_mane_data_from_chr_pos(test_mane_transcript_mappings, braf_select,
                                    braf_plus_clinical):
    """Test that get_mane_data_from_chr_pos method works correctly"""
    resp = test_mane_transcript_mappings.get_mane_data_from_chr_pos(
        "NC_000007.14", 140753336, 140753336)
    assert len(resp) == 2
    assert resp == [braf_select, braf_plus_clinical]

    resp = test_mane_transcript_mappings.get_mane_data_from_chr_pos(
        "NC_000023.11", 37994300, 37994310)
    assert len(resp) == 1
    assert resp == [{
        "#NCBI_GeneID": "GeneID:115482686",
        "Ensembl_Gene": "ENSG00000229674.3",
        "HGNC_ID": "HGNC:53960",
        "symbol": "H2AL3",
        "name": "H2A.L variant histone 3",
        "RefSeq_nuc": "NM_001395555.1",
        "RefSeq_prot": "NP_001382484.1",
        "Ensembl_nuc": "ENST00000448797.3",
        "Ensembl_prot": "ENSP00000498087.1",
        "MANE_status": "MANE Select",
        "GRCh38_chr": "NC_000023.11",
        "chr_start": 37994272,
        "chr_end": 37994904,
        "chr_strand": "+"
    }]

    # Invalid alt_ac (no version)
    resp = test_mane_transcript_mappings.get_mane_data_from_chr_pos(
        "NC_000023", 37994300, 37994310)
    assert resp == []
