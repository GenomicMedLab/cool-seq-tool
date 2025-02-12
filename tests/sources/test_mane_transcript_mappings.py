"""Module for testing MANE Transcript Mapping class."""

from unittest.mock import patch

import polars as pl
import pytest

from cool_seq_tool.schemas import ManeGeneData


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
        "chr_strand": "-",
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
        "chr_strand": "-",
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
        "chr_strand": "-",
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
        "chr_strand": "-",
    }


def test_get_gene_mane_data(
    test_mane_transcript_mappings,
    braf_select,
    braf_plus_clinical,
    ercc6_select,
    ercc6_plus_clinical,
):
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


def test_get_mane_from_transcripts(
    test_mane_transcript_mappings, braf_plus_clinical, braf_select, ercc6_plus_clinical
):
    """Test that get get_mane_from_transcripts method works correctly"""
    transcripts = [
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
    resp = test_mane_transcript_mappings.get_mane_from_transcripts(["NM_012334.34"])
    assert resp == []


def test_get_mane_data_from_chr_pos(
    test_mane_transcript_mappings, braf_select, braf_plus_clinical
):
    """Test that get_mane_data_from_chr_pos method works correctly"""
    resp = test_mane_transcript_mappings.get_mane_data_from_chr_pos(
        "NC_000007.14", 140753336, 140753336
    )
    assert len(resp) == 2
    assert resp == [braf_select, braf_plus_clinical]

    resp = test_mane_transcript_mappings.get_mane_data_from_chr_pos(
        "NC_000023.11", 37994300, 37994310
    )
    assert len(resp) == 1
    assert resp == [
        {
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
            "chr_strand": "+",
        }
    ]

    # Invalid alt_ac (no version)
    resp = test_mane_transcript_mappings.get_mane_data_from_chr_pos(
        "NC_000023", 37994300, 37994310
    )
    assert resp == []


def test_get_genomic_mane_genes(
    test_mane_transcript_mappings, braf_mane_genes, egfr_mane_gene
):
    """Test that get_genomic_mane_genes method works correctly"""
    new_df = pl.DataFrame(
        {
            "#NCBI_GeneID": [
                "GeneID:673",
                "GeneID:673",
                "GeneID:1956",
                "GeneID:1",
                "GeneID:2",
                "GeneID:2",
                "GeneID:3",
            ],
            "Ensembl_Gene": [
                "ENSG00000157764.14",
                "ENSG00000157764.14",
                "ENSG00000146648.21",
                "ENSG1.1",
                "ENSG1.1",
                "ENSG1.1",
                "ENSG1.1",
            ],
            "HGNC_ID": [
                "HGNC:1097",
                "HGNC:1097",
                "HGNC:3236",
                "HGNC:1",
                "HGNC:2",
                "HGNC:2",
                "HGNC:3",
            ],
            "symbol": ["BRAF", "BRAF", "EGFR", "Dummy1", "Dummy2", "Dummy2", "Dummy3"],
            "GRCh38_chr": [
                "NC_000007.14",
                "NC_000007.14",
                "NC_000007.14",
                "NC_000007.14",
                "NC_000007.14",
                "NC_000007.14",
                "NC_000007.14",
            ],
            "chr_start": [
                140719337,
                140730665,
                55019017,
                55019017,
                55019017,
                55019017,
                55019017,
            ],
            "chr_end": [
                140924929,
                140924929,
                55211628,
                55211628,
                55211628,
                55211628,
                55211628,
            ],
            "MANE_status": [
                "MANE Plus Clinical",
                "MANE Select",
                "MANE Select",
                "MANE Plus Clinical",
                "MANE Select",
                "MANE Plus Clinical",
                "MANE Select",
            ],
        }
    )

    with patch.object(test_mane_transcript_mappings, "df", new_df):
        mane_genes = test_mane_transcript_mappings.get_genomic_mane_genes(
            "NC_000007.14", 140753336, 140753336
        )
        assert mane_genes == braf_mane_genes

        mane_genes = test_mane_transcript_mappings.get_genomic_mane_genes(
            "NC_000007.14", 55191822, 55191822
        )
        assert mane_genes == [
            ManeGeneData(
                ncbi_gene_id=2,
                hgnc_id=2,
                symbol="Dummy2",
                status=["mane_select", "mane_plus_clinical"],
            ),
            ManeGeneData(
                ncbi_gene_id=3, hgnc_id=3, symbol="Dummy3", status=["mane_select"]
            ),
            egfr_mane_gene,
            ManeGeneData(
                ncbi_gene_id=1,
                hgnc_id=1,
                symbol="Dummy1",
                status=["mane_plus_clinical"],
            ),
        ]

        # No MANE genes found for given genomic location
        mane_genes = test_mane_transcript_mappings.get_genomic_mane_genes(
            "NC_000007.14", 140718337, 140718337
        )
        assert mane_genes == []
