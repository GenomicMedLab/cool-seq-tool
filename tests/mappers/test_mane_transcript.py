"""Module for testing MANE Transcript class."""
import copy

import polars as pl
import pytest
from mock import patch

from cool_seq_tool.handlers.seqrepo_access import SeqRepoAccess
from cool_seq_tool.schemas import AnnotationLayer, ResidueMode


@pytest.fixture(scope="module")
def test_mane_transcript(test_cool_seq_tool):
    """Build mane transcript test fixture."""
    return test_cool_seq_tool.mane_transcript


@pytest.fixture(scope="module")
def braf_mane_data():
    """Create test fixture for BRAF MANE data."""
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
        "GRCh38_chr": "7",
        "chr_start": 140730665,
        "chr_end": 140924929,
        "chr_strand": "-",
    }


@pytest.fixture(scope="module")
def nm_004333v6_g():
    """Create test fixture for NM_004333.6 genomic data."""
    return {
        "gene": "BRAF",
        "tx_ac": "NM_004333.6",
        "tx_pos_range": (1967, 2086),
        "alt_ac": "NC_000007.14",
        "alt_pos_change_range": (140753336, 140753334),
        "alt_pos_range": (140753274, 140753393),
        "pos_change": (57, 60),
        "strand": "-",
        "alt_aln_method": "splign",
        "tx_exon_id": 7649345,
        "alt_exon_id": 9507338,
        "coding_start_site": 226,
    }


@pytest.fixture(scope="module")
def braf_v600e_mane_p():
    """Create test fixture for BRAF V600E MANE Transcript on p coordinate."""
    return {
        "refseq": "NP_004324.2",
        "ensembl": "ENSP00000493543.1",
        "pos": (599, 599),
        "status": "mane_select",
        "strand": "-",
        "gene": "BRAF",
    }


@pytest.fixture(scope="module")
def egfr_l858r_mane_p():
    """Create test fixture for EGFR L858R MANE Transcript on p coordinate."""
    return {
        "refseq": "NP_005219.2",
        "ensembl": "ENSP00000275493.2",
        "pos": (857, 857),
        "status": "mane_select",
        "strand": "+",
        "gene": "EGFR",
    }


@pytest.fixture(scope="module")
def braf_v600e_mane_c():
    """Create test fixture for BRAF V600E MANE Transcript on c coordinate."""
    return {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_004333.6",
        "ensembl": "ENST00000646891.2",
        "pos": (1798, 1798),
        "status": "mane_select",
        "strand": "-",
        "coding_start_site": 226,
        "coding_end_site": 2527,
        "gene": "BRAF",
    }


@pytest.fixture(scope="module")
def egfr_l858r_mane_c():
    """Create test fixture for EGFR L858R MANE Transcript on c coordinate."""
    return {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_005228.5",
        "ensembl": "ENST00000275493.7",
        "pos": (2572, 2572),
        "status": "mane_select",
        "strand": "+",
        "coding_start_site": 261,
        "coding_end_site": 3894,
        "gene": "EGFR",
    }


@pytest.fixture(scope="module")
def grch38():
    """Create a test fixture for grch38 responses."""
    return {
        "gene": None,
        "refseq": "NC_000007.14",
        "ensembl": None,
        "coding_start_site": None,
        "coding_end_site": None,
        "pos": (55191821, 55191821),
        "strand": None,
        "status": "GRCh38",
        "alt_ac": "NC_000007.14",
    }


@pytest.fixture(scope="module")
def mybpc3_s236g():
    """Create test fixture for MYBPC3 Ser236Gly

    CA1139661942
    https://www.ncbi.nlm.nih.gov/clinvar/variation/922707/?new_evidence=true
    """
    return {
        "refseq": "NP_000247.2",
        "ensembl": "ENSP00000442795.1",
        "pos": (235, 235),
        "status": "mane_select",
        "strand": "-",
        "gene": "MYBPC3",
    }


def test__get_reading_frame(test_mane_transcript):
    """Test that _get_reading_frame works correctly."""
    rf = test_mane_transcript._get_reading_frame(1797)
    assert rf == 3

    rf = test_mane_transcript._get_reading_frame(1798)
    assert rf == 1

    rf = test_mane_transcript._get_reading_frame(1799)
    assert rf == 2

    rf = test_mane_transcript._get_reading_frame(1800)
    assert rf == 3

    rf = test_mane_transcript._get_reading_frame(2573)
    assert rf == 2


def test_p_to_c_pos(test_mane_transcript):
    """Test that _p_to_c_pos method works correctly."""
    # https://civicdb.org/events/genes/5/summary/variants/12/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(600, 600)
    assert c_pos == (1798, 1800)

    # https://civicdb.org/events/genes/19/summary/variants/33/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(858, 858)
    assert c_pos == (2572, 2574)

    # https://civicdb.org/events/genes/29/summary/variants/72/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(576, 576)
    assert c_pos == (1726, 1728)

    # https://civicdb.org/events/genes/38/summary/variants/99/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(842, 842)
    assert c_pos == (2524, 2526)

    # https://civicdb.org/events/genes/30/summary/variants/78/summary#variant
    c_pos = test_mane_transcript._p_to_c_pos(12, 12)
    assert c_pos == (34, 36)

    c_pos = test_mane_transcript._p_to_c_pos(747, 751)
    assert c_pos == (2239, 2253)


@pytest.mark.asyncio
async def test_p_to_c(test_mane_transcript):
    """Test that _p_to_c method works correctly."""
    # Amino Acid Substitution
    expected_pos = 1798, 1800
    ac, pos = await test_mane_transcript._p_to_c("NP_004324.2", 600, 600)
    assert ac == "NM_004333.6"
    assert pos == expected_pos

    ac, pos = await test_mane_transcript._p_to_c("ENSP00000288602.7", 600, 600)
    assert ac == "ENST00000288602.11"
    assert pos == expected_pos

    expected_pos = 2572, 2574
    ac, pos = await test_mane_transcript._p_to_c("NP_005219.2", 858, 858)
    assert ac == "NM_005228.5"
    assert pos == expected_pos

    ac, pos = await test_mane_transcript._p_to_c("ENSP00000275493.2", 858, 858)
    assert ac == "ENST00000275493.7"
    assert pos == expected_pos

    # Polypeptide Truncation
    expected_pos = 553, 555
    ac, pos = await test_mane_transcript._p_to_c("NP_000542.1", 185, 185)
    assert ac == "NM_000551.4"
    assert pos == expected_pos

    ac, pos = await test_mane_transcript._p_to_c("ENSP00000256474.3", 185, 185)
    assert ac == "ENST00000256474.3"
    assert pos == expected_pos

    # Silent Mutation
    expected_pos = 181, 183
    ac, pos = await test_mane_transcript._p_to_c("NP_000542.1", 61, 61)
    assert ac == "NM_000551.4"
    assert pos == expected_pos


@pytest.mark.asyncio
async def test_c_to_g(test_mane_transcript, nm_004333v6_g):
    """Test that _c_to_g method works correctly."""
    tx_ac = "NM_004333.6"
    g = await test_mane_transcript._c_to_g(tx_ac, (1798, 1800))
    assert g == nm_004333v6_g


@pytest.mark.asyncio
async def test__g_to_c(
    test_mane_transcript, braf_mane_data, nm_004333v6_g, braf_v600e_mane_c
):
    """Test that _g_to_c method works correctly."""
    mane_c = await test_mane_transcript._g_to_c(
        g=nm_004333v6_g,
        refseq_c_ac=braf_mane_data["RefSeq_nuc"],
        status="_".join(braf_mane_data["MANE_status"].split()).lower(),
        ensembl_c_ac=braf_mane_data["Ensembl_nuc"],
    )
    expected = copy.deepcopy(braf_v600e_mane_c)
    expected["pos"] = (1798, 1800)
    expected["alt_ac"] = None
    assert mane_c == expected


def test_get_mane_p(test_mane_transcript, braf_mane_data, braf_v600e_mane_p):
    """Test that _get_mane_p method works correctly."""
    mane_p = test_mane_transcript._get_mane_p(braf_mane_data, (1797, 1799))
    assert mane_p == braf_v600e_mane_p


@pytest.mark.asyncio
async def test_p_to_mane_p(test_mane_transcript, braf_v600e_mane_p, egfr_l858r_mane_p):
    """Test that p_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2",
        599,
        AnnotationLayer.PROTEIN,
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2", 600, AnnotationLayer.PROTEIN
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2",
        599,
        AnnotationLayer.PROTEIN,
        residue_mode=ResidueMode.INTER_RESIDUE,
        end_pos=599,
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2", 600, AnnotationLayer.PROTEIN, end_pos=600
    )
    assert mane_p == braf_v600e_mane_p

    # BRAF V600E Ensembl Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7",
        599,
        AnnotationLayer.PROTEIN,
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7", 600, AnnotationLayer.PROTEIN
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7",
        599,
        AnnotationLayer.PROTEIN,
        residue_mode=ResidueMode.INTER_RESIDUE,
        end_pos=599,
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7", 600, AnnotationLayer.PROTEIN, end_pos=600
    )
    assert mane_p == braf_v600e_mane_p

    # EGFR L858R RefSeq Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_005219.2", 858, AnnotationLayer.PROTEIN
    )
    assert mane_p == egfr_l858r_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_005219.2", 858, AnnotationLayer.PROTEIN, end_pos=858
    )
    assert mane_p == egfr_l858r_mane_p

    # EGFR L858R Ensembl Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000275493.2", 858, AnnotationLayer.PROTEIN
    )
    assert mane_p == egfr_l858r_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000275493.2", 858, AnnotationLayer.PROTEIN, end_pos=858
    )
    assert mane_p == egfr_l858r_mane_p

    assert test_mane_transcript.get_mane_transcript(
        "NP_004439.2", 755, AnnotationLayer.PROTEIN, end_pos=759
    )

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000366997.4",
        63,
        AnnotationLayer.PROTEIN,
        gene="DIS3",
        ref="P",
        try_longest_compatible=True,
        end_pos=63,
    )
    assert mane_p == {
        "gene": "DIS3",
        "refseq": "NP_055768.3",
        "ensembl": "ENSP00000366997.4",
        "pos": (62, 62),
        "strand": "-",
        "status": "mane_select",
    }


@pytest.mark.asyncio
async def test_c_to_mane_c(test_mane_transcript, braf_v600e_mane_c, egfr_l858r_mane_c):
    """Test that c_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    cpy_braf_v600e_mane_c = copy.deepcopy(braf_v600e_mane_c)
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.4", 1799, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.4",
        1798,
        AnnotationLayer.CDNA,
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.4",
        1798,
        AnnotationLayer.CDNA,
        residue_mode=ResidueMode.INTER_RESIDUE,
        end_pos=1798,
    )
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.5", 1799, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.6", 1799, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_braf_v600e_mane_c

    # BRAF V600E Ensembl Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000288602.10", 1799, AnnotationLayer.CDNA
    )
    cpy_braf_v600e_mane_c["alt_ac"] = "NC_000007.13"
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000288602.11", 1799, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_braf_v600e_mane_c

    cpy_egfr_l858r_mane_c = copy.deepcopy(egfr_l858r_mane_c)
    # EGFR L858R RefSeq Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_005228.3", 2573, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_005228.4", 2573, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_005228.5", 2573, AnnotationLayer.CDNA, end_pos=2573
    )
    assert mane_c == cpy_egfr_l858r_mane_c

    # EGFR L858R Ensembl Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.7", 2573, AnnotationLayer.CDNA
    )
    cpy_egfr_l858r_mane_c["alt_ac"] = "NC_000007.13"
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.6", 2573, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004448.3", 2264, AnnotationLayer.CDNA, end_pos=2278, ref="TGAGGGAAAACACAT"
    )
    assert mane_c == {
        "alt_ac": "NC_000017.11",
        "refseq": "NM_004448.4",
        "ensembl": "ENST00000269571.10",
        "pos": (2263, 2277),
        "status": "mane_select",
        "strand": "+",
        "coding_start_site": 175,
        "coding_end_site": 3943,
        "gene": "ERBB2",
    }


@pytest.fixture(scope="function")
@patch.object(SeqRepoAccess, "get_reference_sequence")
def test__get_prioritized_transcripts_from_gene(test_get_seqrepo, test_mane_transcript):
    """Test that _get_prioritized_transcripts_from_gene works as expected"""

    def get_reference_sequence(ac):
        """Return test response when getting sequence for a given accession"""
        dummy_tx_ref_seq = {
            "NM_004333.6": ("AC", None),
            "NM_004333.5": ("ACTG", None),
            "NM_001378472.1": ("A", None),
            "NM_001374258.2": ("A", None),
        }
        return dummy_tx_ref_seq[ac]

    test_get_seqrepo.return_value = None
    test_mane_transcript.seqrepo_access.get_reference_sequence = get_reference_sequence

    data = [
        ["NM_004333.6", 2, "NC_000007.13"],
        ["NM_004333.5", 4, "NC_000007.13"],
        ["NM_001378472.1", 1, "NC_000007.13"],
        ["NM_001374258.2", 1, "NC_000007.13"],
        ["NM_004333.6", 2, "NC_000007.14"],
        ["NM_004333.5", 4, "NC_000007.14"],
        ["NM_001378472.1", 1, "NC_000007.14"],
        ["NM_001374258.2", 1, "NC_000007.14"],
    ]
    test_df = pl.DataFrame(data, schema=["tx_ac", "len_of_tx", "alt_ac"])

    resp = test_mane_transcript._get_prioritized_transcripts_from_gene(test_df)
    assert resp == ["NM_004333.6", "NM_001374258.2", "NM_001378472.1"]


@pytest.mark.asyncio
async def test_get_longest_compatible_transcript(test_mane_transcript):
    """Test that get_longest_compatible_transcript method works as expected"""
    mane_transcripts = {
        "ENST00000646891.2",
        "NM_001374258.1",
        "NM_004333.6",
        "ENST00000644969.2",
    }
    expected = {
        "refseq": "NP_001365396.1",
        "ensembl": None,
        "pos": (599, 599),
        "strand": "-",
        "status": "longest_compatible_remaining",
    }
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        599,
        599,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.PROTEIN,
        residue_mode=ResidueMode.INTER_RESIDUE,
        mane_transcripts=mane_transcripts,
    )
    assert resp == expected

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        600,
        600,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.PROTEIN,
        residue_mode=ResidueMode.RESIDUE,
        mane_transcripts=mane_transcripts,
    )
    assert resp == expected

    expected = {
        "refseq": "NM_001378467.1",
        "ensembl": None,
        "pos": (1798, 1798),
        "strand": "-",
        "status": "longest_compatible_remaining",
    }
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        1799,
        1799,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.CDNA,
        mane_transcripts=mane_transcripts,
    )
    assert resp == expected

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        1798,
        1798,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.CDNA,
        residue_mode=ResidueMode.INTER_RESIDUE,
        mane_transcripts=mane_transcripts,
    )
    assert resp == expected

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        140453136,
        140453136,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts,
        alt_ac="NC_000007.13",
    )
    assert resp == {
        "refseq": "NM_001378467.1",
        "ensembl": None,
        "pos": (1807, 1807),
        "strand": "-",
        "status": "longest_compatible_remaining",
    }

    #  CA1139661942 has no other RefSeq accessions
    mane_transcripts = {"NM_000256.3", "ENST00000545968.6"}
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        47348490,
        47348490,
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts,
        alt_ac="NC_000011.10",
    )
    assert resp is None

    mane_transcripts = {"NM_005228.5", "ENST00000275493.7"}
    # CDNA
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        55174776,
        55174793,
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts,
        alt_ac="NC_000007.14",
    )
    assert resp == {
        "refseq": "NM_001346899.2",
        "ensembl": None,
        "pos": (2103, 2120),
        "strand": "+",
        "status": "longest_compatible_remaining",
    }

    # protein
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        55174776,
        55174793,
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts,
        alt_ac="NC_000007.14",
        end_annotation_layer=AnnotationLayer.PROTEIN,
    )
    assert resp == {
        "refseq": "NP_001333828.1",
        "ensembl": None,
        "pos": (701, 706),
        "strand": "+",
        "status": "longest_compatible_remaining",
    }

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        153870419,
        153870476,
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts={"ENST00000370060.7", "NM_001278116.2"},
        alt_ac="NC_000023.11",
        end_annotation_layer=AnnotationLayer.PROTEIN,
    )
    assert resp == {
        "refseq": "NP_000416.1",
        "ensembl": None,
        "pos": (239, 258),
        "strand": "-",
        "status": "longest_compatible_remaining",
    }


@pytest.mark.asyncio
async def test_g_to_mane_c(
    test_mane_transcript, egfr_l858r_mane_c, braf_v600e_mane_c, grch38
):
    """Test that g_to_mane_c method works correctly."""
    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259515, None, gene="EGFR"
    )
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13",
        55259514,
        None,
        gene="EGFR",
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259515, 55259515, gene="EGFR"
    )
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13",
        140453135,
        None,
        gene="BRAF",
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 140453136, AnnotationLayer.GENOMIC, gene="BRAF"
    )
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13",
        140453135,
        AnnotationLayer.GENOMIC,
        gene="BRAF",
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 140453136, None, gene="BRAF"
    )
    assert mane_c == braf_v600e_mane_c

    resp = await test_mane_transcript.g_to_mane_c("NC_000007.13", 55259515, None)
    assert resp == grch38

    resp = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13",
        55259514,
        AnnotationLayer.GENOMIC,
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert resp == grch38

    resp = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 55259515, AnnotationLayer.GENOMIC
    )
    assert resp == grch38

    resp = await test_mane_transcript.g_to_mane_c("NC_000007.13", 140453136, None)
    grch38["pos"] = (140753335, 140753335)
    assert resp == grch38

    resp = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 140453135, None, residue_mode=ResidueMode.INTER_RESIDUE
    )
    assert resp == grch38

    resp = await test_mane_transcript.g_to_mane_c("NC_000007.14", 140753336, None)
    grch38["pos"] = (140753335, 140753335)
    assert resp == grch38

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000012.11", 25398284, None, gene="KRAS"
    )
    assert mane_c == {
        "alt_ac": "NC_000012.12",
        "refseq": "NM_004985.5",
        "ensembl": "ENST00000311936.8",
        "pos": (34, 34),
        "status": "mane_select",
        "strand": "-",
        "coding_start_site": 190,
        "coding_end_site": 757,
        "gene": "KRAS",
    }

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 55249071, AnnotationLayer.GENOMIC, 55249071, "EGFR"
    )
    assert mane_c == {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_005228.5",
        "ensembl": "ENST00000275493.7",
        "pos": (2368, 2368),
        "status": "mane_select",
        "strand": "+",
        "coding_start_site": 261,
        "coding_end_site": 3894,
        "gene": "EGFR",
    }


@pytest.mark.asyncio
async def test_grch38_to_mane_p(test_mane_transcript, mybpc3_s236g):
    """Test that grch38_to_mane_p"""
    # https://www.ncbi.nlm.nih.gov/clinvar/variation/922707/?new_evidence=true
    # Without gene
    resp = await test_mane_transcript.grch38_to_mane_p(
        "NC_000011.10", 47348490, 47348490
    )
    assert resp == mybpc3_s236g

    # With gene
    resp = await test_mane_transcript.grch38_to_mane_p(
        "NC_000011.10", 47348490, 47348491, gene="MYBPC3"
    )
    assert resp == mybpc3_s236g

    # CA645561524 (without gene)
    resp = await test_mane_transcript.grch38_to_mane_p(
        "NC_000007.14", 55174776, 55174793
    )
    resp == {
        "refseq": "NP_005219.2",
        "ensembl": "ENSP00000275493.2",
        "pos": (746, 751),
        "status": "mane_select",
        "strand": "+",
        "gene": " EGFR",
    }

    # CA2499226460 (without gene)
    resp = await test_mane_transcript.grch38_to_mane_p(
        "NC_000023.11", 153870419, 153870476
    )
    assert resp == {
        "refseq": "NP_001265045.1",
        "ensembl": "ENSP00000359077.1",
        "pos": (239, 258),
        "status": "mane_select",
        "strand": "-",
        "gene": "L1CAM",
    }

    # GENE not valid
    resp = await test_mane_transcript.grch38_to_mane_p(
        "NC_000023.11", 153870419, 153870476, gene="FAKE"
    )
    assert resp is None


@pytest.mark.asyncio
async def test_valid(test_mane_transcript):
    """Test that valid queries do not raise any exceptions"""
    resp = await test_mane_transcript.get_mane_transcript(
        "NP_001296812.1",
        201,
        AnnotationLayer.PROTEIN,
        end_pos=201,
        ref="R",
        try_longest_compatible=True,
    )
    assert resp


@pytest.mark.asyncio
async def test_no_matches(test_mane_transcript):
    """Test that invalid queries return None."""
    # Invalid ENST version
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.15645", 2573, AnnotationLayer.CDNA
    )
    assert mane_c is None

    # Invalid residue-mode
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000288602.11", 2573, AnnotationLayer.CDNA, residue_mode="residues"
    )
    assert mane_c is None
