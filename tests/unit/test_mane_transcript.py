"""Module for testing MANE Transcript class."""
import copy

import pytest
from mock import patch
import pandas as pd

from cool_seq_tool.mane_transcript import MANETranscript
from cool_seq_tool.data_sources import SeqRepoAccess, GeneNormalizer
from cool_seq_tool.mane_transcript import MANETranscriptError
from cool_seq_tool.schemas import AnnotationLayer, Assembly, ResidueMode


@pytest.fixture(scope="module")
async def test_mane_transcript(
    test_seqrepo_access, test_transcript_mappings, test_uta_db,
    test_mane_transcript_mappings,
):
    """Build mane transcript test fixture."""
    return MANETranscript(
        test_seqrepo_access, test_transcript_mappings, test_uta_db,
        test_mane_transcript_mappings, GeneNormalizer()
    )


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
        "chr_strand": "-"
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
        "coding_start_site": 226
    }


@pytest.fixture(scope="module")
def braf_v600e_mane_p():
    """Create test fixture for BRAF V600E MANE Transcript on p coordinate."""
    return {
        "refseq": "NP_004324.2",
        "ensembl": "ENSP00000493543.1",
        "pos": (599, 600),
        "status": "mane_select",
        "strand": "-",
        "gene": "BRAF"
    }


@pytest.fixture(scope="module")
def egfr_l858r_mane_p():
    """Create test fixture for EGFR L858R MANE Transcript on p coordinate."""
    return {
        "refseq": "NP_005219.2",
        "ensembl": "ENSP00000275493.2",
        "pos": (857, 858),
        "status": "mane_select",
        "strand": "+",
        "gene": "EGFR"
    }


@pytest.fixture(scope="module")
def braf_v600e_mane_c():
    """Create test fixture for BRAF V600E MANE Transcript on c coordinate."""
    return {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_004333.6",
        "ensembl": "ENST00000646891.2",
        "pos": (1798, 1799),
        "status": "mane_select",
        "strand": "-",
        "coding_start_site": 226,
        "coding_end_site": 2527,
        "gene": "BRAF"
    }


@pytest.fixture(scope="module")
def egfr_l858r_mane_c():
    """Create test fixture for EGFR L858R MANE Transcript on c coordinate."""
    return {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_005228.5",
        "ensembl": "ENST00000275493.7",
        "pos": (2572, 2573),
        "status": "mane_select",
        "strand": "+",
        "coding_start_site": 261,
        "coding_end_site": 3894,
        "gene": "EGFR"
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
        "alt_ac": "NC_000007.14"
    }


@pytest.mark.asyncio
async def test_p_to_mane_p(test_mane_transcript, braf_v600e_mane_p,
                           egfr_l858r_mane_p):
    """Test that p_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2", 599, 600, AnnotationLayer.PROTEIN, residue_mode="inter-residue")
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2", 600, 600, AnnotationLayer.PROTEIN)
    assert mane_p == braf_v600e_mane_p

    # BRAF V600E Ensembl Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7", 599, 600, AnnotationLayer.PROTEIN,
        residue_mode="inter-residue")
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7", 600, 600, AnnotationLayer.PROTEIN)
    assert mane_p == braf_v600e_mane_p

    # EGFR L858R RefSeq Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_005219.2", 858, 858, AnnotationLayer.PROTEIN)
    assert mane_p == egfr_l858r_mane_p

    # EGFR L858R Ensembl Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000275493.2", 858, 858, AnnotationLayer.PROTEIN)
    assert mane_p == egfr_l858r_mane_p

    assert test_mane_transcript.get_mane_transcript(
        "NP_004439.2", 755, 759, AnnotationLayer.PROTEIN)

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000366997.4", 63, 63, AnnotationLayer.PROTEIN, gene="DIS3", ref="P",
        try_longest_compatible=True)
    assert mane_p == {
        "gene": "DIS3",
        "refseq": "NP_055768.3",
        "ensembl": "ENSP00000366997.4",
        "pos": (62, 63),
        "strand": "-",
        "status": "mane_select"
    }


@pytest.mark.asyncio
async def test_c_to_mane_c(test_mane_transcript, braf_v600e_mane_c,
                           egfr_l858r_mane_c):
    """Test that c_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    cpy_braf_v600e_mane_c = copy.deepcopy(braf_v600e_mane_c)
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.4", 1799, 1799, AnnotationLayer.CDNA)
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.4", 1798, 1799, AnnotationLayer.CDNA, residue_mode="inter-residue")
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.5", 1799, 1799, AnnotationLayer.CDNA)
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.6", 1799, 1799, AnnotationLayer.CDNA)
    assert mane_c == cpy_braf_v600e_mane_c

    # BRAF V600E Ensembl Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000288602.10", 1799, 1799, AnnotationLayer.CDNA)
    cpy_braf_v600e_mane_c["alt_ac"] = "NC_000007.13"
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000288602.11", 1799, 1799, AnnotationLayer.CDNA)
    assert mane_c == cpy_braf_v600e_mane_c

    cpy_egfr_l858r_mane_c = copy.deepcopy(egfr_l858r_mane_c)
    # EGFR L858R RefSeq Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_005228.3", 2573, 2573, AnnotationLayer.CDNA)
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_005228.4", 2573, 2573, AnnotationLayer.CDNA)
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_005228.5", 2573, 2573, AnnotationLayer.CDNA)
    assert mane_c == cpy_egfr_l858r_mane_c

    # EGFR L858R Ensembl Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.7", 2573, 2573, AnnotationLayer.CDNA)
    cpy_egfr_l858r_mane_c["alt_ac"] = "NC_000007.13"
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.6", 2573, 2573, AnnotationLayer.CDNA)
    assert mane_c == cpy_egfr_l858r_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004448.3", 2264, 2278, AnnotationLayer.CDNA, ref="TGAGGGAAAACACAT")
    assert mane_c == {
        "alt_ac": "NC_000017.11",
        "refseq": "NM_004448.4",
        "ensembl": "ENST00000269571.10",
        "pos": (2263, 2278),
        "status": "mane_select",
        "strand": "+",
        "coding_start_site": 175,
        "coding_end_site": 3943,
        "gene": "ERBB2"
    }


@pytest.fixture(scope="function")
@patch.object(SeqRepoAccess, "get_reference_sequence")
def test__get_prioritized_transcripts_from_gene(test_get_seqrepo, test_mane_transcript):
    """Test that _get_prioritized_transcripts_from_gene works as expected"""

    def get_reference_sequence(ac):
        """Return test response when getting sequence for a given accession"""
        DUMMY_TX_REF_SEQ = {
            "NM_004333.6": ("AC", None),
            "NM_004333.5": ("ACTG", None),
            "NM_001378472.1": ("A", None),
            "NM_001374258.2": ("A", None)
        }
        return DUMMY_TX_REF_SEQ[ac]

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
    test_df = pd.DataFrame(data, columns=["tx_ac", "len_of_tx", "alt_ac"])

    resp = test_mane_transcript._get_prioritized_transcripts_from_gene(test_df)
    assert resp == ["NM_004333.6", "NM_001374258.2", "NM_001378472.1"]


@pytest.mark.asyncio
async def test_get_longest_compatible_transcript(test_mane_transcript):
    """Test that get_longest_compatible_transcript method works as expected"""
    mane_transcripts = {"ENST00000646891.2", "NM_001374258.1",
                        "NM_004333.6", "ENST00000644969.2"}
    expected = {
        "refseq": "NP_001365396.1",
        "ensembl": None,
        "pos": (599, 600),
        "strand": "-",
        "status": "longest_compatible_remaining"
    }
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        "BRAF", 599, 600, start_annotation_layer=AnnotationLayer.PROTEIN,
        residue_mode="inter-residue", mane_transcripts=mane_transcripts)
    assert resp == expected

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        "BRAF", 600, 600, start_annotation_layer=AnnotationLayer.PROTEIN,
        residue_mode="residue", mane_transcripts=mane_transcripts)
    assert resp == expected

    expected = {
        "refseq": "NM_001378467.1",
        "ensembl": None,
        "pos": (1798, 1799),
        "strand": "-",
        "status": "longest_compatible_remaining"
    }
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        "BRAF", 1799, 1799, start_annotation_layer=AnnotationLayer.CDNA,
        mane_transcripts=mane_transcripts)
    assert resp == expected

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        "BRAF", 1798, 1799, start_annotation_layer=AnnotationLayer.CDNA,
        residue_mode="inter-residue", mane_transcripts=mane_transcripts)
    assert resp == expected

    # NM_001378467.1:c.1808T>A (CA123643)
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        "BRAF", 140453136, 140453136, start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts, alt_ac="NC_000007.13")
    assert resp == {
        "refseq": "NM_001378467.1",
        "ensembl": None,
        "pos": (1807, 1808),
        "strand": "-",
        "status": "longest_compatible_remaining"
    }


@pytest.mark.asyncio
async def test_g_to_mane_c(test_mane_transcript, egfr_l858r_mane_c,
                           braf_v600e_mane_c, grch38):
    """Test that g_to_mane_c method works correctly."""
    # NM_005228.5:c.2573T>G
    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259515, 55259515, gene="EGFR")
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259514, 55259515, gene="EGFR",
        residue_mode="inter-residue")
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259515, 55259515, gene="EGFR")
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 140453135, 140453136, gene="BRAF",
        residue_mode="inter-residue")
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 140453136, 140453136, AnnotationLayer.GENOMIC, gene="BRAF")
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 140453135, 140453136, AnnotationLayer.GENOMIC, gene="BRAF",
        residue_mode="inter-residue")
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 140453136, 140453136, gene="BRAF")
    assert mane_c == braf_v600e_mane_c

    resp = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259515, 55259515)
    assert resp == grch38

    resp = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 55259514, 55259515, AnnotationLayer.GENOMIC,
        residue_mode="inter-residue"
    )
    assert resp == grch38

    resp = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 55259515, 55259515, AnnotationLayer.GENOMIC)
    assert resp == grch38

    resp = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 140453136, 140453136)
    grch38["pos"] = (140753335, 140753335)
    assert resp == grch38

    resp = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 140453135, 140453136,
        residue_mode="inter-residue")
    assert resp == grch38

    resp = await test_mane_transcript.g_to_mane_c(
        "NC_000007.14", 140753336, 140753336)
    grch38["pos"] = (140753335, 140753335)
    assert resp == grch38

    # CA122528
    mane_c = await test_mane_transcript.g_to_mane_c("NC_000012.11", 25398284,
                                                    25398284, gene="KRAS")
    assert mane_c == {
        "alt_ac": "NC_000012.12",
        "refseq": "NM_004985.5",
        "ensembl": "ENST00000311936.8",
        "pos": (34, 35),
        "status": "mane_select",
        "strand": "-",
        "coding_start_site": 190,
        "coding_end_site": 757,
        "gene": "KRAS"
    }

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 55249071, 55249071, AnnotationLayer.GENOMIC, "EGFR")
    assert mane_c == {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_005228.5",
        "ensembl": "ENST00000275493.7",
        "pos": (2368, 2369),
        "status": "mane_select",
        "strand": "+",
        "coding_start_site": 261,
        "coding_end_site": 3894,
        "gene": "EGFR"
    }


@pytest.mark.asyncio
async def test_get_mapped_mane_data(test_mane_transcript):
    """Test that get_mapped_mane_data works correctly"""
    resp = await test_mane_transcript.get_mapped_mane_data(
        "braf", Assembly.GRCH38, 140785808, ResidueMode.INTER_RESIDUE)
    assert resp.dict() == {
        "gene": "BRAF",
        "refseq": "NM_001374258.1",
        "ensembl": "ENST00000644969.2",
        "strand": "-",
        "status": "mane_plus_clinical",
        "alt_ac": "NC_000007.14",
        "assembly": "GRCh38"
    }

    resp = await test_mane_transcript.get_mapped_mane_data(
        "Braf", Assembly.GRCH37, 140485608, ResidueMode.INTER_RESIDUE)
    assert resp.dict() == {
        "gene": "BRAF",
        "refseq": "NM_001374258.1",
        "ensembl": "ENST00000644969.2",
        "strand": "-",
        "status": "mane_plus_clinical",
        "alt_ac": "NC_000007.13",
        "assembly": "GRCh37"
    }

    resp = await test_mane_transcript.get_mapped_mane_data(
        "BRAF", Assembly.GRCH38, 140783157, ResidueMode.INTER_RESIDUE)
    assert resp.dict() == {
        "gene": "BRAF",
        "refseq": "NM_004333.6",
        "ensembl": "ENST00000646891.2",
        "strand": "-",
        "status": "mane_select",
        "alt_ac": "NC_000007.14",
        "assembly": "GRCh38"
    }

    resp = await test_mane_transcript.get_mapped_mane_data(
        "BRAF", Assembly.GRCH37, 140482958, ResidueMode.RESIDUE)
    assert resp.dict() == {
        "gene": "BRAF",
        "refseq": "NM_004333.6",
        "ensembl": "ENST00000646891.2",
        "strand": "-",
        "status": "mane_select",
        "alt_ac": "NC_000007.13",
        "assembly": "GRCh37"
    }

    # Invalid coord given assembly, so no result should be found
    resp = await test_mane_transcript.get_mapped_mane_data(
        "BRAF", Assembly.GRCH38, 140482957, ResidueMode.INTER_RESIDUE)
    assert resp is None

    # Invalid gene
    with pytest.raises(MANETranscriptError) as e:
        await test_mane_transcript.get_mapped_mane_data(
            "dummy", Assembly.GRCH37, 140482958, ResidueMode.RESIDUE)
    assert str(e.value) == "Unable to get HGNC data for gene: dummy"


@pytest.mark.asyncio
async def test_valid(test_mane_transcript):
    """Test that valid queries do not raise any exceptions"""
    # CA126067
    resp = await test_mane_transcript.get_mane_transcript(
        "NP_000507.1", 201, 201, AnnotationLayer.PROTEIN, ref="R",
        try_longest_compatible=True)
    assert resp


@pytest.mark.asyncio
async def test_no_matches(test_mane_transcript):
    """Test that invalid queries return None."""
    # Invalid ENST version
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.15645", 2573, 2573, AnnotationLayer.CDNA
    )
    assert mane_c is None

    # Invalid residue-mode
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000288602.11", 2573, 2573, AnnotationLayer.CDNA, residue_mode="residues"
    )
    assert mane_c is None
