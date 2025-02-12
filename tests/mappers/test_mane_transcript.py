"""Module for testing MANE Transcript class."""

from unittest.mock import patch

import polars as pl
import pytest

from cool_seq_tool.mappers.mane_transcript import (
    CdnaRepresentation,
    DataRepresentation,
    EndAnnotationLayer,
    GenomicRepresentation,
    ProteinAndCdnaRepresentation,
)
from cool_seq_tool.schemas import (
    AnnotationLayer,
    CoordinateType,
    Strand,
    TranscriptPriority,
)
from cool_seq_tool.sources.uta_database import GenomicTxMetadata


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
    params = {
        "gene": "BRAF",
        "tx_ac": "NM_004333.6",
        "tx_pos_range": (1967, 2086),
        "alt_ac": "NC_000007.14",
        "alt_pos_change_range": (140753336, 140753334),
        "alt_pos_range": (140753274, 140753393),
        "pos_change": (57, 60),
        "strand": Strand.NEGATIVE,
        "alt_aln_method": "splign",
        "tx_exon_id": 7649345,
        "alt_exon_id": 9507338,
        "coding_start_site": 226,
    }
    return GenomicTxMetadata(**params)


@pytest.fixture(scope="module")
def braf_v600e_mane_p():
    """Create test fixture for BRAF V600E MANE Transcript on p coordinate (CA123643)."""
    params = {
        "refseq": "NP_004324.2",
        "ensembl": "ENSP00000493543.1",
        "pos": (599, 600),
        "status": TranscriptPriority.MANE_SELECT,
        "strand": Strand.NEGATIVE,
        "gene": "BRAF",
    }
    return DataRepresentation(**params)


@pytest.fixture(scope="module")
def egfr_l858r_mane_p():
    """Create test fixture for EGFR L858R MANE Transcript on p coordinate (CA126713)."""
    params = {
        "refseq": "NP_005219.2",
        "ensembl": "ENSP00000275493.2",
        "pos": (857, 858),
        "status": TranscriptPriority.MANE_SELECT,
        "strand": Strand.POSITIVE,
        "gene": "EGFR",
    }
    return DataRepresentation(**params)


@pytest.fixture(scope="module")
def braf_v600e_mane_c():
    """Create test fixture for BRAF V600E MANE Transcript on c coordinate (CA123643)."""
    params = {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_004333.6",
        "ensembl": "ENST00000646891.2",
        "pos": (1798, 1799),
        "status": TranscriptPriority.MANE_SELECT,
        "strand": Strand.NEGATIVE,
        "coding_start_site": 226,
        "coding_end_site": 2527,
        "gene": "BRAF",
    }
    return CdnaRepresentation(**params)


@pytest.fixture(scope="module")
def egfr_l858r_mane_c():
    """Create test fixture for EGFR L858R MANE Transcript on c coordinate (CA126713)."""
    params = {
        "alt_ac": "NC_000007.14",
        "refseq": "NM_005228.5",
        "ensembl": "ENST00000275493.7",
        "pos": (2572, 2573),
        "status": TranscriptPriority.MANE_SELECT,
        "strand": Strand.POSITIVE,
        "coding_start_site": 261,
        "coding_end_site": 3894,
        "gene": "EGFR",
    }
    return CdnaRepresentation(**params)


@pytest.fixture(scope="module")
def grch38_egfr(egfr_mane_gene):
    """Create a test fixture for grch38 responses (CA126713)."""
    params = {
        "pos": (55191821, 55191822),
        "status": TranscriptPriority.GRCH38.value,
        "ac": "NC_000007.14",
        "mane_genes": [egfr_mane_gene],
    }
    return GenomicRepresentation(**params)


@pytest.fixture(scope="module")
def grch38_braf(braf_mane_genes):
    """Create a test fixture for grch38 responses BRAF V600E (genomic)."""
    params = {
        "pos": (140753335, 140753336),
        "status": TranscriptPriority.GRCH38.value,
        "ac": "NC_000007.14",
        "mane_genes": braf_mane_genes,
    }
    return GenomicRepresentation(**params)


@pytest.fixture(scope="module")
def mybpc3_s236g():
    """Create test fixture for MYBPC3 Ser236Gly

    CA1139661942
    https://www.ncbi.nlm.nih.gov/clinvar/variation/922707/?new_evidence=true

    The reason why the protein position is not (235, 236) is because the cdna position
    spans two codons
    """
    params = {
        "protein": {
            "refseq": "NP_000247.2",
            "ensembl": "ENSP00000442795.1",
            "pos": (234, 236),
            "status": TranscriptPriority.MANE_SELECT,
            "strand": Strand.NEGATIVE,
            "gene": "MYBPC3",
        },
        "cdna": {
            "alt_ac": "NC_000011.10",
            "refseq": "NM_000256.3",
            "ensembl": "ENST00000545968.6",
            "pos": (704, 706),
            "status": TranscriptPriority.MANE_SELECT,
            "strand": Strand.NEGATIVE,
            "coding_start_site": 55,
            "coding_end_site": 3880,
            "gene": "MYBPC3",
        },
    }
    return ProteinAndCdnaRepresentation(**params)


def test_get_reading_frame(test_mane_transcript):
    """Test that get_reading_frame works correctly."""
    rf = test_mane_transcript.get_reading_frame(1797)
    assert rf == 3

    rf = test_mane_transcript.get_reading_frame(1798)
    assert rf == 1

    rf = test_mane_transcript.get_reading_frame(1799)
    assert rf == 2

    rf = test_mane_transcript.get_reading_frame(1800)
    assert rf == 3

    rf = test_mane_transcript.get_reading_frame(2573)
    assert rf == 2


def test_p_to_c_pos(test_mane_transcript):
    """Test that _p_to_c_pos method works correctly."""
    # https://civicdb.org/variants/12/summary
    c_pos = test_mane_transcript._p_to_c_pos(599, 600)
    assert c_pos == (1797, 1799)

    # https://civicdb.org/variants/33/summary
    c_pos = test_mane_transcript._p_to_c_pos(857, 858)
    assert c_pos == (2571, 2573)

    # https://civicdb.org//variants/72/summary
    c_pos = test_mane_transcript._p_to_c_pos(575, 576)
    assert c_pos == (1725, 1727)

    # https://civicdb.org/variants/99/summary
    c_pos = test_mane_transcript._p_to_c_pos(841, 842)
    assert c_pos == (2523, 2525)

    # https://civicdb.org/variants/78/summary
    c_pos = test_mane_transcript._p_to_c_pos(11, 12)
    assert c_pos == (33, 35)

    c_pos = test_mane_transcript._p_to_c_pos(747, 751)
    assert c_pos == (2241, 2252)


@pytest.mark.asyncio()
async def test_p_to_c(test_mane_transcript):
    """Test that _p_to_c method works correctly."""
    # Amino Acid Substitution
    expected_pos = 1797, 1799
    ac, pos = await test_mane_transcript._p_to_c("NP_004324.2", 599, 600)
    assert ac == "NM_004333.6"
    assert pos == expected_pos

    ac, pos = await test_mane_transcript._p_to_c("ENSP00000288602.7", 599, 600)
    assert ac == "ENST00000288602.11"
    assert pos == expected_pos

    expected_pos = 2571, 2573
    ac, pos = await test_mane_transcript._p_to_c("NP_005219.2", 857, 858)
    assert ac == "NM_005228.5"
    assert pos == expected_pos

    ac, pos = await test_mane_transcript._p_to_c("ENSP00000275493.2", 857, 858)
    assert ac == "ENST00000275493.7"
    assert pos == expected_pos

    # Polypeptide Truncation
    expected_pos = 552, 554
    ac, pos = await test_mane_transcript._p_to_c("NP_000542.1", 184, 185)
    assert ac == "NM_000551.4"
    assert pos == expected_pos

    ac, pos = await test_mane_transcript._p_to_c("ENSP00000256474.3", 184, 185)
    assert ac == "ENST00000256474.3"
    assert pos == expected_pos

    # Silent Mutation
    expected_pos = 180, 182
    ac, pos = await test_mane_transcript._p_to_c("NP_000542.1", 60, 61)
    assert ac == "NM_000551.4"
    assert pos == expected_pos


@pytest.mark.asyncio()
async def test_c_to_g(test_mane_transcript, nm_004333v6_g):
    """Test that _c_to_g method works correctly."""
    tx_ac = "NM_004333.6"
    g = await test_mane_transcript._c_to_g(tx_ac, (1798, 1800))
    assert g == nm_004333v6_g


@pytest.mark.asyncio()
async def test_g_to_c(
    test_mane_transcript, braf_mane_data, nm_004333v6_g, braf_v600e_mane_c
):
    """Test that _g_to_c method works correctly."""
    mane_c = await test_mane_transcript._g_to_c(
        g=nm_004333v6_g,
        refseq_c_ac=braf_mane_data["RefSeq_nuc"],
        status="_".join(braf_mane_data["MANE_status"].split()).lower(),
        ensembl_c_ac=braf_mane_data["Ensembl_nuc"],
    )
    expected = braf_v600e_mane_c.model_copy(deep=True)
    expected.pos = (1798, 1800)
    expected.alt_ac = None
    assert mane_c == expected


def test_set_liftover(test_mane_transcript, genomic_tx_data):
    """Test that _set_liftover works correctly."""
    cpy = genomic_tx_data.copy(deep=True)
    expected = genomic_tx_data.copy(deep=True)
    test_mane_transcript._set_liftover(cpy, "alt_pos_range", "chr7", "GRCh38")
    expected.alt_pos_range = (140739811, 140739946)
    assert cpy == expected
    test_mane_transcript._set_liftover(cpy, "alt_pos_change_range", "chr7", "GRCh38")
    expected.alt_pos_change_range = (140739903, 140739903)
    assert cpy == expected


@pytest.mark.asyncio()
async def test_liftover_to_38(test_mane_transcript, genomic_tx_data):
    """Test that liftover_to_38 works correctly."""
    cpy = genomic_tx_data.copy(deep=True)
    expected = genomic_tx_data.copy(deep=True)
    await test_mane_transcript._liftover_to_38(cpy)
    expected.alt_ac = "NC_000007.14"
    expected.alt_pos_change_range = (140739903, 140739903)
    expected.alt_pos_range = (140739811, 140739946)
    assert cpy == expected


def test_get_mane_p(test_mane_transcript, braf_mane_data, braf_v600e_mane_p):
    """Test that _get_mane_p method works correctly."""
    mane_p = test_mane_transcript._get_mane_p(braf_mane_data, (1798, 1799))
    assert mane_p == braf_v600e_mane_p


@pytest.mark.asyncio()
async def test_p_to_mane_p(test_mane_transcript, braf_v600e_mane_p, egfr_l858r_mane_p):
    """Test that p_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2",
        599,
        600,
        AnnotationLayer.PROTEIN,
        coordinate_type=CoordinateType.INTER_RESIDUE,
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_004324.2", 600, 600, AnnotationLayer.PROTEIN
    )
    assert mane_p == braf_v600e_mane_p

    # BRAF V600E Ensembl Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7",
        599,
        600,
        AnnotationLayer.PROTEIN,
        coordinate_type=CoordinateType.INTER_RESIDUE,
    )
    assert mane_p == braf_v600e_mane_p

    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000288602.7", 600, 600, AnnotationLayer.PROTEIN
    )
    assert mane_p == braf_v600e_mane_p

    # EGFR L858R RefSeq Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "NP_005219.2", 858, 858, AnnotationLayer.PROTEIN
    )
    assert mane_p == egfr_l858r_mane_p

    # EGFR L858R Ensembl Accessions
    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000275493.2", 858, 858, AnnotationLayer.PROTEIN
    )
    assert mane_p == egfr_l858r_mane_p

    # CA891842275
    resp = await test_mane_transcript.get_mane_transcript(
        "NP_004439.2", 755, 759, AnnotationLayer.PROTEIN
    )
    assert resp == DataRepresentation(
        gene="ERBB2",
        refseq="NP_004439.2",
        ensembl="ENSP00000269571.4",
        pos=(754, 759),
        strand=Strand.POSITIVE,
        status=TranscriptPriority.MANE_SELECT,
    )

    # CA388294838
    mane_p = await test_mane_transcript.get_mane_transcript(
        "ENSP00000366997.4",
        63,
        63,
        AnnotationLayer.PROTEIN,
        gene="DIS3",
        ref="P",
        try_longest_compatible=True,
    )
    assert mane_p == DataRepresentation(
        gene="DIS3",
        refseq="NP_055768.3",
        ensembl="ENSP00000366997.4",
        pos=(62, 63),
        strand=Strand.NEGATIVE,
        status=TranscriptPriority.MANE_SELECT,
    )


@pytest.mark.asyncio()
async def test_c_to_mane_c(test_mane_transcript, braf_v600e_mane_c, egfr_l858r_mane_c):
    """Test that c_to_mane_p method works correctly."""
    # BRAF V600E RefSeq Accessions
    cpy_braf_v600e_mane_c = braf_v600e_mane_c.model_copy(deep=True)
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.4", 1799, 1799, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.4",
        1798,
        1799,
        AnnotationLayer.CDNA,
        coordinate_type=CoordinateType.INTER_RESIDUE,
    )
    assert mane_c == cpy_braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004333.5", 1799, 1799, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_braf_v600e_mane_c

    # BRAF V600E Ensembl Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000288602.10", 1799, 1799, AnnotationLayer.CDNA
    )
    cpy_braf_v600e_mane_c.alt_ac = "NC_000007.13"
    assert mane_c == cpy_braf_v600e_mane_c

    cpy_egfr_l858r_mane_c = egfr_l858r_mane_c.model_copy(deep=True)
    # EGFR L858R RefSeq Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_005228.3", 2573, 2573, AnnotationLayer.CDNA
    )
    assert mane_c == cpy_egfr_l858r_mane_c

    # EGFR L858R Ensembl Accessions
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.7", 2573, 2573, AnnotationLayer.CDNA
    )
    cpy_egfr_l858r_mane_c.alt_ac = "NC_000007.13"
    assert mane_c == cpy_egfr_l858r_mane_c

    # CA645372623
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NM_004448.3", 2264, 2278, AnnotationLayer.CDNA, ref="TGAGGGAAAACACAT"
    )
    assert mane_c == CdnaRepresentation(
        alt_ac="NC_000017.11",
        refseq="NM_004448.4",
        ensembl="ENST00000269571.10",
        pos=(2263, 2278),
        status=TranscriptPriority.MANE_SELECT,
        strand=Strand.POSITIVE,
        coding_start_site=175,
        coding_end_site=3943,
        gene="ERBB2",
    )


def test_get_prioritized_transcripts_from_gene(
    test_seqrepo_access, test_mane_transcript
):
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

    with patch.object(
        test_seqrepo_access, "get_reference_sequence", get_reference_sequence
    ):
        # cds_start_i does not impact this method
        data = [
            ["NP_004324.2", "NM_004333.6", "NC_000007.13", 0],
            ["NP_004324.2", "NM_004333.5", "NC_000007.13", 0],
            ["NP_001365401.1", "NM_001378472.1", "NC_000007.13", 0],
            ["NP_001365401.1", "NM_001374258.2", "NC_000007.13", 0],
            ["NP_004324.2", "NM_004333.6", "NC_000007.14", 0],
            ["NP_004324.2", "NM_004333.5", "NC_000007.14", 0],
            ["NP_001365401.1", "NM_001378472.1", "NC_000007.14", 0],
            ["NP_001365401.1", "NM_001374258.2", "NC_000007.14", 0],
        ]
        test_df = pl.DataFrame(
            data, schema=["pro_ac", "tx_ac", "alt_ac", "cds_start_i"], orient="row"
        )

        resp = test_mane_transcript._get_prioritized_transcripts_from_gene(test_df)
        assert resp == ["NM_004333.6", "NM_001374258.2", "NM_001378472.1"]


@pytest.mark.asyncio()
async def test_get_longest_compatible_transcript(test_mane_transcript):
    """Test that get_longest_compatible_transcript method works as expected"""
    mane_transcripts = {
        "ENST00000646891.2",
        "NM_001374258.1",
        "NM_004333.6",
        "ENST00000644969.2",
    }

    p_expected = DataRepresentation(
        refseq="NP_001365396.1",
        ensembl=None,
        pos=(599, 600),
        strand=Strand.NEGATIVE,
        gene="BRAF",
        status=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        599,
        600,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.PROTEIN,
        coordinate_type=CoordinateType.INTER_RESIDUE,
        mane_transcripts=mane_transcripts,
    )
    assert resp == p_expected

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        600,
        600,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.PROTEIN,
        coordinate_type=CoordinateType.RESIDUE,
        mane_transcripts=mane_transcripts,
    )
    assert resp == p_expected

    c_expected = CdnaRepresentation(
        refseq="NM_001378467.1",
        ensembl=None,
        pos=(1798, 1799),
        strand=Strand.NEGATIVE,
        gene="BRAF",
        coding_start_site=226,
        coding_end_site=2539,
        status=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        1799,
        1799,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.CDNA,
        mane_transcripts=mane_transcripts,
    )
    assert resp == c_expected

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        1798,
        1799,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.CDNA,
        coordinate_type=CoordinateType.INTER_RESIDUE,
        mane_transcripts=mane_transcripts,
    )
    assert resp == c_expected

    # Both protein and cdna end annotation layer
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        1798,
        1799,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.CDNA,
        coordinate_type=CoordinateType.INTER_RESIDUE,
        mane_transcripts=mane_transcripts,
        end_annotation_layer=EndAnnotationLayer.PROTEIN_AND_CDNA,
    )
    assert resp == ProteinAndCdnaRepresentation(protein=p_expected, cdna=c_expected)

    resp = await test_mane_transcript.get_longest_compatible_transcript(
        140453136,
        140453136,
        gene="BRAF",
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts,
        alt_ac="NC_000007.13",
    )
    assert resp == CdnaRepresentation(
        refseq="NM_001378467.1",
        ensembl=None,
        pos=(1807, 1808),
        strand=Strand.NEGATIVE,
        gene="BRAF",
        coding_start_site=226,
        coding_end_site=2539,
        status=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )

    # CA1139661942 has no other RefSeq accessions
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
    # CDNA (CA645561524, NM_001346899.2:c.2104_2121delinsCAA)
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        55174776,
        55174793,
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts,
        alt_ac="NC_000007.14",
    )
    assert resp == CdnaRepresentation(
        refseq="NM_001346899.2",
        ensembl=None,
        pos=(2103, 2121),
        strand=Strand.POSITIVE,
        gene="EGFR",
        coding_start_site=261,
        coding_end_site=3759,
        status=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )

    # protein (CA645561524, NP_001333828.1:p.Leu702_Ser707delinsGln)
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        55174776,
        55174793,
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts=mane_transcripts,
        alt_ac="NC_000007.14",
        end_annotation_layer=EndAnnotationLayer.PROTEIN,
    )
    assert resp == DataRepresentation(
        refseq="NP_001333828.1",
        ensembl=None,
        pos=(701, 707),
        strand=Strand.POSITIVE,
        gene=None,
        status=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )

    # CA2499226460 (NP_000416.1:p.Pro240_Pro259delinsThrLeuTh...)
    resp = await test_mane_transcript.get_longest_compatible_transcript(
        153870419,
        153870476,
        start_annotation_layer=AnnotationLayer.GENOMIC,
        mane_transcripts={"ENST00000370060.7", "NM_001278116.2"},
        alt_ac="NC_000023.11",
        end_annotation_layer=EndAnnotationLayer.PROTEIN,
    )
    assert resp == DataRepresentation(
        refseq="NP_000416.1",
        ensembl=None,
        pos=(239, 259),
        strand=Strand.NEGATIVE,
        gene=None,
        status=TranscriptPriority.LONGEST_COMPATIBLE_REMAINING,
    )


@pytest.mark.asyncio()
async def test_g_to_grch38(test_mane_transcript, grch38_egfr, grch38_braf):
    """Test that g_to_grch38 method works correctly."""
    resp = await test_mane_transcript.g_to_grch38(
        "NC_000007.13", 55259515, 55259515, get_mane_genes=False
    )
    grch38_egfr_no_genes = grch38_egfr.copy()
    grch38_egfr_no_genes.mane_genes = []
    assert resp == grch38_egfr_no_genes

    resp = await test_mane_transcript.g_to_grch38(
        "NC_000007.13", 55259515, 55259515, get_mane_genes=True
    )
    assert resp == grch38_egfr

    resp = await test_mane_transcript.g_to_grch38(
        "NC_000007.13", 140453136, 140453136, get_mane_genes=True
    )
    assert resp == grch38_braf

    resp = await test_mane_transcript.g_to_grch38(
        "NC_000007.13",
        140453135,
        140453136,
        coordinate_type=CoordinateType.INTER_RESIDUE,
        get_mane_genes=True,
    )
    assert resp == grch38_braf

    # Already on GRCh38
    resp = await test_mane_transcript.g_to_grch38(
        "NC_000007.14", 140753336, 140753336, get_mane_genes=True
    )
    assert resp == grch38_braf


@pytest.mark.asyncio()
async def test_g_to_mane_c(test_mane_transcript, egfr_l858r_mane_c, braf_v600e_mane_c):
    """Test that g_to_mane_c method works correctly."""
    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259515, 55259515, gene="EGFR"
    )
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13",
        55259514,
        55259515,
        gene="EGFR",
        coordinate_type=CoordinateType.INTER_RESIDUE,
    )
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 55259515, 55259515, gene="EGFR"
    )
    assert mane_c == egfr_l858r_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13",
        140453135,
        140453136,
        gene="BRAF",
        coordinate_type=CoordinateType.INTER_RESIDUE,
    )
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 140453136, 140453136, AnnotationLayer.GENOMIC, gene="BRAF"
    )
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13",
        140453135,
        140453136,
        AnnotationLayer.GENOMIC,
        gene="BRAF",
        coordinate_type=CoordinateType.INTER_RESIDUE,
    )
    assert mane_c == braf_v600e_mane_c

    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000007.13", 140453136, 140453136, gene="BRAF"
    )
    assert mane_c == braf_v600e_mane_c

    # CA122528
    mane_c = await test_mane_transcript.g_to_mane_c(
        "NC_000012.11", 25398284, 25398284, gene="KRAS"
    )
    assert mane_c == CdnaRepresentation(
        alt_ac="NC_000012.12",
        refseq="NM_004985.5",
        ensembl="ENST00000311936.8",
        pos=(34, 35),
        status=TranscriptPriority.MANE_SELECT,
        strand=Strand.NEGATIVE,
        coding_start_site=190,
        coding_end_site=757,
        gene="KRAS",
    )

    # CA090928
    mane_c = await test_mane_transcript.get_mane_transcript(
        "NC_000007.13", 55249071, 55249071, AnnotationLayer.GENOMIC, gene="EGFR"
    )
    assert mane_c == CdnaRepresentation(
        alt_ac="NC_000007.14",
        refseq="NM_005228.5",
        ensembl="ENST00000275493.7",
        pos=(2368, 2369),
        status=TranscriptPriority.MANE_SELECT,
        strand=Strand.POSITIVE,
        coding_start_site=261,
        coding_end_site=3894,
        gene="EGFR",
    )


@pytest.mark.asyncio()
async def test_grch38_to_mane_c_p(
    test_mane_transcript,
    braf_v600e_mane_p,
    egfr_l858r_mane_p,
    mybpc3_s236g,
    braf_v600e_mane_c,
    egfr_l858r_mane_c,
):
    """Test that grch38_to_mane_c_p works correctly"""
    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000007.14", 140753336, 140753336, gene="BRAF"
    )
    assert resp == ProteinAndCdnaRepresentation(
        protein=braf_v600e_mane_p, cdna=braf_v600e_mane_c
    )

    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000007.14", 55191822, 55191822
    )
    assert resp == ProteinAndCdnaRepresentation(
        protein=egfr_l858r_mane_p, cdna=egfr_l858r_mane_c
    )

    # https://www.ncbi.nlm.nih.gov/clinvar/variation/922707/?new_evidence=true
    # Without gene
    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000011.10", 47348490, 47348491
    )
    assert resp == mybpc3_s236g

    # With gene
    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000011.10", 47348490, 47348491, gene="MYBPC3"
    )
    assert resp == mybpc3_s236g

    # CA645561524 (without gene)
    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000007.14", 55174776, 55174793
    )
    assert resp == ProteinAndCdnaRepresentation(
        protein=DataRepresentation(
            refseq="NP_005219.2",
            ensembl="ENSP00000275493.2",
            pos=(746, 752),
            status=TranscriptPriority.MANE_SELECT,
            strand=Strand.POSITIVE,
            gene="EGFR",
        ),
        cdna=CdnaRepresentation(
            alt_ac="NC_000007.14",
            refseq="NM_005228.5",
            ensembl="ENST00000275493.7",
            pos=(2238, 2256),
            status=TranscriptPriority.MANE_SELECT,
            strand=Strand.POSITIVE,
            coding_start_site=261,
            coding_end_site=3894,
            gene="EGFR",
        ),
    )

    # CA2499226460 (without gene)
    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000023.11", 153870419, 153870476
    )
    assert resp == ProteinAndCdnaRepresentation(
        protein=DataRepresentation(
            refseq="NP_001265045.1",
            ensembl="ENSP00000359077.1",
            pos=(239, 259),
            status=TranscriptPriority.MANE_SELECT,
            strand=Strand.NEGATIVE,
            gene="L1CAM",
        ),
        cdna=CdnaRepresentation(
            alt_ac="NC_000023.11",
            refseq="NM_001278116.2",
            ensembl="ENST00000370060.7",
            pos=(717, 775),
            status=TranscriptPriority.MANE_SELECT,
            strand=Strand.NEGATIVE,
            coding_start_site=217,
            coding_end_site=3991,
            gene="L1CAM",
        ),
    )

    # CA16042245
    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000002.12",
        74530927,
        74530929,
    )
    assert resp == ProteinAndCdnaRepresentation(
        protein=DataRepresentation(
            refseq="NP_037379.1",
            ensembl="ENSP00000258080.3",
            pos=(242, 244),
            status=TranscriptPriority.MANE_SELECT,
            strand=Strand.POSITIVE,
            gene="HTRA2",
        ),
        cdna=CdnaRepresentation(
            alt_ac="NC_000002.12",
            refseq="NM_013247.5",
            ensembl="ENST00000258080.8",
            pos=(727, 730),
            status=TranscriptPriority.MANE_SELECT,
            strand=Strand.POSITIVE,
            coding_start_site=74,
            coding_end_site=1451,
            gene="HTRA2",
        ),
    )

    # GENE not valid
    resp = await test_mane_transcript.grch38_to_mane_c_p(
        "NC_000023.11", 153870419, 153870476, gene="FAKE"
    )
    assert resp is None


@pytest.mark.asyncio()
async def test_valid(test_mane_transcript):
    """Test that valid queries do not raise any exceptions"""
    resp = await test_mane_transcript.get_mane_transcript(
        "NP_001296812.1",
        201,
        201,
        AnnotationLayer.PROTEIN,
        ref="R",
        try_longest_compatible=True,
    )
    assert resp

    # issue-394
    resp = await test_mane_transcript.get_mane_transcript(
        "ENST00000496384.7",
        1799,
        1799,
        AnnotationLayer.CDNA,
        try_longest_compatible=True,
    )
    assert resp


@pytest.mark.asyncio()
async def test_no_matches(test_mane_transcript):
    """Test that invalid queries return None."""
    # Invalid ENST version
    mane_c = await test_mane_transcript.get_mane_transcript(
        "ENST00000275493.15645", 2573, 2573, AnnotationLayer.CDNA
    )
    assert mane_c is None
