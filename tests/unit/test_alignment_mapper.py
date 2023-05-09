"""Module for testing the Alignment Mapper class"""
import pytest

from cool_seq_tool.alignment_mapper import AlignmentMapper
from cool_seq_tool.data_sources import TranscriptMappings, UTADatabase
from cool_seq_tool.schemas import Assembly, ResidueMode


@pytest.fixture(scope="module")
def test_alignment_mapper(test_seqrepo_access):
    """Build AlignmentMapper test fixture"""
    return AlignmentMapper(test_seqrepo_access, TranscriptMappings(), UTADatabase())


@pytest.fixture(scope="module")
def braf_v600e_c():
    """Create test fixture for BRAF V600E cDNA representation

    Example of residue (top) vs inter-residue (bottom)
    | 1796 |      | 1798 |      | 1799 |      | 1800 |      |
    |      | 1797 |      | 1798 |      | 1799 |      | 1800 |
    """
    return {
        "c_ac": "NM_004333.6",
        "c_start_pos": 1797,
        "c_end_pos": 1800,
        "cds_start": 226,
        "residue_mode": "inter-residue"
    }


@pytest.fixture(scope="module")
def egfr_l858r_c():
    """Create test fixture for EGFR L858R cDNA representation"""
    return {
        "c_ac": "NM_005228.5",
        "c_start_pos": 2571,
        "c_end_pos": 2574,
        "cds_start": 261,
        "residue_mode": "inter-residue"
    }


@pytest.fixture(scope="module")
def braf_v600e_grch37():
    """Create test fixture for BRAF V600E GRCh37 representation"""
    return {
        "g_ac": "NC_000007.13",
        "g_start_pos": 140453134,
        "g_end_pos": 140453137,
        "residue_mode": "inter-residue",
        "gene": "BRAF",
        "strand": "-"
    }


@pytest.fixture(scope="module")
def braf_v600e_grch38():
    """Create test fixture for BRAF V600E GRCh38 representation"""
    return {
        "g_ac": "NC_000007.14",
        "g_start_pos": 140753334,
        "g_end_pos": 140753337,
        "residue_mode": "inter-residue",
        "gene": "BRAF",
        "strand": "-"
    }


@pytest.fixture(scope="module")
def egfr_l858r_grch37():
    """Create test fixture for EGFR L858R GRCh37 representation"""
    return {
        "g_ac": "NC_000007.13",
        "g_start_pos": 55259513,
        "g_end_pos": 55259516,
        "residue_mode": "inter-residue",
        "gene": "EGFR",
        "strand": "+"
    }


@pytest.fixture(scope="module")
def egfr_l858r_grch38():
    """Create test fixture for EGFR L858R GRCh38 representation"""
    return {
        "g_ac": "NC_000007.14",
        "g_start_pos": 55191820,
        "g_end_pos": 55191823,
        "residue_mode": "inter-residue",
        "gene": "EGFR",
        "strand": "+"
    }


@pytest.fixture(scope="module")
def delins_grch37():
    """Create test fixture for CA645544092 delins GRCh37 representation"""
    return {
        "g_ac": "NC_000007.13",
        "g_start_pos": 140453131,
        "g_end_pos": 140453137,
        "residue_mode": "inter-residue",
        "gene": "BRAF",
        "strand": "-"
    }


@pytest.fixture(scope="module")
def hras_t2a():
    """Create test fixture for CA10582926 representation"""

    def _expected(assembly):
        g_ac = "NC_000011.9" if assembly == Assembly.GRCH37 else "NC_000011.10"
        return {
            "g_ac": g_ac,
            "g_start_pos": 534316,
            "g_end_pos": 534319,
            "residue_mode": "inter-residue",
            "gene": "HRAS",
            "strand": "-"
        }

    return _expected


@pytest.mark.asyncio
async def test_get_cds_start(test_alignment_mapper):
    """Test that _get_cds_start method works correctly"""
    # Valid
    resp, w = await test_alignment_mapper._get_cds_start("NM_004333.6")
    assert resp == 226
    assert w is None

    # Invalid (fake accession)
    resp, w = await test_alignment_mapper._get_cds_start("NM_004333.6293702")
    assert resp is None
    assert w == "Accession NM_004333.6293702 not found in UTA db"


@pytest.mark.asyncio
async def test_p_to_c(test_alignment_mapper, braf_v600e_c, egfr_l858r_c):
    """Test that p_to_c works as expected"""
    # BRAF V600E
    for params in [
        ("NP_004324.2", 600, 600, ResidueMode.RESIDUE),
        ("NP_004324.2", 599, 600, ResidueMode.INTER_RESIDUE),
        ("NP_004324.2", 599, 599, ResidueMode.INTER_RESIDUE)

    ]:
        ac, start, end, residue_mode = params
        resp, w = await test_alignment_mapper.p_to_c(ac, start, end, residue_mode)
        assert w is None, params
        assert resp == braf_v600e_c, params

    # EGFR L858R
    for params in [
        ("NP_005219.2", 858, 858, ResidueMode.RESIDUE),
        ("NP_005219.2", 857, 858, ResidueMode.INTER_RESIDUE),
        ("NP_005219.2", 857, 857, ResidueMode.INTER_RESIDUE)
    ]:
        ac, start, end, residue_mode = params
        resp, w = await test_alignment_mapper.p_to_c(ac, start, end, residue_mode)
        assert w is None, params
        assert resp == egfr_l858r_c, params

    # CA16602374
    resp, w = await test_alignment_mapper.p_to_c(
        "NP_005887.2", 132, 132, ResidueMode.RESIDUE)
    assert w is None
    assert resp == {
        "c_ac": "NM_005896.4",
        "c_start_pos": 393,
        "c_end_pos": 396,
        "cds_start": 223,
        "residue_mode": "inter-residue"
    }

    resp, w = await test_alignment_mapper.p_to_c(
        "NP_000542.1", 185, 185, ResidueMode.RESIDUE)
    assert w is None
    assert resp == {
        "c_ac": "NM_000551.4",
        "c_start_pos": 552,
        "c_end_pos": 555,
        "cds_start": 70,
        "residue_mode": "inter-residue"
    }

    resp, w = await test_alignment_mapper.p_to_c(
        "ENSP00000256474.3", 185, 185, ResidueMode.RESIDUE)
    assert w is None
    assert resp == {
        "c_ac": "ENST00000256474.3",
        "c_start_pos": 552,
        "c_end_pos": 555,
        "cds_start": 840,
        "residue_mode": "inter-residue"
    }


@pytest.fixture(scope="module")
async def test_p_to_c_invalid(test_alignment_mapper):
    """Test invalid queries for p_to_c method"""
    # Invalid protein accession
    resp, w = await test_alignment_mapper.p_to_c(
        "NP_005219", 857, 857, ResidueMode.INTER_RESIDUE)
    assert w == "NP_005219 not found in transcript mappings"
    assert resp is None


@pytest.mark.asyncio
async def test_c_to_g(test_alignment_mapper, braf_v600e_grch37, braf_v600e_grch38,
                      egfr_l858r_grch37, egfr_l858r_grch38):
    """Test that c_to_g works as expected"""
    # BRAF V600E
    for params in [
        ("NM_004333.6", 1798, 1800, ResidueMode.RESIDUE, Assembly.GRCH37),
        ("NM_004333.6", 1797, 1800, ResidueMode.INTER_RESIDUE, Assembly.GRCH37),
        ("NM_004333.6", 1798, 1800, ResidueMode.RESIDUE, Assembly.GRCH38),
        ("NM_004333.6", 1797, 1800, ResidueMode.INTER_RESIDUE, Assembly.GRCH38)
    ]:
        ac, start, end, residue_mode, assembly = params
        expected = braf_v600e_grch37 if assembly == Assembly.GRCH37 else braf_v600e_grch38  # noqa: E501
        resp, w = await test_alignment_mapper.c_to_g(
            ac, start, end, residue_mode=residue_mode, target_genome_assembly=assembly)
        assert w is None, params
        assert resp == expected, params

    # EGFR L858R
    for params in [
        ("NM_005228.5", 2572, 2574, ResidueMode.RESIDUE, Assembly.GRCH37),
        ("NM_005228.5", 2571, 2574, ResidueMode.INTER_RESIDUE, Assembly.GRCH37),
        ("NM_005228.5", 2572, 2574, ResidueMode.RESIDUE, Assembly.GRCH38),
        ("NM_005228.5", 2571, 2574, ResidueMode.INTER_RESIDUE, Assembly.GRCH38)
    ]:
        ac, start, end, residue_mode, assembly = params
        expected = egfr_l858r_grch37 if assembly == Assembly.GRCH37 else egfr_l858r_grch38  # noqa: E501
        resp, w = await test_alignment_mapper.c_to_g(
            ac, start, end, residue_mode=residue_mode, target_genome_assembly=assembly)
        assert w is None, params
        assert resp == expected, params


@pytest.mark.asyncio
async def test_c_to_g_invalid(test_alignment_mapper):
    """Test invalid queries for c_to_g method"""
    # Should not expect to find anything given these two positions
    resp, w = await test_alignment_mapper.c_to_g(
        "NM_005228.5", 1, 999999, residue_mode=ResidueMode.RESIDUE)
    assert resp is None
    assert w == "Unable to find genomic and transcript data for NM_005228.5 at position (0, 999999)"  # noqa: E501


@pytest.mark.asyncio
async def test_p_to_g(
    test_alignment_mapper, braf_v600e_grch37, braf_v600e_grch38, egfr_l858r_grch37,
    egfr_l858r_grch38, delins_grch37, hras_t2a
):
    """Test that p_to_g works as expected"""
    # BRAF V600E
    for params in [
        ("NP_004324.2", 600, 600, ResidueMode.RESIDUE, Assembly.GRCH37),
        ("NP_004324.2", 599, 600, ResidueMode.INTER_RESIDUE, Assembly.GRCH37),
        ("NP_004324.2", 600, 600, ResidueMode.RESIDUE, Assembly.GRCH38),
        ("NP_004324.2", 599, 600, ResidueMode.INTER_RESIDUE, Assembly.GRCH38)
    ]:
        ac, start, end, residue_mode, assembly = params
        expected = braf_v600e_grch37 if assembly == Assembly.GRCH37 else braf_v600e_grch38  # noqa: E501
        resp, w = await test_alignment_mapper.p_to_g(
            ac, start, end, residue_mode=residue_mode, target_genome_assembly=assembly)
        assert w is None, params
        assert resp == expected, params

    # EGFR L858R
    for params in [
        ("NP_005219.2", 858, 858, ResidueMode.RESIDUE, Assembly.GRCH37),
        ("NP_005219.2", 857, 858, ResidueMode.INTER_RESIDUE, Assembly.GRCH37),
        ("NP_005219.2", 858, 858, ResidueMode.RESIDUE, Assembly.GRCH38),
        ("NP_005219.2", 857, 858, ResidueMode.INTER_RESIDUE, Assembly.GRCH38)
    ]:
        ac, start, end, residue_mode, assembly = params
        expected = egfr_l858r_grch37 if assembly == Assembly.GRCH37 else egfr_l858r_grch38  # noqa: E501
        resp, w = await test_alignment_mapper.p_to_g(
            ac, start, end, residue_mode=residue_mode, target_genome_assembly=assembly)
        assert w is None, params
        assert resp == expected, params

    # Delins example: CA645544092
    for params in [
        ("NP_004324.2", 600, 601, ResidueMode.RESIDUE, Assembly.GRCH37),
        ("NP_004324.2", 599, 601, ResidueMode.INTER_RESIDUE, Assembly.GRCH37)
    ]:
        ac, start, end, residue_mode, assembly = params
        resp, w = await test_alignment_mapper.p_to_g(
            ac, start, end, residue_mode=residue_mode, target_genome_assembly=assembly)
        assert w is None, params
        assert resp == delins_grch37, params

    # Example not using mane accession: CA10582926
    for params in [
        ("NP_001123914.1", 2, 2, ResidueMode.RESIDUE, Assembly.GRCH37),
        ("NP_001123914.1", 1, 2, ResidueMode.INTER_RESIDUE, Assembly.GRCH37),
        ("NP_001123914.1", 2, 2, ResidueMode.RESIDUE, Assembly.GRCH38),
        ("NP_001123914.1", 1, 2, ResidueMode.INTER_RESIDUE, Assembly.GRCH38)
    ]:
        ac, start, end, residue_mode, assembly = params
        resp, w = await test_alignment_mapper.p_to_g(
            ac, start, end, residue_mode=residue_mode, target_genome_assembly=assembly)
        assert w is None, params
        assert resp == hras_t2a(assembly), params


@pytest.mark.asyncio
async def test_p_to_g_invalid(test_alignment_mapper):
    """Test invalid queries for p_to_g method"""
    # Invalid protein accession
    resp, w = await test_alignment_mapper.p_to_g(
        "NP_0000", 600, 600, ResidueMode.RESIDUE)
    assert resp is None
    assert w == "NP_0000 not found in transcript mappings"


@pytest.mark.asyncio
async def g_to_c(test_alignment_mapper):
    """Test that g_to_c works correctly"""
    # NEGATIVE STRAND (CA123643)
    for params in [
        ("NC_000007.13", 140453136, 140453136, ResidueMode.RESIDUE), # GRCh37
        ("NC_000007.14", 140753336, 140753336, ResidueMode.RESIDUE)  # GRCh38
    ]:
        alt_ac, start_pos, end_pos, residue_mode = params

        # NM_001374258.1:c.1919T>A (MANE Plus Clinical)
        resp, w = await test_alignment_mapper.g_to_c(
            alt_ac, start_pos, end_pos, "NM_001374258.1", residue_mode)
        assert w is None
        assert resp == {
            "c_ac": "NM_001374258.1",
            "c_start_pos": 1918,
            "c_end_pos": 1919,
            "cds_start": 226,
            "residue_mode": "inter-residue"
        }

        # NM_004333.6:c.1799T>A (MANE Select)
        resp, w = await test_alignment_mapper.g_to_c(
            alt_ac, start_pos, end_pos, "NM_004333.6", residue_mode)
        assert w is None
        assert resp == {
            "c_ac": "NM_004333.6",
            "c_start_pos": 1798,
            "c_end_pos": 1799,
            "cds_start": 226,
            "residue_mode": "inter-residue"
        }

        # NM_001378467.1:c.1808T>A
        resp, w = await test_alignment_mapper.g_to_c(
            alt_ac, start_pos, end_pos, "NM_004333.6", residue_mode)
        assert w is None
        assert resp == {
            "c_ac": "NM_004333.6",
            "c_start_pos": 1798,
            "c_end_pos": 1799,
            "cds_start": 226,
            "residue_mode": "inter-residue"
        }

    # POSITIVE STRAND (CA126713)
    for params in [
        ("NC_000007.13", 55259515, 55259515, ResidueMode.RESIDUE),  # GRCh37
        ("NC_000007.14", 55191822, 55191822, ResidueMode.RESIDUE)  # GRCh38
    ]:
        alt_ac, start_pos, end_pos, residue_mode = params

        # NM_005228.5:c.2573T>G (MANE Select)
        resp, w = await test_alignment_mapper.g_to_c(
            alt_ac, start_pos, end_pos, "NM_005228.5", residue_mode)
        assert w is None
        assert resp == {
            "c_ac": "NM_005228.5",
            "c_start_pos": 2572,
            "c_end_pos": 2573,
            "cds_start": 261,
            "residue_mode": "inter-residue"
        }

        # NM_001346941.2:c.1772T>G
        resp, w = await test_alignment_mapper.g_to_c(
            alt_ac, start_pos, end_pos, "NM_001346941.2", residue_mode)
        assert w is None
        assert resp == {
            "c_ac": "NM_001346941.2",
            "c_start_pos": 1771,
            "c_end_pos": 1772,
            "cds_start": 261,
            "residue_mode": "inter-residue"
        }

        # NM_001346900.2:c.2414T>G
        resp, w = await test_alignment_mapper.g_to_c(
            alt_ac, start_pos, end_pos, "NM_001346900.2", residue_mode)
        assert w is None
        assert resp == {
            "c_ac": "NM_001346900.2",
            "c_start_pos": 2413,
            "c_end_pos": 2414,
            "cds_start": 261,
            "residue_mode": "inter-residue"
        }


@pytest.mark.asyncio
async def test_g_to_grch38(test_alignment_mapper):
    """Test that g_to_grch38 works correctly"""
    # GRCh37
    resp = await test_alignment_mapper.g_to_grch38("NC_000007.13", 140453136, 140453136)
    assert resp == {
        "ac": "NC_000007.14",
        "pos": (140753336, 140753336)
    }

    # GRCh38
    resp = await test_alignment_mapper.g_to_grch38("NC_000007.14", 140753336, 140753336)
    assert resp == {
        "ac": "NC_000007.14",
        "pos": (140753336, 140753336)
    }
