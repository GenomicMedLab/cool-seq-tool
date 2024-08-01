"""Module for testing the Alignment Mapper class"""

import pytest

from cool_seq_tool.schemas import Assembly, CoordinateType


@pytest.fixture(scope="module")
def test_alignment_mapper(test_cool_seq_tool):
    """Build AlignmentMapper test fixture"""
    return test_cool_seq_tool.alignment_mapper


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
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }


@pytest.fixture(scope="module")
def egfr_l858r_c():
    """Create test fixture for EGFR L858R cDNA representation"""
    return {
        "c_ac": "NM_005228.5",
        "c_start_pos": 2571,
        "c_end_pos": 2574,
        "cds_start": 261,
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }


@pytest.fixture(scope="module")
def braf_v600e_grch37():
    """Create test fixture for BRAF V600E GRCh37 representation"""
    return {
        "g_ac": "NC_000007.13",
        "g_start_pos": 140453134,
        "g_end_pos": 140453137,
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }


@pytest.fixture(scope="module")
def braf_v600e_grch38():
    """Create test fixture for BRAF V600E GRCh38 representation"""
    return {
        "g_ac": "NC_000007.14",
        "g_start_pos": 140753334,
        "g_end_pos": 140753337,
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }


@pytest.fixture(scope="module")
def egfr_l858r_grch37():
    """Create test fixture for EGFR L858R GRCh37 representation"""
    return {
        "g_ac": "NC_000007.13",
        "g_start_pos": 55259513,
        "g_end_pos": 55259516,
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }


@pytest.fixture(scope="module")
def egfr_l858r_grch38():
    """Create test fixture for EGFR L858R GRCh38 representation"""
    return {
        "g_ac": "NC_000007.14",
        "g_start_pos": 55191820,
        "g_end_pos": 55191823,
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }


@pytest.fixture(scope="module")
def delins_grch37():
    """Create test fixture for CA645544092 delins GRCh37 representation"""
    return {
        "g_ac": "NC_000007.13",
        "g_start_pos": 140453131,
        "g_end_pos": 140453137,
        "coordinate_type": CoordinateType.INTER_RESIDUE,
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
            "coordinate_type": CoordinateType.INTER_RESIDUE,
        }

    return _expected


@pytest.mark.asyncio()
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


@pytest.mark.asyncio()
async def test_p_to_c(test_alignment_mapper, braf_v600e_c, egfr_l858r_c):
    """Test that p_to_c works as expected"""
    # BRAF V600E
    for params in [
        ("NP_004324.2", 600, 600, CoordinateType.RESIDUE),
        ("NP_004324.2", 599, 600, CoordinateType.INTER_RESIDUE),
        ("NP_004324.2", 599, 599, CoordinateType.INTER_RESIDUE),
    ]:
        ac, start, end, coordinate_type = params
        resp, w = await test_alignment_mapper.p_to_c(ac, start, end, coordinate_type)
        assert w is None, params
        assert resp == braf_v600e_c, params

    # EGFR L858R
    for params in [
        ("NP_005219.2", 858, 858, CoordinateType.RESIDUE),
        ("NP_005219.2", 857, 858, CoordinateType.INTER_RESIDUE),
        ("NP_005219.2", 857, 857, CoordinateType.INTER_RESIDUE),
    ]:
        ac, start, end, coordinate_type = params
        resp, w = await test_alignment_mapper.p_to_c(ac, start, end, coordinate_type)
        assert w is None, params
        assert resp == egfr_l858r_c, params

    # CA16602374
    resp, w = await test_alignment_mapper.p_to_c(
        "NP_005887.2", 132, 132, CoordinateType.RESIDUE
    )
    assert w is None
    assert resp == {
        "c_ac": "NM_005896.4",
        "c_start_pos": 393,
        "c_end_pos": 396,
        "cds_start": 223,
        "coordinate_type": CoordinateType.INTER_RESIDUE,
    }


@pytest.mark.asyncio()
async def test_p_to_c_invalid(test_alignment_mapper):
    """Test invalid queries for p_to_c method"""
    # Invalid protein accession
    resp, w = await test_alignment_mapper.p_to_c(
        "NP_005219", 857, 857, CoordinateType.INTER_RESIDUE
    )
    assert w == "NP_005219 not found in transcript mappings"
    assert resp is None


@pytest.mark.asyncio()
async def test_c_to_g(
    test_alignment_mapper,
    braf_v600e_grch37,
    braf_v600e_grch38,
    egfr_l858r_grch37,
    egfr_l858r_grch38,
):
    """Test that c_to_g works as expected"""
    # BRAF V600E
    for params in [
        ("NM_004333.6", 1798, 1800, CoordinateType.RESIDUE, Assembly.GRCH37),
        ("NM_004333.6", 1797, 1800, CoordinateType.INTER_RESIDUE, Assembly.GRCH37),
        ("NM_004333.6", 1798, 1800, CoordinateType.RESIDUE, Assembly.GRCH38),
        ("NM_004333.6", 1797, 1800, CoordinateType.INTER_RESIDUE, Assembly.GRCH38),
    ]:
        ac, start, end, coordinate_type, assembly = params
        expected = (
            braf_v600e_grch37 if assembly == Assembly.GRCH37 else braf_v600e_grch38
        )
        resp, w = await test_alignment_mapper.c_to_g(
            ac,
            start,
            end,
            coordinate_type=coordinate_type,
            target_genome_assembly=assembly,
        )
        assert w is None, params
        assert resp == expected, params

    # EGFR L858R
    for params in [
        ("NM_005228.5", 2572, 2574, CoordinateType.RESIDUE, Assembly.GRCH37),
        ("NM_005228.5", 2571, 2574, CoordinateType.INTER_RESIDUE, Assembly.GRCH37),
        ("NM_005228.5", 2572, 2574, CoordinateType.RESIDUE, Assembly.GRCH38),
        ("NM_005228.5", 2571, 2574, CoordinateType.INTER_RESIDUE, Assembly.GRCH38),
    ]:
        ac, start, end, coordinate_type, assembly = params
        expected = (
            egfr_l858r_grch37 if assembly == Assembly.GRCH37 else egfr_l858r_grch38
        )
        resp, w = await test_alignment_mapper.c_to_g(
            ac,
            start,
            end,
            coordinate_type=coordinate_type,
            target_genome_assembly=assembly,
        )
        assert w is None, params
        assert resp == expected, params


@pytest.mark.asyncio()
async def test_c_to_g_invalid(test_alignment_mapper):
    """Test invalid queries for c_to_g method"""
    # Should not expect to find anything given these two positions
    resp, w = await test_alignment_mapper.c_to_g(
        "NM_005228.5", 1, 999999, coordinate_type=CoordinateType.RESIDUE
    )
    assert resp is None
    assert (
        w
        == "Unable to find genomic and transcript data for NM_005228.5 at position (0, 999999)"
    )

    # c_start_pos and c_end_pos cannot be the same
    resp, w = await test_alignment_mapper.c_to_g(
        "NM_005228.5", 1, 1, coordinate_type=CoordinateType.RESIDUE
    )
    assert resp is None
    assert w == "c_start_pos and c_end_pos are not a valid range for the codon(s)"

    # c_start_pos and c_end_pos range is not a factor of 3
    resp, w = await test_alignment_mapper.c_to_g(
        "NM_005228.5", 1, 2, coordinate_type=CoordinateType.RESIDUE
    )
    assert resp is None
    assert w == "c_start_pos and c_end_pos are not a valid range for the codon(s)"

    # c_start_pos and c_end_pos range is not a factor of 3
    resp, w = await test_alignment_mapper.c_to_g(
        "NM_005228.5", 1, 3, coordinate_type=CoordinateType.INTER_RESIDUE
    )
    assert resp is None
    assert w == "c_start_pos and c_end_pos are not a valid range for the codon(s)"


@pytest.mark.asyncio()
async def test_p_to_g(
    test_alignment_mapper,
    braf_v600e_grch37,
    braf_v600e_grch38,
    egfr_l858r_grch37,
    egfr_l858r_grch38,
    delins_grch37,
    hras_t2a,
):
    """Test that p_to_g works as expected"""
    # BRAF V600E
    for params in [
        ("NP_004324.2", 600, 600, CoordinateType.RESIDUE, Assembly.GRCH37),
        ("NP_004324.2", 599, 600, CoordinateType.INTER_RESIDUE, Assembly.GRCH37),
        ("NP_004324.2", 600, 600, CoordinateType.RESIDUE, Assembly.GRCH38),
        ("NP_004324.2", 599, 600, CoordinateType.INTER_RESIDUE, Assembly.GRCH38),
    ]:
        ac, start, end, coordinate_type, assembly = params
        expected = (
            braf_v600e_grch37 if assembly == Assembly.GRCH37 else braf_v600e_grch38
        )
        resp, w = await test_alignment_mapper.p_to_g(
            ac,
            start,
            end,
            coordinate_type=coordinate_type,
            target_genome_assembly=assembly,
        )
        assert w is None, params
        assert resp == expected, params

    # EGFR L858R
    for params in [
        ("NP_005219.2", 858, 858, CoordinateType.RESIDUE, Assembly.GRCH37),
        ("NP_005219.2", 857, 858, CoordinateType.INTER_RESIDUE, Assembly.GRCH37),
        ("NP_005219.2", 858, 858, CoordinateType.RESIDUE, Assembly.GRCH38),
        ("NP_005219.2", 857, 858, CoordinateType.INTER_RESIDUE, Assembly.GRCH38),
    ]:
        ac, start, end, coordinate_type, assembly = params
        expected = (
            egfr_l858r_grch37 if assembly == Assembly.GRCH37 else egfr_l858r_grch38
        )
        resp, w = await test_alignment_mapper.p_to_g(
            ac,
            start,
            end,
            coordinate_type=coordinate_type,
            target_genome_assembly=assembly,
        )
        assert w is None, params
        assert resp == expected, params

    # Delins example: CA645544092
    for params in [
        ("NP_004324.2", 600, 601, CoordinateType.RESIDUE, Assembly.GRCH37),
        ("NP_004324.2", 599, 601, CoordinateType.INTER_RESIDUE, Assembly.GRCH37),
    ]:
        ac, start, end, coordinate_type, assembly = params
        resp, w = await test_alignment_mapper.p_to_g(
            ac,
            start,
            end,
            coordinate_type=coordinate_type,
            target_genome_assembly=assembly,
        )
        assert w is None, params
        assert resp == delins_grch37, params

    # Example not using mane accession: CA10582926
    for params in [
        ("NP_001123914.1", 2, 2, CoordinateType.RESIDUE, Assembly.GRCH37),
        ("NP_001123914.1", 1, 2, CoordinateType.INTER_RESIDUE, Assembly.GRCH37),
        ("NP_001123914.1", 2, 2, CoordinateType.RESIDUE, Assembly.GRCH38),
        ("NP_001123914.1", 1, 2, CoordinateType.INTER_RESIDUE, Assembly.GRCH38),
    ]:
        ac, start, end, coordinate_type, assembly = params
        resp, w = await test_alignment_mapper.p_to_g(
            ac,
            start,
            end,
            coordinate_type=coordinate_type,
            target_genome_assembly=assembly,
        )
        assert w is None, params
        assert resp == hras_t2a(assembly), params


@pytest.mark.asyncio()
async def test_p_to_g_invalid(test_alignment_mapper):
    """Test invalid queries for p_to_g method"""
    # Invalid protein accession
    resp, w = await test_alignment_mapper.p_to_g(
        "NP_0000", 600, 600, CoordinateType.RESIDUE
    )
    assert resp is None
    assert w == "NP_0000 not found in transcript mappings"
