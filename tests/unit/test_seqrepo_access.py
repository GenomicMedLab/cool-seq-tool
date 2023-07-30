"""Module for testing seqrepo access class"""
import pytest
from biocommons.seqrepo import SeqRepo

from cool_seq_tool.data_sources import SeqRepoAccess
from cool_seq_tool.paths import SEQREPO_ROOT_DIR
from cool_seq_tool.schemas import ResidueMode


@pytest.fixture(scope="module")
def test_seqrepo_access():
    """Create SeqRepoAccess test fixture"""
    return SeqRepoAccess(SeqRepo(root_dir=SEQREPO_ROOT_DIR))


def test_get_reference_sequence(test_seqrepo_access):
    """Test that get_reference_sequence method works correctly"""
    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", start=600, end=600, residue_mode=ResidueMode.RESIDUE
    )
    assert resp == ("V", None)

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", start=600, end=601, residue_mode=ResidueMode.RESIDUE
    )
    assert resp == ("VK", None)

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", start=599, end=600, residue_mode=ResidueMode.INTER_RESIDUE
    )
    assert resp == ("V", None)

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", start=599, end=599, residue_mode=ResidueMode.INTER_RESIDUE
    )
    assert resp == ("V", None)

    # Test getting entire sequence
    resp = test_seqrepo_access.get_reference_sequence("NP_004324.2")
    assert len(resp[0]) == 766

    resp = test_seqrepo_access.get_reference_sequence("NP_0043241311412", start=600)
    assert resp == ("", "Accession, NP_0043241311412, not found in SeqRepo")

    resp = test_seqrepo_access.get_reference_sequence("NP_004324.2", start=600, end=800)
    assert resp == ("", "End inter-residue coordinate (800) "
                    "is out of index on NP_004324.2")

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", start=4654645645654, end=1)
    assert resp == ("", "Start inter-residue coordinate (4654645645653) is "
                    "out of index on NP_004324.2")

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", start=600, end=4654645645654)
    assert resp == ("", "End inter-residue coordinate (4654645645654) "
                    "is out of index on NP_004324.2")


def test_translate_identifier(test_seqrepo_access):
    """Test that translate_identifier method works correctly"""
    expected = (["ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"], None)
    resp = test_seqrepo_access.translate_identifier(
        "NM_152263.3", target_namespaces="ga4gh")
    assert resp == expected

    resp = test_seqrepo_access.translate_identifier(
        "refseq:NM_152263.3", target_namespaces="ga4gh")
    assert resp == expected

    resp = test_seqrepo_access.translate_identifier("refseq:NM_152263.3")
    assert len(resp[0]) > 0
    assert resp[1] is None
    assert expected[0][0] in resp[0]

    resp = test_seqrepo_access.translate_identifier("GRCh38:2")
    assert len(resp[0]) > 0
    assert resp[1] is None
    assert "refseq:NC_000002.12" in resp[0]

    resp = test_seqrepo_access.translate_identifier("NC_000002.12")
    assert len(resp[0]) > 0
    assert resp[1] is None
    assert "refseq:NC_000002.12" in resp[0]

    resp = test_seqrepo_access.translate_identifier("refseq_152263.3")
    assert resp == ([], "SeqRepo unable to get translated identifiers for"
                        " refseq_152263.3")


def test_aliases(test_seqrepo_access):
    """Test that aliases method works correctly"""
    expected = (["ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"], None)
    resp = test_seqrepo_access.translate_alias("NM_152263.3")
    assert len(resp[0]) > 0
    assert resp[1] is None
    assert expected[0][0] in resp[0]

    resp = test_seqrepo_access.translate_alias("NC_000002.12")
    assert len(resp[0]) > 0
    assert resp[1] is None
    assert "GRCh38:2" in resp[0]

    resp = test_seqrepo_access.translate_alias("refseq_152263.3")
    assert resp == ([], "SeqRepo could not translate alias refseq_152263.3")

    resp = test_seqrepo_access.translate_alias("GRCh38:2")
    assert resp == ([], "SeqRepo could not translate alias GRCh38:2")


def test_chromosome_to_acs(test_seqrepo_access):
    """Test that chromosome_to_acs method works correctly"""
    resp = test_seqrepo_access.chromosome_to_acs("7")
    assert resp == (["NC_000007.14", "NC_000007.13"], None)

    resp = test_seqrepo_access.chromosome_to_acs("X")
    assert resp == (["NC_000023.11", "NC_000023.10"], None)

    resp = test_seqrepo_access.chromosome_to_acs("Y")
    assert resp == (["NC_000024.10", "NC_000024.9"], None)

    resp = test_seqrepo_access.chromosome_to_acs("117")
    assert resp == (None, "117 is not a valid chromosome")


def test_ac_to_chromosome(test_seqrepo_access):
    """Test that ac_to_chromosome method works correctly"""
    resp = test_seqrepo_access.ac_to_chromosome("NC_000007.13")
    assert resp == ("7", None)

    resp = test_seqrepo_access.ac_to_chromosome("NC_000007.1323")
    assert resp == (None, "Unable to get chromosome for NC_000007.1323")
