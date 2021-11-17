"""Module for testing seqrepo access class"""
import pytest
from uta_tools.data_sources import SeqRepoAccess


@pytest.fixture(scope='module')
def test_seqrepo_access():
    """Create SeqRepoAccess test fixture"""
    return SeqRepoAccess()


def test__get_start_end(test_seqrepo_access):
    """Test that _get_start_end method works correctly."""
    resp = test_seqrepo_access._get_start_end(600)
    assert resp == ((599, 600), None)

    resp = test_seqrepo_access._get_start_end(600, 600)
    assert resp == ((599, 600), None)

    resp = test_seqrepo_access._get_start_end(
        600, end=600, residue_mode="residue")
    assert resp == ((599, 600), None)

    resp = test_seqrepo_access._get_start_end(
        600, residue_mode="inter-residue")
    assert resp == ((600, 601), None)

    resp = test_seqrepo_access._get_start_end(
        600, end=600, residue_mode="inter-residue")
    assert resp == ((600, 601), None)

    resp = test_seqrepo_access._get_start_end(600, residue_mode="mode")
    assert resp == (None, "residue_mode must be either `inter-residue` or "
                          "`residue`, not `mode`")


def test_is_valid_input_sequence(test_seqrepo_access):
    """Test that is_valid_input_sequence method works correctly"""
    resp = test_seqrepo_access.is_valid_input_sequence("NP_004324.2", 600)
    assert resp == (True, None)

    resp = test_seqrepo_access.is_valid_input_sequence("NP_004324.2", 600, 601)
    assert resp == (True, None)

    resp = test_seqrepo_access.is_valid_input_sequence("NP_004324.2", 600, 600)
    assert resp == (True, None)

    resp = test_seqrepo_access.is_valid_input_sequence("NP_0043241311412", 600)
    assert resp == (False, "Accession, NP_0043241311412, not found in SeqRepo")

    resp = test_seqrepo_access.is_valid_input_sequence("NP_004324.2", 601, 600)
    assert resp == (False, "Invalid inter-residue coordinates: start (600)"
                           " cannot be greater than end (599)")

    resp = test_seqrepo_access.is_valid_input_sequence(
        "NP_004324.2", 4654645645654, 1)
    assert resp == (False, "Start inter-residue coordinate (4654645645653) is "
                           "out of index on NP_004324.2")

    resp = test_seqrepo_access.is_valid_input_sequence(
        "NP_004324.2", 600, 4654645645654)
    assert resp == (False, "End inter-residue coordinate (4654645645653) is"
                           " out of index on NP_004324.2")


def test_get_reference_sequence(test_seqrepo_access):
    """Test that get_reference_sequence method works correctly"""
    resp = test_seqrepo_access.get_reference_sequence("NP_004324.2", 600)
    assert resp == ("V", None)

    resp = test_seqrepo_access.get_reference_sequence("NP_004324.2", 600, 600)
    assert resp == ("V", None)

    resp = test_seqrepo_access.get_reference_sequence("NP_004324.2", 600, 601)
    assert resp == ("V", None)

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", 599, 600, residue_mode="inter-residue")
    assert resp == ("V", None)

    resp = test_seqrepo_access.get_reference_sequence("NP_004324.2", 601, 600)
    assert resp == (None, "Invalid inter-residue coordinates: start (600)"
                          " cannot be greater than end (599)")

    resp = test_seqrepo_access.get_reference_sequence("NP_0043241311412", 600)
    assert resp == (None, "Accession, NP_0043241311412, not found in SeqRepo")

    resp = test_seqrepo_access.get_reference_sequence("NP_004324.2", 600, 800)
    assert resp == (None, "End inter-residue coordinate (799) "
                          "is out of index on NP_004324.2")

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", 4654645645654, 1)
    assert resp == (None, "Start inter-residue coordinate (4654645645653) is "
                          "out of index on NP_004324.2")

    resp = test_seqrepo_access.get_reference_sequence(
        "NP_004324.2", 600, 4654645645654)
    assert resp == (None, "End inter-residue coordinate (4654645645653) "
                          "is out of index on NP_004324.2")


def test_translate_identifier(test_seqrepo_access):
    """Test that translate_identifier method works correctly"""
    expected = (["ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"], None)
    resp = test_seqrepo_access.translate_identifier(
        "NM_152263.3", target_namespace="ga4gh")
    assert resp == expected

    resp = test_seqrepo_access.translate_identifier(
        "refseq:NM_152263.3", target_namespace="ga4gh")
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
    assert resp == ([], 'SeqRepo unable to get translated identifiers for'
                        ' refseq_152263.3')


def test_aliases(test_seqrepo_access):
    """Test that aliases method works correctly"""
    expected = (["ga4gh:SQ.ijXOSP3XSsuLWZhXQ7_TJ5JXu4RJO6VT"], None)
    resp = test_seqrepo_access.aliases("NM_152263.3")
    assert len(resp[0]) > 0
    assert resp[1] is None
    assert expected[0][0] in resp[0]

    resp = test_seqrepo_access.aliases("NC_000002.12")
    assert len(resp[0]) > 0
    assert resp[1] is None
    assert "GRCh38:2" in resp[0]

    resp = test_seqrepo_access.aliases("refseq_152263.3")
    assert resp == ([], 'SeqRepo could not translate alias refseq_152263.3')

    resp = test_seqrepo_access.aliases("GRCh38:2")
    assert resp == ([], 'SeqRepo could not translate alias GRCh38:2')
