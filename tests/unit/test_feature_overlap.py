"""Module for testing Feature Overlap class"""
import pytest

from cool_seq_tool.data_sources.feature_overlap import (
    FeatureOverlap,
    FeatureOverlapError,
)
from cool_seq_tool.schemas import ResidueMode


@pytest.fixture(scope="module")
def test_feature_overlap(test_seqrepo_access):
    """Build Feature Overlap test fixture"""
    return FeatureOverlap(test_seqrepo_access)


def test_df(test_feature_overlap):
    """Test that the dataframe contains correct data"""
    # We only store CDS data
    assert list(test_feature_overlap.df["type"].unique()) == ["CDS"]

    assert set(test_feature_overlap.df.columns) == {
        "type",
        "chromosome",
        "cds_start",
        "cds_stop",
        "info_name",
        "gene",
    }

    assert test_feature_overlap.df["cds_start"].dtype == "int64"
    assert test_feature_overlap.df["cds_stop"].dtype == "int64"

    assert set(test_feature_overlap.df["chromosome"].unique()) == {
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "X",
        "Y",
    }


def test_get_chr_from_alt_ac(test_feature_overlap):
    """Test that _get_chr_from_alt_ac works correctly"""
    resp = test_feature_overlap._get_chr_from_alt_ac("NC_000001.11")
    assert resp == "1"

    resp = test_feature_overlap._get_chr_from_alt_ac("NC_000023.11")
    assert resp == "X"

    # identifier is invalid (no version)
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap._get_chr_from_alt_ac("NC_000001")
    assert str(e.value) == "SeqRepo unable to get translated identifiers for NC_000001"

    # identifier is grch37
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap._get_chr_from_alt_ac("NC_000001.10")
    assert str(e.value) == "Unable to find GRCh38 aliases for: NC_000001.10"


def test_get_grch38_cds_overlap(test_feature_overlap):
    """Test that get_grch38_mane_gene_cds_overlap works correctly"""
    # Variant fully contains exon (negative strand)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        140726490, 140726520, identifier="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
    )
    assert resp == {
        "BRAF": [
            {
                "cds": {
                    "_id": "ga4gh:VSL._H2ST69A4RkWCSRHOoMv-edt-R45fPdq",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726493, "type": "Number"},
                        "end": {"value": 140726516, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL._H2ST69A4RkWCSRHOoMv-edt-R45fPdq",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726493, "type": "Number"},
                        "end": {"value": 140726516, "type": "Number"},
                    },
                },
            }
        ]
    }

    # Using inter-residue (start == stop)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        140726500,
        140726500,
        identifier="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert resp == {
        "BRAF": [
            {
                "cds": {
                    "_id": "ga4gh:VSL._H2ST69A4RkWCSRHOoMv-edt-R45fPdq",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726493, "type": "Number"},
                        "end": {"value": 140726516, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.EqiyoLjrKnKg5F56bjRlBFBUihSkgX5w",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726500, "type": "Number"},
                        "end": {"value": 140726500, "type": "Number"},
                    },
                },
            }
        ]
    }

    # Variant is fully contained within exon (positive strand)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        55019308, 55019341, chromosome="7"
    )
    assert resp == {
        "EGFR": [
            {
                "cds": {
                    "_id": "ga4gh:VSL.fLukxPP69_vJU-1EYdNgM2waELFsJ0gI",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 55019277, "type": "Number"},
                        "end": {"value": 55019365, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.7CpeGak8icEAK4Eobs3ChFo9jqKdSidh",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 55019307, "type": "Number"},
                        "end": {"value": 55019341, "type": "Number"},
                    },
                },
            }
        ]
    }

    # Variant partially overlaps with exon, from the exon's start side (negative strand)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        140726503, 140726520, chromosome="7"
    )
    assert resp == {
        "BRAF": [
            {
                "cds": {
                    "_id": "ga4gh:VSL._H2ST69A4RkWCSRHOoMv-edt-R45fPdq",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726493, "type": "Number"},
                        "end": {"value": 140726516, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.l48V673TzeyfKeLfm3JUN-VTdEPqv80p",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726502, "type": "Number"},
                        "end": {"value": 140726516, "type": "Number"},
                    },
                },
            }
        ]
    }

    # Variant partially overlaps with exon, from the exon's stop side (negative strand)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        140726490, 140726505, identifier="NC_000007.14"
    )
    assert resp == {
        "BRAF": [
            {
                "cds": {
                    "_id": "ga4gh:VSL._H2ST69A4RkWCSRHOoMv-edt-R45fPdq",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726493, "type": "Number"},
                        "end": {"value": 140726516, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.QWkbljVU4ijJQsDS3vTVVfslmuuGbtEd",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 140726493, "type": "Number"},
                        "end": {"value": 140726505, "type": "Number"},
                    },
                },
            }
        ]
    }

    # Variant overlaps with multiple exons (positive strand)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        21522390, 21523491, chromosome="Y"
    )
    assert resp == {
        "RBMY1B": [
            {
                "cds": {
                    "_id": "ga4gh:VSL.b3uUu78bDJGPE2XnzWB4BHaIhP2PZRS1",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 21522382, "type": "Number"},
                        "end": {"value": 21522493, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.-NoMsiuXihe5cwF2ANgmFMDlji3k3qwC",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 21522389, "type": "Number"},
                        "end": {"value": 21522493, "type": "Number"},
                    },
                },
            },
            {
                "cds": {
                    "_id": "ga4gh:VSL.07koDLghr5iKDUK2WtK4R-MFonLIrAAa",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 21522934, "type": "Number"},
                        "end": {"value": 21523045, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.07koDLghr5iKDUK2WtK4R-MFonLIrAAa",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 21522934, "type": "Number"},
                        "end": {"value": 21523045, "type": "Number"},
                    },
                },
            },
            {
                "cds": {
                    "_id": "ga4gh:VSL.QtREnzEQbynxC-nvQLf9hm_3b3kPjBNI",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 21523479, "type": "Number"},
                        "end": {"value": 21523590, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.rdFImJ11G9vMjIeu_lMoeeBtDWuw5ssk",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 21523479, "type": "Number"},
                        "end": {"value": 21523491, "type": "Number"},
                    },
                },
            },
        ]
    }

    # Variant overlaps with multiple exons (negative strand)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        154779177, 154781317, chromosome="X"
    )
    assert resp == {
        "MPP1": [
            {
                "cds": {
                    "_id": "ga4gh:VSL.xwoKqCpJnxAMMx0BNDrsG-T05fFs3vzJ",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 154781238, "type": "Number"},
                        "end": {"value": 154781313, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.xwoKqCpJnxAMMx0BNDrsG-T05fFs3vzJ",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 154781238, "type": "Number"},
                        "end": {"value": 154781313, "type": "Number"},
                    },
                },
            },
            {
                "cds": {
                    "_id": "ga4gh:VSL.KfGMgFOaQRFXpG8X3mK75XbIEFMH6Sg9",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 154779176, "type": "Number"},
                        "end": {"value": 154779353, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.KfGMgFOaQRFXpG8X3mK75XbIEFMH6Sg9",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 154779176, "type": "Number"},
                        "end": {"value": 154779353, "type": "Number"},
                    },
                },
            },
        ]
    }

    # Variant overlap with cds in multiple genes and alt chromosome accession
    # chr19_KI270930v1_alt with exact start/stop CDS
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        135329, 135381, chromosome="19"
    )
    expected = {
        "KIR2DL5B": [
            {
                "cds": {
                    "_id": "ga4gh:VSL.SflvFED1mwAVOHbdM4FE2phzNPRSAi_J",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 135328, "type": "Number"},
                        "end": {"value": 135381, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.SflvFED1mwAVOHbdM4FE2phzNPRSAi_J",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 135328, "type": "Number"},
                        "end": {"value": 135381, "type": "Number"},
                    },
                },
            }
        ],
        "FCGBP": [
            {
                "cds": {
                    "_id": "ga4gh:VSL.6SSrOjGBl-txA6rGNXVfMFvsnNksxG4K",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 135263, "type": "Number"},
                        "end": {"value": 135807, "type": "Number"},
                    },
                },
                "overlap": {
                    "_id": "ga4gh:VSL.SflvFED1mwAVOHbdM4FE2phzNPRSAi_J",
                    "type": "SequenceLocation",
                    "sequence_id": "ga4gh:SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                    "interval": {
                        "type": "SequenceInterval",
                        "start": {"value": 135328, "type": "Number"},
                        "end": {"value": 135381, "type": "Number"},
                    },
                },
            }
        ],
    }
    assert resp == expected

    # Using inter-residue (start != stop)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        135328, 135381, chromosome="19", residue_mode=ResidueMode.INTER_RESIDUE
    )
    assert resp == expected

    # No overlap found
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(1, 2, chromosome="19")
    assert resp is None

    # Testing invalid

    # chromosome does not match regex pattern
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap.get_grch38_mane_gene_cds_overlap(
            154779177, 154781317, chromosome="chrX"
        )
    assert str(e.value) == "`chromosome` must be 1, ..., 22, X, or Y"

    # identifier is GRCh37
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap.get_grch38_mane_gene_cds_overlap(
            154779177, 154781317, identifier="NC_000023.10"
        )
    assert str(e.value) == "Unable to find GRCh38 aliases for: NC_000023.10"

    # no identifier or chromosome provided
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap.get_grch38_mane_gene_cds_overlap(154779177, 154781317)
    assert str(e.value) == "Must provide either `chromosome` or `identifier`"
