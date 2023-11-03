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
        "chrom_normalized",
        "cds_start",
        "cds_stop",
        "info_name",
        "gene",
    }

    assert test_feature_overlap.df["cds_start"].dtype == "int64"
    assert test_feature_overlap.df["cds_stop"].dtype == "int64"

    assert set(test_feature_overlap.df["chrom_normalized"].unique()) == {
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
    """Test that get_grch38_cds_overlap works correctly"""
    # Variant fully contains exon (negative strand)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        140726490, 140726520, identifier="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul"
    )
    assert resp == {
        "BRAF": [
            {
                "cds_start": 140726494,
                "cds_stop": 140726516,
                "overlap_start": 140726494,
                "overlap_stop": 140726516,
            }
        ]
    }

    # Using inter-residue (start == stop)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        140726500,
        140726500,
        identifier="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        residue_mode=ResidueMode.INTER_RESIDUE,
    )
    assert resp == {
        "BRAF": [
            {
                "cds_start": 140726494,
                "cds_stop": 140726516,
                "overlap_start": 140726500,
                "overlap_stop": 140726501,
            }
        ]
    }

    # Variant is fully contained within exon (positive strand)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        55019308, 55019341, chromosome="7"
    )
    assert resp == {
        "EGFR": [
            {
                "cds_start": 55019278,
                "cds_stop": 55019365,
                "overlap_start": 55019308,
                "overlap_stop": 55019341,
            }
        ]
    }

    # Variant partially overlaps with exon, from the exon's start side (negative strand)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        140726503, 140726520, chromosome="7"
    )
    assert resp == {
        "BRAF": [
            {
                "cds_start": 140726494,
                "cds_stop": 140726516,
                "overlap_start": 140726503,
                "overlap_stop": 140726516,
            }
        ]
    }

    # Variant partially overlaps with exon, from the exon's stop side (negative strand)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        140726490, 140726505, identifier="NC_000007.14"
    )
    assert resp == {
        "BRAF": [
            {
                "cds_start": 140726494,
                "cds_stop": 140726516,
                "overlap_start": 140726494,
                "overlap_stop": 140726505,
            }
        ]
    }

    # Variant overlaps with multiple exons (positive strand)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        21522390, 21523491, chromosome="Y"
    )
    assert resp == {
        "RBMY1B": [
            {
                "cds_start": 21522383,
                "cds_stop": 21522493,
                "overlap_start": 21522390,
                "overlap_stop": 21522493,
            },
            {
                "cds_start": 21522935,
                "cds_stop": 21523045,
                "overlap_start": 21522935,
                "overlap_stop": 21523045,
            },
            {
                "cds_start": 21523480,
                "cds_stop": 21523590,
                "overlap_start": 21523480,
                "overlap_stop": 21523491,
            },
        ]
    }

    # Variant overlaps with multiple exons (negative strand)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        154779177, 154781317, chromosome="X"
    )
    assert resp == {
        "MPP1": [
            {
                "cds_start": 154781239,
                "cds_stop": 154781313,
                "overlap_start": 154781239,
                "overlap_stop": 154781313,
            },
            {
                "cds_start": 154779177,
                "cds_stop": 154779353,
                "overlap_start": 154779177,
                "overlap_stop": 154779353,
            },
        ]
    }

    # Variant overlap with cds in multiple genes and alt chromosome accession
    # chr19_KI270930v1_alt with exact start/stop CDS
    resp = test_feature_overlap.get_grch38_cds_overlap(135329, 135381, chromosome="19")
    expected = {
        "KIR2DL5B": [
            {
                "cds_start": 135329,
                "cds_stop": 135381,
                "overlap_start": 135329,
                "overlap_stop": 135381,
            }
        ],
        "FCGBP": [
            {
                "cds_start": 135264,
                "cds_stop": 135807,
                "overlap_start": 135329,
                "overlap_stop": 135381,
            }
        ],
    }
    assert resp == expected

    # Using inter-residue (start != stop)
    resp = test_feature_overlap.get_grch38_cds_overlap(
        135328, 135381, chromosome="19", residue_mode=ResidueMode.INTER_RESIDUE
    )
    assert resp == expected

    # Testing invalid

    # chromosome does not match regex pattern
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap.get_grch38_cds_overlap(
            154779177, 154781317, chromosome="chrX"
        )
    assert str(e.value) == "`chromosome` must be 1, ..., 22, X, or Y"

    # identifier is GRCh37
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap.get_grch38_cds_overlap(
            154779177, 154781317, identifier="NC_000023.10"
        )
    assert str(e.value) == "Unable to find GRCh38 aliases for: NC_000023.10"

    # no identifier or chromosome provided
    with pytest.raises(FeatureOverlapError) as e:
        test_feature_overlap.get_grch38_cds_overlap(154779177, 154781317)
    assert str(e.value) == "Must provide either `chromosome` or `identifier`"
