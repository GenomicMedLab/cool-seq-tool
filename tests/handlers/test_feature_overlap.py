"""Module for testing Feature Overlap class"""

import polars as pl
import pytest

from cool_seq_tool.mappers.feature_overlap import (
    FeatureOverlap,
    FeatureOverlapError,
)
from cool_seq_tool.schemas import CoordinateType


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

    assert test_feature_overlap.df["cds_start"].dtype == pl.Int64
    assert test_feature_overlap.df["cds_stop"].dtype == pl.Int64

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
                    "id": "ga4gh:SL.fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "digest": "fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726493,
                    "end": 140726516,
                },
                "overlap": {
                    "id": "ga4gh:SL.fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "digest": "fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726493,
                    "end": 140726516,
                },
            }
        ]
    }

    expected = {
        "BRAF": [
            {
                "cds": {
                    "id": "ga4gh:SL.fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "digest": "fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726493,
                    "end": 140726516,
                },
                "overlap": {
                    "id": "ga4gh:SL.rMJlP7STVHBdvcCMgHkA4XJXXIdXnsix",
                    "digest": "rMJlP7STVHBdvcCMgHkA4XJXXIdXnsix",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726500,
                    "end": 140726501,
                },
            }
        ]
    }

    # Using residue (start == stop)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        140726501,
        140726501,
        identifier="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        coordinate_type=CoordinateType.RESIDUE,
    )
    assert resp == expected

    # Using inter-residue
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        140726500,
        140726501,
        identifier="ga4gh:SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
        coordinate_type=CoordinateType.INTER_RESIDUE,
    )
    assert resp == expected

    # Variant is fully contained within exon (positive strand)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        55019308, 55019341, chromosome="7"
    )
    assert resp == {
        "EGFR": [
            {
                "cds": {
                    "id": "ga4gh:SL.vjxcgicBFEkN8b8AXhagvUDC7FZgZgCp",
                    "digest": "vjxcgicBFEkN8b8AXhagvUDC7FZgZgCp",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 55019277,
                    "end": 55019365,
                },
                "overlap": {
                    "id": "ga4gh:SL.a_MHSA9TJ5zMkxd52eBuRUNb5ZIXHH7T",
                    "digest": "a_MHSA9TJ5zMkxd52eBuRUNb5ZIXHH7T",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 55019307,
                    "end": 55019341,
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
                    "id": "ga4gh:SL.fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "digest": "fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726493,
                    "end": 140726516,
                },
                "overlap": {
                    "id": "ga4gh:SL.MdSOBEGp0l8wT3y1taeRvVIEi_XDBIGK",
                    "digest": "MdSOBEGp0l8wT3y1taeRvVIEi_XDBIGK",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726502,
                    "end": 140726516,
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
                    "id": "ga4gh:SL.fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "digest": "fYRYzNIAoe6UQF9MT1XaYsFscoU68ZJv",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726493,
                    "end": 140726516,
                },
                "overlap": {
                    "id": "ga4gh:SL.Rjvup1y8hPgveXiYnj7dipqYkt3BaFZE",
                    "digest": "Rjvup1y8hPgveXiYnj7dipqYkt3BaFZE",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.F-LrLMe1SRpfUZHkQmvkVKFEGaoDeHul",
                        "type": "SequenceReference",
                    },
                    "start": 140726493,
                    "end": 140726505,
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
                    "id": "ga4gh:SL.3fbgdG4Z2a1fqkqkY8M2bQMfBhJsVr_i",
                    "digest": "3fbgdG4Z2a1fqkqkY8M2bQMfBhJsVr_i",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "type": "SequenceReference",
                    },
                    "start": 21522382,
                    "end": 21522493,
                },
                "overlap": {
                    "id": "ga4gh:SL.XSqmOKSXECFtfjhcvVDmga72xQ0jVGO7",
                    "digest": "XSqmOKSXECFtfjhcvVDmga72xQ0jVGO7",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "type": "SequenceReference",
                    },
                    "start": 21522389,
                    "end": 21522493,
                },
            },
            {
                "cds": {
                    "id": "ga4gh:SL.wi_fCVQHmZCOUf--3UPKm6johAu3zQYJ",
                    "digest": "wi_fCVQHmZCOUf--3UPKm6johAu3zQYJ",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "type": "SequenceReference",
                    },
                    "start": 21522934,
                    "end": 21523045,
                },
                "overlap": {
                    "id": "ga4gh:SL.wi_fCVQHmZCOUf--3UPKm6johAu3zQYJ",
                    "digest": "wi_fCVQHmZCOUf--3UPKm6johAu3zQYJ",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "type": "SequenceReference",
                    },
                    "start": 21522934,
                    "end": 21523045,
                },
            },
            {
                "cds": {
                    "id": "ga4gh:SL.cnULzfdHZPiY6rSQP2Prfzocf1_YOBGA",
                    "digest": "cnULzfdHZPiY6rSQP2Prfzocf1_YOBGA",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "type": "SequenceReference",
                    },
                    "start": 21523479,
                    "end": 21523590,
                },
                "overlap": {
                    "id": "ga4gh:SL.kHLoQs5cjOqjztK9yZyyVR9ScKCM1S8f",
                    "digest": "kHLoQs5cjOqjztK9yZyyVR9ScKCM1S8f",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.8_liLu1aycC0tPQPFmUaGXJLDs5SbPZ5",
                        "type": "SequenceReference",
                    },
                    "start": 21523479,
                    "end": 21523491,
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
                    "id": "ga4gh:SL.dnHzxh-VwjVdanLcvKkI1otKhZeY223-",
                    "digest": "dnHzxh-VwjVdanLcvKkI1otKhZeY223-",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "type": "SequenceReference",
                    },
                    "start": 154781238,
                    "end": 154781313,
                },
                "overlap": {
                    "id": "ga4gh:SL.dnHzxh-VwjVdanLcvKkI1otKhZeY223-",
                    "digest": "dnHzxh-VwjVdanLcvKkI1otKhZeY223-",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "type": "SequenceReference",
                    },
                    "start": 154781238,
                    "end": 154781313,
                },
            },
            {
                "cds": {
                    "id": "ga4gh:SL.Z4jQtiT0-FplZWGVA1wNhdCaCMKGQ17D",
                    "digest": "Z4jQtiT0-FplZWGVA1wNhdCaCMKGQ17D",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "type": "SequenceReference",
                    },
                    "start": 154779176,
                    "end": 154779353,
                },
                "overlap": {
                    "id": "ga4gh:SL.Z4jQtiT0-FplZWGVA1wNhdCaCMKGQ17D",
                    "digest": "Z4jQtiT0-FplZWGVA1wNhdCaCMKGQ17D",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.w0WZEvgJF0zf_P4yyTzjjv9oW1z61HHP",
                        "type": "SequenceReference",
                    },
                    "start": 154779176,
                    "end": 154779353,
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
                    "id": "ga4gh:SL.tR0TL0hHD3udyK9at0snGQ3zNSmhCz6K",
                    "digest": "tR0TL0hHD3udyK9at0snGQ3zNSmhCz6K",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                        "type": "SequenceReference",
                    },
                    "start": 135328,
                    "end": 135381,
                },
                "overlap": {
                    "id": "ga4gh:SL.tR0TL0hHD3udyK9at0snGQ3zNSmhCz6K",
                    "digest": "tR0TL0hHD3udyK9at0snGQ3zNSmhCz6K",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                        "type": "SequenceReference",
                    },
                    "start": 135328,
                    "end": 135381,
                },
            }
        ],
        "FCGBP": [
            {
                "cds": {
                    "id": "ga4gh:SL.3G3gZfvJ56-y-TRWNSAUHUmyPi_8X3qK",
                    "digest": "3G3gZfvJ56-y-TRWNSAUHUmyPi_8X3qK",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                        "type": "SequenceReference",
                    },
                    "start": 135263,
                    "end": 135807,
                },
                "overlap": {
                    "id": "ga4gh:SL.tR0TL0hHD3udyK9at0snGQ3zNSmhCz6K",
                    "digest": "tR0TL0hHD3udyK9at0snGQ3zNSmhCz6K",
                    "type": "SequenceLocation",
                    "sequenceReference": {
                        "refgetAccession": "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl",
                        "type": "SequenceReference",
                    },
                    "start": 135328,
                    "end": 135381,
                },
            }
        ],
    }
    assert resp == expected

    # Using inter-residue (start != stop)
    resp = test_feature_overlap.get_grch38_mane_gene_cds_overlap(
        135328, 135381, chromosome="19", coordinate_type=CoordinateType.INTER_RESIDUE
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
