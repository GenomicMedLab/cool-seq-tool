"""Module for testing the Liftover class"""
import pytest

from cool_seq_tool.mappers import Liftover
from cool_seq_tool.schemas import Assembly


@pytest.fixture(scope="module")
def liftover():
    """Create test fixture for Liftover"""
    return Liftover()


def test_convert_pos(liftover):
    """Test that the convert_pos method works correctly"""
    grch37_pos = 140453136
    grch38_pos = 140753336

    # GRCh37 -> GRCh38
    pos = liftover.convert_pos("chr7", grch37_pos, Assembly.GRCH38)
    assert pos == grch38_pos

    # GRCh38 -> GRCh37
    pos = liftover.convert_pos("chr7", grch38_pos, Assembly.GRCH37)
    assert pos == grch37_pos

    # invalid chr
    pos = liftover.convert_pos("7", grch37_pos)
    assert pos is None
