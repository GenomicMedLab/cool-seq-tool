"""Module for testing LiftOver class"""

import pytest

from cool_seq_tool.schemas import Assembly


@pytest.fixture(scope="module")
def test_liftover(test_cool_seq_tool):
    """Build LiftOver test fixture"""
    return test_cool_seq_tool.liftover


def test_get_liftover(test_liftover):
    """Test that get_liftover works correctly."""
    resp = test_liftover.get_liftover("chr7", 140453136, Assembly.GRCH38)
    assert resp == ("chr7", 140753336)
    resp = test_liftover.get_liftover("7", 140453136, Assembly.GRCH38)
    assert resp == ("chr7", 140753336)

    resp = test_liftover.get_liftover("chr17", 140453136, Assembly.GRCH38)
    assert resp is None
