"""Integrity of the committed processed dataset (skipped if it is absent)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from supplygraphs import data, leontief as lf

PROC = Path(__file__).resolve().parents[1] / "data" / "processed"
pytestmark = pytest.mark.skipif(not (PROC / "A_detail_2017.csv").exists(), reason="no processed data")


def test_shapes_and_ordering():
    for y in data.DETAIL_YEARS:
        A, ind = data.load(y)
        assert A.shape == (400, 400)
    assert list(ind["code"]) == sorted(ind["code"])
    assert list(ind["index"]) == list(range(1, 401))
    assert ind["sector"].notna().all()


def test_requirements_are_valid():
    for y in data.DETAIL_YEARS:
        A, _ = data.load(y)
        assert np.isfinite(A).all() and (A >= 0).all()
        assert lf.spectral_radius(A) < 1
    for y in data.SUMMARY_YEARS:
        As = data.load_summary(y)
        assert As.shape == (71, 71) and (As >= 0).all()
        assert As.sum(axis=0).max() < 1


def test_two_levels_agree_on_spectral_radius():
    # Not a theorem: a consistency check that the detail and summary tables
    # were built the same way. The two 2017 values agree to about 0.1%.
    A, _ = data.load(2017)
    assert lf.spectral_radius(A) == pytest.approx(lf.spectral_radius(data.load_summary(2017)), rel=0.01)
