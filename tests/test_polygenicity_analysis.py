"""Polygenicity (#147): the unit, the tertile rule, the pre-comment check (report.polygenicity, first half)."""
import numpy as np
import pandas as pd

from report.polygenicity import cutpoints, tertile, unit


def test_cutpoints_are_the_33_and_67_percentiles():
    s = pd.Series(np.arange(1.0, 10.0))     # 1..9
    assert cutpoints(s) == (np.percentile(s, 100 / 3), np.percentile(s, 200 / 3))


def test_tertile_labels_and_none_for_unscored():
    s = pd.Series([1.0, 5.0, 9.0, np.nan])
    assert list(tertile(s, (3.0, 7.0))) == ["bottom", "middle", "top", "none"]
    assert list(tertile(pd.Series([3.0, 7.0]), (3.0, 7.0))) == ["bottom", "middle"]   # a cut value goes down


def test_unit_flags_a_disease_missing_from_the_profile():   # Review Focus 5
    d = pd.DataFrame({"mondo_id": ["M1", "M1", "M2", "M3", "M4"], "phase": [2, 3, 2, 2, 2],
                      "bio_area": ["a", "a", "a", "oncology", "a"], "us_site": [True, True, True, True, False]})
    prof = pd.DataFrame({"mondo_id": ["M1"], "score": [2.0], "ldsc_status": ["ok"]})
    u = unit(d, prof)
    assert list(u.mondo_id) == ["M1", "M2"] and list(u.entries) == [1, 1]
    assert list(u.ldsc_status) == ["ok", "not_profiled"] and list(u.scored) == [True, False]
