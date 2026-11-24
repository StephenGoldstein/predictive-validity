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


def test_tests_table_has_one_row_per_hypothesis(tmp_path, monkeypatch):
    from report.render import polygenicity_tests_table
    t = pd.DataFrame({"family": ["primary"] * 3, "test": ["P1", "P2", "P3"], "rr": [1.2, 0.8, 0.9], "lo": [1, .6, .5],
                      "hi": [1.5, 1, 1.4], "p": [.04, .08, .6], "p_one_sided": [.02, .04, .3], "p_holm": [.06, .08, .3],
                      "reading": ["same direction, not confirmed"] * 3, "exposed_entries": [100] * 3,
                      "exposed_events": [20] * 3, "reference_entries": [100] * 3, "reference_events": [10] * 3})
    (tmp_path / "polygenicity").mkdir()
    t.to_csv(tmp_path / "polygenicity" / "tests.csv", index=False)
    monkeypatch.setattr("report.render.OUT", str(tmp_path))
    out = polygenicity_tests_table()
    assert out.count("\n| P") == 3 and "1.20 (1.00–1.50)" in out


def test_area_and_no_score_tables_read_their_csvs(tmp_path, monkeypatch):
    from report.render import polygenicity_area_table, polygenicity_no_score_table
    (tmp_path / "polygenicity").mkdir()
    pd.DataFrame({"tertile": ["top", "top", "bottom"], "bio_area": ["a", "b", "a"], "entries": [5, 6, 7],
                  "diseases": [1, 2, 3]}).to_csv(tmp_path / "polygenicity" / "area.csv", index=False)
    pd.DataFrame({"ldsc_status": ["no_study", "low_z"], "diseases": [10, 2], "entries": [40, 9]}) \
        .to_csv(tmp_path / "polygenicity" / "no_score.csv", index=False)
    monkeypatch.setattr("report.render.OUT", str(tmp_path))
    area = polygenicity_area_table()
    assert "| a | 1 / 5 | 3 / 7 |" in area and "| b | 2 / 6 | 0 / 0 |" in area   # one row per area, a column per tertile
    ns = polygenicity_no_score_table()
    assert "| No GWAS study of the disease or a subtype | 10 | 40 |" in ns and ns.index("No GWAS") < ns.index("z = 4")


def test_membership_carries_the_study_heritability_prevalence_and_loci():
    from report.polygenicity import membership
    u = pd.DataFrame({"mondo_id": ["M1", "M2"], "label": ["a", "b"], "bio_area": ["x", "x"], "entries": [3, 4],
                      "ldsc_status": ["ok", "ok"], "gcst": ["GCST1", "GCST2"], "trait_level": ["exact", "exact"],
                      "k_source": ["gbd", "panukb"], "score": [10.0, 40.0], "concentration": [0.2, 0.5],
                      "mendelian_form": [True, False], "h2_obs": [0.1, 0.2], "h2_z": [5.0, 6.0], "k": [0.01, 0.02],
                      "h2_liab": [0.3, 0.4], "loci": [3, 16], "n_cases": [1000, 2000], "n_controls": [5000, 6000]})
    m = membership(u, (20.0, 30.0), (0.3, 0.4))
    assert list(m.tertile) == ["bottom", "top"] and list(m.conc_tertile) == ["bottom", "top"]
    assert {"h2_obs", "h2_z", "k", "h2_liab", "loci", "n_cases", "n_controls"} <= set(m.columns)
