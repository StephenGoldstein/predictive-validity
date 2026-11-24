"""Polygenicity (#147): P1-P3, one-sided Holm, sensitivities and descriptives (report.polygenicity_tests) on
synthetic entries. Nothing here reads the database."""
import numpy as np
import pandas as pd

from report.polygenicity_tests import SIDES, cohort, one_sided, p3, registered, results_table, sensitivities


def synthetic(n=600, seed=0):
    rng = np.random.default_rng(seed)
    t = rng.choice(["bottom", "middle", "top"], n)
    g = rng.integers(0, 2, n)
    s4 = (rng.random(n) < np.where(t == "top", 0.15, 0.3) + 0.1 * g).astype(float)
    mode = (rng.random(n) < np.where(t == "top", 0.3, 0.15)).astype(float)
    score = np.where(t == "top", 3.0, np.where(t == "middle", 2.0, 1.0)) + rng.random(n) * 0.5
    return pd.DataFrame({"tertile": t, "genetics": g.astype(float), "s4": s4, "mode_not_causal": mode,
                         "phase": rng.choice([2, 3], n), "mondo_id": rng.integers(0, 60, n).astype(str),
                         "target": rng.integers(0, 120, n).astype(str), "has_target_disease": True,
                         "bio_area": rng.choice(["a", "b"], n), "modality": "m", "era": rng.choice(["e1", "e2"], n),
                         "novel": 0.0, "precedent": 0.0, "us_site": True, "log10_score": score,
                         "conc_tertile": rng.choice(["bottom", "middle", "top"], n),
                         "mendelian_form": rng.integers(0, 2, n).astype(bool), "log10_n_eff": 4 + rng.random(n),
                         "scored": True, "stop_no_benefit": mode})


def test_cohort_keeps_the_two_levels_with_a_resolved_target_and_the_right_phases():
    d = synthetic()
    c = cohort(d, "s4")
    assert set(c.tertile) == {"top", "bottom"} and c.has_target_disease.all() and set(c.phase) == {2}
    assert set(cohort(d, "mode_not_causal").phase) == {2, 3}
    assert list(c.top.unique()) in ([0.0, 1.0], [1.0, 0.0])


def test_one_sided_halves_p_on_the_registered_side():
    assert one_sided({"rr": 1.5, "p": 0.04}, 1)["p_one_sided"] == 0.02
    assert one_sided({"rr": 0.8, "p": 0.04}, 1)["p_one_sided"] == 0.98
    assert one_sided({}, 1) == {}


def test_registered_runs_three_tests_with_holm_and_readings():
    r = registered(synthetic())
    assert set(r) == set(SIDES) == {"P1", "P2", "P3"}
    assert all("p_holm" in r[k] and "reading" in r[k] for k in r)
    assert r["P1"]["rr"] > 1 and r["P2"]["rr"] < 1
    assert r["P1"]["reading"] in ("confirmed", "same direction, not confirmed")


def test_p3_is_not_estimable_when_a_cell_has_no_event():
    d = synthetic()
    d.loc[(d.tertile == "top") & (d.genetics == 1), "s4"] = 0.0
    r = p3(d)
    assert "rr" not in r and r["estimable"] is False and r["cells"]["top_genetics"]["events"] == 0


def test_sensitivities_have_the_four_registered_families():
    s = sensitivities(synthetic())
    assert set(s) == {"gwas_size", "continuous", "concentration", "mendelian_form"}
    assert set(s["continuous"]) == {"P1", "P2", "P3"} and "p_holm" in s["gwas_size"]["P1"]


def test_results_table_rows_and_columns():
    t = results_table({"primary": registered(synthetic())})
    assert list(t.columns) == ["family", "test", "rr", "lo", "hi", "p", "p_one_sided", "p_holm", "reading",
                               "exposed_entries", "exposed_events", "reference_entries", "reference_events"]
    assert list(t.test) == ["P1", "P2", "P3"]
