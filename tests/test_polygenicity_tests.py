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
                         "scored": True, "stop_no_benefit": mode,
                         "stop_unrelated": (rng.random(n) < 0.1).astype(float), "stop_safety": (rng.random(n) < 0.05).astype(float)})


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


def test_sensitivities_have_the_three_amended_families():   # Melissa's 2026-11-24 decision
    s = sensitivities(synthetic())
    assert set(s) == {"gwas_size", "mendelian_form", "negative_control"}
    assert set(s["mendelian_form"]) == {"P1", "P2", "P3"} and "p_holm" in s["gwas_size"]["P1"]
    assert set(s["negative_control"]) == {"P1"}


def test_negative_control_leaves_safety_stops_out_and_is_two_sided():
    from report.polygenicity_tests import negative_control
    d = synthetic()
    r = negative_control(d)
    c = cohort(d, "mode_not_causal")
    assert r["exposed_entries"] + r["reference_entries"] == int((c.stop_safety <= 0).sum())
    assert r["exposed_events"] == int(c.stop_unrelated[(c.top > 0) & (c.stop_safety <= 0)].sum())
    assert "rr" in r and "p" in r and "p_one_sided" not in r and "p_holm" not in r
    assert r["reading"] in ("no association", "association")


def test_power_simulation_reports_power_per_test_and_effect():
    from report.polygenicity_tests import power_simulation
    out = power_simulation(synthetic(), reps=3, seed=1)
    assert set(out["power"]) == {"P1", "P2", "P3"} and out["reps"] == 3
    for by_effect in out["power"].values():
        assert by_effect and all(0.0 <= v["power"] <= 1.0 for v in by_effect.values())
    assert out["base_rates"]["mode_not_causal"] > 0 and out["base_rates"]["s4"] > 0
    assert {"top", "bottom"} <= set(out["entries"])


def test_results_table_rows_and_columns():
    t = results_table({"primary": registered(synthetic())})
    assert list(t.columns) == ["family", "test", "rr", "lo", "hi", "p", "p_one_sided", "p_holm", "reading",
                               "exposed_entries", "exposed_events", "reference_entries", "reference_events"]
    assert list(t.test) == ["P1", "P2", "P3"]


def test_numbers_key_the_tertile_rates_and_cutpoints_and_carry_the_extras():
    from report.polygenicity_tests import numbers
    d = synthetic()
    desc = {"by_tertile": [{"tertile": "bottom", "entries": 10}, {"tertile": "top", "entries": 12}], "no_score": {"diseases": 2}}
    n = numbers(d, {"P1": {}}, {}, desc, cuts=(1.5, 2.5), extra={"daly": {"flags": {"diseases": 3}}})
    assert n["tertiles"]["top"]["entries"] == 12 and n["tertiles"]["bottom"]["entries"] == 10
    assert n["cutpoints"]["score"] == {"lower": 1.5, "upper": 2.5} and "concentration" not in n["cutpoints"]
    assert n["daly"]["flags"]["diseases"] == 3 and n["cohort"]["entries"] == len(d)
    assert n["cohort"]["scored_diseases"] == d.mondo_id.nunique()


def test_descriptives_always_report_each_group_size(monkeypatch):
    from report.polygenicity_tests import descriptives
    monkeypatch.setattr("report.polygenicity_tests.engaged_modes", lambda: pd.DataFrame({"key": [], "not_causal_basis": []}))
    d = synthetic().assign(k_source="panukb", not_causal_basis="history")
    ids = sorted(d.mondo_id.unique())
    u = pd.DataFrame({"mondo_id": ids, "label": ids, "scored": True, "entries": 10, "ldsc_status": "ok", "h2_z": 1.0,
                      "h2_liab": 0.2, "loci": 5, "loci_catalog": None, "gwas_genes_per_10k_cases": None, "me_log10": None,
                      "score": 25.0})
    out = descriptives(d, u)
    assert out["z_at_least_2"]["diseases"] == 0 and out["z_at_least_2"]["tests"] == {}
    assert out["gbd_k_only"]["diseases"] == 0 and out["catalog_loci"]["diseases"] == 0
    assert out["gwas_genes"]["diseases"] == 0 and out["oconnor"]["n"] == 0


def test_holm_keeps_the_registered_family_of_three():
    from report.polygenicity_tests import holm_family
    adj = holm_family({"P1": 0.01, "P2": 0.04, "P3": None})
    assert adj["P1"] == 0.03 and adj["P2"] == 0.08 and "P3" not in adj
