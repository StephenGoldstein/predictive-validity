"""Polygenicity (#147): the share of 2021 WHO DALYs in diseases with a Mendelian form or a drugged large-effect gene."""
import importlib.util
import os

import pandas as pd

from report import frame
from report.polygenicity import daly_share, leaves, who_map

WHO = os.path.join(frame.REPO, "analyses", "model-system-atlas", "lists", "who_ghe_2021_global_dalys.csv")


def burden_map():
    spec = importlib.util.spec_from_file_location("burden", os.path.join(frame.REPO, "analyses", "model-system-atlas", "burden.py"))
    burden = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(burden)
    return burden.MAP


def test_map_covers_every_who_cause():
    who = pd.read_csv(WHO)
    m = who_map()
    assert set(m.code) == set(who.code) and set(m.relation) <= {"exact", "broader", "none"}
    assert m[m.relation.ne("none")].mondo_id.notna().all()


def test_map_agrees_with_burden_py():
    pairs = set(zip(who_map().code.astype(int), who_map().mondo_id.fillna("")))
    for mondo, (grade, codes, _) in burden_map().items():
        if grade in ("exact", "broader"):
            for c in codes:
                assert (c, mondo) in pairs, (c, mondo)


def test_leaves_are_rows_without_children():
    who = pd.DataFrame({"code": [0, 10, 11, 12, 20], "depth": [1, 2, 3, 3, 2], "cause": ["All", "NCD", "A", "B", "Inj"],
                        "dalys_2021": [100, 60, 40, 20, 40]})
    assert list(leaves(who).code) == [11, 12, 20]


def test_daly_share_sums_leaf_causes_with_a_flag():
    who = pd.DataFrame({"code": [0, 10, 11, 12, 20], "depth": [1, 2, 3, 3, 2],
                        "cause": ["All Causes", "Noncommunicable diseases", "A", "B", "Injuries"],
                        "dalys_2021": [100, 60, 40, 20, 40]})
    m = pd.DataFrame({"code": [11, 12, 20], "relation": ["exact", "broader", "none"], "mondo_id": ["M1", "M2", None]})
    flags = pd.DataFrame({"mondo_id": ["M1", "M2"], "mendelian_form": [True, False], "drugged_large_effect": [False, True]})
    s = daly_share(who, m, flags)
    assert s["all"]["mendelian_form"] == 0.4 and s["all"]["drugged_large_effect"] == 0.2 and s["all"]["either"] == 0.6
    assert abs(s["ncd"]["either"] - 1.0) < 1e-9 and s["all"]["total"] == 100


def test_drugged_sets_join_entered_and_approved_targets():
    from report.polygenicity import drugged_sets
    td = pd.DataFrame({"ensembl_id": ["E1", "E2"], "mondo_id": ["M1", "M1"]})
    approved = pd.DataFrame({"ensembl_id": ["E3"], "mondo_id": ["M2"]})
    assert drugged_sets(td, approved) == {"M1": {"E1", "E2"}, "M2": {"E3"}}
