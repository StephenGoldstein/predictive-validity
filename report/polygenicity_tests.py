"""Polygenicity (#147), the half that joins outcomes: P1-P3, the registered sensitivities and the descriptives.

Imported by report.polygenicity.main only once CUTPOINTS is set (the tertiles fixed in the dated comment). The
definitions: entries with a resolved primary target (has_target_disease); P1 on Phase 2 and 3 entries with
COVARS_23, P2 and P3 on Phase 2 entries with COVARS; the top tertile against the bottom, the middle reported; errors
clustered by disease; p one-sided in the registered direction (SIDES) and Holm across the three, read as the #140
confirmation reads its tests (confirmed, same direction not confirmed, not replicated). P3 is the interaction of
the top tertile with human genetic support (g3 > 0) in one model, estimable when each of its four cells holds at
least MIN_TARGETS target sets and one S4 event.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from pipeline.stops import keys
from report.engaged_modes import engaged_modes
from report.evidence_analysis import COVARS_23, traced_stops
from report.failure_modes import with_modes
from report.section3 import COVARS, MIN_TARGETS, coef, describe, estimable, estimate, holm, jsonable, poisson

SIDES = {"P1": 1, "P2": -1, "P3": -1}
OUTCOMES = {"P1": "mode_not_causal", "P2": "s4", "P3": "s4"}
ALPHA = 0.05
CLUSTER = "mondo_id"
PHASES = {"mode_not_causal": (2, 3), "s4": (2, 2)}


def frame_with_outcomes(d: pd.DataFrame, u: pd.DataFrame, cuts: tuple[float, float],
                        conc_cuts: tuple[float, float]) -> pd.DataFrame:
    """The unit's US Phase 2 and 3 entries with their traced stops, failure modes and the profile joined, the
    tertiles by the fixed cutpoints, `genetics` (g3 > 0), `top`, and the sensitivity columns."""
    from report.polygenicity import tertile  # noqa: PLC0415
    e = d[d.us_site & d.phase.between(2, 3) & d.mondo_id.isin(u.mondo_id)].copy()
    e = with_modes(traced_stops(e))
    cols = ["mondo_id", "score", "log10_score", "concentration", "mendelian_form", "n_eff", "ldsc_status", "scored",
            "k_source", "h2_z", "h2_liab", "loci", "loci_catalog", "gwas_genes_per_10k_cases", "me_log10", "label"]
    e = e.merge(u[[c for c in cols if c in u.columns]], on="mondo_id", how="left")
    e["tertile"] = tertile(e.score.where(e.scored.fillna(False)), cuts)
    e["conc_tertile"] = tertile(e.concentration.where(e.scored.fillna(False)), conc_cuts)
    e["genetics"] = (e.g3 > 0).astype(float)
    e["log10_n_eff"] = np.log10(e.n_eff.astype(float))
    e["mendelian_level"] = np.where(e.mendelian_form.fillna(False).astype(bool), "yes", "no")
    e["target"] = e.target.fillna("entity " + e.moiety_id.astype(str))
    h = engaged_modes()
    e["not_causal_basis"] = pd.Series(keys(e), index=e.index).map(dict(zip(h.key, h.not_causal_basis)))
    return e


def cohort(d: pd.DataFrame, outcome: str, levels: tuple[str, str] = ("top", "bottom"), col: str = "tertile") -> pd.DataFrame:
    """Entries with a resolved target at the exposed and reference levels of `col`, in the outcome's phases, with
    `top` = 1 at the exposed level."""
    lo, hi = PHASES[outcome]
    c = d[d.has_target_disease.fillna(False).astype(bool) & d[col].isin(levels) & d.phase.between(lo, hi)].copy()
    c["top"] = (c[col] == levels[0]).astype(float)
    return c


def one_sided(t: dict, side: int) -> dict:
    if not t or "rr" not in t:
        return t
    on_side = side * np.log(t["rr"]) > 0
    return {**t, "p_one_sided": t["p"] / 2 if on_side else 1 - t["p"] / 2}


def group_counts(c: pd.DataFrame, outcome: str) -> dict:
    return {"exposed_entries": int((c.top > 0).sum()), "exposed_events": int(c.loc[c.top > 0, outcome].sum()),
            "reference_entries": int((c.top <= 0).sum()), "reference_events": int(c.loc[c.top <= 0, outcome].sum())}


def p1(d: pd.DataFrame, covars=COVARS_23, extra: tuple[str, ...] = (), levels=("top", "bottom"), col="tertile") -> dict:
    c = cohort(d, "mode_not_causal", levels, col)
    return {**estimate(c, "mode_not_causal", "top", covars + extra, cluster=CLUSTER), **group_counts(c, "mode_not_causal")}


def p2(d: pd.DataFrame, covars=COVARS, extra: tuple[str, ...] = (), levels=("top", "bottom"), col="tertile") -> dict:
    c = cohort(d, "s4", levels, col)
    return {**estimate(c, "s4", "top", covars + extra, cluster=CLUSTER), **group_counts(c, "s4")}


def p3(d: pd.DataFrame, covars=COVARS, extra: tuple[str, ...] = (), levels=("top", "bottom"), col="tertile") -> dict:
    """The interaction of the top tertile with genetic support on S4, and the genetics ratio within each level."""
    c = cohort(d, "s4", levels, col).assign(top_x_genetics=lambda x: x.top * x.genetics)
    cells = {}
    for name, mask in (("top_genetics", (c.top > 0) & (c.genetics > 0)), ("top_no_genetics", (c.top > 0) & (c.genetics <= 0)),
                       ("bottom_genetics", (c.top <= 0) & (c.genetics > 0)), ("bottom_no_genetics", (c.top <= 0) & (c.genetics <= 0))):
        cells[name] = {"entries": int(mask.sum()), "targets": int(c.target[mask].nunique()), "events": int(c.s4[mask].sum())}
    ok = all(v["targets"] >= MIN_TARGETS and v["events"] >= 1 for v in cells.values())
    out = {"cells": cells, "estimable": ok, **group_counts(c, "s4"), "n": int(len(c)), "n_targets": int(c.target.nunique())}
    if not ok:
        return out
    fit = poisson(c, "s4", ["top", "genetics", "top_x_genetics"], covars + extra, cluster=CLUSTER)
    within = {lev: estimate(c[c.top == (1.0 if lev == "top" else 0.0)], "s4", "genetics", covars + extra, cluster=CLUSTER)
              for lev in ("top", "bottom")}
    return {**out, **coef(fit, "top_x_genetics"), "genetics_within": within}


def registered(d: pd.DataFrame, extra: tuple[str, ...] = (), col: str = "tertile",
               levels: tuple[str, str] = ("top", "bottom")) -> dict:
    """P1-P3 with one-sided p in the registered direction, Holm across the three, and the reading."""
    tests = {"P1": p1(d, extra=extra, levels=levels, col=col), "P2": p2(d, extra=extra, levels=levels, col=col),
             "P3": p3(d, extra=extra, levels=levels, col=col)}
    for k, t in tests.items():
        tests[k] = one_sided(t, SIDES[k])
    adj = holm({k: t.get("p_one_sided") for k, t in tests.items()})
    for k, t in tests.items():
        if "rr" in t:
            t["p_holm"] = adj[k]
            t["reading"] = ("confirmed" if adj[k] < ALPHA else
                            "same direction, not confirmed" if SIDES[k] * np.log(t["rr"]) > 0 else "not replicated")
        else:
            t["reading"] = "not estimable"
    return tests


def continuous(d: pd.DataFrame, covars_by_test=None) -> dict:
    """P1-P3 with log10 score in place of the tertile term, on scored entries; P3's product term."""
    s = d[d.scored.fillna(False).astype(bool) & d.log10_score.notna()]
    tests = {}
    for k, outcome in OUTCOMES.items():
        covars = COVARS_23 if k == "P1" else COVARS
        lo, hi = PHASES[outcome]
        c = s[s.has_target_disease.fillna(False).astype(bool) & s.phase.between(lo, hi)].copy()
        desc = describe(c.assign(any=1.0), outcome, "any")
        if k == "P3":
            c = c.assign(log10_score_x_genetics=c.log10_score * c.genetics)
            t = coef(poisson(c, outcome, ["log10_score", "genetics", "log10_score_x_genetics"], covars, cluster=CLUSTER),
                     "log10_score_x_genetics") if estimable(desc) else {}
        else:
            t = coef(poisson(c, outcome, ["log10_score"], covars, cluster=CLUSTER), "log10_score") if estimable(desc) else {}
        tests[k] = {**one_sided(t, SIDES[k]), "n": int(len(c)), "n_targets": int(c.target.nunique()), "per": "log10 unit"}
    adj = holm({k: t.get("p_one_sided") for k, t in tests.items()})
    for k, t in tests.items():
        if "rr" in t:
            t["p_holm"] = adj[k]
            t["reading"] = ("confirmed" if adj[k] < ALPHA else
                            "same direction, not confirmed" if SIDES[k] * np.log(t["rr"]) > 0 else "not replicated")
        else:
            t["reading"] = "not estimable"
    return tests


def sensitivities(d: pd.DataFrame) -> dict:
    """Melissa's four registered sensitivities, each P1-P3 with Holm within itself."""
    if "mendelian_level" not in d.columns:
        d = d.assign(mendelian_level=np.where(d.mendelian_form.fillna(False).astype(bool), "yes", "no"))
    return {"gwas_size": registered(d, extra=("log10_n_eff",)),
            "continuous": continuous(d),
            "concentration": registered(d, col="conc_tertile", levels=("bottom", "top")),
            "mendelian_form": registered(d, col="mendelian_level", levels=("no", "yes"))}


def rates(d: pd.DataFrame, col: str) -> pd.DataFrame:
    """Per level of `col`: diseases, entries, not-causal stops and rate (Phase 2 and 3), S4 and rate (Phase 2)."""
    rows = []
    for level, g in d.groupby(col):
        p2e = g[g.phase.eq(2)]
        rows.append({col: level, "diseases": int(g.mondo_id.nunique()), "entries": int(len(g)),
                     "phase2_entries": int(len(p2e)), "not_causal": int(g.mode_not_causal.sum()),
                     "not_causal_rate": float(g.mode_not_causal.mean()) if len(g) else None,
                     "s4": int(p2e.s4.sum()), "s4_rate": float(p2e.s4.mean()) if len(p2e) else None,
                     "without_target": int((~g.has_target_disease.fillna(False).astype(bool)).sum())})
    return pd.DataFrame(rows)


def descriptives(d: pd.DataFrame, u: pd.DataFrame) -> dict:
    from report.polygenicity import cutpoints, tertile  # noqa: PLC0415
    out = {"by_tertile": rates(d, "tertile").to_dict("records"),
           "area": d.groupby(["tertile", "bio_area"]).agg(entries=("mondo_id", "size"), diseases=("mondo_id", "nunique"))
                    .reset_index().to_dict("records"),
           "no_score": {"diseases": int((~u.scored).sum()), "entries": int(u.entries[~u.scored].sum()),
                        "status": u[~u.scored].groupby("ldsc_status").agg(diseases=("entries", "size"), entries=("entries", "sum"))
                                   .reset_index().to_dict("records"),
                        **{k: v for k, v in rates(d[d.tertile.eq("none")].assign(g="no_score"), "g").to_dict("records")[0].items()
                           if k != "g"}} if (d.tertile == "none").any() else {"diseases": 0, "entries": 0}}
    # inclusion at z >= 2
    z2 = u[u.h2_z.notna() & (u.h2_z >= 2) & u.h2_liab.notna() & (u.h2_liab > 0) & (u.h2_liab <= 1)
           & u.ldsc_status.isin(["ok", "low_z"])]
    if len(z2) >= 6:
        sc = (z2.loci / z2.h2_liab)
        dd = d.merge(z2[["mondo_id"]].assign(score_z2=sc.values), on="mondo_id", how="left")
        dd["tertile_z2"] = tertile(dd.score_z2, cutpoints(sc))
        out["z_at_least_2"] = {"diseases": int(len(z2)), "tests": registered(dd, col="tertile_z2")}
    # GBD-K diseases only
    gk = d[d.k_source.eq("gbd")]
    out["gbd_k_only"] = {"diseases": int(gk.mondo_id.nunique()), "tests": registered(gk)} if gk.mondo_id.nunique() >= 6 else {}
    # the Catalog's all-study loci in place of the clumped loci
    cat = u[u.scored & u.loci_catalog.notna() & u.h2_liab.notna()]
    if len(cat) >= 6:
        sc = cat.loci_catalog / cat.h2_liab
        dd = d.merge(cat[["mondo_id"]].assign(score_cat=sc.values), on="mondo_id", how="left")
        dd["tertile_cat"] = tertile(dd.score_cat, cutpoints(sc))
        out["catalog_loci"] = {"diseases": int(len(cat)), "tests": registered(dd, col="tertile_cat")}
    # gene counts as the exposure
    gc = u[u.gwas_genes_per_10k_cases.notna()]
    if len(gc) >= 6:
        dd = d.copy()
        dd["tertile_genes"] = tertile(dd.gwas_genes_per_10k_cases, cutpoints(gc.gwas_genes_per_10k_cases))
        out["gwas_genes"] = {"diseases": int(len(gc)), "tests": registered(dd, col="tertile_genes")}
    # O'Connor
    me = u[u.me_log10.notna() & u.score.notna()]
    if len(me) >= 3:
        rho, p = spearmanr(me.score, me.me_log10)
        out["oconnor"] = {"n": int(len(me)), "spearman": float(rho), "p": float(p),
                          "diseases": me[["mondo_id", "label", "score", "me_log10"]].to_dict("records")}
    # the basis of P1's not-causal calls
    c = cohort(d, "mode_not_causal")
    nc = c[c.mode_not_causal > 0]
    out["not_causal_basis"] = nc.groupby(["tertile", "not_causal_basis"]).size().reset_index(name="events").to_dict("records")
    return out


def results_table(families: dict) -> pd.DataFrame:
    rows = []
    for family, tests in families.items():
        for k, t in tests.items():
            rows.append({"family": family, "test": k, **{c: t.get(c) for c in ("rr", "lo", "hi", "p", "p_one_sided", "p_holm", "reading",
                                                                            "exposed_entries", "exposed_events",
                                                                            "reference_entries", "reference_events")}})
    return pd.DataFrame(rows, columns=["family", "test", "rr", "lo", "hi", "p", "p_one_sided", "p_holm", "reading",
                                       "exposed_entries", "exposed_events", "reference_entries", "reference_events"])


def run_tests(d: pd.DataFrame, u: pd.DataFrame) -> dict:
    from report import polygenicity as pg  # noqa: PLC0415
    e = frame_with_outcomes(d, u, pg.CUTPOINTS, pg.CONC_CUTPOINTS)
    primary = registered(e)
    sens = sensitivities(e)
    desc = descriptives(e, u)
    results_table({"primary": primary}).to_csv(os.path.join(pg.OUT, "tests.csv"), index=False)
    results_table(sens).to_csv(os.path.join(pg.OUT, "sensitivities.csv"), index=False)
    rates(e, "tertile").to_csv(os.path.join(pg.OUT, "tertiles.csv"), index=False)
    pd.DataFrame(desc["area"]).to_csv(os.path.join(pg.OUT, "area.csv"), index=False)
    pd.DataFrame(desc["no_score"].get("status", [])).to_csv(os.path.join(pg.OUT, "no_score.csv"), index=False)
    n = {"cohort": {"entries": int(len(e)), "phase2": int(e.phase.eq(2).sum()), "with_target": int(e.has_target_disease.fillna(False).sum()),
                    "scored_entries": int(e.scored.fillna(False).sum()), "diseases": int(e.mondo_id.nunique())},
         "cutpoints": {"score": pg.CUTPOINTS, "concentration": pg.CONC_CUTPOINTS}, "build_id": pg.BUILD_ID,
         "comment_id": pg.COMMENT_ID, "covariates": {"P1": list(COVARS_23), "P2": list(COVARS), "P3": list(COVARS)},
         "tests": primary, "sensitivities": sens, "descriptives": desc}
    with open(os.path.join(pg.OUT, "numbers.json"), "w") as f:
        json.dump(jsonable(n), f, indent=1)
    for k, t in primary.items():
        print(k, t.get("rr"), t.get("p_one_sided"), t.get("p_holm"), t.get("reading"))
    for fam, tests in sens.items():
        print(fam, {k: t.get("reading") for k, t in tests.items()})
    return n
