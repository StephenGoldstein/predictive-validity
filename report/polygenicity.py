"""Polygenicity and failure (issue #147): do targets fail more often in more polygenic diseases.

Registered by Melissa on 2026-11-19 (#147, the registration comment), with her dated decisions of 2026-11-23: diseases
without a heritability estimate form a descriptive group outside the tertiles (IC_...qDb0FA); infectious diseases are
out of the unit (IC_...qFfZeg); heritability A, our LDSC on the largest European GWAS with public summary statistics,
loci and heritability from the same GWAS; K from GBD with the Pan-UKB case fraction as fallback; four registered
sensitivities (issuecomment-7123250021). The definitions used, the unit counts, the build id and the tertile
cutpoints are fixed in a dated comment before any outcome is joined (report/audits/polygenicity_definitions.md); the
cutpoints are the constants CUTPOINTS and CONC_CUTPOINTS below, committed with that comment.

Unit (unit): non-oncology, non-infectious diseases with a US-developed Phase 2 entry in 2008-2019, each disease one
observation. Scored diseases (an estimate: preclin.disease_polygenicity.ldsc_status ok) form the tertiles of the
score (loci per unit of liability-scale SNP heritability) and of the heritability concentration; the rest are the
no-score group, reported, never tested.

Registered tests (registered, second half of this module, run only once CUTPOINTS is set): P1, more "target not
causal" stops per entry in the top tertile than the bottom (Phase 2 and 3 entries, COVARS_23); P2, lower Phase 2
advance (S4, COVARS); P3, a smaller S4 risk ratio of human genetic support in the top tertile (interaction, COVARS).
Entries with a resolved primary target; modified Poisson clustered by disease; one-sided p in the registered
direction, Holm across the three, as the #140 confirmation reads its tests. Registered sensitivities, each repeating
P1-P3 with Holm within: GWAS size (log10 effective N as a covariate), the score as a continuous log10 term, the
concentration tertiles (low concentration in place of high polygenicity), and a Mendelian form over all unit diseases
(no Mendelian form in place of the top tertile). Descriptive, not registered: tertile by disease area, the no-score
group, entries without a target, and the checks named in the appendix.

Writes report/output/polygenicity/: unit.csv, membership.csv, then tests.csv, sensitivities.csv, tertiles.csv,
area.csv, no_score.csv, daly.csv and numbers.json.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from pipeline.stops import keys
from report import frame
from report.engaged_modes import engaged_modes
from report.evidence_analysis import traced_stops
from report.evidence_lines import build, us
from report.failure_modes import with_modes
from report.section3 import jsonable

OUT = os.path.join(frame.REPO, "report", "output", "polygenicity")
AUDITS = os.path.join(frame.REPO, "report", "audits")
EXCLUDED_AREAS = ("oncology", "infectious_disease")
LEVELS = ("bottom", "middle", "top")
# Fixed before any outcome was joined (the dated comment on #147 named in polygenicity_definitions.md); None until then.
CUTPOINTS: tuple[float, float] | None = None       # score: 33.3rd and 66.7th percentiles over scored unit diseases
CONC_CUTPOINTS: tuple[float, float] | None = None  # heritability concentration, the same way
BUILD_ID: str | None = None                        # preclin.dev_program_build.build_id the unit was counted on
COMMENT_ID: str | None = None                      # the dated comment


def profile() -> pd.DataFrame:
    """preclin.disease_polygenicity, one row per candidate disease."""
    return frame.query("SELECT * FROM preclin.disease_polygenicity ORDER BY mondo_id")


def unit(d: pd.DataFrame, prof: pd.DataFrame) -> pd.DataFrame:
    """One row per disease in the unit: its Phase 2 entries, the profile's status, `scored`; a disease the profile
    lacks is `not_profiled` (the pre-comment check fails on it)."""
    p2 = d[d.us_site & d.phase.eq(2) & d.mondo_id.notna() & ~d.bio_area.isin(EXCLUDED_AREAS)]
    u = p2.groupby("mondo_id").size().reset_index(name="entries")
    u = u.merge(prof.drop(columns=["entries"], errors="ignore"), on="mondo_id", how="left")
    u["ldsc_status"] = u.ldsc_status.fillna("not_profiled")
    u["scored"] = u.ldsc_status.eq("ok")
    return u.sort_values("mondo_id").reset_index(drop=True)


def cutpoints(score: pd.Series) -> tuple[float, float]:
    s = score.dropna()
    return float(np.percentile(s, 100 / 3)), float(np.percentile(s, 200 / 3))


def tertile(score: pd.Series, cuts: tuple[float, float]) -> pd.Series:
    """bottom, middle, top by the two cutpoints (a value on a cut goes down); none where the score is missing."""
    out = pd.Series("none", index=score.index, dtype=object)
    out[score <= cuts[0]] = "bottom"
    out[(score > cuts[0]) & (score <= cuts[1])] = "middle"
    out[score > cuts[1]] = "top"
    return out


def membership(u: pd.DataFrame, cuts: tuple[float, float] | None, conc_cuts: tuple[float, float] | None) -> pd.DataFrame:
    """The unit with its tertiles: the audit file and the comment's membership list."""
    m = u.copy()
    m["tertile"] = tertile(m.score, cuts) if cuts else "none"
    m["conc_tertile"] = tertile(m.concentration, conc_cuts) if conc_cuts else "none"
    cols = ["mondo_id", "label", "bio_area", "entries", "ldsc_status", "gcst", "trait_level", "k_source", "score",
            "tertile", "concentration", "conc_tertile", "mendelian_form"]
    return m[[c for c in cols if c in m.columns]]


def current_build() -> dict:
    b = frame.query("SELECT build_id::text, built_at, snapshot FROM preclin.dev_program_build ORDER BY built_at DESC LIMIT 1")
    return {"build_id": b.build_id.iloc[0], "built_at": str(b.built_at.iloc[0]), "snapshot": str(b.snapshot.iloc[0])}


def pre_comment_check(d: pd.DataFrame, u: pd.DataFrame) -> dict:
    """What the definitions comment needs and what must be clean before it: the build, the unit's counts, diseases
    not profiled, and engaged-without-benefit stops among the unit's Phase 2 and 3 entries with no target-history
    row. Reads no outcome by tertile."""
    scored = u[u.scored]
    out = {**current_build(), "diseases": int(len(u)), "entries": int(u.entries.sum()),
           "scored_diseases": int(len(scored)), "scored_entries": int(scored.entries.sum()),
           "not_profiled": sorted(u.mondo_id[u.ldsc_status.eq("not_profiled")])}
    e = d[d.us_site & d.phase.between(2, 3) & d.mondo_id.isin(u.mondo_id)]
    e = with_modes(traced_stops(e))
    history = set(engaged_modes().key)
    k = pd.Series(keys(e), index=e.index)
    engaged = e[e.stop_no_benefit > 0]
    out["engaged_stops"] = int(len(engaged))
    out["engaged_stops_without_history"] = int((~k[engaged.index].isin(history)).sum())
    return out


WHO_FILE = os.path.join(frame.REPO, "analyses", "model-system-atlas", "lists", "who_ghe_2021_global_dalys.csv")
WHO_MAP = os.path.join(frame.REPO, "pipeline", "data", "who_ghe_cause_mondo.csv")


def who_map() -> pd.DataFrame:
    """The hand map of WHO causes to MONDO (code, cause, mondo_id, relation, note), extending burden.py's."""
    return pd.read_csv(WHO_MAP, dtype={"code": int, "mondo_id": str})


def leaves(who: pd.DataFrame) -> pd.DataFrame:
    """The causes with no children: the file lists each cause's children right after it, one level deeper."""
    depth = who.depth.tolist()
    keep = [i == len(depth) - 1 or depth[i + 1] <= depth[i] for i in range(len(depth))]
    return who[keep]


def daly_share(who: pd.DataFrame, m: pd.DataFrame, flags: pd.DataFrame) -> dict:
    """Shares of DALYs (all causes; non-communicable diseases) in leaf causes mapped to a disease with a Mendelian
    form, with a drugged large-effect gene, or either; and the largest flagged causes."""
    f = flags.set_index("mondo_id")
    out = {}
    cause_flags = {}
    for code, g in m[m.relation.isin(["exact", "broader"])].groupby("code"):
        ids = [i for i in g.mondo_id if i in f.index]
        cause_flags[code] = {"mendelian_form": bool(f.mendelian_form[ids].any()) if ids else False,
                             "drugged_large_effect": bool(f.drugged_large_effect[ids].any()) if ids else False}
    lf = leaves(who)
    for name, total in (("all", float(who.dalys_2021[who.code == 0].iloc[0])),
                        ("ncd", float(who.dalys_2021[who.cause.str.startswith("Noncommunicable")].sum()))):
        if name == "ncd":
            start = who.index[who.cause.str.startswith("Noncommunicable")]
            if not len(start):
                continue
            i = start[0]
            after = who.loc[i + 1:]
            stop = after.index[after.depth <= who.depth[i]]
            sub = who.loc[i + 1:(stop[0] - 1) if len(stop) else None]
            rows = lf[lf.code.isin(sub.code)]
        else:
            rows = lf
        shares = {}
        for flag in ("mendelian_form", "drugged_large_effect", "either"):
            hit = rows.code.map(lambda c: (cause_flags.get(c, {}).get("mendelian_form", False) or cause_flags.get(c, {}).get("drugged_large_effect", False))
                                if flag == "either" else cause_flags.get(c, {}).get(flag, False))
            shares[flag] = float(rows.dalys_2021[hit].sum() / total) if total else None
        shares["total"] = total
        shares["largest"] = rows[rows.code.map(lambda c: cause_flags.get(c, {}).get("mendelian_form", False)
                                               or cause_flags.get(c, {}).get("drugged_large_effect", False))] \
            .sort_values("dalys_2021", ascending=False).head(10)[["code", "cause", "dalys_2021"]].to_dict("records")
        out[name] = shares
    return out


def drugged_large_effect(prof: pd.DataFrame) -> pd.DataFrame:
    """Per profiled disease: whether a gene with rare or Mendelian evidence for it (the same match as the profile's
    mendelian genes) is a primary target of a programme in the disease that entered Phase 2 or 3, or was approved on
    that target for the disease."""
    from pipeline.polygenicity import MENDELIAN_SOURCES, descendants_table, evidence_frame, members  # noqa: PLC0415
    from report.disease_match import Mondo, mondo_tables  # noqa: PLC0415
    from report.section3 import FAMILY  # noqa: PLC0415
    ev = evidence_frame()
    ev = ev[ev.datasource.map(FAMILY).isin(["burden", "mendelian"])]
    td = frame.query("SELECT ensembl_id, mondo_id FROM preclin.target_disease")
    approved = frame.query("""SELECT DISTINCT t.ensembl_id, p.mondo_id FROM preclin.dev_program p
                              JOIN preclin.moiety_target t USING (moiety_id) WHERE p.mondo_id IS NOT NULL
                              AND EXISTS (SELECT 1 FROM preclin.program_approval a WHERE a.program_id = p.program_id)""") \
        if frame.query("SELECT to_regclass('preclin.program_approval') IS NOT NULL AS ok").ok.iloc[0] else td.iloc[0:0]
    drugged = td.groupby("mondo_id").ensembl_id.apply(set).to_dict()
    for r in approved.itertuples():
        drugged.setdefault(r.mondo_id, set()).add(r.ensembl_id)
    desc, mondo = descendants_table(mondo_tables()[0]), Mondo()
    rows = []
    for mondo_id in prof.mondo_id:
        mem = members(mondo_id, desc)
        close = {d for d in ev.mondo_id.unique() if d not in mem and mondo.matches(d, mondo_id)[1]}
        genes = set(ev.target_ensembl[ev.mondo_id.isin(mem | close)])
        rows.append({"mondo_id": mondo_id, "drugged_large_effect": bool(genes & drugged.get(mondo_id, set()))})
    return pd.DataFrame(rows)


def daly(prof: pd.DataFrame) -> dict:
    who = pd.read_csv(WHO_FILE)
    flags = prof[["mondo_id", "mendelian_form"]].merge(drugged_large_effect(prof), on="mondo_id", how="left")
    flags["mendelian_form"] = flags.mendelian_form.fillna(False).astype(bool)
    flags["drugged_large_effect"] = flags.drugged_large_effect.fillna(False).astype(bool)
    s = daly_share(who, who_map(), flags)
    pd.DataFrame([{"denominator": k, **{kk: vv for kk, vv in v.items() if kk != "largest"}} for k, v in s.items()]) \
        .to_csv(os.path.join(OUT, "daly.csv"), index=False)
    return {"shares": s, "flags": {"diseases": int(len(flags)), "mendelian_form": int(flags.mendelian_form.sum()),
                                   "drugged_large_effect": int(flags.drugged_large_effect.sum())}}


def main() -> None:
    os.makedirs(OUT, exist_ok=True)
    d = build()
    prof = profile()
    u = unit(d, prof)
    u.to_csv(os.path.join(OUT, "unit.csv"), index=False)
    membership(u, CUTPOINTS, CONC_CUTPOINTS).to_csv(os.path.join(OUT, "membership.csv"), index=False)
    check = pre_comment_check(d, u)
    print(json.dumps(jsonable(check), indent=1))
    with open(os.path.join(OUT, "check.json"), "w") as f:
        json.dump(jsonable(check), f, indent=1)
    if CUTPOINTS is None:
        print("CUTPOINTS not set: the tertiles are not fixed; no outcome is joined.")
        return
    from report.polygenicity_tests import run_tests  # noqa: PLC0415  (the second half, after the dated comment)
    run_tests(d, u, extra={"daly": daly(prof)})


if __name__ == "__main__":
    main()
