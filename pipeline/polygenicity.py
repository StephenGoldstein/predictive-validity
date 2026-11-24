"""The per-disease genetic-architecture profile of issue #147 (preclin.disease_polygenicity).

    python -m pipeline.polygenicity load-gbd     # pipeline/data/gbd_2023_prevalence_2021.csv -> preclin.gbd_prevalence
    python -m pipeline.polygenicity choose       # candidate diseases and their LDSC study -> preclin.disease_polygenicity
    python -m pipeline.polygenicity build        # every column: the runs, K, the score, Open Targets, O'Connor, Pan-UKB
    python -m pipeline.polygenicity counts       # scored diseases and entries against the candidates and the 531/2,275 base
    python -m pipeline.polygenicity check-unit   # every disease in today's unit has a row (exit 1 otherwise)

Candidates: non-oncology, non-infectious diseases with a US-developed Phase 2 entry in 2008-2019 (the unit of #147,
counted on the programme build of the day; a disease the unit gains later is added by rerunning choose, ldsc and
build). For each, the LDSC study (choose_study): among the Catalog studies of the disease and its MONDO
descendants, those with public harmonised summary statistics, parsed cases and controls, a European initial sample
(at least 90%), a genome-wide platform (200,000 SNPs or more, or unstated) and a single mapped trait; the disease
itself before a descendant; the largest effective sample size 4 / (1/cases + 1/controls). The status names the
first filter that emptied the set.

K, the population prevalence for the liability scale (k_for): the GBD cause the disease maps to by hand
(pipeline/data/gbd_cause_mondo.csv, exact or broader), rate per 100,000 over 100,000, when under 0.5; else the
Pan-UK Biobank European case fraction of the disease's ICD-10 code (its own or a descendant's, from MONDO's
cross-references), the code with most cases, at least 100; else none, and the disease has no score. Liability
h2 = h2_obs * K^2 (1-K)^2 / (P (1-P) phi(t)^2), t = Phi^-1(1-K), P the study's case fraction (Lee et al. 2011).
An estimate exists when the run is ok on a genome-wide set of SNPs, h2_z = h2_obs / se >= 4 (Pan-UKB's
significant_z threshold) and 0 < h2_liab <= 1; the score is loci / h2_liab.

Secondary measures: distinct genes with GWAS evidence and with rare or Mendelian evidence in preclin.ot_evidence
(Open Targets 26.09, the targets we hold) for the disease, a subtype or a disease at Lin >= 0.8
(report.disease_match.Mondo, the evidence lines' rule), the GWAS count per 10,000 cases of the largest GWAS;
whether a Mendelian form exists (a gene from Orphanet, ClinGen, Genomics England or Gene2Phenotype by the same
match: Melissa's registered sensitivity 4); O'Connor 2019's log10 effective number of associated SNPs (Table 1,
pipeline/data/oconnor2019_me.csv); Pan-UKB's liability h2 and z for the ICD-10 code (a check on our LDSC).

GBD file: Global Burden of Disease Study 2023 (GBD 2023) Results, prevalence, global, both sexes, all ages, 2021,
exported from https://vizhub.healthdata.org/gbd-results/ by Melissa and received 2026-11-23 (pipeline/data/
gbd_citation.txt). IHME publishes it under the free-of-charge non-commercial user agreement.
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from psycopg2.extras import execute_values
from scipy.stats import norm, spearmanr

from pipeline import llm
from pipeline.gwas_catalog import RELEASE, loci_catalog, read_associations, trait_lookup
from report import frame, section3
from report.disease_match import Mondo, mondo_tables
from report.evidence_lines import us
from report.section3 import FAMILY
from report.windows import SECTION1

DATA = os.path.join(llm.REPO, "pipeline", "data")
GBD_FILE = os.path.join(DATA, "gbd_2023_prevalence_2021.csv")
GBD_CITATION = os.path.join(DATA, "gbd_citation.txt")
GBD_MAP = os.path.join(DATA, "gbd_cause_mondo.csv")
ME_FILE = os.path.join(DATA, "oconnor2019_me.csv")
PANUKB = os.path.join(llm.REPO, ".cache", "panukb")
EXCLUDED_AREAS = ("oncology", "infectious_disease")
BASE = {"diseases": 531, "entries": 2275}  # Melissa's unit of 2026-11-23, before the #185 rebuild
EUR_MIN = 0.9
MIN_PLATFORM_SNPS = 200_000
Z_MIN = 4.0
K_MAX = 0.5
MIN_UKB_CASES = 100
MENDELIAN_SOURCES = ("orphanet", "clingen", "genomics_england", "gene2phenotype")
STATUS_ORDER = ("no_study", "no_sumstats", "no_case_control", "not_european", "targeted_array")


def candidates() -> pd.DataFrame:
    """The unit's diseases on the current programme build: mondo_id, label, bio_area, entries."""
    progs, trials, snap = frame.programs(), frame.program_trials(), frame.snapshot()
    e = section3.entries(progs, trials, snap, SECTION1).merge(progs[["program_id", "mondo_id"]], on="program_id", how="left")
    p2 = us(e)
    p2 = p2[p2.phase.eq(2) & p2.mondo_id.notna() & ~p2.bio_area.isin(EXCLUDED_AREAS)]
    out = p2.groupby(["mondo_id", "bio_area"]).size().reset_index(name="entries")
    labels = frame.query("SELECT term_id, label FROM preclin.ontology_term WHERE ontology = 'MONDO'")
    out["label"] = out.mondo_id.map(dict(zip(labels.term_id, labels.label)))
    return out[["mondo_id", "label", "bio_area", "entries"]].sort_values("mondo_id").reset_index(drop=True)


def descendants_table(ancestors: dict) -> dict[str, set[str]]:
    d: dict[str, set[str]] = defaultdict(set)
    for term, anc in ancestors.items():
        for a in anc:
            d[a].add(term)
    return d


def members(mondo_id: str, descendants: dict[str, set[str]]) -> set[str]:
    return {mondo_id} | set(descendants.get(mondo_id, ()))


def n_eff(cases, controls) -> float:
    return 4 / (1 / float(cases) + 1 / float(controls))


def choose_study(studies: pd.DataFrame, mondo_id: str, members: set[str]) -> tuple[str | None, str | None, str]:
    """(gcst, exact or descendant, status) for a disease; None and the reason when no study qualifies."""
    s = studies[studies.mondo_ids.map(lambda m: bool(set(m) & members))]
    if s.empty:
        return None, None, "no_study"
    filters = (("no_sumstats", lambda x: x.sumstats_path.notna()),
               ("no_case_control", lambda x: x.cases.notna() & x.controls.notna()),
               ("not_european", lambda x: x.eur_share.notna() & (pd.to_numeric(x.eur_share) >= EUR_MIN)),
               ("targeted_array", lambda x: x.platform_snps.isna()
                | (pd.to_numeric(x.platform_snps, errors="coerce") >= MIN_PLATFORM_SNPS)),
               ("no_study", lambda x: x.mapped_trait_uris.map(len) == 1))
    for status, keep in filters:
        s = s[keep(s).fillna(False).astype(bool)]
        if s.empty:
            return None, None, status
    exact = s[s.mondo_ids.map(lambda m: mondo_id in m)]
    pool, level = (exact, "exact") if len(exact) else (s, "descendant")
    best = pool.assign(n_eff=[n_eff(c, k) for c, k in zip(pool.cases, pool.controls)]).sort_values("n_eff", ascending=False)
    return best.gcst.iloc[0], level, "ok"


def liability_h2(h2_obs: float, p: float, k: float) -> float:
    """Observed-scale h2 of a case-control GWAS with case fraction p to the liability scale at prevalence k."""
    t = norm.isf(k)
    return h2_obs * k ** 2 * (1 - k) ** 2 / (p * (1 - p) * norm.pdf(t) ** 2)


def k_for(mondo_id: str, gbd_map: dict, gbd: dict, icd_codes: set[str], panukb: pd.DataFrame) -> tuple[float | None, str]:
    """(K, source): the GBD cause's rate when the disease maps exactly or to a broader cause and K < 0.5; else the
    Pan-UKB European case fraction of the ICD-10 code with most cases (at least MIN_UKB_CASES); else none."""
    relation, cause = gbd_map.get(mondo_id, ("none", None))
    if relation in ("exact", "broader") and cause in gbd:
        k = gbd[cause] / 1e5
        if 0 < k < K_MAX:
            return k, "gbd"
    rows = panukb[panukb.phenocode.isin(icd_codes)]
    if len(rows):
        r = rows.sort_values("n_cases_EUR", ascending=False).iloc[0]
        if r.n_cases_EUR >= MIN_UKB_CASES:
            return float(r.n_cases_EUR / (r.n_cases_EUR + r.n_controls_EUR)), "panukb"
    return None, "none"


def gene_counts(mondo_id: str, mondo, members: set[str], evidence: pd.DataFrame) -> tuple[int, int, bool]:
    """(GWAS genes, rare or Mendelian genes, Mendelian form) from evidence rows (target_ensembl, mondo_id,
    datasource) for the disease, its members or a disease at Lin >= 0.8."""
    close = {d for d in evidence.mondo_id.unique() if d not in members and mondo.matches(d, mondo_id)[1]}
    ev = evidence[evidence.mondo_id.isin(members | close)]
    family = ev.datasource.map(FAMILY)
    gwas = int(ev.target_ensembl[family == "gwas"].nunique())
    mendelian = int(ev.target_ensembl[family.isin(("burden", "mendelian"))].nunique())
    form = bool(ev.datasource.isin(MENDELIAN_SOURCES).any())
    return gwas, mendelian, form


def gbd_map() -> dict[str, tuple[str, int]]:
    """Disease -> (relation, GBD cause id): exact before broader, then the larger rate."""
    m = pd.read_csv(GBD_MAP, dtype=str).dropna(subset=["mondo_id"])
    m = m[m.relation.isin(["exact", "broader"])]
    rates = gbd_rates()
    m = m.assign(rank=m.relation.map({"exact": 0, "broader": 1}), rate=m.cause_id.astype(int).map(rates))
    m = m.sort_values(["mondo_id", "rank", "rate"], ascending=[True, True, False]).drop_duplicates("mondo_id")
    return {r.mondo_id: (r.relation, int(r.cause_id)) for r in m.itertuples()}


def gbd_rates() -> dict[int, float]:
    g = pd.read_csv(GBD_FILE)
    g = g[(g.measure_name == "Prevalence") & (g.metric_name == "Rate") & (g.year == 2021)]
    return dict(zip(g.cause_id.astype(int), g.val.astype(float)))


def panukb_manifest() -> pd.DataFrame:
    ph = pd.read_csv(os.path.join(PANUKB, "phenotype_manifest.tsv.bgz"), sep="\t", compression="gzip", low_memory=False)
    ph = ph[ph.trait_type.eq("icd10") & ph.n_cases_EUR.notna() & ph.n_controls_EUR.notna()]
    return ph[["phenocode", "n_cases_EUR", "n_controls_EUR"]].assign(phenocode=lambda x: x.phenocode.astype(str))


def panukb_h2() -> pd.DataFrame:
    h = pd.read_csv(os.path.join(PANUKB, "h2_manifest.tsv.bgz"), sep="\t", compression="gzip", low_memory=False)
    h = h[h.trait_type.eq("icd10") & h.ancestry.eq("EUR") & h["qcflags.pass_all"].astype(str).str.lower().eq("true")]
    return pd.DataFrame({"phenocode": h.phenocode.astype(str), "h2_liab": h["estimates.final.h2_liability"],
                         "z": h["estimates.final.h2_z"]})


def icd_codes(members: set[str], xrefs: pd.DataFrame) -> set[str]:
    codes = xrefs[xrefs.term_id.isin(members)].code
    return {c for code in codes for c in (code, code[:3])}


def evidence_frame() -> pd.DataFrame:
    return frame.query("""SELECT DISTINCT e.target_ensembl, m.mondo AS mondo_id, e.datasource
                          FROM preclin.ot_evidence e JOIN preclin.ot_disease d USING (disease_id),
                               unnest(d.mondo_ids) AS m(mondo)""")


def load_gbd(_args) -> None:
    g = pd.read_csv(GBD_FILE)
    g = g[(g.measure_name == "Prevalence") & (g.metric_name == "Rate")]
    citation = open(GBD_CITATION).read().strip()
    with llm.connect() as conn, conn.cursor() as cur:
        rid = llm.register_release(cur, "gbd_prevalence", "GBD 2023, year 2021", "https://vizhub.healthdata.org/gbd-results/",
                                   llm.sha256_file(GBD_FILE), {"received": "2026-11-23", "from": "Melissa's export",
                                                              "terms": "IHME free-of-charge non-commercial user agreement",
                                                              "citation": citation, "causes": int(len(g))})
        cur.execute("DELETE FROM preclin.gbd_prevalence")
        execute_values(cur, "INSERT INTO preclin.gbd_prevalence (cause_id, cause_name, year, rate_per_100k, source_release_id) VALUES %s",
                       [(int(r.cause_id), r.cause_name, int(r.year), float(r.val), rid) for r in g.itertuples()])
        conn.commit()
    print(f"{len(g)} causes loaded")


def release_id(name: str) -> int:
    r = frame.query("SELECT source_release_id FROM preclin.evidence_source_release WHERE source_name = %s "
                    "ORDER BY retrieved_at DESC LIMIT 1", (name,))
    if r.empty:
        raise SystemExit(f"no release registered for {name}")
    return int(r.source_release_id.iloc[0])


def studies_frame() -> pd.DataFrame:
    return frame.query("""SELECT gcst, mondo_ids, mapped_trait_uris, sumstats_path, cases, controls, eur_share, platform_snps
                          FROM preclin.gwas_catalog_study WHERE mondo_ids <> '{}'""")


def choose(_args) -> None:
    cands = candidates()
    studies = studies_frame()
    desc = descendants_table(mondo_tables()[0])
    rid = release_id("gwas_catalog")
    rows = []
    for c in cands.itertuples():
        gcst, level, status = choose_study(studies, c.mondo_id, members(c.mondo_id, desc))
        rows.append((c.mondo_id, c.label, c.bio_area, int(c.entries), gcst, level, status, [rid]))
    with llm.connect() as conn, conn.cursor() as cur:
        execute_values(cur, """INSERT INTO preclin.disease_polygenicity
                                   (mondo_id, label, bio_area, entries, gcst, trait_level, ldsc_status, source_release_ids)
                               VALUES %s ON CONFLICT (mondo_id) DO UPDATE SET label = EXCLUDED.label, bio_area = EXCLUDED.bio_area,
                                   entries = EXCLUDED.entries, gcst = EXCLUDED.gcst, trait_level = EXCLUDED.trait_level,
                                   ldsc_status = EXCLUDED.ldsc_status, built_at = now()""", rows)
        conn.commit()
    chosen = pd.Series([r[6] for r in rows]).value_counts()
    print(f"{len(rows)} candidate diseases; study chosen for {int(chosen.get('ok', 0))}; " +
          ", ".join(f"{k} {v}" for k, v in chosen.items()))


def build(_args) -> None:
    prof = frame.query("SELECT mondo_id, label, bio_area, entries, gcst, trait_level, ldsc_status FROM preclin.disease_polygenicity")
    runs = frame.query("SELECT gcst, status, n, n_snps, h2_obs, h2_obs_se, intercept, loci, concentration FROM preclin.ldsc_run")
    runs = runs.set_index("gcst")
    studies = frame.query("SELECT gcst, cases, controls, mondo_ids FROM preclin.gwas_catalog_study WHERE mondo_ids <> '{}'")
    studies = studies.set_index("gcst")
    xrefs = frame.query("""SELECT term_id, replace(split_part(xref_id, ':', 2), '.', '') AS code FROM preclin.ontology_xref
                           WHERE ontology = 'MONDO' AND (xref_id LIKE 'ICD10CM:%%' OR xref_id LIKE 'ICD10WHO:%%')""")
    desc = descendants_table(mondo_tables()[0])
    mondo = Mondo()
    gmap, rates = gbd_map(), gbd_rates()
    ukb, ukb_h2 = panukb_manifest(), panukb_h2()
    me = pd.read_csv(ME_FILE, dtype=str).dropna(subset=["mondo_ids"])
    me_map = {m: float(r.log10_me) for r in me.itertuples() for m in r.mondo_ids.split(";")}
    evidence = evidence_frame()
    assoc, lookup = read_associations(), trait_lookup()
    releases = [release_id("gwas_catalog"), release_id("gbd_prevalence"), release_id("open_targets_platform_evidence")]
    out = []
    for d in prof.itertuples():
        mem = members(d.mondo_id, desc)
        row = {"mondo_id": d.mondo_id, "ldsc_status": d.ldsc_status}
        row["gwas_studies"], row["gws_associations"], row["loci_catalog"] = loci_catalog(assoc, mem, lookup)
        codes = icd_codes(mem, xrefs)
        k, source = k_for(d.mondo_id, gmap, rates, codes, ukb)
        row.update({"k": k, "k_source": source})
        gwas, mend, form = gene_counts(d.mondo_id, mondo, mem, evidence)
        row.update({"gwas_genes": gwas, "mendelian_genes": mend, "mendelian_form": form, "me_log10": me_map.get(d.mondo_id)})
        h = ukb_h2[ukb_h2.phenocode.isin(codes)].sort_values("z", ascending=False)
        row.update({"panukb_h2_liab": float(h.h2_liab.iloc[0]) if len(h) else None, "panukb_h2_z": float(h.z.iloc[0]) if len(h) else None})
        cases = None
        if d.gcst is not None and d.gcst in studies.index:
            cases, controls = int(studies.cases[d.gcst]), int(studies.controls[d.gcst])
            row.update({"n_cases": cases, "n_controls": controls, "n_eff": n_eff(cases, controls)})
            r = runs.loc[d.gcst] if d.gcst in runs.index else None
            if r is None or r.status != "ok":
                row["ldsc_status"] = "run_failed"
            else:
                z = float(r.h2_obs / r.h2_obs_se) if r.h2_obs_se and r.h2_obs_se > 0 else None
                row.update({"h2_obs": float(r.h2_obs), "h2_obs_se": float(r.h2_obs_se), "intercept": r.intercept, "h2_z": z,
                            "loci": int(r.loci), "concentration": r.concentration})
                if r.n_snps is None or r.n_snps < MIN_PLATFORM_SNPS:
                    row["ldsc_status"] = "targeted_array"
                elif z is None or z < Z_MIN:
                    row["ldsc_status"] = "low_z"
                elif k is None:
                    row["ldsc_status"] = "no_k"
                else:
                    h2l = liability_h2(float(r.h2_obs), cases / (cases + controls), k)
                    row["h2_liab"] = h2l
                    if not 0 < h2l <= 1:
                        row["ldsc_status"] = "out_of_bounds"
                    else:
                        row["ldsc_status"] = "ok"
                        row["score"] = int(r.loci) / h2l
                        row["log10_score"] = float(np.log10(row["score"])) if r.loci > 0 else None
        if cases is None:
            own = studies[studies.mondo_ids.map(lambda m: bool(set(m) & mem))].cases.dropna()
            cases = int(own.max()) if len(own) else None
        row["gwas_genes_per_10k_cases"] = gwas / cases * 1e4 if cases else None
        row["source_release_ids"] = releases
        out.append(row)
    cols = ["gwas_studies", "gws_associations", "loci_catalog", "n_cases", "n_controls", "n_eff", "h2_obs", "h2_obs_se",
            "intercept", "h2_z", "k", "k_source", "h2_liab", "ldsc_status", "loci", "concentration", "score", "log10_score",
            "gwas_genes", "gwas_genes_per_10k_cases", "mendelian_genes", "mendelian_form", "me_log10", "panukb_h2_liab",
            "panukb_h2_z", "source_release_ids"]
    with llm.connect() as conn, conn.cursor() as cur:
        for row in out:
            cur.execute(f"UPDATE preclin.disease_polygenicity SET {', '.join(f'{c} = %s' for c in cols)}, built_at = now() "
                        f"WHERE mondo_id = %s", [_py(row.get(c)) for c in cols] + [row["mondo_id"]])
        conn.commit()
    table = pd.DataFrame(out)
    print(table.ldsc_status.value_counts().to_string())
    both = table[table.h2_liab.notna() & table.panukb_h2_liab.notna()]
    if len(both) >= 3:
        rho, p = spearmanr(both.h2_liab, both.panukb_h2_liab)
        print(f"Pan-UKB check: Spearman {rho:.2f} (p {p:.3g}) over {len(both)} diseases with both estimates")
    print("K sources:", table.k_source.value_counts().to_dict())


def _py(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.bool_,)):
        return bool(v)
    return v


def counts(_args) -> None:
    t = frame.query("SELECT ldsc_status, entries FROM preclin.disease_polygenicity")
    ok = t[t.ldsc_status.eq("ok")]
    print(f"scored: {len(ok)} diseases, {int(ok.entries.sum())} entries, of {len(t)} candidate diseases with "
          f"{int(t.entries.sum())} entries on this build (Melissa's base: {BASE['diseases']} diseases, {BASE['entries']} entries)")
    print(t.groupby("ldsc_status").agg(diseases=("entries", "size"), entries=("entries", "sum")).to_string())


def check_unit(_args) -> None:
    cands = candidates()
    have = set(frame.query("SELECT mondo_id FROM preclin.disease_polygenicity").mondo_id)
    missing = cands[~cands.mondo_id.isin(have)]
    print(f"unit today: {len(cands)} diseases, {int(cands.entries.sum())} entries; profiled: {len(have)}")
    if len(missing):
        print("not profiled:", ", ".join(f"{r.mondo_id} ({r.label}, {r.entries})" for r in missing.itertuples()))
        sys.exit(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    for name, fn in (("load-gbd", load_gbd), ("choose", choose), ("build", build), ("counts", counts), ("check-unit", check_unit)):
        sub.add_parser(name).set_defaults(fn=fn)
    args = parser.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
