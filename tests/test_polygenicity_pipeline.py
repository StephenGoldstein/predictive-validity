"""Polygenicity (#147): the disease-profile pipeline (pipeline.polygenicity) and its tables."""
import pytest

from pipeline import llm


@pytest.mark.db
def test_migration_156_tables_exist():
    with llm.connect() as conn, conn.cursor() as cur:
        cur.execute("""SELECT table_name FROM information_schema.tables WHERE table_schema = 'preclin'
                       AND table_name IN ('gwas_catalog_study', 'ldsc_run', 'gbd_prevalence', 'disease_polygenicity')""")
        assert {r[0] for r in cur.fetchall()} == {"gwas_catalog_study", "ldsc_run", "gbd_prevalence",
                                                   "disease_polygenicity"}


import numpy as np
import pandas as pd

from pipeline.polygenicity import choose_study, gene_counts, k_for, liability_h2


def test_liability_conversion_matches_the_standard_formula():
    assert abs(liability_h2(0.05, 0.5, 0.01) - 0.0276) < 5e-4    # K^2(1-K)^2 / (P(1-P) phi(t)^2), t = Phi^-1(1-K)


def studies(rows):
    return pd.DataFrame(rows, columns=["gcst", "mondo_ids", "mapped_trait_uris", "sumstats_path", "cases", "controls",
                                       "eur_share", "platform_snps"], dtype=object)


def test_choose_largest_european_exact_single_trait_study():
    s = studies([("G1", ["M1"], ["u1"], "p", 1000, 1000, 1.0, 2_500_000), ("G2", ["M1"], ["u1"], "p", 5000, 5000, 1.0, 2_500_000),
                 ("G3", ["M1"], ["u1"], "p", 9000, 9000, 0.5, 2_500_000),
                 ("G4", ["M1", "M9"], ["u1", "u2"], "p", 9000, 9000, 1.0, 2_500_000)])
    assert choose_study(s, "M1", {"M1", "M2"}) == ("G2", "exact", "ok")   # G3 not European, G4 two traits (Review Focus 3)


def test_descendant_fallback_and_statuses():
    s = studies([("G5", ["M2"], ["u2"], "p", 100, 100, 1.0, 2_500_000)])
    assert choose_study(s, "M1", {"M1", "M2"}) == ("G5", "descendant", "ok")
    assert choose_study(s.assign(sumstats_path=None), "M1", {"M1", "M2"}) == (None, None, "no_sumstats")
    assert choose_study(s.assign(cases=None), "M1", {"M1", "M2"}) == (None, None, "no_case_control")
    assert choose_study(s.assign(eur_share=0.8), "M1", {"M1", "M2"}) == (None, None, "not_european")
    assert choose_study(s.assign(platform_snps=130_000), "M1", {"M1", "M2"}) == (None, None, "targeted_array")
    assert choose_study(s.assign(platform_snps=None), "M1", {"M1", "M2"}) == ("G5", "descendant", "ok")  # unknown platform kept
    assert choose_study(s, "M7", {"M7"}) == (None, None, "no_study")


def test_k_prefers_gbd_then_panukb_then_none():   # Review Focus 4
    gbd_map = {"M1": ("exact", 627), "M2": ("none", None), "M3": ("exact", 1)}
    gbd = {627: 226.6, 1: 60000.0}
    ukb = pd.DataFrame({"phenocode": ["M06", "L40"], "n_cases_EUR": [5000, 50], "n_controls_EUR": [415000, 420000]})
    assert k_for("M1", gbd_map, gbd, {"M06"}, ukb) == (0.002266, "gbd")
    assert k_for("M2", gbd_map, gbd, {"M06"}, ukb) == (5000 / 420000, "panukb")
    assert k_for("M2", gbd_map, gbd, {"L40"}, ukb) == (None, "none")        # fewer than 100 cases
    assert k_for("M3", gbd_map, gbd, set(), ukb) == (None, "none")          # K >= 0.5


class StubMondo:
    def __init__(self, close):
        self.close = close

    def matches(self, evidence, program):
        c = (evidence, program) in self.close
        return c or evidence == program, c


def test_gene_counts_use_members_and_lin_neighbours_and_family():
    ev = pd.DataFrame({"target_ensembl": ["g1", "g2", "g2", "g3", "g4"], "mondo_id": ["M1", "M2", "M1", "M8", "M9"],
                       "datasource": ["gwas_credible_sets", "gwas_credible_sets", "orphanet", "clingen", "gene_burden"]})
    gwas, mend, form = gene_counts("M1", StubMondo({("M8", "M1")}), {"M1", "M2"}, ev)
    assert (gwas, mend, form) == (2, 2, True)      # g3 via Lin neighbour M8; M9 unrelated


def test_score_fields_keep_h2_liab_for_low_z_and_score_only_when_ok():
    from pipeline.polygenicity import score_fields
    low = score_fields(h2_obs=0.10, h2_obs_se=0.05, n_snps=1_000_000, loci=3, p=0.5, k=0.01)
    assert low["ldsc_status"] == "low_z" and low["h2_liab"] is not None and low["score"] is None
    ok = score_fields(h2_obs=0.20, h2_obs_se=0.02, n_snps=1_000_000, loci=10, p=0.5, k=0.01)
    assert ok["ldsc_status"] == "ok" and ok["score"] == pytest.approx(10 / ok["h2_liab"])
    assert ok["log10_score"] == pytest.approx(np.log10(ok["score"]))
    none = score_fields(h2_obs=0.20, h2_obs_se=0.02, n_snps=1_000_000, loci=0, p=0.5, k=0.01)
    assert none["ldsc_status"] == "ok" and none["score"] == 0 and none["log10_score"] is None
    assert score_fields(0.20, 0.02, 50_000, 10, 0.5, 0.01)["ldsc_status"] == "targeted_array"
    nok = score_fields(0.20, 0.02, 1_000_000, 10, 0.5, None)
    assert nok["ldsc_status"] == "no_k" and nok["h2_liab"] is None
    assert score_fields(0.90, 0.02, 1_000_000, 10, 0.5, 0.2)["ldsc_status"] == "out_of_bounds"
