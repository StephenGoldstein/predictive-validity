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
    assert k_for("M1", gbd_map, gbd, {"M06"}, ukb) == (0.002266, "gbd", None, False)
    assert k_for("M2", gbd_map, gbd, {"M06"}, ukb) == (5000 / 420000, "panukb", "M06", False)
    assert k_for("M2", gbd_map, gbd, {"M060"}, ukb) == (5000 / 420000, "panukb", "M06", True)   # through the parent code
    assert k_for("M2", gbd_map, gbd, {"L40"}, ukb) == (None, "none", None, False)        # fewer than 100 cases
    assert k_for("M3", gbd_map, gbd, set(), ukb) == (None, "none", None, False)          # K >= 0.5


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
    assert none["ldsc_status"] == "zero_loci" and none["score"] is None and none["log10_score"] is None
    assert none["h2_liab"] == pytest.approx(ok["h2_liab"])   # the no-score group keeps its heritability
    assert score_fields(0.20, 0.02, 50_000, 10, 0.5, 0.01)["ldsc_status"] == "targeted_array"
    nok = score_fields(0.20, 0.02, 1_000_000, 10, 0.5, None)
    assert nok["ldsc_status"] == "no_k" and nok["h2_liab"] is None
    assert score_fields(0.90, 0.02, 1_000_000, 10, 0.5, 0.2)["ldsc_status"] == "out_of_bounds"


def test_choose_study_skips_excluded_studies_and_breaks_ties_on_accession():
    from pipeline.polygenicity import exome_studies
    s = studies([("G1", ["M1"], ["u1"], "p", 1000, 1000, 1.0, 2_500_000), ("G2", ["M1"], ["u1"], "p", 5000, 5000, 1.0, 2_500_000),
                 ("G0", ["M1"], ["u1"], "p", 5000, 5000, 1.0, 2_500_000)])
    assert choose_study(s, "M1", {"M1"}) == ("G0", "exact", "ok")                      # G0 and G2 tie on n_eff: the accession decides
    assert choose_study(s, "M1", {"M1"}, exclude={"G0", "G2"}) == ("G1", "exact", "ok")
    assert choose_study(s, "M1", {"M1"}, exclude={"G0", "G1", "G2"}) == (None, None, "targeted_array")
    cat = pd.DataFrame({"STUDY ACCESSION": ["G1", "G2", "G3", "G4"],
                        "GENOTYPING TECHNOLOGY": ["Genome-wide genotyping array", "Genome-wide genotyping array, Exome-wide sequencing",
                                                  "Genome-wide sequencing, Exome-wide sequencing", "Exome genotyping array"]})
    assert exome_studies(cat) == {"G2", "G4"}


def test_panukb_k_is_bounded_and_failed_runs_keep_their_reason():
    from pipeline.polygenicity import failed_status
    ukb = pd.DataFrame({"phenocode": ["X1"], "n_cases_EUR": [600], "n_controls_EUR": [400]})
    assert k_for("M1", {}, {}, {"X1"}, ukb) == (None, "none", None, False)
    assert failed_status("no_effect") == "no_effect" and failed_status("fetch_failed") == "run_failed" and failed_status(None) == "run_failed"


def test_evidence_frame_reads_every_gene_of_the_full_association_file(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    from pipeline.polygenicity import evidence_frame
    pq.write_table(pa.table({"diseaseId": ["EFO_1", "EFO_1", "EFO_2", "EFO_3"], "targetId": ["g1", "g2", "g2", "g9"],
                             "aggregationValue": ["gwas_credible_sets", "orphanet", "gwas_credible_sets", "chembl"]}),
                   tmp_path / "part-0.parquet")
    diseases = pd.DataFrame({"disease_id": ["EFO:1", "EFO:2", "EFO:3"], "mondo_id": ["M1", "M2", "M3"]})
    ev = evidence_frame(str(tmp_path), diseases)
    assert set(ev.columns) == {"target_ensembl", "mondo_id", "datasource"}
    assert sorted(map(tuple, ev.values.tolist())) == [("g1", "M1", "gwas_credible_sets"), ("g2", "M1", "orphanet"),
                                                      ("g2", "M2", "gwas_credible_sets")]   # chembl is not a genetic source


def test_association_digest_is_deterministic_and_sensitive_to_the_files(tmp_path):
    from pipeline.polygenicity import association_digest
    (tmp_path / "part-0.parquet").write_bytes(b"a")
    (tmp_path / "part-1.parquet").write_bytes(b"b")
    d1 = association_digest(str(tmp_path))
    assert d1 == association_digest(str(tmp_path)) and len(d1) == 64
    (tmp_path / "part-1.parquet").write_bytes(b"c")
    assert association_digest(str(tmp_path)) != d1


def test_panukb_h2_keeps_european_icd10_rows_that_pass_qc(tmp_path, monkeypatch):
    """The real h2 manifest names the ancestry column `pop` (and pandas has a .pop method, so it is read by key)."""
    import gzip
    from pipeline import polygenicity
    rows = ["trait_type\tphenocode\tpop\tqcflags.pass_all\testimates.final.h2_liability\testimates.final.h2_z",
            "icd10\tM06\tEUR\tTrue\t0.12\t6.5",
            "icd10\tM06\tAFR\tTrue\t0.30\t1.2",
            "icd10\tL40\tEUR\tFalse\t0.20\t2.0",
            "phecode\t714\tEUR\tTrue\t0.15\t7.0"]
    with gzip.open(tmp_path / "h2_manifest.tsv.bgz", "wt") as f:
        f.write("\n".join(rows) + "\n")
    monkeypatch.setattr(polygenicity, "PANUKB", str(tmp_path))
    h = polygenicity.panukb_h2()
    assert h.values.tolist() == [["M06", 0.12, 6.5]]


def test_excluded_runs_skips_studies_with_too_few_merged_snps_or_no_usable_effect():
    """A run that merged fewer than MIN_PLATFORM_SNPS SNPs (a targeted array the Catalog lists as genome-wide) or that
    found no usable autosomal effect (an X-only study, a p-value-only file) excludes its study so the next one is chosen."""
    from pipeline import polygenicity
    runs = pd.DataFrame({"gcst": ["GCST003044", "GCST90444373", "GCST90027158", "GCST000001"],
                         "n_snps": [48_002, None, 1_168_849, None],
                         "status": ["ok", "no_effect", "ok", "fetch_failed"]})
    assert polygenicity.excluded_runs(runs) == {"GCST003044", "GCST90444373"}
