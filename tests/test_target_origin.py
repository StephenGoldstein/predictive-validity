"""pipeline.target_origin and the small generalisations it needs in pipeline.literature. No database."""
from __future__ import annotations

import pandas as pd
import pytest

from pipeline import literature
from pipeline import target_origin as to


def test_pubmed_ids_uses_the_callers_term_builder(monkeypatch):
    seen = []
    monkeypatch.setattr(literature, "pubmed_count", lambda term: seen.append(term) or 3)
    monkeypatch.setattr(literature, "eutils", lambda endpoint, **kw: "<r><IdList><Id>1</Id><Id>2</Id></IdList></r>")
    term, total, ids = literature.pubmed_ids(["x"], 2007, term=lambda lo, hi: f"CUSTOM {lo}..{hi}")
    assert term == "CUSTOM 1900..2007" and total == 3 and ids == ["1", "2"]
    assert seen == ["CUSTOM 1900..2007"]


def test_pubmed_ids_default_term_is_unchanged(monkeypatch):
    monkeypatch.setattr(literature, "pubmed_count", lambda term: 0)
    monkeypatch.setattr(literature, "eutils", lambda endpoint, **kw: "<r><IdList></IdList></r>")
    term, _, _ = literature.pubmed_ids(["x"], 2007)
    assert term == literature.pubmed_term(["x"], "1900", "2007")


def test_epmc_records_uses_the_callers_query(monkeypatch):
    posted = []
    monkeypatch.setattr(literature, "epmc_post",
                        lambda q, cursor="*", log=print, **kw: posted.append(q) or {"hitCount": 0, "resultList": {"result": []}})
    q, total, out = literature.epmc_records(["x"], 2007, query="MY QUERY")
    assert posted == ["MY QUERY"] and q == "MY QUERY" and total == 0 and out == []


# --- names, queries, records, bundle (pipeline.target_origin) ---

def hgnc(**kw) -> pd.DataFrame:
    row = {"symbol": "MET", "name": "MET proto-oncogene, receptor tyrosine kinase", "prev_symbols": [], "alias_symbols": []}
    row.update(kw)
    return pd.DataFrame([row])


def test_gene_names_keeps_the_symbol_even_when_short():
    names = to.gene_names("MET", hgnc(alias_symbols=["HGFR", "c-Met", "RCCP2"]), "hepatocyte growth factor receptor")
    assert names[0] == "MET"
    assert "hepatocyte growth factor receptor" in names and "HGFR" in names and "c-Met" in names


def test_gene_names_drops_identifiers_and_caps():
    aliases = [f"ALIAS{i:02d}X" for i in range(15)] + ["(2S)-2-amino", "Q9UQ80"]
    names = to.gene_names("TNF", hgnc(symbol="TNF", name="tumor necrosis factor", alias_symbols=aliases), None)
    assert names[0] == "TNF" and "(2S)-2-amino" not in names and "Q9UQ80" not in names
    assert len(names) == to.MAX_GENE_NAMES
    assert len({n.lower() for n in names}) == len(names)


def test_disease_names_shortest_first_and_capped():
    synonyms = ["RA", "rheumatoid arthritis", "arthritis, rheumatoid"] + [f"synonym number {i}" for i in range(30)]
    names = to.disease_names("Rheumatoid Arthritis", synonyms)
    assert names[0] == "Rheumatoid Arthritis"
    assert "RA" not in names  # shorter than four characters
    assert "rheumatoid arthritis" not in names[1:]  # case-insensitive duplicate of the label
    assert len(names) == to.MAX_DISEASE_NAMES
    rest = names[1:]
    assert rest == sorted(rest, key=len)


def test_pair_term_has_both_clauses_and_no_review_exclusion():
    term = to.pair_term(["TNF"], ["rheumatoid arthritis"])("1900", "2007")
    assert term == '("TNF"[tiab]) AND ("rheumatoid arthritis"[tiab]) AND ("1900"[dp] : "2007"[dp])'


def test_pair_epmc_query_ends_the_year_before_the_cutoff():
    q = to.pair_epmc_query(["TNF", "TNF-alpha"], ["rheumatoid arthritis"], 2008)
    assert 'TITLE_ABS:"TNF" OR TITLE_ABS:"TNF-alpha"' in q and 'TITLE_ABS:"rheumatoid arthritis"' in q
    assert "PUB_YEAR:[1900 TO 2007]" in q and "NOT SRC:MED" in q and "HAS_ABSTRACT:Y" in q and "review" not in q


def test_records_without_a_year_or_on_or_after_the_cutoff_are_dropped():
    fetched = [("1", "T1", "A", "J", 1999, "10.1/x", "abs"), ("2", "T2", "A", "J", None, None, "abs"),
               ("3", "T3", "A", "J", 2008, None, "abs"), ("4", "T4", "A", "J", 2007, None, None)]
    rows = to.pubmed_rows(fetched, 2008)
    assert [r["record_id"] for r in rows] == ["PMID:1", "PMID:4"]
    assert rows[0] == {"record_id": "PMID:1", "service": "pubmed", "title": "T1", "abstract": "abs", "year": 1999,
                       "doi": "10.1/x"}
    epmc = [("PPR:1", "PPR", "T", "abs", 2001, None), ("PPR:2", "PPR", "T", "abs", 2009, None),
            ("PPR:3", "PPR", "T", "abs", None, None)]
    assert [r["record_id"] for r in to.epmc_rows(epmc, 2008)] == ["PPR:1"]
    assert to.epmc_rows(epmc, 2008)[0]["service"] == "europepmc"


def rows_for(n: int) -> list[dict]:
    return [{"record_id": f"PMID:{i}", "service": "pubmed", "title": f"T{i}", "abstract": "x" * 2000,
             "year": 1950 + i, "doi": None} for i in range(1, n + 1)]


def test_bundle_takes_earliest_25_and_latest_15_without_repeats():
    ranked = to.ranked(rows_for(60))
    assert [r["rank"] for r in ranked][:3] == [1, 2, 3]
    b = to.bundle(ranked)
    assert len(b) == 40
    assert [r["rank"] for r in b] == list(range(1, 26)) + list(range(46, 61))
    assert len(to.bundle(to.ranked(rows_for(30)))) == 30


def test_bundle_cuts_abstracts():
    b = to.bundle(to.ranked(rows_for(3)))
    assert all(len(r["abstract"]) == to.ABSTRACT_CHARS for r in b)


# --- command line (pipeline.target_disease registers the origin commands) ---

def test_cli_registers_every_command():
    from pipeline import target_disease
    ap = target_disease.parser()
    a = ap.parse_args(["search", "--pairs", "window", "--redo", "--max-records", "3000"])
    assert a.fn.__name__ == "search" and a.pairs == "window" and a.redo and a.max_records == 3000
    assert ap.parse_args(["search"]).pairs == "pilot" and ap.parse_args(["search"]).max_records is None
    a = ap.parse_args(["gate", "--hand-check", "x.csv", "--reader", "Stephen"])
    assert a.fn.__name__ == "gate" and a.hand_check == "x.csv" and a.reader == "Stephen" and a.set == "a"
    a = ap.parse_args(["package", "--set", "b", "--out", "p.html"])
    assert a.fn.__name__ == "package" and a.set == "b" and a.out == "p.html"
    assert ap.parse_args(["label", "--hide-provisional"]).hide_provisional
    assert ap.parse_args(["blind"]).pairs == "set_a"
    for cmd in ["check", "sample", "summary", "build", "provisional"]:
        assert ap.parse_args([cmd]).fn.__name__ == cmd


# --- labelling: validation, rendering, rows ---

URL = "https://patents.google.com/patent/US123"


def label_payload() -> dict:
    return {"ensembl_id": "ENSG1", "mondo_id": "MONDO:1", "symbol": "TNF", "gene_name": "tumor necrosis factor",
            "disease": "rheumatoid arthritis", "cutoff_year": 2008, "n_records_total": 2,
            "records": [{"record_id": "PMID:1", "year": 1999, "title": "T1", "abstract": "a", "rank": 1},
                        {"record_id": "MED:2", "year": 2003, "title": "T2", "abstract": "b", "rank": 2}],
            "first_program": {"entity": "drug x", "sponsor": "Sponsor", "phase": 2, "entry": "2008-03-01",
                              "nct_ids": ["NCT00000001"]},
            "provisional_origin": "mouse_genetics", "provisional_year": 1995}


def good_label() -> dict:
    return {"origin": "diseased_tissue", "origin_year": 1999, "secondary_origin": "mouse_genetics", "grade": "stated",
            "quote": "TNF is raised in synovium", "sources": ["PMID:1"], "source_years": [1999], "rationale": "r",
            "discovery_route": "derived_from_active", "route_quote": "derived from", "route_sources": ["MED:2"],
            "route_source_years": [2003], "retrieved_urls": [URL]}


UNCLEAR = {"origin": "unclear", "origin_year": None, "secondary_origin": None, "grade": "none", "quote": "",
           "sources": [], "source_years": [], "discovery_route": "unclear", "route_quote": "", "route_sources": [],
           "route_source_years": []}


@pytest.mark.parametrize("change, ok", [
    ({}, True),
    (UNCLEAR, True),
    ({"sources": [URL], "source_years": [2005], "origin_year": 2005}, True),
    ({"secondary_origin": "diseased_tissue"}, False),
    ({"origin": "unclear"}, False),
    ({"grade": "none"}, False),
    ({"quote": "  "}, False),
    ({"sources": ["PMID:9"], "source_years": [1999]}, False),
    ({"source_years": [2000]}, False),
    ({"sources": [URL], "source_years": [2008], "origin_year": 2008}, False),
    ({"sources": ["PMID:1", URL], "source_years": [1999]}, False),
    ({"origin_year": 1997}, False),
    ({**UNCLEAR, "origin_year": 1999}, False),
    ({"discovery_route": "unclear"}, False),
    ({"route_sources": [], "route_source_years": []}, False),
    ({"route_sources": ["MED:2"], "route_source_years": [2004]}, False),
    ({"origin": "not_a_category"}, False),
])
def test_validate_label(change, ok):
    out = {**good_label(), **change}
    if ok:
        to.validate_label(label_payload(), out)
    else:
        with pytest.raises(ValueError):
            to.validate_label(label_payload(), out)


def test_render_label_omits_the_provisional_line_when_hidden():
    shown = to.render_label(label_payload())
    assert "mouse_genetics" in shown and "1995" in shown and "PMID:1" in shown and "NCT00000001" in shown
    assert "2007 or earlier" in shown
    hidden = to.render_label({**label_payload(), "provisional_origin": None, "provisional_year": None})
    assert "PROVISIONAL" not in hidden and "mouse_genetics" not in hidden


def test_undecided_is_the_items_not_in_decided():
    items = {"a": {"x": 1}, "b": {"x": 2}}
    assert to.undecided(items, [("a", 1, {})]) == {"b": {"x": 2}}


def test_origin_rows_follow_the_insert_order():
    rows = to.origin_rows([("ENSG1|MONDO:1", 77, good_label())])
    assert rows == [("ENSG1", "MONDO:1", "labelled", "diseased_tissue", 1999, "mouse_genetics", "stated", ["PMID:1"],
                     "TNF is raised in synovium", 77, to.LABEL_METHOD, "not_checked", "derived_from_active",
                     "derived from", ["MED:2"])]
    unclear = to.origin_rows([("ENSG1|MONDO:1", 78, {**good_label(), **UNCLEAR})])[0]
    assert unclear[3:9] == ("unclear", None, None, "none", [], None) and unclear[12:15] == ("unclear", None, [])


def test_replace_rows_only_touches_the_keys_it_decided():
    keys = [("ENSG1", "MONDO:1"), ("ENSG2", "MONDO:2")]
    sql, params = to.replace_sql(keys)
    assert "WHERE basis = 'labelled' AND (ensembl_id, mondo_id) IN %s" in sql
    assert params == (tuple(keys),)


# --- check and blind coder ---

@pytest.mark.parametrize("out, outcome", [
    ({"quote_found": "yes", "year_ok": "yes", "origin_follows": "yes"}, "supported"),
    ({"quote_found": "yes", "year_ok": "yes", "origin_follows": "partly"}, "supported"),
    ({"quote_found": "no", "year_ok": "yes", "origin_follows": "yes"}, "unsupported"),
    ({"quote_found": "yes", "year_ok": "no", "origin_follows": "yes"}, "unsupported"),
    ({"quote_found": "yes", "year_ok": "yes", "origin_follows": "no"}, "unsupported"),
    ({"quote_found": "could_not_open", "year_ok": "unknown", "origin_follows": "yes"}, "unverified"),
    ({"quote_found": "could_not_open", "year_ok": "no", "origin_follows": "yes"}, "unsupported"),
])
def test_check_outcome(out, outcome):
    assert to.check_outcome({**out, "note": ""}) == outcome


def blind_payload() -> dict:
    p = label_payload()
    return {k: p[k] for k in ["ensembl_id", "mondo_id", "symbol", "gene_name", "disease", "cutoff_year",
                              "n_records_total", "records"]}


def test_validate_blind_rejects_a_source_outside_the_bundle():
    good = {"origin": "diseased_tissue", "secondary_origin": None, "grade": "pattern", "quote": "q",
            "source": "PMID:1", "source_year": 1999, "note": ""}
    to.validate_blind(blind_payload(), good)
    for bad in [{"source": URL}, {"source": "PMID:9"}, {"source_year": 2000}, {"origin": "unclear"},
                {"grade": "none"}, {"secondary_origin": "diseased_tissue"}]:
        with pytest.raises(ValueError):
            to.validate_blind(blind_payload(), {**good, **bad})
    to.validate_blind(blind_payload(), {"origin": "unclear", "secondary_origin": None, "grade": "none", "quote": "",
                                        "source": None, "source_year": None, "note": ""})


def test_render_blind_shows_no_provisional_origin_and_no_program():
    text = to.render_blind(blind_payload())
    assert "PMID:1" in text and "2007 or earlier" in text
    assert "PROVISIONAL" not in text and "PROGRAM" not in text and "mouse_genetics" not in text


# --- pilot sample and reader package ---

def cohort_and_modes(n: int = 300):
    cohort = pd.DataFrame({"ensembl_id": [f"ENSG{i:04d}" for i in range(n)], "mondo_id": ["MONDO:1"] * n})
    modes = []
    for i in range(n):
        if i % 3 == 0:
            modes.append((f"ENSG{i:04d}", "MONDO:1", "not_causal"))
        if i % 5 == 0:
            modes.append((f"ENSG{i:04d}", "MONDO:1", "unresolved"))
        if i % 7 == 0:
            modes.append((f"ENSG{i:04d}", "MONDO:1", "drug_fell_short"))
    return cohort, pd.DataFrame(modes, columns=["ensembl_id", "mondo_id", "mode"])


def test_sample_excludes_unresolved_and_splits_sets():
    cohort, modes = cohort_and_modes()
    s = to.sample_pairs(cohort, modes, seed=to.PILOT_SEED)
    assert list(s.columns) == ["ensembl_id", "mondo_id", "group", "hand_check_set"]
    assert s.group.value_counts().to_dict() == {"not_causal": 50, "other": 50}
    assert s.groupby(["group", "hand_check_set"]).size().to_dict() == {("not_causal", "a"): 25, ("not_causal", "b"): 25,
                                                                       ("other", "a"): 25, ("other", "b"): 25}
    assert s.ensembl_id.is_unique
    unresolved = set(modes[modes["mode"] == "unresolved"].ensembl_id)
    assert not set(s.ensembl_id) & unresolved
    nc = set(modes[modes["mode"] == "not_causal"].ensembl_id)
    assert set(s[s.group == "not_causal"].ensembl_id) <= nc
    assert not set(s[s.group == "other"].ensembl_id) & nc
    again = to.sample_pairs(cohort, modes, seed=to.PILOT_SEED)
    assert again.equals(s)
    assert not to.sample_pairs(cohort, modes, seed=1).equals(s)


def test_sample_takes_fewer_when_a_group_is_small():
    cohort, modes = cohort_and_modes(40)
    s = to.sample_pairs(cohort, modes, seed=to.PILOT_SEED)
    nc = s[s.group == "not_causal"]
    assert 0 < len(nc) < 50 and set(nc.hand_check_set) <= {"a", "b"}


def package_rows() -> list[dict]:
    p = label_payload()
    return [{**p, "hand_check_set": "a", "group": "other"},
            {**p, "ensembl_id": "ENSG2", "symbol": "IL6", "gene_name": "interleukin 6", "mondo_id": "MONDO:2",
             "disease": "asthma", "hand_check_set": "a", "group": "not_causal",
             "records": [{"record_id": "PMID:77", "year": 1990, "title": "Unique title seven", "abstract": "z",
                          "rank": 1}], "n_records_total": 1}]


def test_package_holds_every_pair_and_no_label():
    html = to.render_package(package_rows(), "# Codebook\n\nDefinitions here.", "a")
    assert "ENSG1|MONDO:1" in html and "ENSG2|MONDO:2" in html
    assert "Unique title seven" in html and "pubmed.ncbi.nlm.nih.gov/77" in html
    assert "patents.google.com" in html and "before=priority:20080101" in html
    assert "Definitions here." in html and "mouse_genetics" in html  # the provisional origin the model also saw
    assert "retrieved_urls" not in html and "labelled" not in html and "diseased_tissue</option>" in html
    assert "<title>Target origin hand check" in html and len(html.encode()) < 16_000_000


# --- gate ---

def test_every_origin_has_a_layer():
    assert set(to.LAYER) == set(to.ORIGINS)
    assert to.LAYER["human_genetics"] == to.LAYER["tumour_genetics"] == "genetics"
    assert to.LAYER["diseased_tissue"] == to.LAYER["observational_biomarker"]
    assert to.LAYER["unclear"] == "unclear"


def gate_frames(disagree: int):
    keys = [(f"ENSG{i:03d}", "MONDO:1") for i in range(50)]
    origins = [to.ORIGINS[i % 12] for i in range(50)]  # every category but unclear, in rotation
    reader = pd.DataFrame([(e, m, o) for (e, m), o in zip(keys, origins)], columns=["ensembl_id", "mondo_id", "origin"])
    stored_origins = list(origins)
    for i in range(disagree):
        stored_origins[i] = to.ORIGINS[(i + 1) % 12]
    stored = pd.DataFrame([(e, m, o) for (e, m), o in zip(keys, stored_origins)],
                          columns=["ensembl_id", "mondo_id", "origin"])
    raw = pd.DataFrame([(e, m, o, None) for (e, m), o in zip(keys, stored_origins)],
                       columns=["ensembl_id", "mondo_id", "origin", "secondary_origin"])
    return reader, stored, raw


def test_gate_passes_at_kappa_point_six_and_gap_ten():
    from sklearn.metrics import cohen_kappa_score
    reader, stored, raw = gate_frames(disagree=10)
    expected = cohen_kappa_score(reader.origin, stored.origin)
    r = to.gate_result(reader, stored, raw, {"not_causal": 0.80, "other": 0.90}, None)
    assert r["n"] == 50 and abs(r["kappa_stored"] - expected) < 1e-9 and r["kappa_raw"] == r["kappa_stored"]
    assert r["gap_points"] == 10.0 and r["pass"] is True and r["kappa_blind_vs_reader"] is None
    assert 0 < r["any_origin_agree"] <= 1
    assert to.gate_result(reader, stored, raw, {"not_causal": 0.78, "other": 0.90}, None)["pass"] is False
    reader, stored, raw = gate_frames(disagree=30)
    assert to.gate_result(reader, stored, raw, {"not_causal": 0.80, "other": 0.90}, None)["pass"] is False


def test_gate_reports_blind_kappas_only_when_given():
    reader, stored, raw = gate_frames(disagree=0)
    blind = reader.rename(columns={"origin": "blind_origin"})
    r = to.gate_result(reader, stored, raw, {"not_causal": 0.9, "other": 0.9}, blind)
    assert r["kappa_blind_vs_reader"] == 1.0 and r["kappa_blind_vs_raw"] == 1.0
    assert r["kappa_stored_layer"] == 1.0


# --- model reference: majority of three blind readings ---

def test_majority_reference_and_consistency():
    readings = [("ENSG1|MONDO:1|r1", "human_genetics"), ("ENSG1|MONDO:1|r2", "human_genetics"),
                ("ENSG1|MONDO:1|r3", "mouse_genetics"),
                ("ENSG2|MONDO:1|r1", "unclear"), ("ENSG2|MONDO:1|r2", "unclear"), ("ENSG2|MONDO:1|r3", "unclear"),
                ("ENSG3|MONDO:1|r1", "cell_screen"), ("ENSG3|MONDO:1|r2", "diseased_tissue"),
                ("ENSG3|MONDO:1|r3", "unclear"),
                ("ENSG4|MONDO:1|r1", "cell_screen"), ("ENSG4|MONDO:1|r2", "cell_screen")]
    ref, stats = to.majority_reference(readings)
    assert ref.set_index("ensembl_id").origin.to_dict() == {"ENSG1": "human_genetics", "ENSG2": "unclear",
                                                            "ENSG4": "cell_screen"}
    assert stats["pairs"] == 4 and stats["unanimous"] == 1 and stats["majority"] == 2 and stats["no_majority"] == 1
    assert stats["replicates"] == 3 and 0 <= stats["pairwise_agreement"] <= 1


def test_blind_items_carry_the_replicate_in_key_and_payload():
    base = {"ENSG1|MONDO:1": {"ensembl_id": "ENSG1", "records": []}}
    items = to.replicate_items(base, 3)
    assert set(items) == {"ENSG1|MONDO:1|r1", "ENSG1|MONDO:1|r2", "ENSG1|MONDO:1|r3"}
    assert items["ENSG1|MONDO:1|r2"]["replicate"] == 2 and items["ENSG1|MONDO:1|r2"]["ensembl_id"] == "ENSG1"
    assert to.split_key("ENSG1|MONDO:1|r2") == ("ENSG1", "MONDO:1")


def test_bundle_tolerates_missing_abstracts_from_a_dataframe():
    import numpy as np
    rows = to.ranked([{"record_id": "PMID:1", "service": "pubmed", "title": "T", "abstract": None, "year": 1999, "doi": None},
                      {"record_id": "PMID:2", "service": "pubmed", "title": "T", "abstract": np.nan, "year": 2000, "doi": np.nan}])
    assert [r["abstract"] for r in to.bundle(rows)] == ["", ""]
    df = pd.DataFrame(rows)
    assert all(r["abstract"] is None for r in to.records_from(df))
