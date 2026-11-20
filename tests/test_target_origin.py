"""pipeline.target_origin and the small generalisations it needs in pipeline.literature. No database."""
from __future__ import annotations

import pandas as pd

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
    close = [f"related disease number {i}" for i in range(30)]
    names = to.disease_names("Rheumatoid Arthritis", ["RA", "rheumatoid arthritis", "arthritis, rheumatoid"], close)
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
