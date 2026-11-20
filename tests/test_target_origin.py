"""pipeline.target_origin and the small generalisations it needs in pipeline.literature. No database."""
from __future__ import annotations

from pipeline import literature


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
