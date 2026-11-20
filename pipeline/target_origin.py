"""What first proposed each target-disease hypothesis, read from sources dated before its first entry (issue #146,
steps 2 to 5). Commands, registered on pipeline.target_disease:

    search   retrieve PubMed and Europe PMC records naming the gene and the disease, dated before the pair's
             cutoff year (the year of its first Phase 2 or 3 entry), into preclin.target_disease_search and
             preclin.target_disease_source; shared with the target-literature screen (#154)
    label    an LLM reads the earliest and latest records (and the web, under the same date rule) and labels the
             origin, its grade and the first drug's discovery route (preclin.target_disease_origin, basis labelled)
    check    a second model verifies each graded label's quote, date and reading; a failed label becomes unclear
    blind    a second model labels the hand-check sets from the records alone, for comparison with the reader
    sample   draw the pilot of 100 pairs (50 with a settled "not causal" stop, 50 without) and the two hand-check sets
    package  write the reader's page for one hand-check set
    gate     score the pilot against the reader's sheet: agreement and coverage gates
    summary  counts of labels, grades, check outcomes and coverage

Every source must be dated before the cutoff year. Every model pass goes through pipeline.llm.
"""
from __future__ import annotations

import os
from collections.abc import Callable

import pandas as pd
from psycopg2.extras import execute_values

from pipeline import attribution, literature, llm
from pipeline.literature import names_clause
from report import frame
from report.section3 import CLOSE, neighbourhood
from report.windows import SECTION1

ORIGINS = ["human_genetics", "tumour_genetics", "diseased_tissue", "mouse_genetics", "curated_pathway", "cell_screen",
           "earlier_drug_other_disease", "observational_biomarker", "human_cell_perturbation",
           "animal_model_pharmacology", "clinical_observation", "phenotypic_screen", "unclear"]
GRADES = ["stated", "pattern", "none"]
ROUTES = ["phenotypic_screen", "target_screen", "de_novo_design", "derived_from_active", "repurposed_drug", "unclear"]

MAX_GENE_NAMES = 12
MAX_DISEASE_NAMES = 20
BUNDLE_EARLIEST = 25  # records shown to the model and the reader: the earliest by (year, id) ...
BUNDLE_LATEST = 15  # ... and the latest before the cutoff
ABSTRACT_CHARS = 1500

PILOT_CSV = os.path.join(frame.REPO, "report", "audits", "target_origin_pilot.csv")
PAIR_SETS = ["pilot", "set_a", "set_b", "window", "all"]

SEARCH_METHOD = "target_origin_search@v1"
SEARCH_DESCRIPTION = f"""For each target-disease pair (preclin.target_disease): PubMed (esearch, efetch) and Europe PMC
(records outside MEDLINE with an abstract) for titles or abstracts naming the gene (HGNC symbol first, then its name,
the target name and aliases passing pipeline.literature.usable, at most {MAX_GENE_NAMES}) and the disease (MONDO label,
synonyms and the labels of terms within Lin {CLOSE}, shortest first, at most {MAX_DISEASE_NAMES}), dated 1900 to the year
before the pair's first entry. Reviews are kept. A record without a year, or dated in or after the cutoff year, is
dropped. Every remaining record is stored with its rank by (year, record_id)."""


def key(ensembl_id: str, mondo_id: str) -> str:
    return f"{ensembl_id}|{mondo_id}"


def split_key(k: str) -> tuple[str, str]:
    e, m = k.split("|")[:2]
    return e, m


# --- names and queries ---

def disease_usable(name: str) -> bool:
    """A disease name a paper could use: 4 to 60 characters, at most six words, no brackets or punctuation."""
    return name_usable(name)


def _dedupe(names: list[str]) -> list[str]:
    seen: set[str] = set()
    out = []
    for n in names:
        n = n.strip()
        if n and n.lower() not in seen:
            seen.add(n.lower())
            out.append(n)
    return out


def name_usable(name: str) -> bool:
    """A gene or protein name a paper could use: 4 to 60 characters, at most six words, no brackets or punctuation
    (HGNC names of the form "X, receptor" are dropped; a paper writes neither half alone)."""
    return 4 <= len(name) <= 60 and len(name.split()) <= 6 and not any(c in name for c in "()[]{},;^")


def gene_names(symbol: str, hgnc: pd.DataFrame, target_name: str | None) -> list[str]:
    """The symbol first, always; then the HGNC name and the target name that pass name_usable, and previous and alias
    symbols that pass pipeline.literature.usable; at most MAX_GENE_NAMES."""
    names = [n for n in [*(hgnc.name if not hgnc.empty else []), target_name] if n and name_usable(n)]
    aliases: list[str] = []
    for r in hgnc.itertuples():
        aliases += list(r.prev_symbols or []) + list(r.alias_symbols or [])
    return _dedupe([symbol] + names + [a for a in aliases if a and literature.usable(a)])[:MAX_GENE_NAMES]


def disease_names(label: str, synonyms: list[str], close_labels: list[str]) -> list[str]:
    """The label first; then synonyms and the labels of close terms that pass disease_usable, shortest first; at most
    MAX_DISEASE_NAMES."""
    rest = sorted((n for n in _dedupe(list(synonyms) + list(close_labels)) if disease_usable(n)
                   and n.lower() != label.lower()), key=len)
    return [label] + rest[:MAX_DISEASE_NAMES - 1]


def pair_term(genes: list[str], diseases: list[str]) -> Callable[[str, str], str]:
    """An esearch term builder for a date span: gene names AND disease names in title or abstract; reviews kept."""
    return lambda lo, hi: f'{names_clause(genes)} AND {names_clause(diseases)} AND ("{lo}"[dp] : "{hi}"[dp])'


def pair_epmc_query(genes: list[str], diseases: list[str], cutoff_year: int) -> str:
    """Europe PMC records outside MEDLINE with an abstract, naming the gene and the disease, dated before the cutoff."""
    def clause(names: list[str]) -> str:
        return "(" + " OR ".join(f'TITLE_ABS:"{n}"' for n in names) + ")"
    return (f"{clause(genes)} AND {clause(diseases)} AND PUB_YEAR:[1900 TO {cutoff_year - 1}] "
            "AND NOT SRC:MED AND HAS_ABSTRACT:Y")


# --- records ---

def pubmed_rows(fetched: list[tuple], cutoff_year: int) -> list[dict]:
    """Rows from pipeline.attribution.efetch tuples (pmid, title, first author, journal, year, doi, abstract), dropping
    records without a year or dated in or after the cutoff year."""
    return [{"record_id": f"PMID:{pmid}", "service": "pubmed", "title": title or "", "abstract": abstract,
             "year": int(year), "doi": doi}
            for pmid, title, _author, _journal, year, doi, abstract in fetched
            if year is not None and int(year) < cutoff_year]


def epmc_rows(records: list[tuple], cutoff_year: int) -> list[dict]:
    """Rows from pipeline.literature.epmc_record tuples (record_id, source, title, abstract, year, doi), same rule."""
    return [{"record_id": rid, "service": "europepmc", "title": title or "", "abstract": abstract, "year": int(year),
             "doi": doi}
            for rid, _source, title, abstract, year, doi in records if year is not None and int(year) < cutoff_year]


def ranked(rows: list[dict]) -> list[dict]:
    """Rows ordered by (year, record_id) with a 1-based rank; a record id appears once."""
    seen: set[str] = set()
    unique = []
    for r in rows:
        if r["record_id"] not in seen:
            seen.add(r["record_id"])
            unique.append(r)
    out = sorted(unique, key=lambda r: (r["year"], r["record_id"]))
    return [{**r, "rank": i + 1} for i, r in enumerate(out)]


def bundle(rows: list[dict]) -> list[dict]:
    """What the model and the reader see: the earliest BUNDLE_EARLIEST and the latest BUNDLE_LATEST ranked rows, with
    abstracts cut to ABSTRACT_CHARS."""
    rows = sorted(rows, key=lambda r: r["rank"])
    chosen = rows[:BUNDLE_EARLIEST] + rows[max(BUNDLE_EARLIEST, len(rows) - BUNDLE_LATEST):]
    return [{**r, "abstract": (r["abstract"] or "")[:ABSTRACT_CHARS]} for r in chosen]


# --- pair sets ---

def pair_frame(which: str) -> pd.DataFrame:
    """The pairs a command works on: ensembl_id, mondo_id, symbol, first_entry, cutoff_year, validated_at_entry
    (and group, hand_check_set for the pilot sets)."""
    hub = frame.query("SELECT ensembl_id, mondo_id, symbol, first_entry, validated_at_entry FROM preclin.target_disease")
    hub["first_entry"] = pd.to_datetime(hub.first_entry)
    hub["cutoff_year"] = hub.first_entry.dt.year
    if which == "all":
        return hub
    if which == "window":
        return hub[hub.first_entry.between(pd.Timestamp(SECTION1[0]), pd.Timestamp(SECTION1[1]))].reset_index(drop=True)
    pilot = pd.read_csv(PILOT_CSV, dtype=str)
    if which == "set_a":
        pilot = pilot[pilot.hand_check_set == "a"]
    elif which == "set_b":
        pilot = pilot[pilot.hand_check_set == "b"]
    elif which != "pilot":
        raise ValueError(f"unknown pair set {which!r}; one of {PAIR_SETS}")
    return hub.merge(pilot[["ensembl_id", "mondo_id", "group", "hand_check_set"]], on=["ensembl_id", "mondo_id"])


# --- search ---

def log(msg: str) -> None:
    print(msg, flush=True)


def mondo_labels() -> dict[str, str]:
    t = frame.query("SELECT term_id, label FROM preclin.ontology_term WHERE ontology = 'MONDO'")
    return dict(zip(t.term_id, t.label))


def pair_names(ensembl_id: str, symbol: str, mondo_id: str, nb, labels: dict[str, str]) -> tuple[list[str], list[str]]:
    hgnc = frame.query("SELECT symbol, name, prev_symbols, alias_symbols FROM preclin.hgnc_gene WHERE ensembl_id = %s",
                       (ensembl_id,))
    tname = frame.query("SELECT min(name) AS name FROM public.targets WHERE ensembl_id = %s", (ensembl_id,))
    target_name = tname.name[0] if len(tname) and isinstance(tname.name[0], str) else None
    term = frame.query("SELECT label, synonyms FROM preclin.ontology_term WHERE term_id = %s", (mondo_id,))
    label = term.label[0] if len(term) else labels.get(mondo_id, mondo_id)
    synonyms = list(term.synonyms[0] or []) if len(term) else []
    close = [labels[t] for t in nb.terms(mondo_id, CLOSE) if t in labels and t != mondo_id]
    return gene_names(symbol, hgnc, target_name), disease_names(label, synonyms, close)


def search_pair(cur, row, nb, labels: dict[str, str], max_records: int | None) -> tuple[int, int, int]:
    """Search both services for one pair and replace its rows; returns (PubMed hits, Europe PMC hits, rows stored)."""
    cutoff = int(row.cutoff_year)
    genes, diseases = pair_names(row.ensembl_id, row.symbol, row.mondo_id, nb, labels)
    term = pair_term(genes, diseases)
    pm_query, pm_hits, pmids = literature.pubmed_ids(genes, cutoff - 1, term=term)
    pm_rows = pubmed_rows(attribution.efetch(pmids), cutoff) if pmids else []
    ep_query = pair_epmc_query(genes, diseases, cutoff)
    _, ep_hits, records = literature.epmc_records(genes, cutoff - 1, log, query=ep_query)
    ep_rows = epmc_rows(records, cutoff)
    if max_records is not None:
        pm_rows = sorted(pm_rows, key=lambda r: (r["year"], r["record_id"]))[:max_records]
        ep_rows = sorted(ep_rows, key=lambda r: (r["year"], r["record_id"]))[:max_records]
    rows = ranked(pm_rows + ep_rows)
    pair = (row.ensembl_id, row.mondo_id)
    cur.execute("DELETE FROM preclin.target_disease_source WHERE (ensembl_id, mondo_id) = (%s, %s)", pair)
    cur.execute("DELETE FROM preclin.target_disease_search WHERE (ensembl_id, mondo_id) = (%s, %s)", pair)
    execute_values(cur, """INSERT INTO preclin.target_disease_search
                             (ensembl_id, mondo_id, service, cutoff_year, gene_names, disease_names, query, n_hits)
                           VALUES %s""",
                   [(*pair, "pubmed", cutoff, genes, diseases, pm_query, int(pm_hits)),
                    (*pair, "europepmc", cutoff, genes, diseases, ep_query, int(ep_hits))])
    if rows:
        execute_values(cur, """INSERT INTO preclin.target_disease_source
                                 (ensembl_id, mondo_id, record_id, service, title, abstract, year, doi, rank)
                               VALUES %s""",
                       [(*pair, r["record_id"], r["service"], r["title"], r["abstract"], r["year"], r["doi"], r["rank"])
                        for r in rows], page_size=2000)
    return int(pm_hits), int(ep_hits), len(rows)


def searched() -> set[tuple[str, str]]:
    s = frame.query("SELECT DISTINCT ensembl_id, mondo_id FROM preclin.target_disease_search")
    return set(zip(s.ensembl_id, s.mondo_id))


def search(args) -> None:
    pairs = pair_frame(args.pairs)
    done = set() if args.redo else searched()
    todo = pairs[[(e, m) not in done for e, m in zip(pairs.ensembl_id, pairs.mondo_id)]]
    log(f"{len(pairs)} pairs in {args.pairs}; {len(todo)} to search")
    nb = neighbourhood()
    labels = mondo_labels()
    with llm.connect() as conn, conn.cursor() as cur:
        llm.record_method(cur, SEARCH_METHOD, "Pre-entry sources per target-disease pair", SEARCH_DESCRIPTION,
                          ["PubMed", "Europe PMC", "preclin.target_disease", "preclin.hgnc_gene", "public.targets",
                           "preclin.ontology_term", "MONDO"],
                          {"max_gene_names": MAX_GENE_NAMES, "max_disease_names": MAX_DISEASE_NAMES, "close": CLOSE,
                           "date_rule": "publication year < year of first_entry", "reviews": "kept",
                           "max_records": args.max_records})
        conn.commit()
    failed, stored = 0, 0
    with llm.connect() as conn, conn.cursor() as cur:
        for i, row in enumerate(todo.itertuples(), 1):
            try:
                pm, ep, n = search_pair(cur, row, nb, labels, args.max_records)
                conn.commit()
            except (literature.SourceDown, RuntimeError) as e:
                conn.rollback()
                failed += 1
                log(f"[{i}/{len(todo)}] {row.symbol} x {labels.get(row.mondo_id, row.mondo_id)} FAILED: {e}")
                continue
            stored += n
            log(f"[{i}/{len(todo)}] {row.symbol} x {labels.get(row.mondo_id, row.mondo_id)} cutoff {row.cutoff_year}: "
                f"pubmed {pm} europepmc {ep} stored {n}")
    log(f"searched {len(todo) - failed} pairs, {stored:,} records stored, {failed} failed")




# --- command line ---

def _pending(name: str):
    """A command registered before its task lands; it carries the command's name so the CLI test reads it."""
    def fn(args) -> None:
        raise NotImplementedError(f"{name} is not built yet")
    fn.__name__ = name
    return fn


def add_commands(sub) -> None:
    """Register the origin commands on pipeline.target_disease's subparsers."""
    def pairs_option(p, default="pilot"):
        p.add_argument("--pairs", choices=PAIR_SETS, default=default, help="which pairs to work on")

    p = sub.add_parser("search", help="retrieve pre-entry PubMed and Europe PMC records per pair")
    pairs_option(p)
    p.add_argument("--redo", action="store_true", help="search pairs already searched")
    p.add_argument("--max-records", type=int, default=None, help="keep only the earliest N records per service")
    p.set_defaults(fn=search)
    p = sub.add_parser("label", help="label the origin and discovery route of each pair from its sources")
    pairs_option(p)
    p.add_argument("--hide-provisional", action="store_true", help="do not show the model the provisional origin")
    p.set_defaults(fn=_pending("label"))
    p = sub.add_parser("check", help="second-model check of each graded label")
    pairs_option(p)
    p.set_defaults(fn=_pending("check"))
    p = sub.add_parser("blind", help="second model labels a hand-check set from the records alone")
    pairs_option(p, default="set_a")
    p.set_defaults(fn=_pending("blind"))
    sub.add_parser("sample", help="draw the pilot of 100 pairs and the two hand-check sets").set_defaults(fn=_pending("sample"))
    p = sub.add_parser("package", help="write the reader's page for a hand-check set")
    p.add_argument("--set", choices=["a", "b"], default="a")
    p.add_argument("--out", required=True, help="path of the HTML file to write")
    p.set_defaults(fn=_pending("package"))
    p = sub.add_parser("gate", help="score the pilot against the reader's sheet")
    p.add_argument("--hand-check", required=True, help="the reader's CSV")
    p.add_argument("--reader", required=True, help="who read the set")
    p.add_argument("--set", choices=["a", "b"], default="a")
    p.set_defaults(fn=_pending("gate"))
    p = sub.add_parser("summary", help="labels, grades, check outcomes and coverage")
    pairs_option(p)
    p.set_defaults(fn=_pending("summary"))
