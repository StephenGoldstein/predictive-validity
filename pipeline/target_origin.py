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
from pipeline.stop_causes import BLOCKED_DOMAINS
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




# --- labelling ---

DEFINITIONS = """- human_genetics: germline variants link the target to the disease (GWAS, gene burden, Mendelian)
- tumour_genetics: the gene is recurrently mutated in the cancer
- diseased_tissue: the gene's RNA level differs in diseased tissue
- mouse_genetics: knocked-out or mutant mice show the disease
- curated_pathway: a curated pathway shows a disease mutation disrupting the target
- cell_screen: a gene-knockout screen in cells
- earlier_drug_other_disease: a drug on the target was approved for another disease
- observational_biomarker: the target or its product tracks the disease in patients
- human_cell_perturbation: perturbing the target in human cells or tissue changes a disease phenotype
- animal_model_pharmacology: modulating the target in an animal disease model changes it
- clinical_observation: an effect seen in people, such as a side effect
- phenotypic_screen: a molecule found without a target hypothesis came first
- unclear: no source settles it"""

GRADE_RULE = """stated: a source says why the target was chosen for the disease (a rationale in an introduction or discussion, a
review written before the cutoff, a patent's background). pattern: no source says so, but a source shows the defining
observation of the category, published before any other kind of evidence for the pair. none: the sources would leave
the origin to inference; then the origin is unclear, the quote is empty and the sources are empty. A label graded
stated or pattern carries a verbatim quote from one of its cited sources."""

ROUTE_DEFINITIONS = """- phenotypic_screen: a screen in cells or animals with no target hypothesis
- target_screen: a high-throughput, fragment or virtual screen against the target; immunization or display for a biologic
- de_novo_design: designed from the sequence or structure (antisense, siRNA, structure-based design)
- derived_from_active: derived from an existing active: a known drug, an endogenous ligand, a natural product
- repurposed_drug: an existing drug tried in this disease
- unclear: no source settles it"""

LABEL_SYSTEM = f"""You work out what first proposed a drug target for a disease: the evidence that made someone think the target was worth a drug in that disease, before any drug on it was tested in patients for it. You are given the gene, the disease, the cutoff year (the year the hypothesis was first tested in Phase 2 or 3), the first program that tested it, a provisional origin read from dated evidence databases, and the records that name the gene and the disease and were published before the cutoff year: the earliest and the latest of them, with the total count. You may search and read the web. Only a source dated before the cutoff year counts: a later source, including a later history or review, is not evidence of what proposed the hypothesis at the time and must not be cited.

ORIGIN. One of these categories, with the database's own definitions:
{DEFINITIONS}
The origin is the kind of evidence that first proposed the pair, not the strongest evidence and not the latest. Where the sources show two kinds, the earlier is the origin and the other is secondary_origin. The provisional origin, when shown, is the earliest dated item in the evidence databases: confirm or correct it from the sources, and never repeat it without a source of your own.

GRADE. {GRADE_RULE}

SOURCES. Cite only record ids from the records given (PMID:... or the Europe PMC id shown) or URLs that you fetched or that a search returned in this conversation. Give one publication year per source, the year the source itself shows; every year must be before the cutoff year. origin_year is the year of the source that carries the origin, and null when the origin is unclear.

DISCOVERY ROUTE. Separately, where did the first program's drug come from? One of:
{ROUTE_DEFINITIONS}
It carries its own verbatim quote and sources under the same rules; unclear carries neither. It describes the drug, not the target hypothesis.

Answer in the JSON schema given. rationale is two or three sentences on how the sources settle the origin and the route, or why they do not."""

LABEL_SCHEMA = {
    "type": "object", "additionalProperties": False,
    "properties": {
        "origin": {"type": "string", "enum": ORIGINS},
        "origin_year": {"type": ["integer", "null"]},
        "secondary_origin": {"type": ["string", "null"], "enum": [*ORIGINS, None]},
        "grade": {"type": "string", "enum": GRADES},
        "quote": {"type": "string"},
        "sources": {"type": "array", "items": {"type": "string"}},
        "source_years": {"type": "array", "items": {"type": "integer"}},
        "rationale": {"type": "string"},
        "discovery_route": {"type": "string", "enum": ROUTES},
        "route_quote": {"type": "string"},
        "route_sources": {"type": "array", "items": {"type": "string"}},
        "route_source_years": {"type": "array", "items": {"type": "integer"}},
    },
    "required": ["origin", "origin_year", "secondary_origin", "grade", "quote", "sources", "source_years",
                 "rationale", "discovery_route", "route_quote", "route_sources", "route_source_years"],
}

TOOLS = [{"type": "web_search_20260209", "name": "web_search", "max_uses": 8, "blocked_domains": BLOCKED_DOMAINS},
         {"type": "web_fetch_20260209", "name": "web_fetch", "max_uses": 4, "blocked_domains": BLOCKED_DOMAINS}]
LABEL_METHOD = "target_origin_label@v1"
CONTEXT = ("provisional_origin", "provisional_year")  # shown beside the input; a change alone does not relabel
LABEL_PASSES = 3  # a pair whose answers fail validation is asked again, up to this many times


def render_records(p: dict) -> list[str]:
    n = p["n_records_total"]
    shown = len(p["records"])
    earliest = min(BUNDLE_EARLIEST, shown)
    which = f"the earliest {earliest} and the latest {shown - earliest} are shown" if n > shown else "all are shown"
    out = [f"RECORDS: {n} records name the gene and the disease before the cutoff year; {which}."]
    for r in p["records"]:
        out.append(f"[{r['record_id']}] ({r['year']}) {r['title']}")
        if r.get("abstract"):
            out.append(r["abstract"])
    if not p["records"]:
        out.append("(none)")
    return out


def render_label(p: dict) -> str:
    fp = p["first_program"]
    out = [f"TARGET: {p['symbol']} ({p['gene_name']}), {p['ensembl_id']}",
           f"DISEASE: {p['disease']} ({p['mondo_id']})",
           f"CUTOFF YEAR: {p['cutoff_year']}. Only sources dated {p['cutoff_year'] - 1} or earlier count.",
           f"FIRST PROGRAM: {fp['entity']}; sponsor {fp['sponsor'] or 'not recorded'}; Phase {fp['phase']} entry "
           f"{fp['entry']}; trials {', '.join(fp['nct_ids']) or 'not recorded'}"]
    if p.get("provisional_origin") is not None:
        year = f" ({p['provisional_year']})" if p.get("provisional_year") else ""
        out.append(f"PROVISIONAL ORIGIN (earliest dated item in the evidence databases): {p['provisional_origin']}{year}")
    return "\n".join(out + render_records(p))


def _check_sources(p: dict, sources: list, years: list, retrieved: set[str], what: str) -> None:
    bundle_years = {r["record_id"]: r["year"] for r in p["records"]}
    if len(years) != len(sources):
        raise ValueError(f"{what}: {len(sources)} sources but {len(years)} years")
    for s, y in zip(sources, years):
        if s in bundle_years:
            if y != bundle_years[s]:
                raise ValueError(f"{what}: {s} is dated {bundle_years[s]}, not {y}")
        elif llm.norm_url(s) not in retrieved:
            raise ValueError(f"{what}: cites a source it did not retrieve: {s}")
        if y >= p["cutoff_year"]:
            raise ValueError(f"{what}: {s} dated {y}, not before the cutoff year {p['cutoff_year']}")


def validate_label(p: dict, out: dict) -> None:
    """Raise ValueError when a label breaks a rule; pipeline.llm then stores the decision as invalid."""
    if out["origin"] not in ORIGINS or out["grade"] not in GRADES or out["discovery_route"] not in ROUTES:
        raise ValueError("origin, grade or discovery_route outside the lists")
    if out["secondary_origin"] is not None and (out["secondary_origin"] not in ORIGINS
                                                or out["secondary_origin"] == out["origin"]):
        raise ValueError(f"secondary_origin {out['secondary_origin']} with origin {out['origin']}")
    if (out["origin"] == "unclear") != (out["grade"] == "none"):
        raise ValueError(f"origin {out['origin']} with grade {out['grade']}")
    retrieved = {llm.norm_url(u) for u in out.get("retrieved_urls", [])}
    if out["grade"] == "none":
        if out["quote"].strip() or out["sources"] or out["origin_year"] is not None:
            raise ValueError("an unclear origin carries no quote, source or year")
    else:
        if not out["quote"].strip() or not out["sources"]:
            raise ValueError("a graded origin needs a quote and a source")
        _check_sources(p, out["sources"], out["source_years"], retrieved, "origin")
        if out["origin_year"] not in out["source_years"]:
            raise ValueError(f"origin_year {out['origin_year']} is not the year of a cited source")
    if (out["discovery_route"] == "unclear") != (not out["route_quote"].strip()):
        raise ValueError(f"discovery_route {out['discovery_route']} with route_quote {out['route_quote']!r}")
    if out["discovery_route"] == "unclear":
        if out["route_sources"]:
            raise ValueError("an unclear route carries no source")
    else:
        if not out["route_sources"]:
            raise ValueError("a graded route needs a source")
        _check_sources(p, out["route_sources"], out["route_source_years"], retrieved, "route")


LABEL_TASK = llm.Task(name="target_origin", prompt_version="v1", schema_version="v1", system=LABEL_SYSTEM,
                      schema=LABEL_SCHEMA, render=render_label, validate=validate_label, subject_type="target_disease",
                      model="claude-opus-5-5", effort="medium", max_tokens=16000, tools=TOOLS)

LABEL_DESCRIPTION = f"""For each pair: task {LABEL_TASK.name} {LABEL_TASK.prompt_version} ({LABEL_TASK.model}, effort
{LABEL_TASK.effort}, web search and fetch with {', '.join(BLOCKED_DOMAINS)} blocked) reads the gene, the disease, the
cutoff year, the first program, the provisional origin and the earliest {BUNDLE_EARLIEST} and latest {BUNDLE_LATEST}
pre-entry records (abstracts cut to {ABSTRACT_CHARS} characters) and returns the origin (13 categories), its year,
a secondary origin, the grade (stated, pattern, none), a verbatim quote, the sources with their years, and the first
drug's discovery route with its own quote and sources. Every cited source is a shown record or a URL retrieved in the
turn, dated before the cutoff year; a label breaking a rule is stored invalid and the pair is asked again, up to
{LABEL_PASSES} times. Rows replace earlier labelled rows for the decided pairs only."""


def first_programs(pairs: pd.DataFrame) -> dict[tuple[str, str], dict]:
    """The first program to test each pair: entity, sponsor, phase, entry date, NCT ids of the trials marking it."""
    from pipeline.target_disease import (
        EARLY,  # the hub module registers this one; import here, not at the top
    )
    from pipeline.target_history import targets
    from report.section1 import phase_rows
    snap = frame.snapshot()
    rows = phase_rows(frame.programs(), frame.program_trials(), snap, window=(EARLY, snap))
    rows = rows[rows.phase.isin([2, 3]) & ~rows.credited][["program_id", "phase", "entry"]]
    progs = frame.query("SELECT program_id, moiety_id, mondo_id, lead_sponsors FROM preclin.dev_program")
    t = targets()[["moiety_id", "ens"]].drop_duplicates()
    d = rows.merge(progs, on="program_id").merge(t, on="moiety_id")
    d = d.sort_values(["ens", "mondo_id", "entry", "program_id"]).drop_duplicates(["ens", "mondo_id"])
    want = set(zip(pairs.ensembl_id, pairs.mondo_id))
    d = d[[(e, m) in want for e, m in zip(d.ens, d.mondo_id)]]
    names = frame.query("SELECT moiety_id, name FROM preclin.moiety").set_index("moiety_id").name.to_dict()
    trials = frame.program_trials()
    trials = trials[trials.marks_entry]
    out = {}
    for r in d.itertuples():
        ncts = trials[(trials.program_id == r.program_id) & (trials.phase == r.phase)].nct_id.tolist()[:3]
        out[(r.ens, r.mondo_id)] = {"entity": names.get(r.moiety_id, str(r.moiety_id)),
                                   "sponsor": (list(r.lead_sponsors or []) or [None])[0], "phase": int(r.phase),
                                   "entry": str(pd.Timestamp(r.entry).date()), "nct_ids": ncts}
    return out


def sources_for(pairs: pd.DataFrame) -> dict[tuple[str, str], list[dict]]:
    """Every stored record of each pair, ranked."""
    if pairs.empty:
        return {}
    s = frame.query("""SELECT ensembl_id, mondo_id, record_id, service, title, abstract, year, rank
                       FROM preclin.target_disease_source WHERE (ensembl_id, mondo_id) IN %s ORDER BY rank""",
                    (tuple(zip(pairs.ensembl_id, pairs.mondo_id)),))
    return {k: g.drop(columns=["ensembl_id", "mondo_id"]).to_dict("records")
            for k, g in s.groupby(["ensembl_id", "mondo_id"])}


def label_items(pairs: pd.DataFrame, hide_provisional: bool) -> dict[str, dict]:
    labels = mondo_labels()
    gene_name = frame.query("SELECT ensembl_id, min(name) AS name FROM preclin.hgnc_gene WHERE ensembl_id IS NOT NULL "
                            "GROUP BY 1").set_index("ensembl_id").name.to_dict()
    prov = frame.query("SELECT ensembl_id, mondo_id, origin, origin_year FROM preclin.target_disease_origin "
                       "WHERE basis = 'provisional'")
    prov = {(e, m): (o, None if pd.isna(y) else int(y)) for e, m, o, y in
            zip(prov.ensembl_id, prov.mondo_id, prov.origin, prov.origin_year)}
    sources = sources_for(pairs)
    programs = first_programs(pairs)
    items = {}
    for r in pairs.itertuples():
        k = (r.ensembl_id, r.mondo_id)
        rows = sources.get(k, [])
        po, py = (None, None) if hide_provisional else prov.get(k, (None, None))
        items[key(*k)] = {"ensembl_id": r.ensembl_id, "mondo_id": r.mondo_id, "symbol": r.symbol,
                          "gene_name": gene_name.get(r.ensembl_id, r.symbol),
                          "disease": labels.get(r.mondo_id, r.mondo_id), "cutoff_year": int(r.cutoff_year),
                          "n_records_total": len(rows), "records": bundle(rows),
                          "first_program": programs.get(k, {"entity": "not recorded", "sponsor": None, "phase": 2,
                                                            "entry": str(r.first_entry.date()), "nct_ids": []}),
                          "provisional_origin": po, "provisional_year": py}
    return items


def undecided(items: dict[str, dict], decided: list[tuple]) -> dict[str, dict]:
    done = {k for k, _, _ in decided}
    return {k: v for k, v in items.items() if k not in done}


def origin_rows(decided: list[tuple[str, int, dict]]) -> list[tuple]:
    """Rows for preclin.target_disease_origin in the INSERT's column order."""
    out = []
    for k, decision_id, o in decided:
        e, m = split_key(k)
        out.append((e, m, "labelled", o["origin"], o["origin_year"], o["secondary_origin"], o["grade"],
                    list(o["sources"]), o["quote"].strip() or None, decision_id, LABEL_METHOD, "not_checked",
                    o["discovery_route"], o["route_quote"].strip() or None, list(o["route_sources"])))
    return out


def replace_sql(keys: list[tuple[str, str]]) -> tuple[str, tuple]:
    """The DELETE that clears the labelled rows of exactly these pairs before their new rows go in."""
    return ("DELETE FROM preclin.target_disease_origin WHERE basis = 'labelled' AND (ensembl_id, mondo_id) IN %s",
            (tuple(keys),))


INSERT_ORIGIN = """INSERT INTO preclin.target_disease_origin
                     (ensembl_id, mondo_id, basis, origin, origin_year, secondary_origin, grade, sources, quote,
                      llm_decision_id, method_id, check_outcome, discovery_route, route_quote, route_sources)
                   VALUES %s"""


def label(args) -> None:
    pairs = pair_frame(args.pairs)
    items = label_items(pairs, args.hide_provisional)
    params = {"task": LABEL_TASK.name, "prompt_version": LABEL_TASK.prompt_version, "model": LABEL_TASK.model,
              "effort": LABEL_TASK.effort, "tools": TOOLS, "origins": ORIGINS, "routes": ROUTES,
              "bundle": [BUNDLE_EARLIEST, BUNDLE_LATEST], "abstract_chars": ABSTRACT_CHARS,
              "hide_provisional": bool(args.hide_provisional), "passes": LABEL_PASSES, "runs": []}
    inputs = ["preclin.target_disease", "preclin.target_disease_source", "preclin.target_disease_origin",
              "preclin.dev_program", "preclin.dev_program_trial", "web search"]
    with llm.connect() as conn, conn.cursor() as cur:
        llm.record_method(cur, LABEL_METHOD, "What first proposed each target-disease hypothesis", LABEL_DESCRIPTION,
                          inputs, params)
        conn.commit()
    decided: list[tuple] = []
    todo = items
    for _ in range(LABEL_PASSES):
        if not todo:
            break
        run_id = llm.run(LABEL_TASK, todo, context=CONTEXT)
        params["runs"].append(str(run_id))
        new = llm.decisions(run_id)
        decided += new
        todo = undecided(todo, new)
        log(f"run {run_id}: {len(new)} decided, {len(todo)} undecided")
    rows = origin_rows(decided)
    with llm.connect() as conn, conn.cursor() as cur:
        llm.record_method(cur, LABEL_METHOD, "What first proposed each target-disease hypothesis", LABEL_DESCRIPTION,
                          inputs, params)
        if rows:
            cur.execute(*replace_sql([(r[0], r[1]) for r in rows]))
            execute_values(cur, INSERT_ORIGIN, rows, page_size=1000)
        conn.commit()
    log(f"pairs={len(items)} decided={len(rows)} undecided={len(todo)} runs={params['runs']}")


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
    p.set_defaults(fn=label)
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
