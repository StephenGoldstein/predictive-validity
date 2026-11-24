"""GWAS Catalog studies mapped to MONDO, with cases, controls and European ancestry share (issue #147).

    python -m pipeline.gwas_catalog fetch   # the release files -> .cache/gwas_catalog/<RELEASE>/
    python -m pipeline.gwas_catalog load    # one row per study -> preclin.gwas_catalog_study

Release RELEASE of the Catalog: the studies file (v1.0.3.1, with the summary-statistics columns), the ancestry
file (one row per sample group with its broad ancestral category), the EFO trait mappings, the associations file
(all associations with their ontology annotations) and the list of harmonised summary-statistics files on the
Catalog's FTP site. A study's traits (MAPPED_TRAIT_URI, comma-separated) map to MONDO through the Open Targets
disease index (preclin.ot_disease), equivalent MONDO cross-references (preclin.ontology_xref) and the id itself
when it is a MONDO id. Cases and controls are read from the initial sample text ("2,080 European ancestry cases,
6,435 European ancestry controls"), summed over groups; a text with no cases gives none. The European share is the
initial sample's "European" category over its total in the ancestry file.

loci_catalog counts, for a disease and its MONDO descendants, the studies, the genome-wide significant
associations (P <= 5e-8, one chromosome and position) and the loci they merge into within 500 kb, the extended MHC
(chr6:25-34 Mb) as one locus. That count is descriptive: the registered loci come from LD clumping of one study's
summary statistics (pipeline.ldsc).
"""
from __future__ import annotations

import argparse
import os
import re
import time
import zipfile
from collections import defaultdict

import pandas as pd
import requests
from psycopg2.extras import execute_values

from pipeline import llm

RELEASE = "2026-11-19"
CACHE = os.path.join(llm.REPO, ".cache", "gwas_catalog", RELEASE)
RELEASES = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
SUMSTATS_BASE = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
STUDIES = "gwas-catalog-download-studies-v1.0.3.1.tsv"
ANCESTRY = "gwas-catalog-download-ancestries-v1.0.3.1.tsv"
TRAIT_MAPPINGS = "gwas-efo-trait-mappings.tsv"
ASSOCIATIONS_ZIP = "gwas-catalog-associations_ontology-annotated-full.zip"
ASSOCIATIONS = "gwas-catalog-download-associations-alt-full.tsv"
HARMONISED_LIST = "harmonised_list.txt"
GWS = 5e-8
LOCUS_BP = 500_000
MHC = ("6", 25_000_000, 34_000_000)
CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X"}

PLATFORM = re.compile(r"\[[^\]]*?~?\s*(\d[\d,]*)")
CASES = re.compile(r"(\d[\d,]*)\s+[^,\d]*?\bcases\b", re.I)
CONTROLS = re.compile(r"(\d[\d,]*)\s+[^,\d]*?\bcontrols\b", re.I)


def download(url: str, dest: str, attempts: int = 6) -> str:
    """Stream url to dest through dest.part, resuming a partial file, with attempts tries; returns dest."""
    if os.path.exists(dest):
        return dest
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    part = dest + ".part"
    for attempt in range(1, attempts + 1):
        try:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            headers = {"Range": f"bytes={have}-"} if have else {}
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                if r.status_code == 416:  # nothing left to fetch
                    break
                r.raise_for_status()
                mode = "ab" if have and r.status_code == 206 else "wb"
                with open(part, mode) as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            break
        except (requests.RequestException, OSError) as e:
            if attempt == attempts:
                raise
            time.sleep(15 * attempt)
            print(f"retry {attempt} for {url}: {e}")
    os.replace(part, dest)
    return dest


def fetch(_args) -> None:
    for name in (STUDIES, ANCESTRY, TRAIT_MAPPINGS, ASSOCIATIONS_ZIP):
        download(RELEASES + name, os.path.join(CACHE, name))
    download(SUMSTATS_BASE + HARMONISED_LIST, os.path.join(CACHE, HARMONISED_LIST))
    if not os.path.exists(os.path.join(CACHE, ASSOCIATIONS)):
        with zipfile.ZipFile(os.path.join(CACHE, ASSOCIATIONS_ZIP)) as z:
            member = next(n for n in z.namelist() if n.endswith(".tsv"))
            z.extract(member, CACHE)
            os.replace(os.path.join(CACHE, member), os.path.join(CACHE, ASSOCIATIONS))
    for name in (STUDIES, ANCESTRY, TRAIT_MAPPINGS, ASSOCIATIONS, HARMONISED_LIST):
        print(name, os.path.getsize(os.path.join(CACHE, name)), "bytes")


def _count(pattern: re.Pattern, text: str) -> int | None:
    found = [int(m.replace(",", "")) for m in pattern.findall(text)]
    return sum(found) if found else None


def parse_sample_size(text: str | None) -> tuple[int | None, int | None]:
    """(cases, controls) summed over the groups named in an initial-sample text; (None, None) without cases."""
    if not isinstance(text, str):
        return None, None
    cases = _count(CASES, text)
    if cases is None:
        return None, None
    return cases, _count(CONTROLS, text)


def platform_snps(text: str | None) -> int | None:
    """The largest bracketed SNP count in the Catalog's platform text ("Illumina [~2,500,000] (imputed)")."""
    if not isinstance(text, str):
        return None
    found = [int(m.replace(",", "")) for m in PLATFORM.findall(text)]
    return max(found) if found else None


def _ancestry_columns(ancestry: pd.DataFrame) -> pd.DataFrame:
    return ancestry.rename(columns={"NUMBER OF INDIVDUALS": "NUMBER OF INDIVIDUALS"})  # the Catalog's spelling


def eur_shares(ancestry: pd.DataFrame) -> dict[str, float]:
    """Study -> the initial sample's European share: individuals in the "European" category over all; studies with
    no initial-stage rows or no individuals are absent."""
    a = _ancestry_columns(ancestry)
    a = a[a["STAGE"].str.lower() == "initial"]
    n = pd.to_numeric(a["NUMBER OF INDIVIDUALS"], errors="coerce").fillna(0)
    total = n.groupby(a["STUDY ACCESSION"]).sum()
    eur = n.where(a["BROAD ANCESTRAL CATEGORY"] == "European", 0).groupby(a["STUDY ACCESSION"]).sum()
    share = (eur / total)[total > 0]
    return {k: float(v) for k, v in share.items()}


def eur_share(ancestry: pd.DataFrame, gcst: str) -> float | None:
    return eur_shares(ancestry).get(gcst)


def curie(uri: str) -> str:
    """http://www.ebi.ac.uk/efo/EFO_0000685 -> EFO:0000685; .../obo/MONDO_0005148 -> MONDO:0005148."""
    tail = uri.strip().rstrip("/").rsplit("/", 1)[-1]
    return tail.replace("_", ":", 1)


class TraitLookup:
    """Trait URIs to MONDO ids through a CURIE -> MONDO-id map; a MONDO CURIE also maps to itself."""

    def __init__(self, mapping: dict[str, set[str]]):
        self.mapping = mapping

    def to_mondo(self, uri: str) -> set[str]:
        c = curie(uri)
        out = set(self.mapping.get(c, ()))
        if c.startswith("MONDO:"):
            out.add(c)
        return out


def trait_lookup() -> TraitLookup:
    """The CURIE -> MONDO map from the Open Targets disease index and MONDO's equivalent cross-references."""
    mapping: dict[str, set[str]] = defaultdict(set)
    with llm.connect() as conn, conn.cursor() as cur:
        cur.execute("SELECT disease_id, mondo_ids FROM preclin.ot_disease WHERE mondo_ids <> '{}'")
        for disease_id, mondo_ids in cur.fetchall():
            mapping[disease_id].update(mondo_ids)
        cur.execute("SELECT xref_id, term_id FROM preclin.ontology_xref WHERE ontology = 'MONDO' AND equivalent")
        for xref_id, term_id in cur.fetchall():
            mapping[xref_id].add(term_id)
    return TraitLookup(dict(mapping))


def split_uris(text) -> list[str]:
    return [u.strip() for u in str(text).split(",") if isinstance(text, str) and u.strip()]


def harmonised_paths(path: str) -> dict[str, str]:
    """GCST accession -> path of its harmonised file, relative to SUMSTATS_BASE."""
    out = {}
    for line in open(path):
        p = line.strip().lstrip("./")
        parts = p.split("/")
        if len(parts) >= 3 and parts[-2] == "harmonised":
            out[parts[-3]] = p
    return out


def _relative(path: str | None) -> str | None:
    return path.lstrip("./") if isinstance(path, str) else None


def study_frame(studies: pd.DataFrame, ancestry: pd.DataFrame, lookup: TraitLookup,
                harmonised: dict[str, str]) -> pd.DataFrame:
    """One row per study with the preclin.gwas_catalog_study columns (object columns, None for missing)."""
    shares = eur_shares(ancestry)
    rows = []
    for _, s in studies.iterrows():
        uris = split_uris(s["MAPPED_TRAIT_URI"])
        mondo = sorted(set().union(*(lookup.to_mondo(u) for u in uris)) if uris else set())
        cases, controls = parse_sample_size(s["INITIAL SAMPLE SIZE"])
        published = pd.to_datetime(s["DATE"], errors="coerce")
        count = pd.to_numeric(s["ASSOCIATION COUNT"], errors="coerce")
        rows.append({"gcst": s["STUDY ACCESSION"], "pubmed_id": s["PUBMED ID"] if isinstance(s["PUBMED ID"], str) else None,
                     "published": published.date() if not pd.isna(published) else None,
                     "trait": s["DISEASE/TRAIT"], "mapped_trait_uris": uris, "mondo_ids": mondo,
                     "initial_sample": s["INITIAL SAMPLE SIZE"] if isinstance(s["INITIAL SAMPLE SIZE"], str) else None,
                     "cases": cases, "controls": controls, "eur_share": shares.get(s["STUDY ACCESSION"]),
                     "sumstats_path": _relative(harmonised.get(s["STUDY ACCESSION"])),
                     "association_count": int(count) if not pd.isna(count) else None,
                     "platform_snps": platform_snps(s.get("PLATFORM [SNPS PASSING QC]"))})
    return pd.DataFrame(rows, dtype=object)


def loci_catalog(assoc: pd.DataFrame, members: set[str], lookup: TraitLookup) -> tuple[int, int, int]:
    """(studies, genome-wide significant associations, loci) for a disease's members in the associations file."""
    uri_text = assoc["MAPPED_TRAIT_URI"].astype(object)
    hits = {t: any(lookup.to_mondo(u) & members for u in split_uris(t)) for t in uri_text.dropna().unique()}
    a = assoc[uri_text.map(hits).fillna(False).astype(bool)]
    studies = int(a["STUDY ACCESSION"].nunique())
    p = pd.to_numeric(a["P-VALUE"], errors="coerce")
    pos = pd.to_numeric(a["CHR_POS"], errors="coerce")
    g = a[(p <= GWS) & a["CHR_ID"].astype(str).isin(CHROMOSOMES) & pos.notna()].assign(pos=pos)
    loci = 0
    for chrom, group in g.groupby("CHR_ID"):
        last = None
        for position in sorted(group.pos):
            if chrom == MHC[0] and MHC[1] <= position <= MHC[2]:
                position = MHC[1]  # every MHC association joins one locus
            if last is None or position - last > LOCUS_BP or (position == MHC[1] and last != MHC[1] and chrom == MHC[0]):
                loci += 1
            last = position
    return studies, int(len(g)), loci


def read_studies() -> pd.DataFrame:
    return pd.read_csv(os.path.join(CACHE, STUDIES), sep="\t", dtype=str, low_memory=False)


def read_ancestry() -> pd.DataFrame:
    return _ancestry_columns(pd.read_csv(os.path.join(CACHE, ANCESTRY), sep="\t", dtype=str, low_memory=False))


def read_associations() -> pd.DataFrame:
    cols = ["STUDY ACCESSION", "MAPPED_TRAIT_URI", "CHR_ID", "CHR_POS", "P-VALUE"]
    return pd.read_csv(os.path.join(CACHE, ASSOCIATIONS), sep="\t", dtype=str, usecols=cols, low_memory=False)


def load(_args) -> None:
    studies, ancestry = read_studies(), read_ancestry()
    harmonised = harmonised_paths(os.path.join(CACHE, HARMONISED_LIST))
    frame = study_frame(studies, ancestry, trait_lookup(), harmonised)
    with llm.connect() as conn, conn.cursor() as cur:
        rid = llm.register_release(cur, "gwas_catalog", RELEASE, RELEASES, llm.sha256_file(os.path.join(CACHE, STUDIES)),
                                   {"files": [STUDIES, ANCESTRY, TRAIT_MAPPINGS, ASSOCIATIONS, HARMONISED_LIST],
                                    "studies": int(len(frame)), "harmonised": len(harmonised)})
        cur.execute("DELETE FROM preclin.disease_polygenicity")
        cur.execute("DELETE FROM preclin.ldsc_run")
        cur.execute("DELETE FROM preclin.gwas_catalog_study")
        cols = ["gcst", "pubmed_id", "published", "trait", "mapped_trait_uris", "mondo_ids", "initial_sample", "cases",
                "controls", "eur_share", "sumstats_path", "association_count", "platform_snps"]
        rows = [tuple(r[c] for c in cols) + (rid,) for r in frame.to_dict("records")]
        execute_values(cur, f"INSERT INTO preclin.gwas_catalog_study ({', '.join(cols)}, source_release_id) VALUES %s",
                       rows, page_size=5000)
        conn.commit()
    mapped = frame.mondo_ids.map(len).gt(0)
    both = frame.cases.notna() & frame.controls.notna()
    print(f"{len(frame)} studies: {int(mapped.sum())} mapped to MONDO, {frame.sumstats_path.notna().sum()} with "
          f"harmonised summary statistics, {int(both.sum())} with cases and controls, "
          f"{int((mapped & frame.sumstats_path.notna() & both).sum())} with all three")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    sub.add_parser("fetch").set_defaults(fn=fetch)
    sub.add_parser("load").set_defaults(fn=load)
    args = parser.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
