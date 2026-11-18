"""Literature and public-dataset search: Europe PMC (papers with abstracts) and NCBI GEO."""

from __future__ import annotations

import html
import re
from typing import List

from .http import get_json, ncbi_get_json

EUROPEPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def _clip(text: str, n: int) -> str:
    text = re.sub(r"<[^>]+>", "", html.unescape(text or "")).strip()
    return text if len(text) <= n else text[: n - 1].rstrip() + "…"


SORTS = {"cited": "CITED desc", "recent": "P_PDATE_D desc"}


def search_europepmc(query: str, max_results: int = 8, from_year: int = None,
                     abstract_chars: int = 900, sort: str = "cited") -> List[dict]:
    """Papers (preprints flagged) matching `query`, with clipped abstracts.

    sort='cited' surfaces the landmark studies; 'recent' the newest work. (The API's
    default order is effectively by date, and it rejects sort=RELEVANCE with a 503.)
    """
    q = query
    if from_year:
        q = f"({q}) AND PUB_YEAR:[{int(from_year)} TO 3000]"
    data = get_json(EUROPEPMC, {"query": q, "format": "json", "resultType": "core",
                                "pageSize": max(1, min(int(max_results), 25)),
                                "sort": SORTS.get(sort, SORTS["cited"])})
    papers = []
    for r in data.get("resultList", {}).get("result", []):
        pmid = r.get("pmid")
        journal = (r.get("journalInfo") or {}).get("journal", {}).get("title") or r.get("journalTitle", "")
        papers.append({
            "pmid": pmid,
            "id": f"PMID:{pmid}" if pmid else r.get("id"),
            "title": _clip(r.get("title", ""), 300),
            "authors": _clip(r.get("authorString", ""), 120),
            "journal": journal,
            "year": r.get("pubYear"),
            "doi": r.get("doi"),
            "preprint": r.get("source") == "PPR",
            "cited_by": r.get("citedByCount", 0),
            "abstract": _clip(r.get("abstractText", ""), abstract_chars),
        })
    return papers


def search_geo(query: str, organism: str = None, max_results: int = 10) -> List[dict]:
    """GEO Series matching `query` (e.g. 'Tox ChIP-seq CD8 T cell'), with linked PubMed IDs."""
    term = f"({query}) AND gse[ETYP]"
    if organism:
        term += f' AND "{organism}"[Organism]'
    search = ncbi_get_json("esearch", {"db": "gds", "term": term, "retmax": max(1, min(int(max_results), 30)),
                                       "sort": "relevance"})
    ids = search.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []
    summary = ncbi_get_json("esummary", {"db": "gds", "id": ",".join(ids)}).get("result", {})
    out = []
    for uid in summary.get("uids", []):
        s = summary[uid]
        out.append({
            "accession": s.get("accession"),
            "title": _clip(s.get("title", ""), 250),
            "summary": _clip(s.get("summary", ""), 500),
            "type": s.get("gdstype"),
            "organism": s.get("taxon"),
            "n_samples": s.get("n_samples"),
            "date": s.get("pdat"),
            "pubmed_ids": [str(p) for p in s.get("pubmedids", [])],
            "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={s.get('accession')}",
        })
    return out
