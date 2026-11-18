"""Pathway networks from enrichment results.

Enrichment map: significant terms are nodes, edges join terms whose gene sets overlap
(Jaccard), and modularity clustering groups redundant terms into modules.
Gene-pathway network: the lead terms linked to their genes; genes shared by several terms are
the hubs. Node and edge tables are written for Cytoscape (x, y = figure layout, inches).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

from . import figures
from .common import file_tag, palette, signif, write_tsv


def split_genes(text: str) -> List[str]:
    return list(dict.fromkeys(g.strip() for g in text.replace("/", ";").replace(",", ";").split(";") if g.strip()))


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


def enrichment_nodes(rows: List[dict], padj_cutoff: float = 0.05, top_n: Optional[int] = 30) -> List[dict]:
    """One node per term: genes pooled over groups, best padj, and the groups where it is significant."""
    nodes: Dict[str, dict] = {}
    for r in rows:
        if float(r["padj"]) > padj_cutoff:
            continue
        nd = nodes.setdefault(r["term"], {"term": r["term"], "padj": 1.0, "groups": set(), "genes": {}})
        nd["padj"] = min(nd["padj"], float(r["padj"]))
        nd["groups"].add(r["group"])
        nd["genes"].update(dict.fromkeys(split_genes(r["genes"])))
    out = [{"term": n["term"], "padj": n["padj"], "groups": ";".join(sorted(n["groups"])),
            "genes": ";".join(n["genes"]), "count": len(n["genes"])} for n in nodes.values()]
    out.sort(key=lambda n: n["padj"])
    return out[:top_n] if top_n else out


def modularity_modules(n: int, edges: List[tuple]) -> List[int]:
    """Greedy modularity clustering (Clauset-Newman-Moore) of a weighted graph; deterministic."""
    total = sum(w for _, _, w in edges)
    if not total:
        return list(range(n))
    comm = list(range(n))
    degree = [0.0] * n
    between: Dict[tuple, float] = {}
    for i, j, w in edges:
        degree[i] += w
        degree[j] += w
        between[(min(i, j), max(i, j))] = between.get((min(i, j), max(i, j)), 0.0) + w
    while True:
        best, pair = 1e-12, None
        for (a, b), w in sorted(between.items()):
            gain = w / total - degree[a] * degree[b] / (2 * total * total)
            if gain > best:
                best, pair = gain, (a, b)
        if pair is None:
            return comm
        a, b = pair
        comm = [a if c == b else c for c in comm]
        degree[a] += degree[b]
        merged: Dict[tuple, float] = {}
        for (x, y), w in between.items():
            x, y = (a if x == b else x), (a if y == b else y)
            if x != y:
                merged[(min(x, y), max(x, y))] = merged.get((min(x, y), max(x, y)), 0.0) + w
        between = merged


def build_enrichment_map(rows: List[dict], padj_cutoff: float = 0.05, top_n: int = 30, min_jaccard: float = 0.25,
                         k: int = 3) -> Optional[dict]:
    """Edges: gene-set overlap >= min_jaccard, keeping each term's k strongest links, so that
    redundant ontologies do not connect almost every pair."""
    nodes = enrichment_nodes(rows, padj_cutoff, top_n)
    if len(nodes) < 3:
        return None
    sets = [split_genes(n["genes"]) for n in nodes]
    pairs = [(i, j, jaccard(sets[i], sets[j])) for i in range(len(nodes)) for j in range(i + 1, len(nodes))]
    pairs = [p for p in pairs if p[2] >= min_jaccard]
    kept = set()
    for node in range(len(nodes)):
        mine = sorted((p for p in pairs if node in p[:2]), key=lambda p: -p[2])[:k]
        kept.update(mine)
    kept = sorted(kept)
    membership = modularity_modules(len(nodes), kept)
    best: Dict[int, float] = {}
    for nd, m in zip(nodes, membership):
        best[m] = min(best.get(m, 1.0), nd["padj"])
    rank = {m: i + 1 for i, m in enumerate(sorted(best, key=lambda m: best[m]))}   # modules ordered by best padj
    for nd, m in zip(nodes, membership):
        nd["module"] = rank[m]
    lead = {}
    for nd in sorted(nodes, key=lambda n: n["padj"]):
        lead.setdefault(nd["module"], nd["term"])
    for nd in nodes:
        nd["module_label"] = lead[nd["module"]]
    nodes.sort(key=lambda n: (n["module"], n["padj"]))
    index = {i: n["term"] for i, n in enumerate(enrichment_nodes(rows, padj_cutoff, top_n))}
    edges = [{"from": index[i], "to": index[j], "jaccard": round(w, 3)} for i, j, w in kept]
    return {"nodes": nodes, "edges": edges}


def gene_concept(rows: List[dict], padj_cutoff: float, lead_terms: List[str], n_terms: int = 5,
                 max_genes: int = 25) -> Optional[dict]:
    nodes = enrichment_nodes(rows, padj_cutoff, None)
    by_term = {n["term"]: n for n in nodes}
    chosen = list(dict.fromkeys([t for t in lead_terms if t in by_term] + [n["term"] for n in nodes]))[:n_terms]
    if len(chosen) < 2:
        return None
    terms = [by_term[t] for t in chosen]
    links = [(t["term"], g) for t in terms for g in split_genes(t["genes"])[:max_genes]]
    degree: Dict[str, int] = {}
    for _, gene in links:
        degree[gene] = degree.get(gene, 0) + 1
    hubs = sorted(({"gene": g, "n_pathways": k} for g, k in degree.items() if k >= 2), key=lambda h: (-h["n_pathways"], h["gene"]))
    return {"terms": terms, "links": links, "hub_genes": hubs}


def run_pathway_networks(rows: List[dict], out_dir, prefix: str, padj_cutoff: float = 0.05, top_n: int = 30,
                         min_jaccard: float = 0.25) -> dict:
    """Networks for every method of a standardised enrichment table."""
    out_dir = Path(out_dir)
    figs, tables, modules, hubs, warnings = {}, {}, {}, {}, []
    for method in dict.fromkeys(r["method"] for r in rows):
        sub = [r for r in rows if r["method"] == method]
        tag = file_tag(method)
        colours: Dict[str, str] = {}
        leads: List[str] = []
        try:
            em = build_enrichment_map(sub, padj_cutoff, top_n, min_jaccard)
            if em:
                png, pos, colours = figures.enrichment_map_figure(em["nodes"], em["edges"],
                                                                  f"Pathway enrichment map: {method}", out_dir,
                                                                  f"{prefix}_network_{tag}")
                figs[f"map_{method}"] = png
                cols = ("term", "module", "module_label", "padj", "count", "groups", "genes", "x", "y")
                for nd in em["nodes"]:
                    nd["x"], nd["y"] = (round(v, 3) for v in pos[nd["term"]])
                tables[f"network_nodes_{method}"] = write_tsv(out_dir / f"{prefix}_network_{tag}_nodes.tsv", cols,
                                                              [[nd[c] for c in cols] for nd in em["nodes"]])
                tables[f"network_edges_{method}"] = write_tsv(out_dir / f"{prefix}_network_{tag}_edges.tsv",
                                                              ("from", "to", "jaccard"),
                                                              [[e["from"], e["to"], e["jaccard"]] for e in em["edges"]])
                mods = {}
                for nd in em["nodes"]:
                    mods.setdefault(nd["module"], []).append(nd)
                modules[method] = [{"module": m, "lead_term": v[0]["module_label"], "n_terms": len(v),
                                    "best_padj": signif(min(n["padj"] for n in v)),
                                    "other_terms": "; ".join([n["term"] for n in v if n["term"] != v[0]["module_label"]][:3])}
                                   for m, v in sorted(mods.items())]
                leads = [m["lead_term"] for m in modules[method]]
        except Exception as err:
            warnings.append(f"Enrichment map for {method} failed: {err}")
        try:
            cn = gene_concept(sub, padj_cutoff, leads)
            if cn:
                spare = palette().colors(len(cn["terms"]))
                term_colour = {t["term"]: colours.get(t["term"], spare[i]) for i, t in enumerate(cn["terms"])}
                figs[f"cnet_{method}"] = figures.gene_concept_figure(cn["terms"], cn["links"], term_colour,
                                                                     f"Gene-pathway network: {method}", out_dir,
                                                                     f"{prefix}_cnet_{tag}")
                if cn["hub_genes"]:
                    hubs[method] = cn["hub_genes"][:15]
        except Exception as err:
            warnings.append(f"Gene-pathway network for {method} failed: {err}")
    return {"figures": figs, "tables": tables, "modules": modules, "hub_genes": hubs, "warnings": warnings}
