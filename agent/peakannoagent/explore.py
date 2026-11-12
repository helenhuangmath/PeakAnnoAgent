"""Data for the interactive pathway explorer (static/explore.js): the significant terms of an
enrich_pathways result and its enrichment-map networks. The web page serves it as JSON, the
self-contained report inlines it together with the script."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict, List

from .evidence import Evidence

ENRICH_COLS = ("method", "group", "term_id", "term", "gene_ratio", "padj", "count", "genes")
NODE_COLS = ("term", "module", "module_label", "padj", "count", "groups", "genes", "x", "y")
STATIC = Path(__file__).resolve().parent / "static"


def explore_data(ev: Evidence, table: Callable[[str], List[dict]], max_terms_per_set: int = 60) -> dict:
    """`table(name)` returns the rows of one of `ev.tables`. Keeps the top terms per method and
    peak group (by adjusted p) and every network's nodes and edges."""
    enrich: List[dict] = []
    if "enrichment" in ev.tables:
        by_set: Dict[tuple, List[dict]] = {}
        for r in table("enrichment"):
            by_set.setdefault((r.get("method"), r.get("group")), []).append(r)
        for rows in by_set.values():
            rows.sort(key=lambda r: float(r.get("padj") or 1))
            enrich.extend({k: r.get(k) for k in ENRICH_COLS if k in r} for r in rows[:max_terms_per_set])
    networks: Dict[str, dict] = {}
    for name in ev.tables:
        if name.startswith("network_nodes_"):
            key = name[len("network_nodes_"):]
            nodes = [{k: r.get(k) for k in NODE_COLS if k in r} for r in table(name)]
            edges = table(f"network_edges_{key}") if f"network_edges_{key}" in ev.tables else []
            if nodes:
                networks[key] = {"nodes": nodes,
                                 "edges": [{"from": e.get("from"), "to": e.get("to"), "jaccard": e.get("jaccard")} for e in edges]}
    return {"evidence": ev.id, "enrichment": enrich, "networks": networks}


def has_explore(ev: Evidence) -> bool:
    return ev.tool == "enrich_pathways" and ev.status == "ok" and bool(ev.tables)


def explore_js() -> str:
    return (STATIC / "explore.js").read_text()


def inline_json(data: dict) -> str:
    """JSON safe to embed inside a <script type="application/json"> element."""
    return json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
