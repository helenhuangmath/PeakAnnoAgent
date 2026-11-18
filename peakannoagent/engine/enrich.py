"""Pathway over-representation analysis and pathway networks.

Gene sets come from public GMT files (MSigDB for human and mouse, Enrichr for mouse KEGG),
downloaded once and cached, or from the user's own GMT. Each peak-linked gene list is tested
against the genome background with the hypergeometric test and Benjamini-Hochberg correction.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from . import figures
from .common import bh_adjust, cache_dir, download, file_tag, read_tsv, signif, write_tsv
from .reference import normalize_genome

MSIGDB = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1."
ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName="
# method -> (human file, mouse file). A mouse entry of None uses the human sets, matched by symbol.
LIBRARIES: Dict[str, Tuple[Optional[str], Optional[str]]] = {
    "go_bp": ("Hs/c5.go.bp.v2024.1.Hs.symbols.gmt", "Mm/m5.go.bp.v2024.1.Mm.symbols.gmt"),
    "go_mf": ("Hs/c5.go.mf.v2024.1.Hs.symbols.gmt", "Mm/m5.go.mf.v2024.1.Mm.symbols.gmt"),
    "go_cc": ("Hs/c5.go.cc.v2024.1.Hs.symbols.gmt", "Mm/m5.go.cc.v2024.1.Mm.symbols.gmt"),
    "hallmark": ("Hs/h.all.v2024.1.Hs.symbols.gmt", "Mm/mh.all.v2024.1.Mm.symbols.gmt"),
    "reactome": ("Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt", "Mm/m2.cp.reactome.v2024.1.Mm.symbols.gmt"),
    "wikipathways": ("Hs/c2.cp.wikipathways.v2024.1.Hs.symbols.gmt", "Mm/m2.cp.wikipathways.v2024.1.Mm.symbols.gmt"),
    "kegg": ("Hs/c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt", "enrichr:KEGG_2019_Mouse"),
    "tf_targets": ("Hs/c3.tft.gtrd.v2024.1.Hs.symbols.gmt", "Mm/m3.gtrd.v2024.1.Mm.symbols.gmt"),
    "immunesigdb": ("Hs/c7.immunesigdb.v2024.1.Hs.symbols.gmt", None),
    "cell_types": ("Hs/c8.all.v2024.1.Hs.symbols.gmt", "Mm/m8.all.v2024.1.Mm.symbols.gmt"),
}
ALIASES = {"msigdb_h": "hallmark", "msigdb_c2:cp:reactome": "reactome", "msigdb_c2:cp:wikipathways": "wikipathways",
           "msigdb_c3:tft:gtrd": "tf_targets", "msigdb_c7:immunesigdb": "immunesigdb", "msigdb_c8": "cell_types",
           "msigdb_c5:go:bp": "go_bp", "msigdb_c5:go:mf": "go_mf", "msigdb_c5:go:cc": "go_cc"}
METHODS = list(LIBRARIES) + ["custom"]


def canonical_method(method: str) -> str:
    return ALIASES.get(method.lower(), method.lower() if method.lower() in LIBRARIES or method.lower() == "custom" else method)


def read_gmt(path) -> Dict[str, List[str]]:
    sets: Dict[str, List[str]] = {}
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = list(dict.fromkeys(g.split(",")[0].strip() for g in parts[2:] if g.strip()))
    return sets


def library(method: str, genome: str, cache: Optional[str] = None) -> Tuple[Dict[str, List[str]], str, bool]:
    """Gene sets of a public library for a genome. Returns (sets, source, cross_species)."""
    key = canonical_method(method)
    if key not in LIBRARIES:
        raise ValueError(f"Unknown enrichment method '{method}'. Choose from: {', '.join(METHODS)}.")
    mouse = normalize_genome(genome).startswith("mm")
    human_file, mouse_file = LIBRARIES[key]
    name = mouse_file if mouse and mouse_file else human_file
    cross = mouse and not mouse_file
    if name.startswith("enrichr:"):
        lib = name.split(":", 1)[1]
        url, dest, source = ENRICHR + lib, f"enrichr_{lib}.gmt", f"Enrichr {lib}"
    else:
        url, dest, source = MSIGDB + name, name.split("/", 1)[1], f"MSigDB {name.split('/', 1)[1]}"
    path = download(url, cache_dir(cache) / "gene_sets" / dest, f"{key} gene sets")
    return read_gmt(path), source, cross


def hypergeom_sf(k: int, population: int, successes: int, draws: int) -> float:
    """P(X >= k) for a hypergeometric draw."""
    if k <= 0:
        return 1.0
    def logc(n, r):
        return math.lgamma(n + 1) - math.lgamma(r + 1) - math.lgamma(n - r + 1)
    denom = logc(population, draws)
    total = 0.0
    for i in range(k, min(successes, draws) + 1):
        if draws - i > population - successes:
            continue
        total += math.exp(logc(successes, i) + logc(population - successes, draws - i) - denom)
    return min(1.0, total)


def ora(gene_lists: Dict[str, List[str]], sets: Dict[str, List[str]], universe: List[str], method: str,
        min_size: int = 10, max_size: int = 500, annotated_background: bool = True) -> List[dict]:
    """Over-representation of each gene list in each set. Symbols are matched case-insensitively
    to the universe. With `annotated_background` the background is the universe genes that occur
    in the library, as clusterProfiler does; without it, the whole universe."""
    lookup = {g.upper(): g for g in universe}
    mapped = {}
    for term, genes in sets.items():
        hit = list(dict.fromkeys(lookup[g.upper()] for g in genes if g.upper() in lookup))
        if min_size <= len(hit) <= max_size:
            mapped[term] = set(hit)
    background = set().union(*mapped.values()) if (annotated_background and mapped) else set(universe)
    rows: List[dict] = []
    for group, genes in gene_lists.items():
        query = {lookup[g.upper()] for g in genes if g.upper() in lookup} & background
        block = []
        for term, members in mapped.items():
            hit = sorted(query & members)
            if not hit and annotated_background:
                continue
            block.append({"method": method, "group": group, "term_id": term, "term": term,
                          "gene_ratio": f"{len(hit)}/{len(query)}", "bg_ratio": f"{len(members)}/{len(background)}",
                          "pvalue": hypergeom_sf(len(hit), len(background), len(members), len(query)),
                          "count": len(hit), "genes": ";".join(hit)})
        for row, padj in zip(block, bh_adjust([r["pvalue"] for r in block])):
            row["padj"] = padj
        rows.extend(block)
    return rows


ENRICH_COLUMNS = ("method", "group", "term_id", "term", "gene_ratio", "bg_ratio", "pvalue", "padj", "count", "genes")


def read_gene_sets(path: str) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}
    for row in read_tsv(path):
        if row.get("SYMBOL"):
            groups.setdefault(row.get("Cluster") or "All", []).append(row["SYMBOL"])
    return {k: list(dict.fromkeys(v)) for k, v in groups.items()}


def run_enrichment(gene_sets: str, out_dir, prefix: str, genome: str = "mm10", methods=("go_bp",),
                   universe: Optional[str] = None, groups=None, padj: float = 0.05, top: int = 8,
                   custom_gmt: Optional[str] = None, cache: Optional[str] = None) -> dict:
    from .networks import run_pathway_networks
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    methods = [methods] if isinstance(methods, str) else list(methods)
    if not methods:
        raise ValueError("At least one pathway method is required.")
    gene_lists = read_gene_sets(gene_sets)
    if groups:
        gene_lists = {g: gene_lists[g] for g in ([groups] if isinstance(groups, str) else groups) if g in gene_lists}
    gene_lists = {g: v for g, v in gene_lists.items() if len(v) >= 5}
    if not gene_lists:
        raise ValueError("No gene group has at least 5 genes to test.")
    warnings: List[str] = []
    background = [r["SYMBOL"] for r in read_tsv(universe) if r.get("SYMBOL")] if universe else None
    if background is None:
        warnings.append("No gene universe given: each library's own genes are used as the background.")

    combined: List[dict] = []
    failed: Dict[str, str] = {}
    sources: Dict[str, str] = {}
    for method in methods:
        try:
            if canonical_method(method) == "custom":
                if not custom_gmt:
                    raise ValueError("method 'custom' needs a GMT file of your own signatures.")
                sets = read_gmt(custom_gmt)
                bg = background or sorted({g for genes in sets.values() for g in genes})
                rows = ora(gene_lists, sets, bg, method, min_size=3, max_size=10 ** 9, annotated_background=False)
                sources[method] = f"GMT {Path(custom_gmt).name}"
            else:
                sets, sources[method], cross = library(method, genome, cache)
                bg = background or sorted({g for genes in sets.values() for g in genes})
                rows = ora(gene_lists, sets, bg, method)
                if cross:
                    warnings.append(f"{method}: no mouse version exists; the human gene sets were matched to mouse "
                                    f"genes by symbol.")
            combined.extend(rows)
        except Exception as err:
            failed[method] = str(err)
            warnings.append(f"{method} failed: {err}")
    if len(failed) == len(methods):
        raise RuntimeError("All pathway methods failed: " + "; ".join(f"{m}: {e}" for m, e in failed.items()))
    table_path = write_tsv(out_dir / f"{prefix}_enrichment.tsv", ENRICH_COLUMNS,
                           [[r[c] for c in ENRICH_COLUMNS] for r in combined])

    figs: Dict[str, str] = {}
    sig = [r for r in combined if r["padj"] <= padj]
    for method in dict.fromkeys(r["method"] for r in sig):
        rows = [r for r in sig if r["method"] == method]
        try:
            figs[method] = figures.dotplot(rows, f"Pathway enrichment: {method}", out_dir,
                                           f"{prefix}_enrichment_{file_tag(method)}_dotplot", top)
        except Exception as err:
            warnings.append(f"Dot plot for {method} failed: {err}")
    networks = run_pathway_networks(combined, out_dir, prefix, padj) if combined else \
        {"figures": {}, "tables": {}, "modules": {}, "hub_genes": {}, "warnings": []}
    figs.update(networks["figures"])
    warnings.extend(networks["warnings"])

    top_terms, n_sig = [], []
    for method in dict.fromkeys(r["method"] for r in sig):
        for group in dict.fromkeys(r["group"] for r in sig if r["method"] == method):
            rows = sorted((r for r in sig if r["method"] == method and r["group"] == group), key=lambda r: r["padj"])
            n_sig.append({"method": method, "group": group, "N": len(rows)})
            top_terms.extend({"method": method, "group": group, "term": r["term"], "padj": signif(r["padj"]),
                              "count": r["count"], "genes": ";".join(r["genes"].split(";")[:12])} for r in rows[:top])
    return {
        "methods": methods, "failed_methods": failed, "gene_set_sources": sources,
        "groups": {g: len(v) for g, v in gene_lists.items()},
        "universe_size": len(background or []), "n_significant": n_sig, "top_terms": top_terms,
        "pathway_modules": networks["modules"], "hub_genes": networks["hub_genes"], "warnings": warnings,
        "tables": dict({"enrichment": table_path}, **networks["tables"]), "figures": figs,
        "versions": dict(sources),
    }
