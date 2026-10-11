"""Analyses done in Python on the engine's output tables: RNA-seq integration,
positive-control checks, and overlap with known TF target genes."""

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set


def read_table(path: str) -> List[dict]:
    with open(path, newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        delim = "\t" if sample.count("\t") >= sample.count(",") else ","
        return list(csv.DictReader(handle, delimiter=delim))


def write_table(rows: List[dict], path: str) -> str:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else ["empty"]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def hypergeom_sf(k: int, population: int, successes: int, draws: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(population, successes, draws)."""
    if k <= 0:
        return 1.0
    hi = min(successes, draws)
    if k > hi:
        return 0.0

    def logc(n, r):
        return math.lgamma(n + 1) - math.lgamma(r + 1) - math.lgamma(n - r + 1)

    denom = logc(population, draws)
    terms = [logc(successes, i) + logc(population - successes, draws - i) - denom for i in range(k, hi + 1)]
    top = max(terms)
    return min(1.0, math.exp(top) * sum(math.exp(t - top) for t in terms))


def _pick(columns: Iterable[str], candidates: Iterable[str]) -> Optional[str]:
    lower = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


# ---------------------------------------------------------------- peak -> gene links

def peak_gene_links(peak_table: str) -> Dict[str, dict]:
    """gene -> {'n_peaks', 'classes', 'min_abs_distance'} from the ChIPpeakAnno peak table."""
    genes: Dict[str, dict] = defaultdict(lambda: {"n_peaks": 0, "classes": set(), "min_abs_distance": None})
    for row in read_table(peak_table):
        linked = set()
        if row.get("promoter_symbols") and row["promoter_symbols"] not in ("NA", ""):
            for g in row["promoter_symbols"].split(";"):
                linked.add((g, "promoter", 0))
        if row.get("nearest_symbol") and row["nearest_symbol"] not in ("NA", "") and row.get("region_class") != "distal_far":
            try:
                dist = abs(int(float(row.get("distance_to_tss") or 0)))
            except ValueError:
                dist = None
            linked.add((row["nearest_symbol"], row.get("region_class", "distal"), dist))
        for gene, cls, dist in linked:
            info = genes[gene]
            info["n_peaks"] += 1
            info["classes"].add(cls)
            if dist is not None and (info["min_abs_distance"] is None or dist < info["min_abs_distance"]):
                info["min_abs_distance"] = dist
    return dict(genes)


def peak_genes_by_class(peak_table: str) -> Dict[str, Set[str]]:
    links = peak_gene_links(peak_table)
    return {
        "promoter": {g for g, i in links.items() if "promoter" in i["classes"]},
        "any": set(links),
    }


# ---------------------------------------------------------------- RNA-seq integration

def integrate_rnaseq(de_path: str, peak_table: str, out_path: str, padj_cutoff: float = 0.05,
                     lfc_cutoff: float = 0.5) -> dict:
    rows = read_table(de_path)
    if not rows:
        raise ValueError("The RNA-seq table is empty.")
    cols = rows[0].keys()
    gene_col = _pick(cols, ["SYMBOL", "gene", "gene_name", "gene_symbol", "Gene", "symbol", "external_gene_name", ""])
    lfc_col = _pick(cols, ["log2FoldChange", "logFC", "log2FC", "lfc", "avg_log2FC"])
    padj_col = _pick(cols, ["padj", "adj.P.Val", "FDR", "qvalue", "p_val_adj", "q_value"])
    if not (gene_col is not None and lfc_col and padj_col):
        raise ValueError(f"Could not find gene / log2FC / padj columns in {list(cols)}.")

    links = peak_gene_links(peak_table)
    tested, up, down = set(), set(), set()
    lfc, padj = {}, {}
    for r in rows:
        g = (r.get(gene_col) or "").strip()
        try:
            fc, q = float(r[lfc_col]), float(r[padj_col])
        except (TypeError, ValueError):
            continue
        if not g or math.isnan(fc) or math.isnan(q):
            continue
        tested.add(g)
        lfc[g], padj[g] = fc, q
        if q <= padj_cutoff and fc >= lfc_cutoff:
            up.add(g)
        elif q <= padj_cutoff and fc <= -lfc_cutoff:
            down.add(g)

    def stats(de: Set[str], bound: Set[str]) -> dict:
        bound_tested = bound & tested
        k = len(de & bound_tested)
        n = len(de)
        bg = len(bound_tested) / len(tested) if tested else 0
        return {"n_de": n, "n_with_peak": k, "pct_with_peak": round(100 * k / n, 1) if n else 0.0,
                "background_pct": round(100 * bg, 1),
                "fold": round((k / n) / bg, 2) if n and bg else None,
                "pvalue": hypergeom_sf(k, len(tested), len(bound_tested), n) if n else None}

    promoter = {g for g, i in links.items() if "promoter" in i["classes"]}
    anylink = set(links)
    result = {
        "n_tested_genes": len(tested), "n_up": len(up), "n_down": len(down),
        "thresholds": {"padj": padj_cutoff, "abs_log2fc": lfc_cutoff},
        "up_any_peak": stats(up, anylink), "down_any_peak": stats(down, anylink),
        "up_promoter_peak": stats(up, promoter), "down_promoter_peak": stats(down, promoter),
    }
    table = []
    for g in sorted((up | down) & anylink, key=lambda x: padj[x]):
        info = links[g]
        table.append({"gene": g, "direction": "up" if g in up else "down", "log2fc": round(lfc[g], 3),
                      "padj": f"{padj[g]:.3g}", "n_peaks": info["n_peaks"],
                      "classes": ";".join(sorted(info["classes"])), "min_abs_distance": info["min_abs_distance"]})
    write_table(table, out_path)
    result["top_bound_de_genes"] = [f"{t['gene']}({t['direction']},{t['classes']})" for t in table[:30]]
    warnings = []
    if len(tested) and not (tested & anylink):
        warnings.append("No overlap between RNA-seq genes and peak genes: check that both use gene symbols "
                        "of the same species.")
    result["warnings"] = warnings
    result["tables"] = {"rnaseq_integration": out_path}
    return result


# ---------------------------------------------------------------- positive controls

def _tf_matches(query: str, name: str) -> bool:
    q = query.lower()
    parts = [p.strip().lower() for p in name.replace("::", ",").replace("/", ",").split(",")]
    return q in parts or name.lower() == q


def check_expectations(expected_genes: List[str], expected_tfs: List[str], peak_table: Optional[str],
                       motif_table: Optional[str], chipatlas_summary: Optional[List[dict]],
                       padj_cutoff: float = 0.05) -> dict:
    out = {"genes": [], "tfs": []}
    links = peak_gene_links(peak_table) if peak_table else {}
    lower_links = {g.lower(): g for g in links}
    for gene in expected_genes:
        key = lower_links.get(gene.lower())
        if key:
            info = links[key]
            out["genes"].append({"gene": gene, "bound": True, "n_peaks": info["n_peaks"],
                                 "classes": sorted(info["classes"]), "min_abs_distance": info["min_abs_distance"]})
        else:
            out["genes"].append({"gene": gene, "bound": False})

    motifs = read_table(motif_table) if motif_table else []
    sig = [m for m in motifs if float(m.get("padj") or 1) <= padj_cutoff and float(m.get("log2_enrichment") or 0) > 0]
    for tf in expected_tfs:
        entry = {"tf": tf}
        if motif_table:
            hits = [i for i, m in enumerate(sig) if _tf_matches(tf, m.get("tf", ""))]
            entry["motif_rank"] = hits[0] + 1 if hits else None
            entry["motif_significant"] = bool(hits)
            entry["n_significant_motifs"] = len(sig)
        if chipatlas_summary is not None:
            ranks = [i for i, s in enumerate(chipatlas_summary) if s["antigen"].lower() == tf.lower()]
            entry["chipatlas_rank"] = ranks[0] + 1 if ranks else None
        out["tfs"].append(entry)

    n_bound = sum(g["bound"] for g in out["genes"])
    out["genes_recovered"] = f"{n_bound}/{len(out['genes'])}"
    found_tf = sum(1 for t in out["tfs"] if t.get("motif_significant") or t.get("chipatlas_rank"))
    out["tfs_recovered"] = f"{found_tf}/{len(out['tfs'])}"
    warnings = []
    if out["genes"] and n_bound / len(out["genes"]) < 0.3:
        warnings.append("Fewer than 30% of expected genes have a linked peak: check genome build, peak "
                        "calling, or whether the expectations fit this assay.")
    if out["tfs"] and found_tf == 0:
        warnings.append("None of the expected TFs were recovered by motifs or public-data overlap.")
    out["warnings"] = warnings
    return out


def target_overlap(targets: Dict[str, float], peak_table: str, universe: Iterable[str]) -> dict:
    """Do the peak-linked genes overlap known ChIP-Atlas targets of a TF more than expected?"""
    universe = set(universe)
    groups = peak_genes_by_class(peak_table)
    known = set(targets) & universe
    out = {"n_known_targets": len(known), "universe_size": len(universe)}
    for name, genes in groups.items():
        genes = genes & universe
        k = len(genes & known)
        out[name] = {"n_peak_genes": len(genes), "n_known_targets": k,
                     "pct": round(100 * k / len(genes), 1) if genes else 0.0,
                     "background_pct": round(100 * len(known) / len(universe), 1) if universe else 0.0,
                     "pvalue": hypergeom_sf(k, len(universe), len(known), len(genes)) if genes else None}
    top = sorted((g for g in groups["any"] if g in targets), key=lambda g: -targets[g])[:25]
    out["strongest_shared_targets"] = top
    return out
