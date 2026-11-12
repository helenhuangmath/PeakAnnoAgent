"""Optional companion tools, run as subprocesses and folded into the evidence store:

- SimpNetworkViz (`snv`): EnrichmentMap/AutoAnnotate-style pathway and gene networks drawn
  from this run's enrichment table.
- PathwayAgent (`pathwayagent`): several enrichment engines on the peak-linked gene lists,
  combined into cross-tool consensus themes and key genes.

Both are found on PATH, through an environment variable, or in a sibling checkout of the
repository. A tool that is not installed is reported as skipped, never as an error."""

from __future__ import annotations

import csv
import json
import os
import re
import shutil
import subprocess
import shlex
from pathlib import Path
from typing import Dict, List, Optional

from .analysis import read_table

REPO = Path(__file__).resolve().parents[2]
SIBLINGS = [REPO.parent, REPO.parent.parent]        # the checkout may sit one level below its siblings
ORGANISM = {"mm10": "mmusculus", "mm39": "mmusculus", "hg38": "hsapiens", "hg19": "hsapiens"}
SNV_INSTALL = "pip install git+https://github.com/helenhuangmath/SimpNetworkViz"
PA_INSTALL = ("clone https://github.com/helenhuangmath/PathwayAgent next to this repository and run its setup, "
              "or set PEAKANNOAGENT_PATHWAYAGENT to its `pathwayagent` command")


# ------------------------------------------------------------------ discovery

def _split_cmd(value: str) -> List[str]:
    return shlex.split(value)


def find_snv() -> Optional[List[str]]:
    env = os.environ.get("PEAKANNOAGENT_SNV")
    if env:
        return _split_cmd(env)
    exe = shutil.which("snv") or shutil.which("simpnetworkviz")
    if exe:
        return [exe]
    for root in SIBLINGS:
        cand = root / "SimpNetworkViz" / ".venv" / "bin" / "snv"
        if cand.is_file():
            return [str(cand)]
    return None


def find_pathwayagent() -> Optional[List[str]]:
    env = os.environ.get("PEAKANNOAGENT_PATHWAYAGENT")
    if env:
        return _split_cmd(env)
    exe = shutil.which("pathwayagent")
    if exe:
        return [exe]
    for root in SIBLINGS:
        py = root / "PathwayAgent" / ".venv" / "bin" / "python"
        if py.is_file():
            return [str(py), "-m", "pathwayagent"]
    return None


def available() -> Dict[str, Optional[str]]:
    snv, pa = find_snv(), find_pathwayagent()
    return {"snv": " ".join(snv) if snv else None, "pathwayagent": " ".join(pa) if pa else None}


def _run(cmd: List[str], cwd: Optional[str], timeout: int, log_path: Path) -> subprocess.CompletedProcess:
    env = dict(os.environ, MPLBACKEND="Agg")
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout, env=env)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(f"$ {' '.join(cmd)}\n\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout).strip().splitlines()[-6:]
        raise RuntimeError(f"{Path(cmd[0]).name} exited with {proc.returncode}: " + " | ".join(tail))
    return proc


def _num(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# ------------------------------------------------------------------ SimpNetworkViz

def enrichment_to_clusterprofiler(enrichment_tsv: str, method: str, group: str, out_csv: str) -> int:
    """Rows of the engine's enrichment table for one method/group, in clusterProfiler layout."""
    rows = [r for r in read_table(enrichment_tsv) if r.get("method") == method and r.get("group") == group]
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ID", "Description", "pvalue", "p.adjust", "geneID", "Count"])
        for r in rows:
            w.writerow([r.get("term_id") or r.get("term"), r.get("term"), r.get("pvalue"), r.get("padj"),
                        (r.get("genes") or "").replace(";", "/"), r.get("count")])
    return len(rows)


def enrichment_methods(enrichment_tsv: str) -> Dict[str, List[str]]:
    """method -> groups present in the table (significant rows only, as the engine writes them)."""
    out: Dict[str, List[str]] = {}
    for r in read_table(enrichment_tsv):
        out.setdefault(r.get("method", ""), [])
        if r.get("group") not in out[r["method"]]:
            out[r["method"]].append(r.get("group"))
    return out


def run_snv(enrichment_tsv: str, out_dir: str, method: str, group: str, gmt: Optional[str] = None,
            max_pathways: int = 40, gene_network: str = "cnet", cutoff: Optional[float] = None, title: str = "",
            timeout: int = 600) -> dict:
    """Similarity between terms is computed from their hit genes unless a GMT with the full sets is given
    (the user's custom signatures). Hit genes overlap far more than full sets, so the default edge cutoff
    is 0.6 in that case (0.375, the EnrichmentMap default, with full sets); otherwise dense GO results
    turn into a grey mass of edges."""
    cmd = find_snv()
    if not cmd:
        raise RuntimeError(f"SimpNetworkViz is not installed ({SNV_INSTALL}).")
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    table = out / "enrichment_input.csv"
    n = enrichment_to_clusterprofiler(enrichment_tsv, method, group, str(table))
    if n == 0:
        raise RuntimeError(f"No enriched terms for method '{method}' and group '{group}'.")
    args = cmd + ["run", "--enrichment", str(table), "-o", str(out), "--max-pathways", str(int(max_pathways)),
                  "--gene-network", gene_network, "--formats", "png,pdf,html,csv,cyjs",
                  "--title", title or f"{method}, {group} peaks"]
    full = bool(gmt and method == "custom")
    if cutoff is None:
        cutoff = 0.375 if full else 0.6
    args += ["--cutoff", str(cutoff)]
    if full:
        args += ["--gmt", gmt, "--similarity-genes", "full"]
    else:
        args += ["--similarity-genes", "hits"]
    _run(args, None, timeout, out / "snv.log")

    def rows(name):
        p = out / name
        return read_table(str(p)) if p.is_file() else []

    nodes, edges, clusters = rows("pathway_network_nodes.csv"), rows("pathway_network_edges.csv"), rows("pathway_network_clusters.csv")
    gnodes, gedges = rows("gene_network_nodes.csv"), rows("gene_network_edges.csv")
    clusters_out = [{"cluster": c.get("cluster"), "label": c.get("label"), "n_pathways": c.get("n_nodes"),
                     "best_fdr": c.get("best_fdr"),
                     "members": "; ".join((c.get("members") or "").split("; ")[:5])
                     + (" …" if len((c.get("members") or "").split("; ")) > 5 else "")} for c in clusters]
    figures = {k: str(out / f"{k}.png") for k in ("pathway_network", "gene_network", "cluster_barplot")
               if (out / f"{k}.png").is_file()}
    tables = {k: str(out / f) for k, f in (("clusters", "pathway_network_clusters.csv"),
                                          ("pathway_nodes", "pathway_network_nodes.csv"),
                                          ("pathway_edges", "pathway_network_edges.csv"),
                                          ("gene_nodes", "gene_network_nodes.csv")) if (out / f).is_file()}
    interactive = {k: str(out / f"{k}.html") for k in ("pathway_network", "gene_network") if (out / f"{k}.html").is_file()}
    version = None
    params = out / "run_params.json"
    if params.is_file():
        try:
            version = json.loads(params.read_text()).get("simpnetworkviz")
        except ValueError:
            pass
    return {"method": method, "group": group, "n_input_terms": n, "similarity_cutoff": cutoff,
            "similarity_genes": "full gene sets" if full else "hit genes", "n_pathways": len(nodes), "n_edges": len(edges),
            "n_clusters": len(clusters), "clusters": clusters_out,
            "gene_network": {"type": gene_network, "n_nodes": len(gnodes), "n_edges": len(gedges)} if gnodes else None,
            "interactive": interactive, "cytoscape": str(out / "pathway_network.cyjs"),
            "figures": figures, "tables": tables, "versions": {"SimpNetworkViz": version} if version else {},
            "warnings": [] if clusters else ["No pathway clusters were formed: too few similar pathways passed the FDR cut."]}


# ------------------------------------------------------------------ PathwayAgent

def write_gene_lists(gene_sets_tsv: str, out_dir: str, lists: List[str]) -> Dict[str, str]:
    """One file per list (promoter / distal / all) from the engine's gene-set table."""
    groups: Dict[str, List[str]] = {}
    for r in read_table(gene_sets_tsv):
        groups.setdefault(r.get("Cluster", ""), []).append(r.get("SYMBOL", ""))
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    files = {}
    for name in lists:
        genes = [g for g in groups.get(name, []) if g]
        if genes:
            p = out / f"{name}.txt"
            p.write_text("\n".join(dict.fromkeys(genes)) + "\n")
            files[name] = str(p)
    return files


def run_pathwayagent(gene_sets_tsv: str, universe_tsv: Optional[str], out_dir: str, genome: str,
                     lists: List[str] = ("promoter", "distal", "all"), preset: str = "quick",
                     local_only: bool = True, title: str = "", timeout: int = 1800) -> dict:
    cmd = find_pathwayagent()
    if not cmd:
        raise RuntimeError(f"PathwayAgent is not installed ({PA_INSTALL}).")
    out = Path(out_dir)
    files = write_gene_lists(gene_sets_tsv, str(out / "input"), list(lists))
    if not files:
        raise RuntimeError("No gene lists to analyse (run annotate_peaks first).")
    args = cmd + ["run", *files.values(), "--organism", ORGANISM.get(genome, "auto"), "--preset", preset,
                  "--out", str(out), "--title", title or "PeakAnnoAgent peak-linked genes"]
    universe = None
    if universe_tsv and Path(universe_tsv).is_file():
        universe = out / "input" / "universe.txt"
        universe.write_text("\n".join(r.get("SYMBOL", "") for r in read_table(universe_tsv) if r.get("SYMBOL")) + "\n")
        args += ["--universe", str(universe)]
    if local_only:
        args.append("--local-only")
    cwd = str(Path(cmd[0]).resolve().parents[2]) if cmd[0].endswith("python") and "-m" in cmd else None
    _run(args, cwd, timeout, out / "pathwayagent.log")

    result = {"organism": ORGANISM.get(genome, "auto"), "preset": preset, "local_only": local_only,
              "universe_size": sum(1 for _ in open(universe)) if universe else None,
              "lists": {}, "figures": {}, "tables": {}, "warnings": []}
    manifest = out / "manifest.json"
    if manifest.is_file():
        try:
            m = json.loads(manifest.read_text())
            result["versions"] = {"PathwayAgent": m.get("pathwayagent_version")}
            result["warnings"] += [str(w) for w in (m.get("warnings") or [])][:6]
        except ValueError:
            pass
    if (out / "report.html").is_file():
        result["interactive"] = {"pathwayagent_report": str(out / "report.html")}
    tools_ok = set()
    for name in files:
        d = out / name
        if not d.is_dir():
            continue
        themes = read_table(str(d / "themes.tsv")) if (d / "themes.tsv").is_file() else []
        terms = read_table(str(d / "consensus_terms.tsv")) if (d / "consensus_terms.tsv").is_file() else []
        key = read_table(str(d / "key_genes.tsv")) if (d / "key_genes.tsv").is_file() else []
        status = read_table(str(d / "tool_status.tsv")) if (d / "tool_status.tsv").is_file() else []
        tools_ok |= {s["tool"] for s in status if s.get("status") == "ok"}
        tiers = {}
        for t in terms:
            tiers[t.get("tier_label", "")] = tiers.get(t.get("tier_label", ""), 0) + 1
        themes.sort(key=lambda t: -_num(t.get("theme_score"), 0))
        result["lists"][name] = {
            "n_genes": sum(1 for _ in open(files[name])),
            "tools": {s["tool"]: s.get("status") for s in status},
            "n_terms_by_tier": tiers,
            "top_themes": [{"theme": t.get("theme_name"), "tier": t.get("tier_label"), "score": _num(t.get("theme_score")),
                            "n_terms": t.get("n_terms"), "databases": t.get("databases"), "best_padj": _num(t.get("best_padj")),
                            "core_genes": ";".join((t.get("core_genes") or "").replace(";", "|").split("|")[:8])}
                           for t in themes[:8]],
            "key_genes": [{"gene": k.get("gene"), "score": _num(k.get("gene_score")), "n_themes": k.get("n_themes"),
                           "driver_of": (k.get("driver_of") or "")[:80]} for k in key[:15]],
        }
        for fig, label in (("01_consensus_bar", "consensus_terms"), ("02_theme_network", "theme_network"),
                           ("05_key_gene_bubble", "key_genes"), ("06_gene_theme_heatmap", "gene_theme_heatmap")):
            p = d / "figures" / f"{fig}.png"
            if p.is_file():
                result["figures"][f"{name}_{label}"] = str(p)
        for tbl in ("themes", "consensus_terms", "key_genes", "tool_status"):
            p = d / f"{tbl}.tsv"
            if p.is_file():
                result["tables"][f"{name}_{tbl}"] = str(p)
    result["tools_ok"] = sorted(tools_ok)
    if len(tools_ok) < 2:
        result["warnings"].append("Fewer than two engines returned results, so no cross-tool consensus was possible; "
                                  "themes are single-tool and exploratory. Use preset 'standard' with network access "
                                  "for a real consensus.")
    return result
