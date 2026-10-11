"""Tool registry shared by the Claude agent loop, the deterministic pipeline, and
the MCP server. Each tool records its outcome in the EvidenceStore and returns
the Evidence object; the agent sees a compact JSON view of it."""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional

from . import analysis, literature, public_data
from .config import RunConfig
from .evidence import Evidence, EvidenceStore
from .rbridge import RBridge

RATIONALE = {"type": "string", "description": "One sentence: why this step, with these settings, for this "
                                              "dataset. Shown in the report's analysis log."}


@dataclass
class RunContext:
    config: RunConfig
    evidence: EvidenceStore = field(default_factory=EvidenceStore)
    bridge: Optional[RBridge] = None
    state: dict = field(default_factory=dict)
    report_path: Optional[str] = None
    log: Callable[[str], None] = print

    def __post_init__(self):
        Path(self.config.out_dir).mkdir(parents=True, exist_ok=True)
        if self.bridge is None:
            self.bridge = RBridge(self.config.rscript, log_dir=str(self.tables_dir / "logs"))

    @property
    def tables_dir(self) -> Path:
        return Path(self.config.out_dir) / "results"

    @property
    def peaks(self) -> str:
        """Peaks after QC cleaning when QC has run, otherwise the user's file."""
        return self.state.get("clean_peaks") or self.config.peaks

    def base_params(self, prefix_suffix: str = "") -> dict:
        return {"out_dir": str(self.tables_dir), "prefix": self.config.prefix + prefix_suffix,
                "genome": self.config.genome}


@dataclass
class Tool:
    name: str
    description: str
    properties: dict
    handler: Callable[[RunContext, dict], Evidence]
    required: List[str] = field(default_factory=list)
    needs_network: bool = False

    def schema(self) -> dict:
        props = dict(self.properties)
        props["rationale"] = RATIONALE
        return {"name": self.name, "description": self.description,
                "input_schema": {"type": "object", "properties": props, "required": self.required,
                                 "additionalProperties": False}}


def _record(ctx: RunContext, tool: str, title: str, params: dict, outcome: dict) -> Evidence:
    ev = ctx.evidence.add(tool, title, params, outcome, rationale=params.get("rationale", ""))
    status = "ok" if ev.status == "ok" else f"{ev.status}: {(ev.error or '')[:160]}"
    ctx.log(f"  [{ev.id}] {tool} -> {status} ({ev.elapsed_sec:.0f}s)")
    return ev


def _python_outcome(fn: Callable[[], dict]) -> dict:
    start = time.time()
    try:
        result = fn()
        return {"status": "ok", "result": result, "elapsed_sec": round(time.time() - start, 1)}
    except Exception as err:  # tool failures go back to the agent, not up the stack
        return {"status": "error", "error": f"{type(err).__name__}: {err}",
                "elapsed_sec": round(time.time() - start, 1)}


def _skip(reason: str) -> dict:
    return {"status": "skipped", "error": reason}


# ------------------------------------------------------------------ R-backed tools

def peak_qc(ctx: RunContext, args: dict) -> Evidence:
    params = dict(ctx.base_params(), peaks=ctx.config.peaks, blacklist=ctx.config.blacklist)
    outcome = ctx.bridge.run("peak_qc", params)
    if outcome.get("status") == "ok":
        ctx.state["clean_peaks"] = outcome["result"].get("clean_peaks")
    return _record(ctx, "peak_qc", "Peak QC", args, outcome)


def annotate_peaks(ctx: RunContext, args: dict) -> Evidence:
    params = dict(ctx.base_params(), peaks=ctx.peaks,
                  promoter_upstream=int(args.get("promoter_upstream", 2000)),
                  promoter_downstream=int(args.get("promoter_downstream", 500)),
                  distal_max_distance=int(args.get("distal_max_distance", 100000)),
                  chipseeker=bool(args.get("chipseeker_crosscheck", False)))
    outcome = ctx.bridge.run("annotate_peaks", params)
    if outcome.get("status") == "ok":
        res = outcome["result"]
        ctx.state.update(peak_table=res["tables"]["peak_genes"], gene_sets=res.get("gene_sets"),
                         universe=res.get("universe"))
    return _record(ctx, "annotate_peaks", "ChIPpeakAnno annotation", args, outcome)


def enrich_pathways(ctx: RunContext, args: dict) -> Evidence:
    if not ctx.state.get("gene_sets"):
        return _record(ctx, "enrich_pathways", "Pathway enrichment", args, _skip("Run annotate_peaks first."))
    methods = list(args.get("methods") or ["go_bp", "msigdb_H", "kegg"])
    if "custom" in methods and not ctx.config.biology.custom_gmt:
        methods.remove("custom")
    n = ctx.state["n_enrich"] = ctx.state.get("n_enrich", 0) + 1
    params = dict(ctx.base_params(f"_enrich{n}"), gene_sets=ctx.state["gene_sets"],
                  universe=ctx.state.get("universe"), methods=methods,
                  groups=args.get("groups") or ["promoter", "distal", "all"],
                  padj=float(args.get("padj", 0.05)), top=int(args.get("top", 8)),
                  custom_gmt=ctx.config.biology.custom_gmt)
    outcome = ctx.bridge.run("enrich_pathways", params)
    return _record(ctx, "enrich_pathways", "Pathway enrichment (" + ", ".join(methods) + ")", args, outcome)


def great_enrichment(ctx: RunContext, args: dict) -> Evidence:
    params = dict(ctx.base_params(), peaks=ctx.peaks, gene_sets=args.get("gene_sets", "GO:BP"),
                  padj=float(args.get("padj", 0.05)))
    outcome = ctx.bridge.run("great_enrichment", params)
    return _record(ctx, "great_enrichment", f"GREAT region enrichment ({params['gene_sets']})", args, outcome)


def _subset_bed(ctx: RunContext, regions: str) -> Optional[str]:
    if regions == "all":
        return ctx.peaks
    table = ctx.state.get("peak_table")
    if not table:
        return None
    keep = {"promoter": ("promoter",), "distal": ("distal", "distal_far")}[regions]
    rows = [r for r in analysis.read_table(table) if r.get("region_class") in keep]
    path = ctx.tables_dir / f"{ctx.config.prefix}_{regions}_peaks.bed"
    path.write_text("".join(f"{r['Chr']}\t{r['Start']}\t{r['End']}\t{r['PeakID']}\n" for r in rows))
    return str(path)


def find_motifs(ctx: RunContext, args: dict) -> Evidence:
    regions = args.get("regions", "all")
    bed = _subset_bed(ctx, regions)
    if bed is None:
        return _record(ctx, "find_motifs", f"Motif enrichment ({regions})", args,
                       _skip("Run annotate_peaks before motif analysis of a promoter/distal subset."))
    params = dict(ctx.base_params(f"_{regions}"), peaks=bed, backend=args.get("backend", "auto"),
                  size=int(args.get("size", 200)), max_peaks=int(args.get("max_peaks", 5000)),
                  padj=float(args.get("padj", 0.05)))
    outcome = ctx.bridge.run("find_motifs", params)
    if outcome.get("status") == "ok":
        ctx.state.setdefault("motif_tables", {})[regions] = outcome["result"]["tables"]["motifs"]
    return _record(ctx, "find_motifs", f"Known-motif enrichment ({regions} peaks)", args, outcome)


# ------------------------------------------------------------------ public data

def chipatlas_submit(ctx: RunContext, args: dict) -> Evidence:
    antigen = args.get("antigen_class", "TFs and others")
    cell = args.get("cell_class", "All cell types")

    def run():
        regions = args.get("regions", "all")
        bed = _subset_bed(ctx, regions) or ctx.peaks
        rid = public_data.submit_enrichment(bed, ctx.config.genome, antigen_class=antigen, cell_class=cell,
                                            threshold=int(args.get("threshold", 50)),
                                            distance=int(args.get("distance", 5000)),
                                            title=f"{ctx.config.prefix}_{regions}")
        ctx.state.setdefault("chipatlas_jobs", {})[rid] = {"antigen_class": antigen, "cell_class": cell,
                                                          "regions": regions, "submitted": time.time()}
        return {"request_id": rid, "antigen_class": antigen, "cell_class": cell, "regions": regions,
                "note": "Job queued on the DDBJ supercomputer (typically 5-30 min). Continue with other "
                        "analyses and call chipatlas_collect later."}

    return _record(ctx, "chipatlas_submit", f"ChIP-Atlas submission ({antigen}, {cell})", args,
                   _python_outcome(run))


def chipatlas_collect(ctx: RunContext, args: dict) -> Evidence:
    jobs = ctx.state.get("chipatlas_jobs", {})
    rid = args.get("request_id") or next((j for j in reversed(list(jobs)) if not jobs[j].get("collected")), None)
    if not rid:
        return _record(ctx, "chipatlas_collect", "ChIP-Atlas enrichment", args, _skip("No submitted job."))
    job = jobs.get(rid, {})
    wait_min = float(args.get("wait_minutes", 0))

    def run():
        deadline = time.time() + wait_min * 60
        status = public_data.job_status(rid)
        while status not in ("finished", "unavailable") and time.time() < deadline:
            time.sleep(30)
            status = public_data.job_status(rid)
        if status != "finished":
            return {"request_id": rid, "job_status": status, "ready": False,
                    "note": "Not finished yet; call chipatlas_collect again later (optionally with wait_minutes)."}
        rows = public_data.fetch_enrichment(rid)
        job["collected"] = True
        table = ctx.tables_dir / f"{ctx.config.prefix}_chipatlas_{rid[-6:]}.tsv"
        analysis.write_table(rows, str(table))
        by_antigen = public_data.summarize_by_antigen(rows)
        ctx.state["chipatlas_summary"] = by_antigen
        top_exp = sorted((r for r in rows if r["fold_enrichment"] > 1), key=lambda r: r["log10_q"])[:15]
        for r in top_exp:
            ctx.evidence.add_dataset(r["experiment"], {"antigen": r["antigen"], "cell": r["cell"],
                                                        "url": f"https://chip-atlas.org/view?id={r['experiment']}"})
        return {"request_id": rid, "ready": True, "n_experiments_tested": len(rows),
                "antigen_class": job.get("antigen_class"), "cell_class": job.get("cell_class"),
                "regions": job.get("regions"), "top_antigens": by_antigen,
                "top_experiments": [{k: r[k] for k in ("experiment", "antigen", "cell", "log10_q",
                                                        "fold_enrichment", "overlaps_user")} for r in top_exp],
                "tables": {"chipatlas": str(table)}}

    outcome = _python_outcome(run)
    ready = outcome.get("result", {}).get("ready")
    title = "ChIP-Atlas co-binding enrichment" if ready else "ChIP-Atlas status check"
    return _record(ctx, "chipatlas_collect" if ready else "chipatlas_status", title, args, outcome)


def chipatlas_target_genes(ctx: RunContext, args: dict) -> Evidence:
    tf = args["tf"]

    def run():
        if not ctx.state.get("peak_table"):
            raise RuntimeError("Run annotate_peaks first.")
        targets = public_data.target_genes(tf, ctx.config.genome, int(args.get("distance_kb", 5)))
        universe = [r["SYMBOL"] for r in analysis.read_table(ctx.state["universe"])] if ctx.state.get("universe") \
            else list(targets)
        out = analysis.target_overlap(targets, ctx.state["peak_table"], universe)
        out["tf"] = tf
        return out

    return _record(ctx, "chipatlas_target_genes", f"Overlap with public {tf} targets", args, _python_outcome(run))


def search_literature(ctx: RunContext, args: dict) -> Evidence:
    def run():
        papers = literature.search_europepmc(args["query"], int(args.get("max_results", 8)),
                                             args.get("from_year"), sort=args.get("sort", "cited"))
        for p in papers:
            if p.get("pmid"):
                ctx.evidence.add_reference(p["pmid"], p)
        return {"query": args["query"], "papers": papers,
                "note": "Cite papers as [PMID:n]. Only claims supported by these abstracts may be cited."}

    return _record(ctx, "search_literature", f"Literature: {args['query']}", args, _python_outcome(run))


def search_geo(ctx: RunContext, args: dict) -> Evidence:
    def run():
        series = literature.search_geo(args["query"], ctx.config.organism, int(args.get("max_results", 10)))
        for s in series:
            if s.get("accession"):
                ctx.evidence.add_dataset(s["accession"], s)
        return {"query": args["query"], "datasets": series, "note": "Cite datasets by accession, e.g. [GSE12345]."}

    return _record(ctx, "search_geo", f"GEO datasets: {args['query']}", args, _python_outcome(run))


# ------------------------------------------------------------------ integration / checks

def integrate_rnaseq(ctx: RunContext, args: dict) -> Evidence:
    def run():
        if not ctx.config.biology.rnaseq:
            raise RuntimeError("No RNA-seq table was provided for this run.")
        if not ctx.state.get("peak_table"):
            raise RuntimeError("Run annotate_peaks first.")
        out = str(ctx.tables_dir / f"{ctx.config.prefix}_rnaseq_integration.tsv")
        return analysis.integrate_rnaseq(ctx.config.biology.rnaseq, ctx.state["peak_table"], out,
                                         float(args.get("padj", 0.05)), float(args.get("lfc", 0.5)))

    return _record(ctx, "integrate_rnaseq", "Peak-RNA-seq integration", args, _python_outcome(run))


def check_expectations(ctx: RunContext, args: dict) -> Evidence:
    bio = ctx.config.biology
    genes = list(dict.fromkeys(bio.expected_genes + list(args.get("extra_genes") or [])))
    tfs = list(dict.fromkeys(bio.expected_tfs + list(args.get("extra_tfs") or [])))
    if bio.is_tf_assay and bio.target and bio.target not in tfs:
        tfs.insert(0, bio.target)

    def run():
        if not genes and not tfs:
            raise RuntimeError("No expected genes or TFs were given.")
        motif_tables = ctx.state.get("motif_tables", {})
        motif = motif_tables.get("all") or next(iter(motif_tables.values()), None)
        return analysis.check_expectations(genes, tfs, ctx.state.get("peak_table"), motif,
                                           ctx.state.get("chipatlas_summary"))

    return _record(ctx, "check_expectations", "Positive-control check", args, _python_outcome(run))


def inspect_table(ctx: RunContext, args: dict) -> Evidence:
    """Read rows from a table produced earlier in this run (never arbitrary paths)."""
    def run():
        ev = ctx.evidence.items.get(args["evidence_id"])
        if not ev or not ev.tables:
            raise RuntimeError(f"{args['evidence_id']} has no tables.")
        key = args.get("table") or next(iter(ev.tables))
        if key not in ev.tables:
            raise RuntimeError(f"Table '{key}' not in {list(ev.tables)}.")
        rows = analysis.read_table(ev.tables[key])
        col, val = args.get("filter_column"), args.get("filter_value")
        if col:
            rows = [r for r in rows if val is None or str(val).lower() in str(r.get(col, "")).lower()]
        sort_by = args.get("sort_by")
        if sort_by:
            def keyf(r):
                try:
                    return float(r.get(sort_by))
                except (TypeError, ValueError):
                    return float("inf")
            rows.sort(key=keyf, reverse=bool(args.get("descending")))
        n = max(1, min(int(args.get("n", 20)), 60))
        return {"source": ev.id, "table": key, "n_rows_matching": len(rows), "rows": rows[:n]}

    outcome = _python_outcome(run)
    ev = ctx.evidence.add("inspect_table", f"Inspect {args.get('evidence_id')}", args, outcome,
                          rationale=args.get("rationale", ""))
    ev.status = "lookup" if ev.status == "ok" else ev.status   # lookups are not new evidence
    return ev


# ------------------------------------------------------------------ registry

ENRICH_METHODS = ["go_bp", "go_mf", "go_cc", "kegg", "reactome", "msigdb_H", "msigdb_C2:CP:REACTOME",
                  "msigdb_C2:CP:WIKIPATHWAYS", "msigdb_C3:TFT:GTRD", "msigdb_C7:IMMUNESIGDB", "msigdb_C8", "custom"]

TOOLS: Dict[str, Tool] = {t.name: t for t in [
    Tool("peak_qc", "Peak-set QC: counts, widths, chromosomes, duplicates, non-standard contigs and blacklist "
         "overlap. Writes a cleaned peak set that later tools use automatically. Run first.", {}, peak_qc),
    Tool("annotate_peaks", "ChIPpeakAnno annotation: nearest gene by TSS for every peak, all genes whose TSS "
         "window overlaps the peak (bidirectional promoters), genomic-region distribution, TSS binding profile. "
         "Classifies peaks as promoter / distal / distal_far and builds gene sets (promoter, distal, all) plus a "
         "genome background for enrichment. Choose the promoter window to fit the assay (tight for "
         "promoter-bound TFs, standard 2000/500 otherwise).",
         {"promoter_upstream": {"type": "integer", "description": "bp upstream of TSS (default 2000)"},
          "promoter_downstream": {"type": "integer", "description": "bp downstream of TSS (default 500)"},
          "distal_max_distance": {"type": "integer",
                                  "description": "Distal peaks further than this from any TSS are not linked "
                                                 "to a gene (default 100000)"},
          "chipseeker_crosscheck": {"type": "boolean", "description": "Also run ChIPseeker for comparison"}},
         annotate_peaks),
    Tool("enrich_pathways", "Over-representation analysis (clusterProfiler) of peak-linked gene sets against the "
         "genome background. Methods: " + ", ".join(ENRICH_METHODS) + ". msigdb_<collection>[:<subcollection>] "
         "takes any MSigDB collection; 'custom' uses the user's GMT signatures. Pick collections that fit the "
         "biology (e.g. C7:IMMUNESIGDB for immune cells). Also draws pathway networks per method: an enrichment "
         "map (terms linked by shared genes, grouped into modules - returned as pathway_modules) and a "
         "gene-pathway network (hub genes shared by several top pathways - returned as hub_genes). Interpret "
         "modules rather than redundant individual terms.",
         {"methods": {"type": "array", "items": {"type": "string"}},
          "groups": {"type": "array", "items": {"type": "string", "enum": ["promoter", "distal", "all"]}},
          "padj": {"type": "number"}, "top": {"type": "integer", "description": "Top terms per group to return"}},
         enrich_pathways, required=["methods"]),
    Tool("great_enrichment", "Region-based GREAT enrichment (rGREAT, local) - better than gene ORA for distal "
         "enhancer peaks because it models each gene's regulatory domain.",
         {"gene_sets": {"type": "string", "description": "e.g. GO:BP, GO:MF, msigdb:H, msigdb:C2:CP:KEGG"},
          "padj": {"type": "number"}}, great_enrichment),
    Tool("find_motifs", "Known-motif enrichment at peak centres (HOMER if installed, else monaLisa + JASPAR "
         "CORE vs GC-matched genomic background). Use regions='distal' to find enhancer-driving TFs without "
         "promoter GC bias. Not meaningful for broad histone domains.",
         {"regions": {"type": "string", "enum": ["all", "promoter", "distal"]},
          "backend": {"type": "string", "enum": ["auto", "homer", "monalisa"]},
          "size": {"type": "integer", "description": "Window around peak centre, bp (default 200)"},
          "max_peaks": {"type": "integer", "description": "Cap on peaks scanned (strongest kept; default 5000)"},
          "padj": {"type": "number"}}, find_motifs),
    Tool("chipatlas_submit", "Submit the peaks to ChIP-Atlas Enrichment Analysis: overlap with every public "
         "ChIP-seq/ATAC-seq experiment for this genome, vs random permutation. Identifies TFs/marks that co-bind "
         "and the most similar public experiments. Asynchronous: submit early, collect later.",
         {"antigen_class": {"type": "string", "enum": list(public_data.ANTIGEN_CLASSES)},
          "cell_class": {"type": "string", "description": "ChIP-Atlas cell-type class, e.g. 'All cell types', "
                                                            "'Blood', 'Neural', 'Liver', 'Pluripotent stem cell'"},
          "regions": {"type": "string", "enum": ["all", "promoter", "distal"]},
          "threshold": {"type": "integer", "enum": [50, 100, 200, 500],
                        "description": "MACS2 -10log10(q) threshold for public peaks"},
          "distance": {"type": "integer"}}, chipatlas_submit, needs_network=True),
    Tool("chipatlas_collect", "Fetch results of a ChIP-Atlas job (summarised per antigen). Returns ready=false "
         "if still running; set wait_minutes to block.",
         {"request_id": {"type": "string"}, "wait_minutes": {"type": "number"}}, chipatlas_collect,
         needs_network=True),
    Tool("chipatlas_target_genes", "Known targets of a TF from public ChIP-seq (ChIP-Atlas, peaks within "
         "+/- 1/5/10 kb of TSS) and their overlap with this run's peak-linked genes (hypergeometric test).",
         {"tf": {"type": "string", "description": "Official symbol as in ChIP-Atlas, e.g. Tox, CTCF"},
          "distance_kb": {"type": "integer", "enum": [1, 5, 10]}}, chipatlas_target_genes, required=["tf"],
         needs_network=True),
    Tool("integrate_rnaseq", "Integrate the user's differential-expression table: are up/down-regulated genes "
         "enriched for peaks (promoter or any)? Lists bound DE genes as candidate direct targets.",
         {"padj": {"type": "number"}, "lfc": {"type": "number", "description": "abs log2FC cutoff"}},
         integrate_rnaseq),
    Tool("check_expectations", "Positive controls: are the user's expected genes bound, and are expected TFs "
         "(and the ChIPped TF itself) recovered by motifs / ChIP-Atlas? Run after motifs and ChIP-Atlas.",
         {"extra_genes": {"type": "array", "items": {"type": "string"}},
          "extra_tfs": {"type": "array", "items": {"type": "string"}}}, check_expectations),
    Tool("search_literature", "Search Europe PMC for papers (with abstracts). Use specific queries combining "
         "TFs/genes/pathways found here with the cell type and condition. Returns PMIDs you may cite.",
         {"query": {"type": "string"}, "max_results": {"type": "integer"},
          "from_year": {"type": "integer"},
          "sort": {"type": "string", "enum": ["cited", "recent"],
                   "description": "'cited' (default) finds landmark studies; 'recent' the newest work"}},
         search_literature, required=["query"], needs_network=True),
    Tool("search_geo", "Search NCBI GEO for related public datasets (series) for this organism.",
         {"query": {"type": "string"}, "max_results": {"type": "integer"}}, search_geo, required=["query"],
         needs_network=True),
    Tool("inspect_table", "Read rows of a table produced earlier (by evidence id), optionally filtered and "
         "sorted - e.g. check a specific gene's peaks or a motif's statistics.",
         {"evidence_id": {"type": "string"}, "table": {"type": "string"},
          "filter_column": {"type": "string"}, "filter_value": {"type": "string"},
          "sort_by": {"type": "string"}, "descending": {"type": "boolean"}, "n": {"type": "integer"}},
         inspect_table, required=["evidence_id"]),
]}

WRITE_REPORT_SCHEMA = {
    "name": "write_report",
    "description": "Finish the run: write the interpretation into the HTML report. Call exactly once, last. "
                   "Every factual claim must cite evidence ids [E3], papers [PMID:123] or datasets [GSE123] "
                   "returned by tools in this run. Use markdown (paragraphs, '- ' bullets, **bold**).",
    "input_schema": {
        "type": "object",
        "properties": {
            "title": {"type": "string"},
            "executive_summary": {"type": "string", "description": "5-8 sentences: what the peaks are, "
                                  "what they regulate, which TFs drive them, how this fits the known biology."},
            "sections": {"type": "array", "items": {
                "type": "object",
                "properties": {"heading": {"type": "string"}, "body": {"type": "string"}},
                "required": ["heading", "body"], "additionalProperties": False}},
            "hypotheses": {"type": "array", "items": {
                "type": "object",
                "properties": {"statement": {"type": "string"},
                               "evidence": {"type": "array", "items": {"type": "string"}},
                               "confidence": {"type": "string", "enum": ["high", "medium", "low"]},
                               "next_experiment": {"type": "string"}},
                "required": ["statement", "evidence", "confidence", "next_experiment"],
                "additionalProperties": False}},
            "caveats": {"type": "array", "items": {"type": "string"}},
        },
        "required": ["title", "executive_summary", "sections", "hypotheses", "caveats"],
        "additionalProperties": False,
    },
}


def tool_schemas(offline: bool) -> List[dict]:
    return [t.schema() for t in TOOLS.values() if not (offline and t.needs_network)] + [WRITE_REPORT_SCHEMA]


def run_tool(ctx: RunContext, name: str, args: dict) -> Evidence:
    tool = TOOLS.get(name)
    if tool is None:
        raise KeyError(f"Unknown tool {name}")
    if tool.needs_network and ctx.config.offline:
        return _record(ctx, name, name, args, _skip("Offline run: network tools are disabled."))
    missing = [k for k in tool.required if k not in args]
    if missing:
        return _record(ctx, name, name, args, {"status": "error", "error": f"Missing argument(s): {missing}"})
    return tool.handler(ctx, args)
