"""Command line: `peakanno-agent run|check|tools|web`."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import os
import shutil
import sys

from .config import DEFAULT_MODEL, anthropic_credentials, SUPPORTED_GENOMES, RunConfig, build_biology, load_biology_file
from .engine import Engine, REFERENCE_OPTIONS

PALETTE_HELP = ("Figure colours: chippeakanno (default), npg (Nature style), spectral, or a Matplotlib colormap name "
                "(e.g. Set2, Dark2, tab10, viridis, RdYlBu_r)")
from .tools import TOOLS


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="peakanno-agent", description=__doc__)
    from . import __version__
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Annotate peaks, find pathways/motifs/regulators, write an HTML report.")
    run.add_argument("--peaks", required=True, help="BED (or TSV with Chr/Start/End) peak file")
    run.add_argument("--genome", required=True, help=f"{', '.join(SUPPORTED_GENOMES)} (or mouse/human)")
    run.add_argument("--out", required=True, help="Output directory")
    run.add_argument("--prefix", default="run", help="Prefix for output files")
    bio = run.add_argument_group("known biology")
    bio.add_argument("--biology", help="Free-text description of the experiment and what is known")
    bio.add_argument("--biology-file", help="JSON with biology fields, or a text/markdown notes file")
    bio.add_argument("--assay", help='e.g. "TF ChIP-seq", "ATAC-seq", "H3K27ac CUT&Tag"')
    bio.add_argument("--target", help="ChIPped TF or histone mark, e.g. Tox or H3K27ac")
    bio.add_argument("--cell-type")
    bio.add_argument("--condition", help="Treatment or comparison")
    bio.add_argument("--expected-genes", help="Comma-separated genes expected to be bound (positive controls)")
    bio.add_argument("--expected-tfs", help="Comma-separated TFs expected to drive these peaks")
    bio.add_argument("--gmt", help="Custom gene signatures (GMT) to test")
    bio.add_argument("--rnaseq", help="Differential expression table (gene, log2FC, padj columns)")
    run.add_argument("--blacklist", help="Blacklist BED to remove (e.g. ENCODE mm10 blacklist)")
    run.add_argument("--reference", default="gencode",
                     help="Gene models for annotation: gencode (default), gencode:<release> "
                          "(e.g. gencode:v44, gencode:vM25), refseq or ucsc")
    run.add_argument("--gtf", help="Annotate against this GTF file instead (GENCODE, Ensembl, RefSeq or custom)")
    run.add_argument("--palette", default="chippeakanno", help=PALETTE_HELP)
    run.add_argument("--offline", action="store_true", help="No ChIP-Atlas / literature / GEO calls. Gene models and gene sets are still "
                          "downloaded once if they are not in the cache")
    run.add_argument("--no-llm", action="store_true", help="Fixed analysis plan, no Claude interpretation")
    run.add_argument("--model", default=DEFAULT_MODEL)
    run.add_argument("--effort", default="high", choices=["low", "medium", "high", "xhigh", "max"])
    run.add_argument("--max-turns", type=int, default=40)
    run.add_argument("--chipatlas-wait", type=float, default=25.0,
                     help="Minutes the deterministic plan waits for ChIP-Atlas (default 25)")

    ann = sub.add_parser("annotate", help="Annotate peaks only: features, nearest genes, tables and figures. "
                                          "Fast; no agent, no enrichment.")
    ann.add_argument("--peaks", required=True, help="BED (or TSV with Chr/Start/End) peak file")
    ann.add_argument("--genome", required=True, help=f"{', '.join(SUPPORTED_GENOMES)} (or mouse/human)")
    ann.add_argument("--out", required=True, help="Output directory")
    ann.add_argument("--prefix", default="run", help="Prefix for output files")
    ann.add_argument("--reference", default="gencode",
                     help="Gene models: " + "; ".join(f"{k}: {v}" for k, v in REFERENCE_OPTIONS.items()))
    ann.add_argument("--gtf", help="Annotate against this GTF file instead")
    ann.add_argument("--promoter-upstream", type=int, default=2000, help="Promoter window upstream of a TSS, bp")
    ann.add_argument("--promoter-downstream", type=int, default=500, help="Promoter window downstream of a TSS, bp")
    ann.add_argument("--distal-max-distance", type=int, default=100000,
                     help="Peaks further than this from any TSS are not linked to a gene, bp")
    ann.add_argument("--gene-types", default="auto", choices=["auto", "protein_coding", "all"],
                     help="Genes to annotate against: protein-coding (default) or all")
    ann.add_argument("--palette", default="chippeakanno", help=PALETTE_HELP)
    ann.add_argument("--no-figures", action="store_true", help="Write tables only")
    ann.add_argument("--no-download", action="store_true", help="Use cached gene models only; never download")

    chk = sub.add_parser("check", help="Check the installation: packages, HOMER, cache, credentials, network.")
    chk.add_argument("--genome", default="mm10")
    chk.add_argument("--offline", action="store_true", help="Check local dependencies without contacting public services.")

    sub.add_parser("tools", help="List the agent's tools.")

    web = sub.add_parser("web", help="Start the web app: submit runs, watch progress, search results, chat.")
    web.add_argument("--runs-dir", default="runs", help="Where runs are written (default ./runs)")
    web.add_argument("--host", default="127.0.0.1")
    web.add_argument("--port", type=int, default=8766)
    web.add_argument("--model", default=DEFAULT_MODEL, help="Default model for agent runs and the chat")
    web.add_argument("--open", action="store_true", help="Open the browser")
    web.add_argument("--debug", action="store_true")
    return p


def cmd_run(args) -> int:
    file_fields = load_biology_file(args.biology_file) if args.biology_file else {}
    biology = build_biology(file_fields, description=args.biology, assay=args.assay, target=args.target,
                            cell_type=args.cell_type, condition=args.condition,
                            expected_genes=args.expected_genes, expected_tfs=args.expected_tfs,
                            custom_gmt=args.gmt, rnaseq=args.rnaseq)
    config = RunConfig(peaks=args.peaks, genome=args.genome, out_dir=args.out, prefix=args.prefix,
                       biology=biology, blacklist=args.blacklist, reference=args.reference, gtf=args.gtf,
                       palette=args.palette, offline=args.offline, use_llm=not args.no_llm,
                       model=args.model, effort=args.effort, max_turns=args.max_turns,
                       chipatlas_wait_min=args.chipatlas_wait)
    if config.use_llm:
        try:
            import anthropic  # noqa: F401
        except ImportError:
            print("Agent mode needs the Anthropic SDK: pip install 'peakannoagent[agent]' (or use --no-llm).",
                  file=sys.stderr)
            return 2
    from .agent import execute
    manifest = execute(config)
    print(manifest["report"])
    print(f"Analysis status: {manifest['status']}", file=sys.stderr)
    return {"success": 0, "partial": 3, "failed": 1}[manifest["status"]]


def cmd_annotate(args) -> int:
    from . import engine
    from .engine import common
    common.OFFLINE = bool(args.no_download)
    palette = engine.resolve_palette(args.palette)
    result = engine.run_annotation(engine.read_peak_file(args.peaks), args.out, args.prefix, genome=args.genome,
                                   reference=args.reference, gtf=args.gtf,
                                   promoter_upstream=args.promoter_upstream,
                                   promoter_downstream=args.promoter_downstream,
                                   distal_max_distance=args.distal_max_distance, gene_types=args.gene_types,
                                   draw_figures=not args.no_figures, palette=palette)
    summary = Path(args.out) / f"{args.prefix}_annotation.json"
    summary.write_text(json.dumps(result, indent=1))
    ref, t = result["reference"], result["timing_sec"]
    print(f"{result['n_peaks']:,} peaks annotated against {ref['source']} ({ref['genome']}, "
          f"{ref['genes_used']:,} {ref['gene_types'].replace('_', '-')} genes) in {sum(t.values()):.1f} s",
          file=sys.stderr)
    for feature, pct in result["feature_pct"].items():
        print(f"  {pct:6.2f}%  {feature}", file=sys.stderr)
    for warning in result["warnings"]:
        print(f"warning: {warning}", file=sys.stderr)
    print(summary)
    return 0


def cmd_check(args) -> int:
    info = Engine().run("session_info", {"genome": args.genome})["result"]
    ready = not info["required_missing"]
    print(f"Python {info['python']}")
    for pkg, ver in sorted(info["packages"].items()):
        print(f"  ok       {pkg} {ver}")
    for pkg in info["missing"]:
        print(f"  missing  {pkg}")
    print(f"HOMER findMotifsGenome.pl: {'found' if info['homer'] else 'not found (motif enrichment is skipped)'}")
    print(f"Cache: {info['cache_dir']} ({len(info['cached_references'])} gene-model reference(s) built)")
    from .external import available
    for name, cmd in available().items():
        print(f"{ {'snv': 'SimpNetworkViz', 'pathwayagent': 'PathwayAgent'}[name]}: {cmd or 'not found (optional)'}")
    try:
        import anthropic
        print(f"anthropic SDK: {anthropic.__version__}")
    except ImportError:
        print("anthropic SDK: not installed (needed for agent mode; --no-llm works without it)")
    print(f"ant CLI: {'found' if shutil.which('ant') else 'not found'}; "
          f"Anthropic credentials: {anthropic_credentials() or 'none (run `ant auth login` or set ANTHROPIC_API_KEY)'}")
    if args.offline:
        return 0 if ready else 1
    from .http import get_text
    for name, url in [("Europe PMC", "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=CTCF&format=json&pageSize=1"),
                      ("NCBI E-utilities", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?retmode=json"),
                      ("ChIP-Atlas", "https://chip-atlas.org/wabi_endpoint_status")]:
        try:
            get_text(url, timeout=20)
            print(f"{name}: reachable")
        except Exception as err:
            print(f"{name}: unreachable ({err})")
    return 0 if ready else 1


def cmd_tools(_args) -> int:
    for tool in TOOLS.values():
        net = " [network]" if tool.needs_network else ""
        print(f"{tool.name}{net}\n  {tool.description}\n")
    return 0


def cmd_web(args) -> int:
    try:
        from .web import serve
    except ImportError:
        print("The web app needs Flask: pip install 'peakannoagent[web]'", file=sys.stderr)
        return 2
    serve(args.runs_dir, args.host, args.port, model=args.model,
          open_browser=args.open, debug=args.debug)
    return 0


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return {"run": cmd_run, "annotate": cmd_annotate, "check": cmd_check, "tools": cmd_tools, "web": cmd_web}[args.command](args)
    except (FileNotFoundError, ValueError, RuntimeError, json.JSONDecodeError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
