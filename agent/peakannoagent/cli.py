"""Command line: `peakanno-agent run|check|tools`."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys

from .config import DEFAULT_MODEL, SUPPORTED_GENOMES, RunConfig, build_biology, load_biology_file
from .rbridge import RBridge
from .tools import TOOLS


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="peakanno-agent", description=__doc__)
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
    run.add_argument("--offline", action="store_true", help="No ChIP-Atlas / literature / GEO calls")
    run.add_argument("--no-llm", action="store_true", help="Fixed analysis plan, no Claude interpretation")
    run.add_argument("--model", default=DEFAULT_MODEL)
    run.add_argument("--effort", default="high", choices=["low", "medium", "high", "xhigh", "max"])
    run.add_argument("--max-turns", type=int, default=40)
    run.add_argument("--rscript", default=os.environ.get("PEAKANNOAGENT_RSCRIPT", "Rscript"),
                     help="Rscript with the Bioconductor packages installed")
    run.add_argument("--chipatlas-wait", type=float, default=25.0,
                     help="Minutes the deterministic plan waits for ChIP-Atlas (default 25)")

    chk = sub.add_parser("check", help="Check R packages, HOMER, API credentials and network access.")
    chk.add_argument("--genome", default="mm10")
    chk.add_argument("--rscript", default=os.environ.get("PEAKANNOAGENT_RSCRIPT", "Rscript"))

    sub.add_parser("tools", help="List the agent's tools.")
    return p


def cmd_run(args) -> int:
    file_fields = load_biology_file(args.biology_file) if args.biology_file else {}
    biology = build_biology(file_fields, description=args.biology, assay=args.assay, target=args.target,
                            cell_type=args.cell_type, condition=args.condition,
                            expected_genes=args.expected_genes, expected_tfs=args.expected_tfs,
                            custom_gmt=args.gmt, rnaseq=args.rnaseq)
    config = RunConfig(peaks=args.peaks, genome=args.genome, out_dir=args.out, prefix=args.prefix,
                       biology=biology, blacklist=args.blacklist, offline=args.offline, use_llm=not args.no_llm,
                       model=args.model, effort=args.effort, max_turns=args.max_turns, rscript=args.rscript,
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
    return 0


def cmd_check(args) -> int:
    bridge = RBridge(args.rscript)
    print(f"Rscript: {args.rscript} ({'found' if bridge.available() else 'NOT FOUND'})")
    if bridge.available():
        info = bridge.run("session_info", {"genome": args.genome})
        if info.get("status") == "ok":
            r = info["result"]
            print(f"  {r['r_version']}")
            for pkg, ver in sorted(r["packages"].items()):
                print(f"  ok       {pkg} {ver}")
            for pkg in r["missing"]:
                print(f"  missing  {pkg}")
            print(f"  HOMER findMotifsGenome.pl: {'yes' if r['homer'] else 'no (monaLisa will be used)'}")
        else:
            print(f"  R error: {info.get('error')}")
    try:
        import anthropic
        print(f"anthropic SDK: {anthropic.__version__}")
    except ImportError:
        print("anthropic SDK: not installed (needed for agent mode; --no-llm works without it)")
    print(f"ant CLI: {'found' if shutil.which('ant') else 'not found'}; "
          f"ANTHROPIC_API_KEY: {'set' if os.environ.get('ANTHROPIC_API_KEY') else 'not set'}")
    from .http import get_text
    for name, url in [("Europe PMC", "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=CTCF&format=json&pageSize=1"),
                      ("NCBI E-utilities", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi?retmode=json"),
                      ("ChIP-Atlas", "https://chip-atlas.org/wabi_endpoint_status")]:
        try:
            get_text(url, timeout=20)
            print(f"{name}: reachable")
        except Exception as err:
            print(f"{name}: unreachable ({err})")
    return 0


def cmd_tools(_args) -> int:
    for tool in TOOLS.values():
        net = " [network]" if tool.needs_network else ""
        print(f"{tool.name}{net}\n  {tool.description}\n")
    return 0


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return {"run": cmd_run, "check": cmd_check, "tools": cmd_tools}[args.command](args)
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
