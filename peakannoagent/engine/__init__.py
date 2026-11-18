"""The analysis engine: peak QC, fast annotation, pathway enrichment, pathway networks and motif
enrichment, in Python. `Engine.run(tool, params)` is what the agent's tools call."""
from __future__ import annotations

import platform
import shutil
import time
import traceback
from importlib import metadata
from pathlib import Path
from typing import Optional

from .annotate import annotate_peaks, run_annotation
from . import common
from .common import PALETTE_OPTIONS, SkipTool, cache_dir, resolve_palette, set_palette
from .enrich import METHODS as ENRICH_METHODS, run_enrichment
from .motifs import run_motifs
from .peaks import Peaks, read_peak_file
from .qc import run_peak_qc
from .reference import GENOMES, REFERENCE_OPTIONS, Reference, build_reference, normalize_genome, resolve_reference

__all__ = ["Engine", "Peaks", "Reference", "annotate_peaks", "build_reference", "read_peak_file", "run_annotation",
           "run_enrichment", "run_motifs", "run_peak_qc", "resolve_reference", "REFERENCE_OPTIONS", "ENRICH_METHODS",
           "GENOMES", "PALETTE_OPTIONS", "resolve_palette", "set_palette"]


def _peak_qc(p: dict) -> dict:
    return run_peak_qc(read_peak_file(p["peaks"]), p["out_dir"], p.get("prefix", "run"), p.get("blacklist"))


def _annotate(p: dict) -> dict:
    return run_annotation(read_peak_file(p["peaks"]), p["out_dir"], p.get("prefix", "run"),
                          genome=p.get("genome", "mm10"), reference=p.get("reference"), gtf=p.get("gtf"),
                          promoter_upstream=int(p.get("promoter_upstream", 2000)),
                          promoter_downstream=int(p.get("promoter_downstream", 500)),
                          distal_max_distance=int(p.get("distal_max_distance", 100000)),
                          downstream_bp=int(p.get("downstream_bp", 300)), gene_types=p.get("gene_types", "auto"),
                          cache=p.get("cache_dir"), palette=p.get("palette"))


def _enrich(p: dict) -> dict:
    return run_enrichment(p["gene_sets"], p["out_dir"], p.get("prefix", "run"), genome=p.get("genome", "mm10"),
                          methods=p.get("methods") or ["go_bp"], universe=p.get("universe"), groups=p.get("groups"),
                          padj=float(p.get("padj", 0.05)), top=int(p.get("top", 8)), custom_gmt=p.get("custom_gmt"),
                          cache=p.get("cache_dir"))


def _motifs(p: dict) -> dict:
    return run_motifs(read_peak_file(p["peaks"]), p["out_dir"], p.get("prefix", "run"), genome=p.get("genome", "mm10"),
                      size=int(p.get("size", 200)), max_peaks=int(p.get("max_peaks", 5000)),
                      padj=float(p.get("padj", 0.05)), top=int(p.get("top", 25)))


def session_info(p: Optional[dict] = None) -> dict:
    """What is installed: Python packages the engine needs, HOMER, and the cache location."""
    packages, missing = {}, []
    for name in ("numpy", "matplotlib"):
        try:
            packages[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            missing.append(name)
    cache = cache_dir((p or {}).get("cache_dir"))
    return {"python": platform.python_version(), "genome": normalize_genome((p or {}).get("genome", "mm10")),
            "packages": packages, "required_missing": missing, "missing": missing,
            "homer": bool(shutil.which("findMotifsGenome.pl")), "cache_dir": str(cache),
            "cached_references": sorted(f.name for f in cache.glob("ref*.pkl"))}


TOOLS = {"peak_qc": _peak_qc, "annotate_peaks": _annotate, "enrich_pathways": _enrich, "find_motifs": _motifs,
         "session_info": session_info}


class Engine:
    """Runs one analysis tool and returns {"status": "ok" | "error" | "skipped", "result": ...}."""

    def __init__(self, log_dir: Optional[str] = None):
        self.log_dir = Path(log_dir) if log_dir else None

    def available(self) -> bool:
        return True

    def run(self, tool: str, params: dict) -> dict:
        start = time.time()
        if tool not in TOOLS:
            return {"status": "error", "tool": tool, "error": f"Unknown tool '{tool}'. Valid: {', '.join(TOOLS)}."}
        common.OFFLINE = bool(params.get("offline", False))
        try:
            common.set_palette(params.get("palette"))
            common.set_previews(params.get("previews", False))
            outcome = {"status": "ok", "result": TOOLS[tool](params)}
        except SkipTool as skip:
            outcome = {"status": "skipped", "error": str(skip)}
        except Exception as err:
            outcome = {"status": "error", "error": f"{type(err).__name__}: {err}" if not str(err) else str(err)}
            if self.log_dir:
                self.log_dir.mkdir(parents=True, exist_ok=True)
                (self.log_dir / f"{tool}-{time.time_ns()}.log").write_text(traceback.format_exc())
        outcome.update(tool=tool, elapsed_sec=round(time.time() - start, 1))
        return outcome
