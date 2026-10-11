"""ChIP-Atlas: overlap of the user's peaks with every public ChIP-seq/ATAC-seq
experiment (Enrichment Analysis), and per-TF target gene lists."""

from __future__ import annotations

import csv
import io
import math
import re
from collections import defaultdict
from typing import Dict, Iterable, List

from .http import get_text, post_form

SUBMIT_URL = "https://chip-atlas.org/wabi_chipatlas"           # the site's own proxy to DDBJ WABI
STATUS_URL = "https://chip-atlas.org/wabi_chipatlas?id={id}"
RESULT_URL = "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/{id}?info=result&format={fmt}"
TARGET_URL = "https://chip-atlas.dbcls.jp/data/{genome}/target/{tf}.{distance}.tsv"
CHIPATLAS_GENOMES = ("hg38", "hg19", "mm10", "mm9")
ANTIGEN_CLASSES = ("TFs and others", "Histone", "ATAC-Seq", "DNase-seq", "RNA polymerase", "Input control",
                   "Bisulfite-Seq")

RESULT_COLUMNS = ["experiment", "antigen_class", "antigen", "cell_class", "cell", "n_peaks",
                  "overlaps_user", "overlaps_control", "log10_p", "log10_q", "fold_enrichment"]


def bed_payload(bed_path: str, max_peaks: int = 50000) -> str:
    """Chr/start/end lines with only the characters the service accepts (as the web form does)."""
    lines = []
    with open(bed_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[1].isdigit():
                continue
            lines.append("\t".join(fields[:3]))
            if len(lines) >= max_peaks:
                break
    return re.sub(r"[^a-zA-Z0-9\t_\n]", "_", "\n".join(lines))


def submit_enrichment(bed_path: str, genome: str, antigen_class: str = "TFs and others",
                      cell_class: str = "All cell types", threshold: int = 50,
                      distance: int = 5000, title: str = "PeakAnnoAgent") -> str:
    if genome not in CHIPATLAS_GENOMES:
        raise ValueError(f"ChIP-Atlas supports {', '.join(CHIPATLAS_GENOMES)}; got {genome}.")
    if antigen_class not in ANTIGEN_CLASSES:
        raise ValueError(f"antigen_class must be one of {ANTIGEN_CLASSES}.")
    fields = {
        "address": "", "format": "text", "result": "www", "genome": genome,
        "antigenClass": antigen_class, "cellClass": cell_class, "threshold": str(threshold),
        "typeA": "bed", "bedAFile": bed_payload(bed_path), "typeB": "rnd", "bedBFile": "empty",
        "permTime": "1", "title": re.sub(r"[^A-Za-z0-9_. -]", "_", title)[:60],
        "descriptionA": "user_peaks", "descriptionB": "random_permutation",
        "distanceUp": str(distance), "distanceDown": str(distance),
        "sbatchOptions": "-p epyc -t 180",
    }
    text = post_form(SUBMIT_URL, fields)
    match = re.search(r"wabi_chipatlas_[\w-]+", text)
    if not match:
        raise RuntimeError(f"ChIP-Atlas did not return a request id: {text[:200]}")
    return match.group(0)


def job_status(request_id: str) -> str:
    """'running', 'finished', 'unavailable', or the raw text the service returned."""
    return get_text(STATUS_URL.format(id=request_id), timeout=30).strip()


def _num(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def parse_enrichment_tsv(text: str) -> List[dict]:
    """Parse the result TSV (header optional) into dicts keyed by RESULT_COLUMNS."""
    rows = []
    for fields in csv.reader(io.StringIO(text), delimiter="\t"):
        if len(fields) < len(RESULT_COLUMNS) or not re.match(r"^[SDE]RX\d+", fields[0]):
            continue
        row = dict(zip(RESULT_COLUMNS, fields[: len(RESULT_COLUMNS)]))
        for key in ("log10_p", "log10_q", "fold_enrichment"):
            row[key] = _num(row[key])
        row["n_peaks"] = int(_num(row["n_peaks"])) if not math.isnan(_num(row["n_peaks"])) else None
        rows.append(row)
    return rows


def fetch_enrichment(request_id: str) -> List[dict]:
    text = get_text(RESULT_URL.format(id=request_id, fmt="tsv"), timeout=180)
    if "NOT FOUND" in text[:300]:
        raise RuntimeError("Result not ready.")
    return parse_enrichment_tsv(text)


def summarize_by_antigen(rows: Iterable[dict], q_cutoff: float = -2.0, top: int = 25) -> List[dict]:
    """Collapse experiments to one line per antigen (TF/mark): best q, # significant experiments, cells."""
    groups: Dict[str, list] = defaultdict(list)
    for r in rows:
        groups[r["antigen"]].append(r)
    summary = []
    for antigen, exps in groups.items():
        sig = [e for e in exps if e["log10_q"] <= q_cutoff and e["fold_enrichment"] > 1]
        if not sig:
            continue
        best = min(sig, key=lambda e: (e["log10_q"], -e["fold_enrichment"]))
        cells = sorted({e["cell"] for e in sig})
        summary.append({
            "antigen": antigen,
            "best_log10_q": round(best["log10_q"], 2),
            "best_fold_enrichment": round(best["fold_enrichment"], 2),
            "best_experiment": best["experiment"],
            "best_cell": best["cell"],
            "n_significant_experiments": len(sig),
            "n_experiments": len(exps),
            "cells": cells[:6],
        })
    summary.sort(key=lambda s: (s["best_log10_q"], -s["n_significant_experiments"]))
    return summary[:top]


def target_genes(tf: str, genome: str, distance_kb: int = 5) -> Dict[str, float]:
    """Genes with a ChIP-Atlas peak for `tf` within +/- distance_kb of the TSS -> average MACS2 score."""
    if distance_kb not in (1, 5, 10):
        raise ValueError("distance_kb must be 1, 5 or 10.")
    text = get_text(TARGET_URL.format(genome=genome, tf=tf, distance=distance_kb), timeout=120)
    if text.lstrip().startswith("<"):
        raise RuntimeError(f"ChIP-Atlas has no target-gene table for '{tf}' on {genome}.")
    out = {}
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    next(reader, None)
    for fields in reader:
        if len(fields) >= 2:
            score = _num(fields[1])
            if score > 0:
                out[fields[0]] = score
    return out
