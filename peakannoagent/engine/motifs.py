"""Known-motif enrichment with HOMER (findMotifsGenome.pl), if it is installed."""
from __future__ import annotations

import math
import re
import shutil
import subprocess
from pathlib import Path
from typing import List

from . import figures
from .common import SkipTool, read_tsv, write_tsv
from .peaks import Peaks

HOMER_INSTALL = "install HOMER and its genome package (http://homer.ucsd.edu/homer/) so that findMotifsGenome.pl is on PATH"


def motif_tf_name(name: str) -> str:
    """"BATF(bZIP)/Th17-BATF-ChIP-Seq/Homer" -> "BATF"."""
    return re.split(r"[(/]", name, maxsplit=1)[0].strip()


def read_homer_known_results(path) -> List[dict]:
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]

    def pct(text: str) -> float:
        return float(text.replace("%", ""))
    p_col, q_col = header.index("P-value"), header.index("q-value (Benjamini)")
    out = []
    for r in rows:
        family = re.sub(r"^[^(]*\(([^)]*)\).*$", r"\1", r[0]) if "(" in r[0] else ""
        target, background = pct(r[6]), pct(r[8])
        out.append({"motif": r[0], "tf": motif_tf_name(r[0]), "family": family, "pct_target": target,
                    "pct_background": background,
                    "log2_enrichment": round(math.log2((target + 0.1) / (background + 0.1)), 3),
                    "pvalue": float(r[p_col]), "padj": float(r[q_col])})
    return out


def run_motifs(peaks: Peaks, out_dir, prefix: str, genome: str = "mm10", size: int = 200, max_peaks: int = 5000,
               padj: float = 0.05, top: int = 25, threads: int = 4) -> dict:
    exe = shutil.which("findMotifsGenome.pl")
    if not exe:
        raise SkipTool(f"HOMER is not installed; motif enrichment was skipped ({HOMER_INSTALL}).")
    out_dir = Path(out_dir)
    homer_dir = out_dir / f"{prefix}_homer_motifs"
    homer_dir.mkdir(parents=True, exist_ok=True)
    n = min(len(peaks), int(max_peaks))
    bed = homer_dir / "input_peaks.bed"
    with open(bed, "w") as handle:
        for i in range(n):
            handle.write(f"{peaks.chrom[i]}\t{peaks.start[i]}\t{peaks.end[i]}\t{peaks.ids[i]}\t.\t+\n")
    proc = subprocess.run([exe, str(bed), genome, str(homer_dir), "-size", str(size), "-p", str(threads), "-nomotif"],
                          capture_output=True, text=True)
    known = homer_dir / "knownResults.txt"
    if proc.returncode != 0 or not known.is_file():
        tail = "\n".join((proc.stderr or proc.stdout).strip().splitlines()[-15:])
        raise RuntimeError(f"HOMER findMotifsGenome.pl failed:\n{tail}")
    table = read_homer_known_results(known)
    cols = ("motif", "tf", "family", "pct_target", "pct_background", "log2_enrichment", "pvalue", "padj")
    table_path = write_tsv(out_dir / f"{prefix}_motifs_homer.tsv", cols, [[r[c] for c in cols] for r in table])
    sig = [r for r in table if r["padj"] <= padj and r["log2_enrichment"] > 0]
    figs, warnings, nominal = {}, [], None
    if sig:
        figs["motifs"] = figures.motif_figure(sig[:min(top, 20)], "HOMER", out_dir, f"{prefix}_motifs_homer")
    else:
        warnings.append("No motif passed the adjusted p-value cutoff; top_nominal lists the best motifs by raw "
                        "p-value (not significant).")
        nominal = [{k: r[k] for k in ("tf", "log2_enrichment", "pvalue", "padj")}
                   for r in sorted((r for r in table if r["log2_enrichment"] > 0), key=lambda r: r["pvalue"])[:10]]
    return {"backend": "HOMER", "database": "HOMER known motifs", "n_sequences": n, "window_bp": size,
            "n_significant": len(sig),
            "top_motifs": [{k: r[k] for k in ("tf", "family", "pct_target", "pct_background", "log2_enrichment", "padj")}
                           for r in sig[:top]],
            "top_nominal": nominal, "warnings": warnings, "tables": {"motifs": table_path}, "figures": figs,
            "versions": {}}
