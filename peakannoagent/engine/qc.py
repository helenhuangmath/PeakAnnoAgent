"""Peak-set QC: counts, widths, chromosomes, duplicates, non-standard contigs, blacklist."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

import numpy as np

from . import figures
from .common import write_tsv
from .peaks import Peaks, read_peak_file

STANDARD = re.compile(r"chr([0-9]+|X|Y|M|MT)$")


def blacklist_overlaps(peaks: Peaks, blacklist_path: str) -> np.ndarray:
    """Boolean mask of peaks overlapping a blacklist BED (half-open: touching intervals do not overlap)."""
    bl = read_peak_file(blacklist_path)
    if np.any(bl.start < 0) or np.any(bl.end <= bl.start):
        raise ValueError("Blacklist BED must contain valid zero-based, half-open intervals.")
    out = np.zeros(len(peaks), dtype=bool)
    for chrom in set(bl.chrom):
        rows = np.flatnonzero(peaks.chrom == chrom)
        if not len(rows):
            continue
        sel = bl.chrom == chrom
        order = np.argsort(bl.start[sel])
        s, e = bl.start[sel][order], np.maximum.accumulate(bl.end[sel][order])
        k = np.searchsorted(s, peaks.end[rows], side="left")       # blacklist intervals starting before the peak ends
        hit = k > 0
        hit[hit] = e[k[hit] - 1] > peaks.start[rows][hit]
        out[rows] = hit
    return out


def run_peak_qc(peaks: Peaks, out_dir, prefix: str, blacklist: Optional[str] = None) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    n = len(peaks)
    if n == 0:
        raise ValueError("The peak file contains no peaks.")
    width = peaks.end - peaks.start
    warnings = []
    valid = (peaks.start >= 0) & (peaks.end > peaks.start)
    if (~valid).sum():
        warnings.append(f"{int((~valid).sum())} peaks have invalid BED coordinates (require 0 <= Start < End).")
    keys = list(zip(peaks.chrom, peaks.start, peaks.end))
    seen, first = set(), np.zeros(n, dtype=bool)
    for i, key in enumerate(keys):
        if key not in seen:
            seen.add(key)
            first[i] = True
    dups = int(n - first.sum())
    if dups:
        warnings.append(f"{dups} duplicated peak intervals.")
    standard = np.array([bool(STANDARD.match(c)) for c in peaks.chrom])
    nonstd = int((~standard).sum())
    if nonstd:
        warnings.append(f"{nonstd} peaks on non-standard contigs (random/Un/alt).")
    if not any(c.startswith("chr") for c in peaks.chrom):
        warnings.append("Chromosome names lack the 'chr' prefix; annotation expects chr1, chr2, ...")
    q = [float(v) for v in np.quantile(width, [0, 0.25, 0.5, 0.75, 1])]
    if q[2] > 5000:
        warnings.append("Median peak width > 5 kb: looks like broad domains; avoid summit-centred motif scans.")
    if n < 200:
        warnings.append(f"Only {n} peaks: enrichment and motif statistics will have low power.")
    black = np.zeros(n, dtype=bool)
    if blacklist:
        black = blacklist_overlaps(peaks, blacklist) & valid
        if black.sum():
            warnings.append(f"{int(black.sum())} peaks ({100 * black.sum() / n:.1f}%) overlap the blacklist.")

    counts = {}
    for c in peaks.chrom:
        counts[c] = counts.get(c, 0) + 1
    chrom = sorted(counts.items(), key=lambda kv: -kv[1])
    chrom_path = write_tsv(out_dir / f"{prefix}_qc_chromosomes.tsv", ("Chr", "N", "Pct"),
                           [(c, k, round(100 * k / n, 2)) for c, k in chrom])
    clean = peaks.subset(standard & valid & ~black & first)
    if not len(clean):
        raise ValueError("No valid peaks remain after QC and blacklist filtering.")
    clean_path = clean.write_bed(out_dir / f"{prefix}_qc_clean_peaks.bed")
    figs = figures.qc_figures(width, q[2], [(c, k) for c, k in chrom if STANDARD.match(c)] or chrom, out_dir, prefix)
    return {
        "n_peaks": n, "n_clean_peaks": len(clean),
        "width": dict(zip(("min", "q25", "median", "q75", "max"), q)),
        "n_duplicated": dups, "n_nonstandard_chrom": nonstd, "n_blacklisted": int(black.sum()),
        "top_chromosomes": [{"Chr": c, "N": k, "Pct": round(100 * k / n, 2)} for c, k in chrom[:5]],
        "warnings": warnings, "tables": {"chromosomes": chrom_path}, "clean_peaks": clean_path, "figures": figs,
    }
