"""Fast peak annotation on a cached gene-model reference (see reference.py).

For every peak: one genomic feature by priority (promoter, 5' UTR, 3' UTR, exon, intron,
downstream, distal intergenic - the ChIPseeker convention), the nearest transcription start
site and its gene, and every gene whose promoter window the peak overlaps. `run_annotation`
adds the tables, the gene sets for enrichment and the figure set.
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from . import figures
from .common import set_palette, write_tsv
from .peaks import Peaks, check_intervals
from .reference import Reference, build_reference

FEATURE_SETS = ("Promoter", "5' UTR", "3' UTR", "Exon", "Intron", "Downstream", "Distal Intergenic")
REGIONS = ("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns", "Intergenic.Region")
TSS_BINS = ("0-1kb", "1-3kb", "3-5kb", "5-10kb", "10-100kb", ">100kb")
COLUMNS = ("PeakID", "Chr", "Start", "End", "feature", "feature_overlaps", "nearest_symbol", "nearest_gene_id",
           "nearest_gene_type", "nearest_transcript", "distance_to_tss", "promoter_symbols", "promoter_gene_ids",
           "region_class")


def _bounds(chrom: np.ndarray, n_chrom: int):
    """Start and end row of each chromosome in a table sorted by chromosome."""
    return np.searchsorted(chrom, np.arange(n_chrom)), np.searchsorted(chrom, np.arange(n_chrom), side="right")


def _any_overlap(table: Dict[str, np.ndarray], mask: Optional[np.ndarray], n_chrom: int,
                 q_chrom: np.ndarray, q_start: np.ndarray, q_end: np.ndarray) -> np.ndarray:
    """Does each query interval (closed) overlap any interval of the table?"""
    out = np.zeros(len(q_chrom), dtype=bool)
    chrom, start, end = table["chrom"], table["start"], table["end"]
    if mask is not None:
        chrom, start, end = chrom[mask], start[mask], end[mask]
    if not len(chrom):
        return out
    lo, hi = _bounds(chrom, n_chrom)
    for c in np.unique(q_chrom[q_chrom >= 0]):
        if hi[c] == lo[c]:
            continue
        s, e = start[lo[c]:hi[c]], end[lo[c]:hi[c]]
        reach = np.maximum.accumulate(e)                 # furthest end among intervals starting so far
        rows = np.flatnonzero(q_chrom == c)
        k = np.searchsorted(s, q_end[rows], side="right")
        hit = k > 0
        hit[hit] = reach[k[hit] - 1] >= q_start[rows][hit]
        out[rows] = hit
    return out


def promoter_label(gap: np.ndarray, window_bp: int) -> List[str]:
    if window_bp <= 1000:
        return ["Promoter"] * len(gap)
    kb = np.maximum(1, np.ceil(gap / 1000).astype(int))
    return ["Promoter (<=1kb)" if k == 1 else f"Promoter ({k - 1}-{k}kb)" for k in kb]


def annotate_peaks(peaks: Peaks, reference: Reference, promoter_upstream: int = 2000, promoter_downstream: int = 500,
                   downstream_bp: int = 300, distal_max_distance: int = 100000,
                   gene_types: str = "auto") -> Dict[str, list]:
    """Annotate peaks. Returns the annotation as a dict of equal-length columns (see COLUMNS)."""
    if not isinstance(reference, Reference):
        raise TypeError("`reference` must come from build_reference().")
    check_intervals(peaks)
    ref = reference
    only_pc = ref.gene_types(gene_types) == "protein_coding"
    up, down = int(promoter_upstream), int(promoter_downstream)
    n, n_chrom = len(peaks), len(ref.chroms)
    chrom_index = {c: i for i, c in enumerate(ref.chroms)}
    q_chrom = np.array([chrom_index.get(c, -1) for c in peaks.chrom], dtype=np.int32)
    q_start, q_end = peaks.start + 1, peaks.end.copy()     # one-based, closed
    mid = (q_start + q_end) // 2

    tss = ref.tss
    keep = tss["pc"] if only_pc else np.ones(len(tss["pos"]), dtype=bool)
    t_chrom, t_pos, t_strand, t_gene = tss["chrom"][keep], tss["pos"][keep], tss["strand"][keep], tss["gene"][keep]
    t_row = np.flatnonzero(keep)
    lo, hi = _bounds(t_chrom, n_chrom)

    nearest = np.full(n, -1, dtype=np.int64)               # row of the nearest TSS in the kept table
    gap = np.full(n, -1, dtype=np.int64)                   # distance to the nearest promoter TSS, if any
    prom_genes: List[List[int]] = [[] for _ in range(n)]
    reach = max(up, down)
    for c in np.unique(q_chrom[q_chrom >= 0]):
        if hi[c] == lo[c]:
            continue
        pos, strand, gene = t_pos[lo[c]:hi[c]], t_strand[lo[c]:hi[c]], t_gene[lo[c]:hi[c]]
        rows = np.flatnonzero(q_chrom == c)
        # nearest TSS to the peak centre
        k = np.searchsorted(pos, mid[rows])
        left, right = np.clip(k - 1, 0, len(pos) - 1), np.clip(k, 0, len(pos) - 1)
        pick = np.where(np.abs(mid[rows] - pos[left]) <= np.abs(pos[right] - mid[rows]), left, right)
        nearest[rows] = lo[c] + pick
        # every TSS whose promoter window the peak overlaps
        a = np.searchsorted(pos, q_start[rows] - reach)
        b = np.searchsorted(pos, q_end[rows] + reach, side="right")
        counts = b - a
        if counts.sum():
            who = np.repeat(np.arange(len(rows)), counts)
            idx = np.repeat(a, counts) + (np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts))
            w_start = np.where(strand[idx] < 0, pos[idx] - down, pos[idx] - up)
            w_end = np.where(strand[idx] < 0, pos[idx] + up, pos[idx] + down)
            qs, qe = q_start[rows][who], q_end[rows][who]
            ok = (w_start <= qe) & (w_end >= qs)
            who, idx, qs, qe = who[ok], idx[ok], qs[ok], qe[ok]
            dist = np.maximum(0, np.maximum(qs - pos[idx], pos[idx] - qe))
            for peak_i, g, d in zip(rows[who], gene[idx], dist):
                if gap[peak_i] < 0 or d < gap[peak_i]:
                    gap[peak_i] = d
                if g not in prom_genes[peak_i]:
                    prom_genes[peak_i].append(int(g))

    def hits(table, flag=None):
        mask = None
        if only_pc and flag:
            mask = table["first_pc"]
        elif only_pc:
            mask = table["pc"]
        elif flag:
            mask = table["first"]
        return _any_overlap(table, mask, n_chrom, q_chrom, q_start, q_end)

    tes = ref.tes
    tes_table = {"chrom": tes["chrom"],
                 "start": np.where(tes["strand"] < 0, tes["pos"] - int(downstream_bp), tes["pos"] + 1),
                 "end": np.where(tes["strand"] < 0, tes["pos"] - 1, tes["pos"] + int(downstream_bp)), "pc": tes["pc"]}
    order = np.lexsort((tes_table["start"], tes_table["chrom"]))
    tes_table = {k: v[order] for k, v in tes_table.items()}
    promoter = gap >= 0
    utr5, utr3 = hits(ref.utr5), hits(ref.utr3)
    exon1, exon = hits(ref.exons, "first"), hits(ref.exons)
    intron1, intron = hits(ref.introns, "first"), hits(ref.introns)
    downstream = hits(tes_table)

    down_label = f"Downstream (<={downstream_bp // 1000}kb)" if downstream_bp >= 1000 else f"Downstream (<={downstream_bp}bp)"
    feature = np.full(n, "Distal Intergenic", dtype=object)
    for flag, label in ((downstream, down_label), (intron, "Other Intron"), (intron1, "1st Intron"),
                        (exon, "Other Exon"), (exon1, "1st Exon"), (utr3, "3' UTR"), (utr5, "5' UTR")):
        feature[flag] = label                              # written from lowest to highest priority
    if promoter.any():
        feature[promoter] = promoter_label(gap[promoter], max(up, down))
    member = np.column_stack([promoter, utr5, utr3, exon, intron, downstream])
    member = np.column_stack([member, ~member.any(axis=1)])

    genes = ref.genes
    has = nearest >= 0
    gene_row = np.where(has, t_gene[np.clip(nearest, 0, None)] if len(t_gene) else 0, -1)
    near_pos = np.where(has, t_pos[np.clip(nearest, 0, None)] if len(t_pos) else 0, 0)
    near_strand = np.where(has, t_strand[np.clip(nearest, 0, None)] if len(t_strand) else 1, 1)
    distance = np.where(near_strand < 0, near_pos - mid, mid - near_pos)
    out: Dict[str, list] = {c: [] for c in COLUMNS}
    for i in range(n):
        g = int(gene_row[i])
        is_prom = bool(promoter[i])
        d = int(distance[i]) if has[i] else None
        out["PeakID"].append(peaks.ids[i])
        out["Chr"].append(peaks.chrom[i])
        out["Start"].append(int(peaks.start[i]))
        out["End"].append(int(peaks.end[i]))
        out["feature"].append(feature[i])
        out["feature_overlaps"].append(";".join(s for s, m in zip(FEATURE_SETS, member[i]) if m))
        out["nearest_symbol"].append(genes["symbol"][g] if has[i] else None)
        out["nearest_gene_id"].append(genes["gene_id"][g] if has[i] else None)
        out["nearest_gene_type"].append(genes["gene_type"][g] if has[i] else None)
        out["nearest_transcript"].append(tss["tx_id"][t_row[nearest[i]]] if has[i] else None)
        out["distance_to_tss"].append(d)
        symbols = [genes["symbol"][k] for k in prom_genes[i] if genes["symbol"][k]]
        out["promoter_symbols"].append(";".join(dict.fromkeys(symbols)) or None)
        out["promoter_gene_ids"].append(";".join(genes["gene_id"][k] for k in prom_genes[i]) or None)
        out["region_class"].append("promoter" if is_prom else
                                   ("distal" if d is not None and abs(d) <= distal_max_distance else "distal_far"))
    order = sorted(range(n), key=lambda i: (out["Chr"][i], out["Start"][i]))
    return {c: [v[i] for i in order] for c, v in out.items()}


def _feature_rank(label: str):
    for rank, prefix in enumerate(("Promoter", "5' UTR", "3' UTR", "1st Exon", "Other Exon", "1st Intron",
                                   "Other Intron", "Downstream")):
        if label.startswith(prefix):
            return (rank, 0 if "<=" in label else 1, label)       # "Promoter (<=1kb)" before "(1-2kb)"
    return (9, 0, label)


def feature_distribution(anno: Dict[str, list]) -> List[dict]:
    counts: Dict[str, int] = {}
    for f in anno["feature"]:
        counts[f] = counts.get(f, 0) + 1
    total = max(1, len(anno["feature"]))
    return [{"Annotation": f, "Freq": counts[f], "Pct": round(100 * counts[f] / total, 2)}
            for f in sorted(counts, key=_feature_rank)]


def region_distribution(anno: Dict[str, list]) -> List[dict]:
    """The categories of ChIPpeakAnno::assignChromosomeRegion, one region per peak."""
    def region(f: str) -> str:
        if f.startswith("Promoter"):
            return "Promoters"
        if f.startswith("Downstream"):
            return "immediateDownstream"
        return {"5' UTR": "fiveUTRs", "3' UTR": "threeUTRs", "1st Exon": "Exons", "Other Exon": "Exons",
                "1st Intron": "Introns", "Other Intron": "Introns"}.get(f, "Intergenic.Region")
    regions = [region(f) for f in anno["feature"]]
    total = max(1, len(regions))
    return [{"Region": r, "Pct": round(100 * regions.count(r) / total, 2)} for r in REGIONS]


def tss_distance_bins(anno: Dict[str, list]) -> List[dict]:
    d = np.array([v for v in anno["distance_to_tss"] if v is not None], dtype=float)
    edges = [1000, 3000, 5000, 10000, 100000]
    which = np.searchsorted(edges, np.abs(d), side="left")
    rows = []
    for side, sel in (("upstream", d < 0), ("downstream", d >= 0)):
        for k, label in enumerate(TSS_BINS):
            count = int(np.sum(sel & (which == k)))
            rows.append({"Bin": label, "Side": side, "N": count, "Pct": round(100 * count / max(1, len(d)), 2)})
    return rows


def tss_profile(peaks: Peaks, reference: Reference, gene_types: str = "auto", radius: int = 5000,
                nbins: int = 50) -> List[dict]:
    """Peak centres per bin around all TSSs (what ChIPpeakAnno::binOverFeature draws)."""
    tss = reference.tss
    keep = tss["pc"] if reference.gene_types(gene_types) == "protein_coding" else np.ones(len(tss["pos"]), dtype=bool)
    key = np.unique(np.column_stack([tss["chrom"][keep], tss["pos"][keep], tss["strand"][keep]]), axis=0)
    chrom_index = {c: i for i, c in enumerate(reference.chroms)}
    q_chrom = np.array([chrom_index.get(c, -1) for c in peaks.chrom])
    mid = (peaks.start + 1 + peaks.end) // 2
    edges = np.linspace(-radius, radius, nbins + 1)
    counts = np.zeros(nbins, dtype=np.int64)
    for c in np.unique(q_chrom[q_chrom >= 0]):
        sel = key[:, 0] == c
        pos, strand = key[sel, 1], key[sel, 2]
        if not len(pos):
            continue
        m = np.sort(mid[q_chrom == c])
        a, b = np.searchsorted(m, pos - radius), np.searchsorted(m, pos + radius, side="right")
        n = b - a
        if not n.sum():
            continue
        idx = np.repeat(a, n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
        rel = (m[idx] - np.repeat(pos, n)) * np.repeat(strand, n)
        counts += np.histogram(rel, bins=edges)[0]
    return [{"distance": float((edges[i] + edges[i + 1]) / 2), "peaks": int(counts[i])} for i in range(nbins)]


def run_annotation(peaks: Peaks, out_dir, prefix: str, genome: str = "mm10", reference=None, gtf: Optional[str] = None,
                   promoter_upstream: int = 2000, promoter_downstream: int = 500, distal_max_distance: int = 100000,
                   downstream_bp: int = 300, gene_types: str = "auto", draw_figures: bool = True,
                   cache: Optional[str] = None, palette: Optional[str] = None) -> dict:
    """Annotate peaks and write the tables, the gene sets and background for enrichment, and the
    figures. Returns a JSON-serialisable summary with the paths of everything written."""
    started = time.perf_counter()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref = reference if isinstance(reference, Reference) else build_reference(genome, reference, gtf, cache)
    loaded = time.perf_counter()
    anno = annotate_peaks(peaks, ref, promoter_upstream, promoter_downstream, downstream_bp, distal_max_distance,
                          gene_types)
    annotated = time.perf_counter()
    n = len(anno["PeakID"])

    def path(name: str) -> str:
        return str(out_dir / f"{prefix}_{name}.tsv")

    feature_tab, region_tab, bins = feature_distribution(anno), region_distribution(anno), tss_distance_bins(anno)
    profile = tss_profile(peaks, ref, gene_types)
    write_tsv(path("peak_genes"), COLUMNS, zip(*(anno[c] for c in COLUMNS)))
    for name, rows in (("feature_distribution", feature_tab), ("region_distribution", region_tab),
                       ("tss_distance_bins", bins)):
        write_tsv(path(name), list(rows[0]), [list(r.values()) for r in rows])

    kept = ref.gene_types(gene_types)
    allowed = set(ref.symbols(gene_types))
    promoter_genes = list(dict.fromkeys(s for cell in anno["promoter_symbols"] if cell
                                        for s in cell.split(";") if s in allowed))
    distal_genes = list(dict.fromkeys(s for s, rc in zip(anno["nearest_symbol"], anno["region_class"])
                                      if rc == "distal" and s in allowed))
    all_genes = list(dict.fromkeys(promoter_genes + distal_genes))
    sets = [("promoter", g) for g in promoter_genes] + [("distal", g) for g in distal_genes] + [("all", g) for g in all_genes]
    write_tsv(path("peak_gene_sets"), ("Cluster", "SYMBOL"), sets)
    write_tsv(path("gene_universe"), ("SYMBOL",), [(s,) for s in sorted(allowed)])

    warnings: List[str] = []
    classes = anno["region_class"]
    off = sum(d is None for d in anno["distance_to_tss"])
    if off:
        warnings.append(f"{off} peaks lie on sequences absent from the reference and have no nearest gene.")
    if classes.count("promoter") > 0.6 * n:
        warnings.append("Most peaks are promoter-proximal: typical of promoter-bound factors or H3K4me3.")
    if classes.count("distal_far") > 0.3 * n:
        warnings.append(f"Over 30% of peaks are > {distal_max_distance // 1000} kb from any TSS; "
                        f"nearest-gene links are weak.")
    figs: Dict[str, str] = {}
    if draw_figures:
        if palette:
            set_palette(palette)
        figs, failed = figures.annotation_figures(anno, ref, feature_tab, region_tab, bins, profile, out_dir, prefix)
        warnings.extend(failed)
    done = time.perf_counter()
    dist = [abs(d) for d in anno["distance_to_tss"] if d is not None]
    meta = ref.meta
    return {
        "reference": {"genome": meta["genome"], "source": meta["source"], "genes": meta["n_genes"],
                      "transcripts": meta["n_transcripts"], "gene_types": kept, "genes_used": len(allowed),
                      "cached": bool(meta.get("from_cache"))},
        "n_peaks": n,
        "region_class_counts": {k: classes.count(k) for k in ("promoter", "distal", "distal_far") if classes.count(k)},
        "feature_pct": {r["Annotation"]: r["Pct"] for r in feature_tab},
        "region_distribution_pct": {r["Region"]: r["Pct"] for r in region_tab},
        "median_abs_distance_to_tss": float(np.median(dist)) if dist else None,
        "gene_set_sizes": {"promoter": len(promoter_genes), "distal": len(distal_genes), "all": len(all_genes)},
        "example_promoter_genes": promoter_genes[:40],
        "promoter_window": {"upstream": int(promoter_upstream), "downstream": int(promoter_downstream)},
        "timing_sec": {"reference": round(loaded - started, 2), "annotation": round(annotated - loaded, 2),
                       "tables_and_figures": round(done - annotated, 2)},
        "warnings": warnings,
        "tables": {"peak_genes": path("peak_genes"), "feature_distribution": path("feature_distribution"),
                   "region_distribution": path("region_distribution"), "tss_distance_bins": path("tss_distance_bins")},
        "gene_sets": path("peak_gene_sets"),
        "universe": path("gene_universe"),
        "figures": figs,
        "versions": {"annotation": meta["label"]},
    }
