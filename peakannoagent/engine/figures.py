"""Figures of the analysis engine, drawn with matplotlib. Every figure is written as a vector
PDF with editable text (see common.save_figure); colours come from the chosen palette."""
from __future__ import annotations

import math
import re
import textwrap
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .common import palette, save_figure


FEATURE_ORDER = ("5' UTR", "3' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "Downstream",
                 "Distal Intergenic")


def feature_colors(features: Sequence[str]) -> Dict[str, str]:
    """One colour per feature type. With a categorical palette each feature keeps its colour from
    run to run and the promoter bins are shades of the first colour; a continuous palette is
    sampled in the order given."""
    pal = palette()
    if not pal.qualitative:
        return dict(zip(features, pal.colors(len(features))))
    fixed = pal.colors(len(FEATURE_ORDER) + 1)
    promoters = [f for f in features if f.startswith("Promoter")]
    out = {}
    for f in features:
        if f in promoters:
            i = promoters.index(f)
            if len(promoters) == 1:
                out[f] = fixed[0]
            elif pal.promoters:
                out[f] = pal.promoters[min(i, len(pal.promoters) - 1)]
            else:
                out[f] = _lighten(fixed[0], 0.55 * i / (len(promoters) - 1))
        else:
            key = "Downstream" if f.startswith("Downstream") else f
            if key in pal.features:
                out[f] = pal.features[key]
            else:
                out[f] = fixed[1 + FEATURE_ORDER.index(key)] if key in FEATURE_ORDER else "#9E9E9E"
    return out


def _ink(color: str, limit: float = 0.45) -> str:
    """A colour dark enough for text and thin lines on white: pale colours are darkened."""
    from matplotlib.colors import to_hex, to_rgb
    c = np.array(to_rgb(color))
    lum = float(c @ np.array([0.2126, 0.7152, 0.0722]))
    return to_hex(c * limit / lum) if lum > limit else to_hex(c)


def natural_chrom_order(chroms: Sequence[str]) -> List[str]:
    def key(c: str):
        core = c[3:] if c.startswith("chr") else c
        return (0, int(core), "") if core.isdigit() else (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(core, 3), core)
    return sorted(set(chroms), key=key)


# ---- annotation ------------------------------------------------------------------------------

def _pie(fig, tab, title):
    cols = feature_colors([r["Annotation"] for r in tab])
    ax = fig.add_axes([0.02, 0.04, 0.5, 0.82])
    ax.pie([r["Freq"] for r in tab], colors=[cols[r["Annotation"]] for r in tab], startangle=90, counterclock=False,
           wedgeprops={"edgecolor": "white", "linewidth": 0.8})
    ax.set_title(title)
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=cols[r["Annotation"]], edgecolor="black", linewidth=0.5,
                     label=f"{r['Annotation']} ({r['Pct']:.2f}%)") for r in tab]
    fig.legend(handles=handles, title="Feature", loc="center left", bbox_to_anchor=(0.55, 0.5), alignment="left")


def _feature_bar(fig, tab):
    cols = feature_colors([r["Annotation"] for r in tab])
    ax = fig.add_axes([0.05, 0.3, 0.6, 0.48])
    left = 0.0
    for r in tab:
        ax.barh(0, r["Pct"], left=left, color=cols[r["Annotation"]], height=0.8, label=r["Annotation"])
        left += r["Pct"]
    ax.set(xlim=(0, 100), ylim=(-0.5, 0.5), xlabel="Percentage (%)", title="Feature distribution", yticks=[])
    ax.set_xticks(range(0, 101, 25))
    ax.spines["left"].set_visible(False)
    fig.legend(title="Feature", loc="center left", bbox_to_anchor=(0.68, 0.5), alignment="left")


def _tss_bins(fig, bins):
    up = [r["Pct"] for r in bins if r["Side"] == "upstream"]
    down = [r["Pct"] for r in bins if r["Side"] == "downstream"]
    lim = max(10, math.ceil(max(sum(up), sum(down)) / 10) * 10)
    ax = fig.add_axes([0.05, 0.3, 0.68, 0.48])
    for side, sign in ((down, 1), (up, -1)):
        left = 0.0
        for pct, col, label in zip(side, palette().ordered(len(up)), [r["Bin"] for r in bins[:len(up)]]):
            ax.barh(0, sign * pct, left=sign * left, color=col, height=0.8, label=label if sign == 1 else None)
            left += pct
    ax.axvline(0, color="black", linewidth=1.2)
    ticks = list(range(-lim, lim + 1, 10))
    ax.set(xlim=(-lim, lim), ylim=(-0.5, 0.5), yticks=[], xlabel="Binding sites (%) (5' -> 3')",
           title="Distribution of peaks relative to TSS")
    ax.set_xticks(ticks)
    ax.set_xticklabels(["TSS" if t == 0 else str(abs(t)) for t in ticks])
    ax.spines["left"].set_visible(False)
    fig.legend(title="Distance", loc="center left", bbox_to_anchor=(0.76, 0.5), alignment="left")


def _upset(fig, anno, sets, max_sets=12):
    combos: Dict[str, int] = {}
    for cell in anno["feature_overlaps"]:
        combos[cell] = combos.get(cell, 0) + 1
    top = sorted(combos.items(), key=lambda kv: -kv[1])[:max_sets]
    totals = [sum(s in cell.split(";") for cell in anno["feature_overlaps"]) for s in sets]
    x = np.arange(len(top))
    bars = fig.add_axes([0.34, 0.45, 0.63, 0.45])
    bars.bar(x, [c for _, c in top], color=palette().bar, width=0.7)
    for i, (_, c) in enumerate(top):
        bars.text(i, c, f"{c:,}", ha="center", va="bottom", fontsize=13)
    bars.set(xlim=(-0.6, len(top) - 0.4), xticks=[], ylabel="Peaks", title="Peaks by overlapping features")
    dots = fig.add_axes([0.34, 0.06, 0.63, 0.36])
    y = np.arange(len(sets))[::-1]
    dots.scatter(np.repeat(x, len(sets)), np.tile(y, len(top)), s=90, color="#D9D9D9")
    for i, (cell, _) in enumerate(top):
        yy = [y[sets.index(s)] for s in cell.split(";") if s in sets]
        dots.plot([i, i], [min(yy), max(yy)], color="black", linewidth=2)
        dots.scatter([i] * len(yy), yy, s=90, color="black", zorder=3)
    dots.set(xlim=(-0.6, len(top) - 0.4), ylim=(-0.6, len(sets) - 0.4), xticks=[])
    dots.set_yticks(y)
    dots.set_yticklabels([f"{s} ({t:,})" for s, t in zip(sets, totals)], fontsize=14)
    dots.tick_params(left=False)
    for side in ("left", "bottom"):
        dots.spines[side].set_visible(False)


def _coverage(fig, anno, reference):
    chroms = natural_chrom_order(anno["Chr"])
    standard = [c for c in chroms if re.fullmatch(r"chr([0-9]+|[XYM])", c)]
    chroms = standard or chroms
    ends: Dict[str, int] = {}
    for c, e in zip(anno["Chr"], anno["End"]):
        ends[c] = max(ends.get(c, 0), e)
    sizes = [max(ends.get(c, 0), reference.chrom_ends.get(c, 0)) / 1e6 for c in chroms]
    row = {c: len(chroms) - 1 - i for i, c in enumerate(chroms)}
    ax = fig.add_axes([0.1, 0.1, 0.86, 0.8])
    ax.hlines([row[c] for c in chroms], 0, sizes, color="#BFBFBF", linewidth=1)
    xs = [(s + e) / 2e6 for c, s, e in zip(anno["Chr"], anno["Start"], anno["End"]) if c in row]
    ys = [row[c] for c in anno["Chr"] if c in row]
    ax.vlines(xs, np.array(ys) - 0.32, np.array(ys) + 0.32, color=palette().primary, linewidth=0.6,
              alpha=0.35 if len(xs) > 5000 else 0.8)
    ax.set(xlim=(0, max(sizes) * 1.02), ylim=(-0.7, len(chroms) - 0.3), xlabel="Chromosome position (Mb)",
           title="Peaks over chromosomes")
    ax.set_yticks([row[c] for c in chroms])
    ax.set_yticklabels(chroms, fontsize=min(14, 280 / max(1, len(chroms))))
    ax.tick_params(left=False)
    ax.spines["left"].set_visible(False)


def annotation_figures(anno, reference, feature_tab, region_tab, bins, profile, out_dir, prefix) -> Tuple[dict, list]:
    from .annotate import FEATURE_SETS
    title = f"{prefix} (n = {len(anno['PeakID']):,} peaks)"

    def region(fig):
        ax = fig.add_axes([0.36, 0.17, 0.6, 0.71])
        rows = region_tab[::-1]
        ax.barh([r["Region"] for r in rows], [r["Pct"] for r in rows], color=palette().bar)
        ax.set(xlabel="% of peaks", title="Genomic region distribution")

    def tss_prof(fig):
        ax = fig.add_axes([0.14, 0.16, 0.82, 0.72])
        ax.plot([r["distance"] for r in profile], [r["peaks"] for r in profile], color=palette().primary, linewidth=2)
        ax.axvline(0, color=palette().accent, linestyle="--", linewidth=1)
        ax.set(xlabel="Distance from TSS (bp)", ylabel="Peak count", title="Peak distribution around TSS")

    def tss_hist(fig):
        d = np.array([v for v in anno["distance_to_tss"] if v is not None], dtype=float)
        ax = fig.add_axes([0.14, 0.16, 0.82, 0.72])
        ax.hist(np.clip(d, -100000, 100000) / 1000, bins=80, color=palette().bar, edgecolor="white", linewidth=0.4)
        ax.axvline(0, color=palette().accent, linestyle="--", linewidth=1)
        ax.set(xlabel="kb (clipped to +/-100 kb)", ylabel="Frequency", title="Distance to nearest TSS")

    plan = [("feature_pie", lambda f: _pie(f, feature_tab, title), 8.5, 5),
            ("feature_bar", lambda f: _feature_bar(f, feature_tab), 11, 3.8),
            ("tss_distance_bins", lambda f: _tss_bins(f, bins), 11, 3.8),
            ("feature_upset", lambda f: _upset(f, anno, list(FEATURE_SETS)), 10, 7),
            ("region_distribution", region, 9, 5.5), ("tss_profile", tss_prof, 8, 5.5),
            ("tss_distance", tss_hist, 8, 5.5),
            ("chromosome_coverage", lambda f: _coverage(f, anno, reference), 10, 7)]
    figs, failed = {}, []
    for name, draw, w, h in plan:
        try:
            figs[name] = save_figure(draw, out_dir, f"{prefix}_{name}", w, h)
        except Exception as err:                       # a failed figure must not lose the annotation
            failed.append(f"Figure {name} failed: {err}")
    return figs, failed


# ---- QC, motifs, enrichment ------------------------------------------------------------------

def qc_figures(widths, median_width, chrom_counts, out_dir, prefix) -> dict:
    def width(fig):
        ax = fig.add_axes([0.14, 0.16, 0.82, 0.72])
        ax.hist(np.log10(np.maximum(widths, 1)), bins=40, color=palette().bar, edgecolor="white", linewidth=0.4)
        ax.axvline(math.log10(max(median_width, 1)), color=palette().accent, linestyle="--", linewidth=1)
        ax.set(xlabel="log10(peak width, bp)", ylabel="Frequency", title="Peak width distribution")

    def chrom(fig):
        order = natural_chrom_order([c for c, _ in chrom_counts])
        counts = dict(chrom_counts)
        ax = fig.add_axes([0.1, 0.16, 0.87, 0.72])
        ax.bar([c[3:] if c.startswith("chr") else c for c in order], [counts[c] for c in order], color=palette().bar)
        ax.set(ylabel="Peaks", title="Peaks per chromosome")
        ax.tick_params(axis="x", labelrotation=90)

    return {"width": save_figure(width, out_dir, f"{prefix}_qc_peak_width", 8, 5.5),
            "chromosomes": save_figure(chrom, out_dir, f"{prefix}_qc_chromosomes", 10, 5.5)}


def motif_figure(rows, backend, out_dir, name) -> str:
    def draw(fig):
        show = rows[::-1]
        ax = fig.add_axes([0.3, 0.14, 0.66, 0.76])
        ax.barh([r["tf"] for r in show], [min(-math.log10(max(r["padj"], 1e-300)), 300) for r in show],
                color=palette().bar)
        ax.set(xlabel="-log10 adjusted p", title=f"Known motif enrichment ({backend})")
    return save_figure(draw, out_dir, name, 9, max(5.5, len(rows) * 0.45 + 2.5))


def pretty_term(term: str) -> str:
    """MSigDB-style names (HALLMARK_TNFA_SIGNALING_VIA_NFKB) -> "TNFA signaling via NFKB"."""
    if not re.fullmatch(r"[A-Z0-9]+(_[A-Z0-9]+)+", term):
        return term
    words = re.sub(r"^(HALLMARK|KEGG|KEGG_MEDICUS|REACTOME|WP|BIOCARTA|GOBP|GOMF|GOCC)_", "", term).lower().split("_")
    small = {"via", "and", "of", "in", "to", "the", "by", "for", "up", "dn", "vs", "or", "on", "at"}
    # Short words are usually acronyms (T, CD4, DNA, NFKB); common four-letter English words are not.
    plain = {"cell", "beta", "type", "like", "acid", "host", "self", "bone", "stem", "with", "from", "into", "cycle",
             "gene", "fate", "body", "tube", "node", "ion", "late", "early", "wall", "tail", "head", "limb", "lung",
             "skin", "axis", "site", "base", "pore", "ring", "rate", "size", "loss", "gain", "pair", "long", "fold",
             "lipid", "iron", "zinc", "heme", "bile", "salt", "drug", "pain", "heat", "cold", "male", "germ", "egg",
             "eye", "ear", "gut", "fat", "end", "cap", "arm", "leg", "non", "pre", "pro", "anti", "post", "down", "over"}
    words = [w.upper() if (len(w) <= 4 or any(ch.isdigit() for ch in w)) and w not in small and w not in plain
             else w for w in words]
    if words:
        words[0] = words[0][0].upper() + words[0][1:]
    return " ".join(words)


def wrap_label(term: str, width: int = 42, max_lines: int = 0) -> str:
    text = pretty_term(term)
    if text[:1].islower() and re.fullmatch(r"[a-z][a-z-]*", text.split(" ")[0]):
        text = text[0].upper() + text[1:]
    lines = textwrap.wrap(text, width) or [text]
    if max_lines and len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1][:max(1, width - 3)].rstrip(" ,;:-") + "..."
    return "\n".join(lines)


def dotplot(rows: List[dict], title: str, out_dir, name: str, top: int = 8) -> str:
    """Top terms per peak group: dot size = genes, colour = -log10 adjusted p."""
    groups = list(dict.fromkeys(r["group"] for r in rows))
    terms: List[str] = []
    for g in groups:
        for r in sorted((r for r in rows if r["group"] == g), key=lambda r: r["padj"])[:top]:
            if r["term"] not in terms:
                terms.append(r["term"])
    shown = [r for r in rows if r["term"] in terms]
    labels = [wrap_label(t, 40) for t in terms]
    n_lines = sum(label.count("\n") + 1 for label in labels)
    height = max(5.5, n_lines * 0.37 + 3)

    def draw(fig):
        cmap = palette().cmap()
        ax = fig.add_axes([0.44, 1.7 / height, 0.29, 1 - 2.6 / height])
        ypos = {t: len(terms) - 1 - i for i, t in enumerate(terms)}
        counts = np.array([r["count"] for r in shown], dtype=float)
        size = 40 + 360 * counts / counts.max()
        sc = ax.scatter([groups.index(r["group"]) for r in shown], [ypos[r["term"]] for r in shown], s=size,
                        c=[min(-math.log10(max(r["padj"], 1e-300)), 50) for r in shown], cmap=cmap, alpha=0.95,
                        edgecolors="#555555", linewidths=0.4)
        ax.set(xlim=(-0.6, len(groups) - 0.4), ylim=(-0.7, len(terms) - 0.3), title=title)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels(groups, rotation=45, ha="right")
        ax.set_yticks([ypos[t] for t in terms])
        ax.set_yticklabels(labels, fontsize=13)
        ax.grid(axis="y", color="#EBEBEB", linewidth=0.5)
        ax.set_axisbelow(True)
        cax = fig.add_axes([0.8, 0.55, 0.022, 0.25])
        fig.colorbar(sc, cax=cax).set_label("-log10 padj", fontsize=14)
        ticks = sorted({int(v) for v in np.unique(np.round(np.quantile(counts, [0, 0.5, 1])))} - {0}) or [1]
        handles = [ax.scatter([], [], s=40 + 360 * t / counts.max(), color="black") for t in ticks]
        fig.legend(handles, [str(t) for t in ticks], title="# genes", loc="center left", bbox_to_anchor=(0.77, 0.3),
                   labelspacing=1.2)
    return save_figure(draw, out_dir, name, 10.5, height)


# ---- pathway networks ------------------------------------------------------------------------

def _lighten(color: str, amount: float) -> str:
    from matplotlib.colors import to_hex, to_rgb
    c = np.array(to_rgb(color))
    return to_hex(c + (1 - c) * amount)


def _radius(values, rmin, rmax):
    v = np.sqrt(np.asarray(values, dtype=float))
    lo, hi = v.min(), v.max()
    return np.full(len(v), (rmin + rmax) / 2) if hi == lo else rmin + (rmax - rmin) * (v - lo) / (hi - lo)


NODE_FONT, MODULE_FONT = 14.0, 17.0          # enrichment map: term labels, module names (pt)


def _spring(n: int, links: List[Tuple[int, int, float]], iters: int = 300) -> np.ndarray:
    """Deterministic force-directed layout (Fruchterman-Reingold) of one module."""
    if n == 1:
        return np.zeros((1, 2))
    angle = 2 * math.pi * np.arange(n) / n
    pos = np.column_stack([np.cos(angle), np.sin(angle)]) * math.sqrt(n)
    for step in range(iters):
        delta = pos[:, None, :] - pos[None, :, :]
        dist = np.maximum(np.hypot(delta[..., 0], delta[..., 1]), 0.01)
        move = (delta / dist[..., None] * (1.0 / dist)[..., None]).sum(axis=1)          # repulsion
        for i, j, w in links:                                                           # attraction along edges
            pull = delta[i, j] * dist[i, j] * (0.5 + w)
            move[i] -= pull
            move[j] += pull
        move -= 0.05 * pos                                                              # keep the module together
        length = np.maximum(np.hypot(move[:, 0], move[:, 1]), 1e-9)
        temp = 0.3 * (1 - step / iters) + 0.01
        pos += move / length[:, None] * np.minimum(length, temp)[:, None]
    return pos - pos.mean(axis=0)


def _separate(pos: np.ndarray, half_w: np.ndarray, up: np.ndarray, down: np.ndarray, pad: float = 0.12) -> np.ndarray:
    """Push nodes apart until their boxes (circle plus the label under it) no longer overlap."""
    pos = pos.copy()
    n = len(pos)
    for _ in range(400):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                dx = pos[j, 0] - pos[i, 0]
                ci, cj = pos[i, 1] + (up[i] - down[i]) / 2, pos[j, 1] + (up[j] - down[j]) / 2
                dy = cj - ci
                ox = half_w[i] + half_w[j] + pad - abs(dx)
                oy = (up[i] + down[i] + up[j] + down[j]) / 2 + pad - abs(dy)
                if ox > 0 and oy > 0:
                    if ox < oy * 1.6:                      # prefer sideways moves: labels are wide, not tall
                        sgn = 1 if dx >= 0 else -1
                        pos[i, 0] -= sgn * (ox / 2 + 0.005)
                        pos[j, 0] += sgn * (ox / 2 + 0.005)
                    else:
                        sgn = 1 if dy >= 0 else -1
                        pos[i, 1] -= sgn * (oy / 2 + 0.005)
                        pos[j, 1] += sgn * (oy / 2 + 0.005)
                    moved = True
        if not moved:
            break
    return pos


def _rounded_hull(points: np.ndarray, pad: float) -> np.ndarray:
    """Convex hull of the points, grown by `pad` with round corners."""
    pts = sorted(set(map(tuple, np.round(points, 6))))

    def half(seq):
        out: List[Tuple[float, float]] = []
        for p in seq:
            while len(out) >= 2 and ((out[-1][0] - out[-2][0]) * (p[1] - out[-2][1])
                                     - (out[-1][1] - out[-2][1]) * (p[0] - out[-2][0])) <= 0:
                out.pop()
            out.append(p)
        return out[:-1]
    hull = np.array(half(pts) + half(pts[::-1]))                      # counter-clockwise
    n = len(hull)
    edges = np.roll(hull, -1, axis=0) - hull
    normal = np.arctan2(-edges[:, 0], edges[:, 1])                    # outward normal angle of each edge
    out = []
    for i in range(n):
        a0, a1 = normal[i - 1], normal[i]
        while a1 < a0:
            a1 += 2 * math.pi
        for a in np.linspace(a0, a1, 10):
            out.append((hull[i, 0] + pad * math.cos(a), hull[i, 1] + pad * math.sin(a)))
    return np.array(out)


def layout_modules(nodes: List[dict], edges: Sequence[dict] = (), gap: float = 0.5, aspect: float = 1.35):
    """Lay each module out by its own links (force-directed, labels kept apart), outline it with a
    rounded hull and pack the modules in rows (units: inches, y downwards from 0).
    Returns node positions, module outlines {id: polygon}, module label anchors {id: (x, y)},
    node radii, wrapped labels, the canvas size and the wrapped module names."""
    radius = dict(zip([n["term"] for n in nodes], _radius([n["count"] for n in nodes], 0.1, 0.24)))
    text = {n["term"]: wrap_label(n["term"], 18) for n in nodes}       # long terms wrap, never cut
    char, line = 0.55 * NODE_FONT / 72, 1.2 * NODE_FONT / 72
    modules: Dict[int, List[dict]] = {}
    for nd in nodes:
        modules.setdefault(nd["module"], []).append(nd)
    blocks, titles = [], {}
    for mid in sorted(modules, key=lambda m: (-len(modules[m]), m)):
        members = [n["term"] for n in modules[mid]]
        index = {t: i for i, t in enumerate(members)}
        links = [(index[e["from"]], index[e["to"]], float(e["jaccard"])) for e in edges
                 if e["from"] in index and e["to"] in index]
        half_w = np.array([max(max(len(l) for l in text[t].split("\n")) * char / 2 + 0.04, radius[t]) for t in members])
        up = np.array([radius[t] for t in members])
        down = np.array([radius[t] + 0.06 + (text[t].count("\n") + 1) * line for t in members])
        pos = _separate(_spring(len(members), links) * 1.6, half_w, up, down, pad=0.2)
        corners = np.concatenate([np.column_stack([pos[:, 0] + sx * half_w, pos[:, 1] + sy]) for sx in (-1, 1)
                                  for sy in (up, -down)])
        named = len(members) > 1
        hull = _rounded_hull(corners, 0.3) if named else corners
        x0, y0, x1, y1 = hull[:, 0].min(), hull[:, 1].min(), hull[:, 0].max(), hull[:, 1].max()
        mchar = 0.56 * MODULE_FONT / 72
        # The module name wraps onto as many rows as it needs within the width of its outline.
        titles[mid] = wrap_label(f"Module {mid}: {pretty_term(modules[mid][0]['module_label'])}",
                                 max(24, int((x1 - x0) / mchar)))
        rows = titles[mid].split("\n")
        label_h = 0.2 + 1.2 * MODULE_FONT / 72 * len(rows) if named else 0.0
        label_w = max(len(l) for l in rows) * mchar
        width = max(x1 - x0, label_w if named else 0)
        blocks.append((mid, members, pos, hull, named, x0 - (width - (x1 - x0)) / 2, y1, width, y1 - y0 + label_h, label_h))
    total = sum(b[7] * b[8] for b in blocks)
    row_width = max(max(b[7] for b in blocks), math.sqrt(total * aspect))
    positions, outlines, anchors = {}, {}, {}
    x = y = row_h = canvas_w = 0.0
    for mid, members, pos, hull, named, left, top, w, h, label_h in blocks:
        if x > 0 and x + w > row_width + 1e-9:
            x, y, row_h = 0.0, y + row_h + gap, 0.0
        shift = np.array([x - left, -(y + label_h) - top])            # block's top-left corner -> (x, -y)
        for t, p in zip(members, pos + shift):
            positions[t] = (float(p[0]), float(p[1]))
        if named:
            outlines[mid] = hull + shift
            anchors[mid] = (x + w / 2, -(y + label_h) + 0.08)
        x += w + gap
        canvas_w, row_h = max(canvas_w, x - gap), max(row_h, h)
    return positions, outlines, anchors, radius, text, canvas_w, y + row_h, titles


def enrichment_map_figure(nodes: List[dict], edges: List[dict], title: str, out_dir, name: str):
    """Enrichment map: terms linked by shared genes, each module inside a tinted, named outline.
    Returns the figure path, node positions and the colour of each term."""
    pos, outlines, anchors, radii, text, width, height, titles = layout_modules(nodes, edges)
    colour = dict(zip(sorted(outlines), palette().colors(len(outlines))))
    node_colour = {n["term"]: colour.get(n["module"], "#BDBDBD") for n in nodes}
    fig_w, fig_h = max(width, 9.5) + 0.8, height + 2.3

    def draw(fig):
        from matplotlib import patheffects
        from matplotlib.patches import Circle, Polygon
        halo = [patheffects.withStroke(linewidth=3, foreground="white")]
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set(xlim=(-0.4, fig_w - 0.4), ylim=(-height - 1.5, 0.8), aspect="equal")
        ax.axis("off")
        ax.text(-0.2, 0.5, title, fontsize=20, fontweight="bold", va="center")
        for mid, poly in outlines.items():
            ax.add_patch(Polygon(poly, closed=True, facecolor=_lighten(colour[mid], 0.85),
                                 edgecolor=_ink(colour[mid], 0.6), linewidth=1.2, joinstyle="round", zorder=0))
            ax.text(*anchors[mid], titles[mid], fontsize=MODULE_FONT, fontweight="bold", color=_ink(colour[mid]),
                    ha="center", va="bottom", linespacing=1.1)
        for e in edges:
            (x0, y0), (x1, y1) = pos[e["from"]], pos[e["to"]]
            same = node_colour[e["from"]] == node_colour[e["to"]]
            ax.plot([x0, x1], [y0, y1], color="#555555" if same else "#BBBBBB",
                    linewidth=(0.5 + 2.2 * float(e["jaccard"])) * (1 if same else 0.6),
                    alpha=0.45 if same else 0.3, linestyle="-" if same else (0, (4, 3)), zorder=1)
        for n in nodes:
            x, y = pos[n["term"]]
            ax.add_patch(Circle((x, y), radii[n["term"]], facecolor=node_colour[n["term"]],
                                edgecolor=_ink(node_colour[n["term"]], 0.35), linewidth=0.6, zorder=2))
            ax.text(x, y - radii[n["term"]] - 0.04, text[n["term"]], fontsize=NODE_FONT, ha="center",
                    va="top", linespacing=1.05, zorder=3, path_effects=halo)
        ax.text(-0.2, -height - 0.55, "Node size: genes in the term.  Edge width: shared genes (Jaccard); dashed edges "
                "link modules.\nModule: a group of terms that share genes, numbered by significance (Module 1 holds "
                "the most\nsignificant term) and named after its most significant term.",
                fontsize=13, color="#444444", va="top", linespacing=1.3)
    fig_path = save_figure(draw, out_dir, name, fig_w, fig_h)
    return fig_path, pos, {t: c for t, c in node_colour.items() if c != "#BDBDBD"}


def gene_concept_figure(terms: List[dict], links: List[Tuple[str, str]], term_colour: Dict[str, str], title: str,
                        out_dir, name: str) -> str:
    """Gene-pathway network: pathways on a ring, each gene near the pathways it belongs to."""
    n = len(terms)
    ring = 2.6 + 0.3 * n
    centre = {t["term"]: (ring * math.cos(math.pi / 2 - 2 * math.pi * i / n), ring * math.sin(math.pi / 2 - 2 * math.pi * i / n))
              for i, t in enumerate(terms)}
    members: Dict[str, List[str]] = {}
    for term, gene in links:
        members.setdefault(gene, []).append(term)
    own: Dict[str, List[str]] = {t["term"]: [] for t in terms}
    shared = []
    for gene, ts in members.items():
        (own[ts[0]] if len(ts) == 1 else shared).append(gene)
    pos: Dict[str, Tuple[float, float]] = {}
    for term, genes in own.items():
        cx, cy = centre[term]
        out_angle = math.atan2(cy, cx)
        for k, gene in enumerate(sorted(genes)):
            layer, per = divmod(k, 7)
            a = out_angle + (per - 3) * 0.47
            r = 1.05 + 0.55 * layer
            pos[gene] = (cx + r * math.cos(a), cy + r * math.sin(a))
    shared.sort(key=lambda g: (-len(members[g]), g))
    pts = []
    for k, gene in enumerate(shared):
        x = sum(centre[t][0] for t in members[gene]) / len(members[gene]) * 0.4
        y = sum(centre[t][1] for t in members[gene]) / len(members[gene]) * 0.4
        a = 2.399963 * k                                   # golden-angle spiral start
        r = 0.5 * math.sqrt(k + 1)
        pts.append([x + r * math.cos(a), y + r * math.sin(a)])
    for _ in range(60):                                    # push shared genes apart so labels do not collide
        moved = False
        for i in range(len(pts)):
            for j in range(i + 1, len(pts)):
                dx, dy = pts[j][0] - pts[i][0], pts[j][1] - pts[i][1]
                if abs(dx) < 0.8 and abs(dy) < 0.38:
                    push_x, push_y = (0.8 - abs(dx)) / 2 + 0.01, (0.38 - abs(dy)) / 2 + 0.01
                    if push_y <= push_x:
                        sgn = 1 if dy >= 0 else -1
                        pts[i][1] -= sgn * push_y
                        pts[j][1] += sgn * push_y
                    else:
                        sgn = 1 if dx >= 0 else -1
                        pts[i][0] -= sgn * push_x
                        pts[j][0] += sgn * push_x
                    moved = True
        if not moved:
            break
    for gene, pt in zip(shared, pts):
        pos[gene] = (pt[0], pt[1])
    label_pos = {}
    for term, genes in own.items():                        # pathway name beyond its own genes
        cx, cy = centre[term]
        a = math.atan2(cy, cx)
        r = 1.05 + 0.55 * ((len(genes) - 1) // 7 if genes else -1) + 0.85
        label_pos[term] = (cx + r * math.cos(a), cy + r * math.sin(a))
    extent = ring + 2.9
    radii = dict(zip([t["term"] for t in terms], _radius([t["count"] for t in terms], 0.2, 0.34)))

    def draw(fig):
        from matplotlib.patches import Circle
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set(xlim=(-extent - 1.0, extent + 1.0), ylim=(-extent - 0.3, extent + 0.5), aspect="equal")
        ax.axis("off")
        ax.text(-extent - 0.8, extent + 0.15, title, fontsize=17, fontweight="bold", va="center")
        for term, gene in links:
            (x0, y0), (x1, y1) = centre[term], pos[gene]
            ax.plot([x0, x1], [y0, y1], color=_ink(term_colour[term], 0.6), linewidth=0.5, alpha=0.45, zorder=1)
        for gene, (x, y) in pos.items():
            hub = len(members[gene]) > 1
            fill = "#333333" if hub else _lighten(term_colour[members[gene][0]], 0.45)
            ax.add_patch(Circle((x, y), 0.075 if hub else 0.055, facecolor=fill,
                                edgecolor="white" if hub else _ink(fill, 0.5), linewidth=0.4, zorder=2))
            ax.text(x, y - 0.1, gene, fontsize=10, ha="center", va="top", fontweight="bold" if hub else "normal", zorder=3)
        for t in terms:
            x, y = centre[t["term"]]
            ax.add_patch(Circle((x, y), radii[t["term"]], facecolor=term_colour[t["term"]],
                                edgecolor=_ink(term_colour[t["term"]], 0.35), linewidth=0.8, zorder=4))
            lx, ly = label_pos[t["term"]]
            ax.text(lx, ly, wrap_label(t["term"], 20, 3), fontsize=12.5, fontweight="bold", ha="center", va="center",
                    color=_ink(term_colour[t["term"]]), zorder=5,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 1.5})
        ax.text(-extent - 0.8, -extent, "Dark genes belong to several of the pathways shown (hubs).", fontsize=12,
                color="#444444")
    return save_figure(draw, out_dir, name, (2 * extent + 2.0) * 0.75, (2 * extent + 0.8) * 0.75)
