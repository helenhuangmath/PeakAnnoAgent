"""Shared helpers of the analysis engine: tables, the cache, downloads and figure output."""
from __future__ import annotations

import contextvars
import csv
import os
import re
import shutil
import urllib.request
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence

# Figure colours. "chippeakanno" (the default colours of ChIPpeakAnno's genomicElementDistribution)
# is the default; "npg" is the Nature Publishing Group palette; "spectral" is ColorBrewer
# Spectral; any Matplotlib colormap name (viridis, Set2, Dark2, RdYlBu_r...) also works.
PALETTE_OPTIONS = {"chippeakanno": "ChIPpeakAnno (default)", "npg": "NPG, Nature style", "spectral": "Spectral",
                   "Set2": "Set2 (soft, colour-blind friendly)", "Dark2": "Dark2", "tab10": "tab10 (Matplotlib default)",
                   "viridis": "viridis (colour-blind friendly)", "RdYlBu_r": "RdYlBu (red-yellow-blue)"}
NET_PALETTE = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148",
               "#B09C85", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#EE4C97", "#0072B5", "#BC3C29", "#925E9F",
               "#FDAF91", "#AD002A"]

DEFAULT_PALETTE = "chippeakanno"
PREVIEW_DPI = 300      # resolution of the PNG copies embedded in the HTML report

# ChIPpeakAnno's defaults: `labelColors` and the promoter-level colours of genomicElementDistribution().
CHIPPEAKANNO = {"promoter": "#D55E00", "geneDownstream": "#E69F00", "geneBody": "#51C6E6",
                "distalIntergenic": "#AAAAAA", "exon": "#009DDA", "intron": "#666666", "intergenic": "#DDDDDD",
                "utr5": "#0072B2", "utr3": "#56B4E9", "CDS": "#0033BF", "otherExon": "#009E73"}
CHIPPEAKANNO_PROMOTER = ["#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32"]     # far from the TSS -> at the TSS

_palette = contextvars.ContextVar("peakannoagent_palette", default=DEFAULT_PALETTE)
_previews = contextvars.ContextVar("peakannoagent_previews", default=False)


def resolve_palette(name: Optional[str]) -> str:
    """The canonical name of a palette, or a ValueError that lists the choices."""
    key = (name or DEFAULT_PALETTE).strip()
    if key.lower() in ("chippeakanno", "spectral", "npg"):
        return key.lower()
    import matplotlib
    for known in matplotlib.colormaps:
        if known.lower() == key.lower():
            return known
    raise ValueError(f"Unknown palette '{name}'. Choose chippeakanno, npg, spectral, or a Matplotlib colormap name "
                     f"such as Set2, Dark2, tab10, viridis or RdYlBu_r.")


def set_palette(name: Optional[str]) -> str:
    """Choose the colours of the figures drawn from here on (in this thread)."""
    name = resolve_palette(name)
    _palette.set(name)
    return name


def set_previews(on: bool) -> None:
    """Also write a PNG of each figure to `<out_dir>/.preview` (the HTML report embeds them)."""
    _previews.set(bool(on))


class Palette:
    """Colours for one palette: `colors(n)` for categories, `cmap()` for significance, and the
    single colours used for bars, lines and reference lines."""

    def __init__(self, name: str):
        import matplotlib
        from matplotlib.colors import LinearSegmentedColormap, ListedColormap, to_hex
        self.name = name
        self._list: Optional[List[str]] = None
        self._ordered = None
        self.features: Dict[str, str] = {}     # fixed colour per genomic feature, where the palette has one
        self.promoters: List[str] = []         # promoter bins, nearest the TSS first
        if name == "chippeakanno":
            c = CHIPPEAKANNO
            # No greys (grey is kept for "unassigned"): where ChIPpeakAnno uses grey for introns and
            # intergenic regions, tints of its green and blue are used instead.
            self._list = [c["promoter"], c["utr5"], c["otherExon"], c["geneDownstream"], c["utr3"], c["CDS"],
                          "#FF8E32", c["exon"], "#7FCEB9", c["geneBody"]]
            self.features = {"5' UTR": c["utr5"], "3' UTR": c["utr3"], "1st Exon": c["CDS"],
                             "Other Exon": c["otherExon"], "1st Intron": "#7FCEB9", "Other Intron": "#BFE7DC",
                             "Downstream": c["geneDownstream"], "Distal Intergenic": "#CFE8F7"}
            self.promoters = list(reversed(CHIPPEAKANNO_PROMOTER))
            self._seq = LinearSegmentedColormap.from_list("chippeakanno", CHIPPEAKANNO_PROMOTER + [c["promoter"]])
            self._ordered = LinearSegmentedColormap.from_list(
                "chippeakanno_ordered", [c["promoter"], "#FF8E32", "#FFCA99", c["utr3"], c["utr5"], c["CDS"]])
            self.bar = self.primary = c["utr5"]
            self.accent = c["promoter"]
            return
        if name == "spectral":
            self._cmap, self._range = matplotlib.colormaps["Spectral"], (0.0, 1.0)
            self._seq = matplotlib.colormaps["Spectral_r"]
            self.bar = self.primary = "#3288BD"
            self.accent = "#D53E4F"
            return
        if name == "npg":
            self._list = [c for c in NET_PALETTE if c != "#DC0000"] + ["#DC0000"]
            self._seq = LinearSegmentedColormap.from_list("npg", ["#FBE3DE", "#F39B7F", "#E64B35", "#A52A1A"])
            self._ordered = LinearSegmentedColormap.from_list(
                "npg_ordered", ["#E64B35", "#F39B7F", "#91D1C2", "#4DBBD5", "#3C5488"])
            self.bar = self.primary = "#3C5488"
            self.accent = "#E64B35"
            return
        cmap = matplotlib.colormaps[name]
        if isinstance(cmap, ListedColormap) and cmap.N <= 20:              # qualitative: Set2, tab10...
            self._list = [to_hex(c) for c in cmap.colors]
            self.bar = self.primary = self._list[0]
            self.accent = self._list[1 % len(self._list)]
            self._seq = LinearSegmentedColormap.from_list(name, ["#EFEFEF", self._list[0]])
        else:                                                              # continuous: viridis, RdYlBu...
            self._cmap, self._range, self._seq = cmap, (0.08, 0.92), cmap
            self.bar = self.primary = to_hex(cmap(0.75))
            self.accent = to_hex(cmap(0.2))

    def colors(self, n: int) -> List[str]:
        if self._list is not None:
            return [self._list[i % len(self._list)] for i in range(n)]
        from matplotlib.colors import to_hex
        lo, hi = self._range
        return [to_hex(self._cmap(lo + (hi - lo) * i / max(1, n - 1))) for i in range(n)]

    @property
    def qualitative(self) -> bool:
        return self._list is not None

    def ordered(self, n: int) -> List[str]:
        """Colours for ordered categories (distance bins): a ramp where the palette has one."""
        if self._ordered is None:
            return self.colors(n)
        from matplotlib.colors import to_hex
        return [to_hex(self._ordered(i / max(1, n - 1))) for i in range(n)]

    def cmap(self):
        return self._seq


def palette() -> Palette:
    return Palette(_palette.get())


OFFLINE = False        # set per run: no downloads; cached files only


class SkipTool(Exception):
    """Raised by a tool that cannot run here for a reason that is not an error (e.g. HOMER absent)."""


def cache_dir(path: Optional[str] = None) -> Path:
    """Folder for downloaded gene models, gene-set libraries and built references."""
    root = path or os.environ.get("PEAKANNOAGENT_CACHE")
    if not root:
        base = os.environ.get("XDG_CACHE_HOME") or str(Path.home() / ".cache")
        root = str(Path(base) / "peakannoagent")
    out = Path(root).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    return out


def _ssl_contexts():
    """Certificate stores to try, in order: Python's default, certifi, the system bundle. Some
    Python installations (python.org on macOS) ship without root certificates; verification is
    never turned off."""
    import ssl
    yield ssl.create_default_context()
    try:
        import certifi
        yield ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        pass
    for bundle in ("/etc/ssl/cert.pem", "/etc/ssl/certs/ca-certificates.crt", "/etc/pki/tls/certs/ca-bundle.crt"):
        if Path(bundle).is_file():
            yield ssl.create_default_context(cafile=bundle)


def download(url: str, dest: Path, what: str = "file") -> Path:
    """Download once; a partial download never replaces a good file."""
    import ssl
    dest = Path(dest)
    if dest.is_file() and dest.stat().st_size > 0:
        return dest
    if OFFLINE:
        raise RuntimeError(f"The {what} is not in the cache ({dest}) and this run is offline. Run once with "
                           f"network access, or pass a local file.")
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + f".{os.getpid()}.part")
    error: Exception = RuntimeError("no certificate store")
    try:
        for context in _ssl_contexts():
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "peakannoagent"})
                with urllib.request.urlopen(request, timeout=120, context=context) as response, open(tmp, "wb") as handle:
                    shutil.copyfileobj(response, handle, 1 << 20)
                if tmp.stat().st_size == 0:
                    raise OSError("empty response")
                tmp.replace(dest)
                return dest
            except Exception as err:
                error = err
                reason = getattr(err, "reason", err)
                if not isinstance(reason, ssl.SSLCertVerificationError):
                    break                      # only a missing certificate store is worth another try
    finally:
        tmp.unlink(missing_ok=True)
    raise RuntimeError(f"Could not download the {what} from {url} ({error}). "
                       f"Download it yourself and pass the file instead.")


def write_tsv(path, columns: Sequence[str], rows: Iterable[Sequence]) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE, escapechar="\\")
        writer.writerow(columns)
        for row in rows:
            writer.writerow(["NA" if v is None else v for v in row])
    return str(path)


def read_tsv(path) -> List[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def file_tag(text: str) -> str:
    """Method names such as "msigdb_C7:IMMUNESIGDB" -> safe file-name parts."""
    return re.sub(r"[^A-Za-z0-9_.-]+", "-", text)


def signif(value: float, digits: int = 3) -> float:
    return float(f"{value:.{digits}g}")


def bh_adjust(pvalues: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg adjusted p-values, in the input order."""
    n = len(pvalues)
    order = sorted(range(n), key=lambda i: pvalues[i], reverse=True)
    out, running = [0.0] * n, 1.0
    for rank, i in zip(range(n, 0, -1), order):
        running = min(running, pvalues[i] * n / rank)
        out[i] = running
    return out


def _style():
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",       # text stays editable
        "font.size": 15, "axes.titlesize": 18, "axes.titleweight": "bold", "axes.labelsize": 16,
        "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "legend.title_fontsize": 15,
        "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.8,
        "axes.edgecolor": "black", "text.color": "black", "axes.labelcolor": "black",
        "xtick.color": "black", "ytick.color": "black", "xtick.direction": "out", "ytick.direction": "out",
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    import matplotlib.pyplot as plt
    return plt


def save_figure(draw: Callable, out_dir, name: str, width: float = 7, height: float = 4.5) -> str:
    """Draw once and write a vector PDF in `out_dir`; `draw(fig)` adds axes to the figure. Returns
    the PDF path. With previews on (full runs), a PNG is also written to `out_dir/.preview`."""
    plt = _style()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(width, height))
    try:
        draw(fig)
        fig.savefig(out_dir / f"{name}.pdf")
        if _previews.get():
            (out_dir / ".preview").mkdir(exist_ok=True)
            fig.savefig(preview_path(out_dir / f"{name}.pdf"), dpi=PREVIEW_DPI)
    finally:
        plt.close(fig)
    return str(out_dir / f"{name}.pdf")


def preview_path(pdf) -> Path:
    """Where the PNG preview of a figure PDF is written during a full run."""
    pdf = Path(pdf)
    return pdf.parent / ".preview" / f"{pdf.stem}.png"
