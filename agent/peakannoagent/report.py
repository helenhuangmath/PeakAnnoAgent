"""Self-contained HTML report. Figures are embedded; tables are previewed and linked.
Narrative citations are checked against the EvidenceStore: unknown ones are
visibly flagged as unverified rather than silently trusted."""

from __future__ import annotations

import base64
import html
import json
import os
import re
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from . import __version__
from .analysis import read_table
from .evidence import CITATION_RE, Evidence, EvidenceStore
from .explore import explore_data, explore_js, has_explore, inline_json

IGV_JS = "https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"
HIDDEN_TOOLS = ("inspect_table", "chipatlas_status")
TOOL_SECTIONS = [
    ("Peak QC", ("peak_qc",)),
    ("Annotation", ("annotate_peaks",)),
    ("Pathways", ("enrich_pathways", "network_viz", "pathway_consensus", "great_enrichment")),
    ("Motifs and regulators", ("find_motifs", "chipatlas_submit", "chipatlas_collect", "chipatlas_target_genes")),
    ("Integration and controls", ("integrate_rnaseq", "check_expectations")),
    ("Literature and public data", ("search_literature", "search_geo")),
]


def esc(x) -> str:
    return html.escape("" if x is None else str(x))


# ------------------------------------------------------------------ narrative markdown

def _cite_html(token: str, store: EvidenceStore) -> str:
    t = token.replace(" ", "")
    ok = store.known_citation(t)
    if t.startswith("PMID:"):
        href = f"https://pubmed.ncbi.nlm.nih.gov/{t[5:]}/"
    elif re.fullmatch(r"E\d+", t):
        href = f"#{t}"
    elif t.startswith(("GSE", "GSM", "GDS")):
        href = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={t}"
    else:
        href = f"https://chip-atlas.org/view?id={t}"
    cls = "cite" if ok else "cite unverified"
    title = "" if ok else ' title="Not produced by any tool in this run - unverified"'
    target = "" if href.startswith("#") else ' target="_blank" rel="noopener"'
    return f'<a class="{cls}" href="{href}"{target}{title}>{esc(t)}</a>'


def md_inline(s: str, store: EvidenceStore) -> str:
    """Inline markdown (**bold**, *italic*, `code`) plus linked, verified citations."""
    s = esc(s)
    s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?!\w)", r"<em>\1</em>", s)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    return CITATION_RE.sub(
        lambda m: "[" + ", ".join(_cite_html(t, store) for t in re.split(r"\s*[,;]\s*", m.group(1))) + "]", s)


def md(text: str, store: EvidenceStore) -> str:
    """Minimal block markdown: paragraphs, '- ' / '1.' lists and '#' headings."""
    if not text:
        return ""
    out = []
    bullet = re.compile(r"^\s*([-*]|\d+\.)\s+")
    for block in re.split(r"\n\s*\n", text.strip()):
        lines = [ln for ln in block.splitlines() if ln.strip()]
        if all(bullet.match(ln) for ln in lines):
            out.append("<ul>" + "".join(f"<li>{md_inline(bullet.sub('', ln), store)}</li>" for ln in lines) + "</ul>")
        elif lines[0].startswith("#"):
            level = min(len(lines[0]) - len(lines[0].lstrip("#")) + 2, 5)
            out.append(f"<h{level}>{md_inline(lines[0].lstrip('#').strip(), store)}</h{level}>")
            if len(lines) > 1:
                out.append(f"<p>{md_inline(' '.join(lines[1:]), store)}</p>")
        else:
            out.append(f"<p>{md_inline(' '.join(lines), store)}</p>")
    return "\n".join(out)


# ------------------------------------------------------------------ evidence rendering

def _png_width(raw: bytes) -> int:
    """Pixel width from a PNG header (0 if not a PNG)."""
    return int.from_bytes(raw[16:20], "big") if raw[:8] == b"\x89PNG\r\n\x1a\n" else 0


def _img(path: str, report_dir: Path | None = None) -> str:
    try:
        raw = Path(path).read_bytes()
    except OSError:
        return ""
    data = base64.b64encode(raw).decode()
    # Previews are 120 dpi; show them at ~80 px per inch so a figure's type size stays constant
    # whatever its canvas size (wide figures are still capped at the column width).
    width = _png_width(raw) * 2 // 3
    pdf = Path(path).with_suffix(".pdf")
    if Path(path).parent.name == ".preview":          # previews sit one level below their PDFs
        pdf = Path(path).parent.parent / pdf.name
    links = [f'<a href="data:image/png;base64,{data}" download="{esc(Path(path).stem)}.png">PNG</a>']
    if report_dir is not None and pdf.is_file():
        links.append(f'<a href="{esc(_rel(str(pdf), report_dir))}" download>PDF (vector)</a>')
    caption = f'<figcaption><span>{esc(Path(path).stem.replace("_", " "))}</span><span class="dl">Download: {" · ".join(links)}</span></figcaption>'
    style = f' style="width:{width}px"' if width else ""
    return (f'<figure><img alt="{esc(Path(path).stem)}"{style} src="data:image/png;base64,{data}">'
            f'{caption}</figure>')


def _rel(path: str, report_dir: Path) -> str:
    try:
        return os.path.relpath(path, report_dir)
    except ValueError:
        return path


def _fmt(v) -> str:
    if isinstance(v, float):
        if v != 0 and (abs(v) < 1e-3 or abs(v) >= 1e6):
            return f"{v:.2e}"
        return f"{v:.3g}" if abs(v) < 100 else f"{v:,.0f}"
    if isinstance(v, list):
        if all(not isinstance(x, (dict, list)) for x in v):
            return esc(", ".join(str(x) for x in v[:40]) + (" …" if len(v) > 40 else ""))
    return esc(v)


def _table(rows: List[dict], max_rows: int = 15, cls: str = "data") -> str:
    if not rows:
        return ""
    cols = list(rows[0].keys())
    head = "".join(f"<th>{esc(c)}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td>{_fmt(r.get(c))}</td>" for c in cols) + "</tr>" for r in rows[:max_rows])
    more = f'<p class="muted">Showing {min(max_rows, len(rows))} of {len(rows)} rows.</p>' if len(rows) > max_rows else ""
    return f'<div class="scroll"><table class="{cls}"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>{more}'


def _summary_html(summary: dict) -> str:
    scalars, blocks = [], []
    for key, value in summary.items():
        label = key.replace("_", " ")
        if isinstance(value, list) and value and isinstance(value[0], dict):
            blocks.append(f"<h4>{esc(label)}</h4>" + _table(value, 12))
        elif isinstance(value, dict) and value and all(not isinstance(v, (dict, list)) for v in value.values()):
            scalars.append((label, ", ".join(f"{k}: {_fmt(v)}" for k, v in value.items())))
        elif isinstance(value, dict) and value:
            blocks.append(f"<h4>{esc(label)}</h4>" + _summary_html(value))
        elif isinstance(value, str) and os.path.isabs(value):
            continue                      # file paths belong in the table links, not the summary
        elif value not in (None, "", [], {}):
            scalars.append((label, _fmt(value)))
    kv = "".join(f"<tr><th>{esc(k)}</th><td>{v}</td></tr>" for k, v in scalars)
    return (f'<table class="kv">{kv}</table>' if kv else "") + "".join(blocks)


def _paper_list(papers: List[dict]) -> str:
    items = []
    for p in papers:
        link = f"https://pubmed.ncbi.nlm.nih.gov/{p['pmid']}/" if p.get("pmid") else (
            f"https://doi.org/{p['doi']}" if p.get("doi") else "#")
        tag = ' <span class="tag">preprint</span>' if p.get("preprint") else ""
        items.append(f'<li><a href="{link}" target="_blank" rel="noopener">{esc(p.get("title"))}</a>{tag}'
                     f'<br><span class="muted">{esc(p.get("authors"))} · <em>{esc(p.get("journal"))}</em> '
                     f'{esc(p.get("year"))} · {esc(p.get("id"))}</span></li>')
    return "<ol class='refs'>" + "".join(items) + "</ol>"


def _dataset_list(datasets: List[dict]) -> str:
    items = []
    for d in datasets:
        pm = ", ".join(f"PMID:{x}" for x in d.get("pubmed_ids", [])[:3])
        items.append(f'<li><a href="{esc(d.get("url"))}" target="_blank" rel="noopener">{esc(d.get("accession"))}'
                     f'</a> - {esc(d.get("title"))}<br><span class="muted">{esc(d.get("type"))} · '
                     f'{esc(d.get("n_samples"))} samples · {esc(d.get("date"))} {esc(pm)}</span></li>')
    return "<ul class='refs'>" + "".join(items) + "</ul>"


def evidence_card(ev: Evidence, report_dir: Path) -> str:
    status_cls = {"ok": "ok", "error": "err", "skipped": "skip"}.get(ev.status, "skip")
    parts = [f'<article class="card" id="{ev.id}"><header><span class="eid">{ev.id}</span> '
             f'<h3>{esc(ev.title)}</h3><span class="status {status_cls}">{esc(ev.status)}</span></header>']
    if ev.rationale:
        parts.append(f'<p class="why">Why: {esc(ev.rationale)}</p>')
    if ev.error:
        parts.append(f'<pre class="error">{esc(ev.error)}</pre>')
    summary = dict(ev.summary)
    papers = summary.pop("papers", None)
    datasets = summary.pop("datasets", None)
    summary.pop("note", None)
    interactive = summary.pop("interactive", None) or {}
    summary.pop("cytoscape", None)
    if has_explore(ev):
        def table(name: str) -> List[dict]:
            path = Path(ev.tables[name])
            for cand in (path, report_dir / "results" / path.name, report_dir / path.name):
                if cand.is_file():
                    try:
                        return read_table(str(cand))
                    except (OSError, UnicodeDecodeError):
                        return []
            return []
        parts.append(f'<section class="explore-box"><h4>Explore the pathways</h4><div class="explore" id="explore-{ev.id}">'
                     f'<p class="muted">Interactive view needs JavaScript.</p></div></section>'
                     f'<script type="application/json" id="explore-data-{ev.id}">{inline_json(explore_data(ev, table))}</script>')
    if summary:
        parts.append(_summary_html(summary))
    if interactive:
        parts.append("<p>" + " · ".join(f'<a href="{esc(_rel(v, report_dir))}" target="_blank" rel="noopener">'
                                        f'{esc(k.replace("_", " "))} (interactive) ↗</a>' for k, v in interactive.items()) + "</p>")
    if papers:
        parts.append(_paper_list(papers))
    if datasets:
        parts.append(_dataset_list(datasets))
    for fig in ev.figures.values():
        parts.append(_img(fig, report_dir))
    if ev.warnings:
        parts.append("<ul class='warn'>" + "".join(f"<li>{esc(w)}</li>" for w in ev.warnings) + "</ul>")
    for name, path in ev.tables.items():
        preview = ""
        try:
            preview = _table(read_table(path), 10)
        except (OSError, UnicodeDecodeError):
            pass
        parts.append(f'<details><summary>Table: {esc(name)} - <a href="{esc(_rel(path, report_dir))}">'
                     f'{esc(Path(path).name)}</a></summary>{preview}</details>')
    parts.append("</article>")
    return "\n".join(parts)


# ------------------------------------------------------------------ methods text

def methods_text(store: EvidenceStore, config) -> str:
    versions = {}
    for ev in store.items.values():
        versions.update(ev.versions or {})
    v = lambda p: f" (v{versions[p]})" if p in versions else ""
    s = []
    ok = lambda tool: any(e.tool == tool and e.status == "ok" for e in store.items.values())
    if ok("peak_qc"):
        s.append("Peaks were filtered to standard chromosomes, de-duplicated" +
                 (" and cleared of blacklisted regions" if config.blacklist else "") + ".")
    ann = store.latest("annotate_peaks")
    if ann:
        w = ann.summary.get("promoter_window", {})
        s.append(f"Peaks were annotated with ChIPpeakAnno{v('ChIPpeakAnno')} against "
                 f"{next((p for p in versions if p.startswith('TxDb')), 'the UCSC knownGene')} gene models: each "
                 f"peak was assigned its nearest gene by TSS distance, and all genes whose TSS window "
                 f"(-{w.get('upstream', 2000)}/+{w.get('downstream', 500)} bp) overlapped the peak were recorded as "
                 f"promoter targets. Genomic-region distribution used assignChromosomeRegion and the TSS profile "
                 f"used binOverFeature.")
    if ok("enrich_pathways"):
        s.append(f"Over-representation of peak-linked genes was tested with clusterProfiler{v('clusterProfiler')} "
                 f"(GO via org.*.eg.db, KEGG, ReactomePA{v('ReactomePA')}, MSigDB via msigdbr{v('msigdbr')}) "
                 "against all annotated genes as background, with Benjamini-Hochberg correction.")
    nv = store.latest("network_viz")
    if nv:
        s.append(f"Pathway networks were drawn with SimpNetworkViz{v('SimpNetworkViz')}: enriched terms were linked "
                 "by shared hit genes (combined Jaccard/overlap similarity ≥ 0.375), clustered with Markov clustering "
                 "and labelled automatically, following EnrichmentMap and AutoAnnotate.")
    pc = store.latest("pathway_consensus")
    if pc:
        s.append(f"Peak-linked gene lists were also analysed with PathwayAgent{v('PathwayAgent')} "
                 f"({', '.join(pc.summary.get('tools_ok') or [])}); terms were matched across engines and combined "
                 "into a consensus score, grouped into themes by hit-gene similarity, and key genes were ranked by "
                 "their membership in supported themes.")
    if ok("great_enrichment"):
        s.append(f"Region-based enrichment used rGREAT{v('rGREAT')} (basal plus extension rule).")
    mot = store.latest("find_motifs")
    if mot:
        if mot.summary.get("backend") == "HOMER":
            s.append(f"Known-motif enrichment used HOMER findMotifsGenome.pl on {mot.summary.get('window_bp')} bp "
                     "windows centred on peaks.")
        else:
            s.append(f"Known-motif enrichment used monaLisa{v('monaLisa')} with {mot.summary.get('database')} CORE "
                     f"vertebrate PWMs on {mot.summary.get('window_bp')} bp windows centred on peaks, against a "
                     "GC-matched random genomic background.")
    if ok("chipatlas_collect"):
        s.append("Overlap with public ChIP-seq/ATAC-seq experiments was tested with ChIP-Atlas Enrichment "
                 "Analysis (random-permutation background, Fisher's exact test, BH-adjusted).")
    if ok("chipatlas_target_genes"):
        s.append("Known TF target genes were taken from ChIP-Atlas target-gene tables; overlap was assessed "
                 "with a hypergeometric test.")
    if ok("integrate_rnaseq"):
        s.append("Differentially expressed genes were intersected with peak-linked genes; enrichment was "
                 "assessed with a hypergeometric test over all tested genes.")
    if ok("search_literature") or ok("search_geo"):
        s.append("Related literature was retrieved from Europe PMC and related datasets from NCBI GEO.")
    if config.use_llm:
        s.append(f"Analysis steps were selected and results interpreted by an LLM agent ({config.model}); all "
                 "statistics come from the tools above and every claim in the narrative is linked to its "
                 "source result.")
    return " ".join(s), versions


# ------------------------------------------------------------------ auto narrative (no LLM)

def auto_narrative(store: EvidenceStore) -> dict:
    facts = []
    qc = store.latest("peak_qc")
    if qc:
        facts.append(f"- {qc.summary.get('n_peaks')} peaks; {qc.summary.get('n_clean_peaks')} after QC "
                     f"(median width {qc.summary.get('width', {}).get('median')} bp) [{qc.id}].")
    ann = store.latest("annotate_peaks")
    if ann:
        rc = ann.summary.get("region_class_counts", {})
        total = sum(rc.values()) or 1
        facts.append(f"- {100 * rc.get('promoter', 0) / total:.0f}% of peaks are promoter-proximal and "
                     f"{100 * (rc.get('distal', 0) + rc.get('distal_far', 0)) / total:.0f}% distal [{ann.id}].")
    for ev in store.items.values():
        if ev.tool == "enrich_pathways" and ev.status == "ok":
            for t in (ev.summary.get("top_terms") or [])[:5]:
                facts.append(f"- {t['method']} ({t['group']}): {t['term']} (padj {t['padj']}) [{ev.id}].")
            for method, mods in (ev.summary.get("pathway_modules") or {}).items():
                if isinstance(mods, list) and len(mods) > 1:
                    names = "; ".join(f"M{m['module']} {m['lead_term']} ({m['n_terms']})" for m in mods[:6])
                    facts.append(f"- {method} pathway modules: {names} [{ev.id}].")
            for method, hubs in (ev.summary.get("hub_genes") or {}).items():
                if isinstance(hubs, list) and hubs:
                    facts.append(f"- {method} hub genes shared across top pathways: "
                                 f"{', '.join(h['gene'] for h in hubs[:10])} [{ev.id}].")
    mot = store.latest("find_motifs")
    if mot and mot.summary.get("top_motifs"):
        names = ", ".join(m["tf"] for m in mot.summary["top_motifs"][:8])
        facts.append(f"- Most enriched known motifs: {names} [{mot.id}].")
    ca = store.latest("chipatlas_collect")
    if ca and ca.summary.get("top_antigens"):
        names = ", ".join(a["antigen"] for a in ca.summary["top_antigens"][:8])
        facts.append(f"- Public ChIP-seq most overlapping these peaks: {names} [{ca.id}].")
    exp = store.latest("check_expectations")
    if exp:
        facts.append(f"- Positive controls: expected genes bound {exp.summary.get('genes_recovered')}, "
                     f"expected TFs recovered {exp.summary.get('tfs_recovered')} [{exp.id}].")
    return {
        "title": "Peak annotation report",
        "executive_summary": "Automated summary (no LLM interpretation was run). Key results:\n\n" + "\n".join(facts),
        "sections": [], "hypotheses": [],
        "caveats": ["This report was generated without the LLM agent: steps followed a fixed plan and no "
                    "biological interpretation is provided."],
    }


# ------------------------------------------------------------------ page

CSS = """
:root{--bg:#ffffff;--fg:#2c3e50;--muted:#7b8a8b;--card:#fff;--line:#dce4ec;--grey:#ecf0f1;--brand:#2c3e50;--accent:#18bc9c;
--accent2:#3498db;--ok:#18bc9c;--err:#e74c3c;--warn-bg:#fef5e7;--warn:#b9770e;--code:#ecf0f1;--bad:#e74c3c}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#1b2631;--fg:#ecf0f1;--muted:#95a5a6;
--card:#212f3d;--line:#2e4053;--grey:#2e4053;--brand:#2c3e50;--accent:#1dd1a1;--accent2:#5dade2;--ok:#1dd1a1;--err:#ec7063;
--warn-bg:#3a3020;--warn:#f5b041;--code:#2e4053;--bad:#ec7063}}
:root[data-theme="dark"]{--bg:#1b2631;--fg:#ecf0f1;--muted:#95a5a6;--card:#212f3d;--line:#2e4053;--grey:#2e4053;--brand:#2c3e50;
--accent:#1dd1a1;--accent2:#5dade2;--ok:#1dd1a1;--err:#ec7063;--warn-bg:#3a3020;--warn:#f5b041;--code:#2e4053;--bad:#ec7063}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 "Lato","Helvetica Neue",
Helvetica,Arial,sans-serif}
.topbar{background:var(--brand);color:#fff;min-height:60px;display:flex;align-items:center;gap:14px;padding:0 22px;font-size:.95em}
.topbar .brand{font-size:1.2rem;font-weight:700;color:#fff;text-decoration:none}.topbar .brand span{color:var(--accent)}
.topbar .right{margin-left:auto;opacity:.85}main{max-width:1080px;margin:0 auto;padding:24px 16px 80px}
h1,h2{color:var(--fg);font-weight:700}
h1{font-size:1.9rem;line-height:1.25;margin:.2em 0}
h2{font-size:1.3rem;margin:2.2em 0 .6em;padding:8px 12px;background:var(--grey);border-left:4px solid var(--accent);border-radius:4px}
h3{font-size:1.02rem;margin:0}h4{margin:1em 0 .3em;font-size:.92rem;color:var(--muted);text-transform:capitalize}
a{color:var(--accent)}a:hover{color:var(--brand)}.muted{color:var(--muted);font-size:.88em}nav{display:flex;flex-wrap:wrap;gap:6px 8px;
margin:12px 0 0;font-size:.9em}nav a{background:var(--grey);color:var(--fg);padding:4px 10px;border-radius:4px;text-decoration:none}nav a:hover{background:var(--brand);color:#fff}.meta{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:8px;
margin:14px 0}.meta div{background:var(--card);border:1px solid var(--line);border-top:3px solid var(--accent);border-radius:4px;padding:8px 10px}
.meta b{display:block;font-size:.75em;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
.card{background:var(--card);border:1px solid var(--line);border-radius:4px;padding:14px 16px;margin:12px 0;box-shadow:0 1px 1px rgba(0,0,0,.05)}
.card header{display:flex;align-items:center;gap:10px;flex-wrap:wrap}.eid{font:600 .8em ui-monospace,monospace;color:#fff;
background:var(--brand);padding:2px 6px;border-radius:4px}.status{margin-left:auto;font-size:.78em;font-weight:600}
.status.ok{color:var(--ok)}.status.err{color:var(--err)}.status.skip{color:var(--muted)}
.why{color:var(--muted);font-style:italic;margin:.5em 0}pre.error{white-space:pre-wrap;color:var(--err);
background:var(--code);padding:8px;border-radius:6px;font-size:.85em}
.scroll{overflow-x:auto}table{border-collapse:collapse;font-size:.85em;margin:.4em 0}
th,td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;vertical-align:top}
th{color:var(--brand);border-bottom:2px solid var(--line)}table.data th{position:sticky;top:0;background:var(--card)}table.kv th{color:var(--muted);font-weight:500;
white-space:nowrap;text-transform:capitalize}figure{margin:10px 0}figcaption{font-size:.8em}figure img{max-width:100%;height:auto;
border-radius:6px;background:#fff}ul.warn{background:var(--warn-bg);color:var(--warn);border-radius:6px;
padding:8px 8px 8px 28px;font-size:.9em}details{margin:.5em 0}summary{cursor:pointer;color:var(--muted)}
.summary{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:4px;
padding:4px 18px}a.cite{text-decoration:none;font-size:.85em;font-family:ui-monospace,monospace;background:var(--grey);padding:0 4px;border-radius:3px}
a.cite.unverified{color:var(--bad);text-decoration:line-through wavy}.tag{font-size:.72em;border:1px solid
var(--line);border-radius:4px;padding:0 4px;color:var(--muted)}ol.refs li,ul.refs li{margin:.35em 0}
.conf{font-size:.75em;font-weight:700;text-transform:uppercase}.conf.high{color:var(--ok)}.conf.medium{color:var(--warn)}
.conf.low{color:var(--muted)}#igv{border:1px solid var(--line);border-radius:8px;min-height:120px;background:#fff}
select{font:inherit;padding:4px;margin:6px 0;max-width:100%}code{background:var(--code);padding:0 4px;border-radius:4px}
.flag{background:var(--warn-bg);color:var(--warn);padding:8px 12px;border-radius:6px;font-size:.9em}
.search{display:flex;gap:8px;align-items:center;margin:14px 0 0}.search input{flex:1;font:inherit;padding:6px 10px;
border:1px solid var(--line);border-radius:6px;background:var(--card);color:var(--fg);max-width:480px}
.search .n{font-size:.85em;color:var(--muted)}mark{background:#ffe97a;color:#111;padding:0 1px;border-radius:2px}
.card.dim{opacity:.35}tr.dim{display:none}
figcaption{display:flex;justify-content:space-between;gap:10px;color:var(--muted);margin-top:4px}figcaption .dl a{margin-left:4px}
.explore-box{margin:14px 0;border:1px solid var(--line);border-radius:4px;background:var(--card)}
.explore-box>h4{margin:0;padding:8px 12px;background:var(--grey);color:var(--brand);text-transform:none;font-size:.9em}
.explore{padding:10px 12px}.ex-head{display:flex;flex-wrap:wrap;align-items:center;gap:8px 16px;margin-bottom:8px}
.ex-tabs{display:flex;gap:4px}.ex-tab{font:inherit;font-size:.88em;padding:6px 12px;border:1px solid var(--line);background:var(--card);color:var(--fg);border-radius:4px;cursor:pointer}
.ex-tab.on{background:var(--brand);color:#fff;border-color:var(--brand)}
.ex-controls{display:flex;flex-wrap:wrap;gap:10px;font-size:.86em;color:var(--muted)}.ex-controls label{display:inline-flex;align-items:center;gap:6px}
.ex-controls select{margin:0;padding:3px 6px}.ex-dl{margin-left:auto;display:flex;align-items:center;gap:8px}
.ex-link{font:inherit;font-size:.86em;background:none;border:none;color:var(--accent);cursor:pointer;padding:0;text-decoration:underline}
.ex-body{display:grid;grid-template-columns:minmax(0,1fr) 280px;gap:14px}@media (max-width:900px){.ex-body{grid-template-columns:1fr}}
.ex-stage{min-width:0;overflow-x:auto}.ex-svg{max-width:100%;height:auto;display:block;background:#fff;border:1px solid var(--line);border-radius:4px}
.ex-svg .node.sel circle,.ex-svg .bar.sel rect{stroke:var(--brand);stroke-width:3}.ex-svg .bar.sel text{font-weight:700}
.ex-panel{border:1px solid var(--line);border-radius:4px;padding:10px 12px;background:var(--bg);font-size:.9em;align-self:start;max-height:640px;overflow:auto}
.ex-panel h4{margin:0 0 4px;text-transform:none;font-size:1em;color:var(--brand)}.ex-panel h5{margin:12px 0 4px;font-size:.82em;color:var(--muted);text-transform:uppercase}
.ex-genes{display:flex;flex-wrap:wrap;gap:4px;margin:8px 0}
.ex-gene{font:inherit;font-size:.82em;padding:2px 8px;border-radius:4px;border:1px solid var(--line);background:var(--card);color:var(--fg);cursor:pointer}.ex-gene:hover{border-color:var(--accent)}
.ex-actions{display:flex;gap:12px;margin:4px 0 8px}.ex-geneinfo{background:var(--grey);border-radius:4px;padding:8px 10px;margin:6px 0;font-size:.92em}
.ex-geneinfo .ex-term{display:inline-block;margin:4px 6px 0 0;font-size:.9em}.ex-list{margin:4px 0;padding-left:18px}
.ex-legend{display:flex;flex-wrap:wrap;gap:4px 14px;margin:8px 0 0;font-size:.82em;color:var(--muted)}.ex-legend i{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:4px;vertical-align:-1px}
"""

REPORT_EXPLORE_JS = """
(function(){
  document.querySelectorAll('script[id^="explore-data-"]').forEach(function(sc){
    var eid = sc.id.replace('explore-data-', ''), box = document.getElementById('explore-' + eid);
    if (!box || !window.PeakExplore) return;
    try {
      var data = JSON.parse(sc.textContent);
      box.innerHTML = '';
      PeakExplore.mount(box, data, {title: document.title + '_' + eid, onGene: function(g){
        var s = document.getElementById('rsearch'); if (!s) return;
        s.value = g; s.dispatchEvent(new Event('input', {bubbles: true})); s.scrollIntoView({behavior: 'smooth', block: 'start'});
      }});
    } catch (e) { box.textContent = 'Interactive view unavailable: ' + e; }
  });
})();
"""

REPORT_SEARCH_JS = """
(function(){
  const box = document.getElementById('rsearch'), out = document.getElementById('rsearch-n');
  const cards = [...document.querySelectorAll('article.card')];
  let marks = [];
  function clearMarks(){ marks.forEach(m => m.replaceWith(document.createTextNode(m.textContent))); marks = [];
    document.querySelectorAll('.card.dim, tr.dim').forEach(e => e.classList.remove('dim')); }
  function markIn(el, q){
    const walker = document.createTreeWalker(el, NodeFilter.SHOW_TEXT, null);
    const nodes = []; let n; while ((n = walker.nextNode())) if (n.nodeValue.toLowerCase().includes(q)) nodes.push(n);
    nodes.forEach(node => { const i = node.nodeValue.toLowerCase().indexOf(q); const m = document.createElement('mark');
      const after = node.splitText(i); after.nodeValue = after.nodeValue.slice(q.length);
      m.textContent = node.nodeValue.slice(i); node.nodeValue = node.nodeValue.slice(0, i);
      node.parentNode.insertBefore(m, after); marks.push(m); });
    return nodes.length;
  }
  function run(){
    clearMarks(); const q = box.value.trim().toLowerCase();
    if (!q) { out.textContent = ''; return; }
    let hits = 0, first = null;
    cards.forEach(card => {
      card.querySelectorAll('details').forEach(d => {
        d.querySelectorAll('tbody tr').forEach(tr => { if (!tr.textContent.toLowerCase().includes(q)) tr.classList.add('dim'); });
        if ([...d.querySelectorAll('tbody tr')].some(tr => !tr.classList.contains('dim')) && d.textContent.toLowerCase().includes(q)) d.open = true;
      });
      const n = markIn(card, q); hits += n;
      if (!n) card.classList.add('dim'); else if (!first) first = card;
    });
    out.textContent = hits ? hits + ' match(es) in ' + cards.filter(c => !c.classList.contains('dim')).length + ' result(s)' : 'no match in results';
    if (first) first.scrollIntoView({behavior: 'smooth', block: 'start'});
  }
  let t; box.addEventListener('input', () => { clearTimeout(t); t = setTimeout(run, 250); });
  box.addEventListener('keydown', e => { if (e.key === 'Escape') { box.value = ''; run(); } });
})();
"""


def _igv_section(config, store: EvidenceStore, ctx_state: dict) -> str:
    try:
        peaks = []
        with open(ctx_state.get("clean_peaks") or config.peaks) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 3 and f[1].isdigit():
                    peaks.append({"chr": f[0], "start": int(f[1]), "end": int(f[2]),
                                  "name": f[3] if len(f) > 3 else ""})
                if len(peaks) >= 5000:
                    break
    except OSError:
        return ""
    loci = list(dict.fromkeys(config.biology.expected_genes +
                              list((store.latest("annotate_peaks").summary.get("example_promoter_genes") or [])
                                   if store.latest("annotate_peaks") else [])))[:60]
    if not loci and peaks:
        loci = [f"{p['chr']}:{p['start']}-{p['end']}" for p in peaks[:30]]
    options = "".join(f"<option>{esc(l)}</option>" for l in loci)
    data = json.dumps(peaks).replace("</", "<\\/")
    return f"""
<h2 id="browser">Genome browser</h2>
<p class="muted">Peaks from this run on {esc(config.genome)} (igv.js; needs internet). Jump to a gene:</p>
<select id="locus">{options}</select>
<div id="igv"></div>
<script src="{IGV_JS}"></script>
<script>
(function(){{
  const peaks = {data};
  const sel = document.getElementById('locus');
  if (!window.igv) {{ document.getElementById('igv').textContent = 'igv.js could not be loaded (offline?).'; return; }}
  igv.createBrowser(document.getElementById('igv'), {{
    genome: {json.dumps(config.genome)}, locus: sel.value || 'all',
    tracks: [{{name: 'Peaks', type: 'annotation', format: 'bed', features: peaks, displayMode: 'EXPANDED',
               color: '#2c3e50'}}]
  }}).then(b => sel.addEventListener('change', () => b.search(sel.value)));
}})();
</script>"""


def build_report(store: EvidenceStore, config, narrative: Optional[dict], out_path: str,
                 state: Optional[dict] = None) -> dict:
    """Write the HTML report. Returns {'path', 'unverified_citations'}."""
    out = Path(out_path)
    report_dir = out.parent
    narrative = narrative or auto_narrative(store)
    all_text = " ".join([narrative.get("executive_summary", "")] +
                        [s.get("body", "") for s in narrative.get("sections", [])] +
                        [h.get("statement", "") + " " + " ".join(f"[{e}]" for e in h.get("evidence", []))
                         for h in narrative.get("hypotheses", [])] + list(narrative.get("caveats", [])))
    unverified = store.check_citations(all_text)
    bio = config.biology
    visible = [e for e in store.items.values() if e.tool not in HIDDEN_TOOLS]
    n_peaks = (store.latest("peak_qc").summary.get("n_peaks") if store.latest("peak_qc") else "")

    p = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' "
         f"content='width=device-width,initial-scale=1'><title>{esc(narrative.get('title') or 'Peak report')}</title>"
         f"<link rel='preconnect' href='https://fonts.googleapis.com'><link rel='preconnect' href='https://fonts.gstatic.com' crossorigin>"
         f"<link href='https://fonts.googleapis.com/css2?family=Lato:wght@400;700&display=swap' rel='stylesheet'>"
         f"<style>{CSS}</style></head><body>"
         f"<div class='topbar'><a class='brand' href='https://github.com/helenhuangmath/PeakAnnoAgent'><span>◆</span> PeakAnnoAgent</a>"
         f"<span class='right'>Report · {datetime.now():%Y-%m-%d %H:%M} · v{__version__}</span></div><main>"]
    p.append(f"<h1>{esc(narrative.get('title') or 'Peak annotation report')}</h1>")
    from .outcome import assess_run
    completion = (state or {}).get("completion") or assess_run(store, config)
    p.append(f"<p class='{'flag' if completion['status'] != 'success' else 'muted'}' role='status'>"
             f"Analysis status: <strong>{esc(completion['status'])}</strong>"
             + (": " + esc("; ".join(completion["reasons"])) if completion["reasons"] else "") + "</p>")
    meta = [("Peaks", f"{Path(config.peaks).name} ({n_peaks})"), ("Genome", config.genome),
            ("Assay", bio.assay or "-"), ("Target", bio.target or "-"), ("Cell type", bio.cell_type or "-"),
            ("Condition", bio.condition or "-"),
            ("Mode", f"agent ({config.model})" if config.use_llm else "deterministic (no LLM)")]
    p.append("<div class='meta'>" + "".join(f"<div><b>{esc(k)}</b>{esc(v)}</div>" for k, v in meta) + "</div>")
    nav = ["summary", "hypotheses", "interpretation", "results", "browser", "references", "methods", "log"]
    p.append("<nav>" + "".join(f"<a href='#{n}'>{n.capitalize()}</a>" for n in nav) + "</nav>")
    p.append("<div class='search'><input id='rsearch' type='search' placeholder='Search results: gene, peak id, "
             "pathway, motif, TF, PMID…' aria-label='Search this report'><span class='n' id='rsearch-n'></span></div>")
    if bio.description or bio.expected_genes or bio.expected_tfs:
        p.append("<details><summary>Known biology provided</summary>")
        if bio.description:
            p.append(f"<p>{esc(bio.description)}</p>")
        if bio.expected_genes:
            p.append(f"<p><b>Expected genes:</b> {esc(', '.join(bio.expected_genes))}</p>")
        if bio.expected_tfs:
            p.append(f"<p><b>Expected TFs:</b> {esc(', '.join(bio.expected_tfs))}</p>")
        p.append("</details>")

    p.append("<h2 id='summary'>Summary</h2>")
    if unverified:
        p.append(f"<p class='flag'>{len(unverified)} citation(s) in the narrative do not match any result from "
                 f"this run and are struck through: {esc(', '.join(unverified))}.</p>")
    p.append(f"<div class='summary'>{md(narrative.get('executive_summary', ''), store)}</div>")

    hyps = narrative.get("hypotheses") or []
    p.append("<h2 id='hypotheses'>Hypotheses and next experiments</h2>")
    if hyps:
        for i, h in enumerate(hyps, 1):
            ev = ", ".join(_cite_html(e, store) for e in h.get("evidence", []))
            p.append(f"<div class='card'><header><h3>H{i}. {md_inline(h.get('statement', ''), store)}</h3>"
                     f"<span class='conf {esc(h.get('confidence'))}'>{esc(h.get('confidence'))}</span></header>"
                     f"<p class='muted'>Evidence: {ev}</p><p><b>Test:</b> {esc(h.get('next_experiment'))}</p></div>")
    else:
        p.append("<p class='muted'>None generated.</p>")

    p.append("<h2 id='interpretation'>Interpretation</h2>")
    for s in narrative.get("sections") or []:
        p.append(f"<h3>{esc(s.get('heading'))}</h3>{md(s.get('body', ''), store)}")
    if narrative.get("caveats"):
        p.append("<h3>Caveats</h3><ul>" + "".join(f"<li>{md_inline(c, store)}</li>" for c in narrative["caveats"])
                 + "</ul>")

    p.append("<h2 id='results'>Results</h2><p class='muted'>Each result has an id (E1, E2, …) that the "
             "narrative cites. Tables are linked relative to this report.</p>")
    for heading, tools in TOOL_SECTIONS:
        cards = [e for e in visible if e.tool in tools]
        if cards:
            p.append(f"<h3 style='margin-top:1.4em'>{esc(heading)}</h3>")
            p.extend(evidence_card(e, report_dir) for e in cards)

    p.append(_igv_section(config, store, state or {}))

    p.append("<h2 id='references'>References and datasets</h2>")
    p.append(_paper_list(list(store.references.values())) if store.references
             else "<p class='muted'>No literature searched.</p>")
    geo = [d for d in store.datasets.values() if d.get("accession")]
    if geo:
        p.append("<h4>GEO series</h4>" + _dataset_list(geo))
    exps = [(k, d) for k, d in store.datasets.items() if not d.get("accession")]
    if exps:
        p.append("<h4>Public experiments (ChIP-Atlas)</h4><ul class='refs'>" + "".join(
            f"<li><a href='{esc(d.get('url'))}' target='_blank' rel='noopener'>{esc(k)}</a> - "
            f"{esc(d.get('antigen'))}, {esc(d.get('cell'))}</li>" for k, d in exps) + "</ul>")

    text, versions = methods_text(store, config)
    p.append(f"<h2 id='methods'>Methods</h2><p>{esc(text)}</p>")
    if versions:
        p.append("<details><summary>Software versions</summary>" +
                 _table([{"package": k, "version": v} for k, v in sorted(versions.items())], 100) + "</details>")

    p.append("<h2 id='log'>Analysis log</h2>" + _table(
        [{"step": i + 1, "time": t["time"], "id": t["evidence"], "tool": t["tool"], "status": t["status"],
          "rationale": t["rationale"]} for i, t in enumerate(store.trace)], 500))
    p.append(f"<script>{REPORT_SEARCH_JS}</script>")
    if any(has_explore(e) for e in visible):
        p.append(f"<script>{explore_js()}</script><script>{REPORT_EXPLORE_JS}</script>")
    p.append("</main></body></html>")
    out.write_text("\n".join(p))
    return {"path": str(out), "unverified_citations": unverified}
