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

from .analysis import read_table
from .evidence import CITATION_RE, Evidence, EvidenceStore

IGV_JS = "https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"
HIDDEN_TOOLS = ("inspect_table", "chipatlas_status")
TOOL_SECTIONS = [
    ("Peak QC", ("peak_qc",)),
    ("Annotation", ("annotate_peaks",)),
    ("Pathways", ("enrich_pathways", "great_enrichment")),
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
    caption = ""
    pdf = Path(path).with_suffix(".pdf")
    if Path(path).parent.name == ".preview":          # previews sit one level below their PDFs
        pdf = Path(path).parent.parent / pdf.name
    if report_dir is not None and pdf.is_file():
        caption = f'<figcaption><a href="{esc(_rel(str(pdf), report_dir))}">PDF (vector)</a></figcaption>'
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
    if summary:
        parts.append(_summary_html(summary))
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
:root{--bg:#ffffff;--fg:#111111;--muted:#5c6370;--card:#fff;--line:#e2e5ea;--brand:#000000;--accent:#000000;
--accent2:#000000;--ok:#1a1a1a;--err:#000000;--warn-bg:#f3f3f3;--warn:#333333;--code:#f1f3f7;--bad:#000000}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#15181d;--fg:#e6e8eb;--muted:#9aa3ad;
--card:#1d2127;--line:#2e343c;--brand:#f2f2f2;--accent:#f2f2f2;--accent2:#f2f2f2;--ok:#e6e8eb;--err:#f2f2f2;
--warn-bg:#262626;--warn:#d9d9d9;--code:#262b32;--bad:#f2f2f2}}
:root[data-theme="dark"]{--bg:#15181d;--fg:#e6e8eb;--muted:#9aa3ad;--card:#1d2127;--line:#2e343c;--brand:#f2f2f2;
--accent:#f2f2f2;--accent2:#f2f2f2;--ok:#e6e8eb;--err:#f2f2f2;--warn-bg:#262626;--warn:#d9d9d9;--code:#262b32;--bad:#f2f2f2}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.55 -apple-system,
BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
body{border-top:5px solid var(--brand)}main{max-width:1080px;margin:0 auto;padding:24px 16px 80px}
h1,h2{font-family:Georgia,"Times New Roman",serif;color:var(--brand);font-weight:700}
h1{font-size:1.8rem;line-height:1.25;margin:.2em 0}
h2{font-size:1.3rem;margin:2.2em 0 .6em;padding-bottom:.3em;border-bottom:2px solid var(--brand)}
h3{font-size:1.02rem;margin:0}h4{margin:1em 0 .3em;font-size:.92rem;color:var(--muted);text-transform:capitalize}
a{color:var(--accent)}.muted{color:var(--muted);font-size:.88em}nav{display:flex;flex-wrap:wrap;gap:6px 14px;
margin:12px 0 0;font-size:.9em}.meta{display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:8px;
margin:14px 0}.meta div{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px 10px}
.meta b{display:block;font-size:.75em;color:var(--muted);text-transform:uppercase;letter-spacing:.04em}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px 16px;margin:12px 0}
.card header{display:flex;align-items:center;gap:10px;flex-wrap:wrap}.eid{font:600 .8em ui-monospace,monospace;color:#fff;
background:var(--brand);padding:2px 6px;border-radius:5px}.status{margin-left:auto;font-size:.78em;font-weight:600}
.status.ok{color:var(--ok)}.status.err{color:var(--err)}.status.skip{color:var(--muted)}
.why{color:var(--muted);font-style:italic;margin:.5em 0}pre.error{white-space:pre-wrap;color:var(--err);
background:var(--code);padding:8px;border-radius:6px;font-size:.85em}
.scroll{overflow-x:auto}table{border-collapse:collapse;font-size:.85em;margin:.4em 0}
th,td{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;vertical-align:top}
table.data th{position:sticky;top:0;background:var(--card)}table.kv th{color:var(--muted);font-weight:500;
white-space:nowrap;text-transform:capitalize}figure{margin:10px 0}figcaption{font-size:.8em}figure img{max-width:100%;height:auto;
border-radius:6px;background:#fff}ul.warn{background:var(--warn-bg);color:var(--warn);border-radius:6px;
padding:8px 8px 8px 28px;font-size:.9em}details{margin:.5em 0}summary{cursor:pointer;color:var(--muted)}
.summary{background:var(--card);border:1px solid var(--line);border-left:4px solid var(--accent);border-radius:8px;
padding:4px 18px}a.cite{text-decoration:none;font-size:.85em;font-family:ui-monospace,monospace}
a.cite.unverified{color:var(--bad);text-decoration:line-through wavy}.tag{font-size:.72em;border:1px solid
var(--line);border-radius:4px;padding:0 4px;color:var(--muted)}ol.refs li,ul.refs li{margin:.35em 0}
.conf{font-size:.75em;font-weight:700;text-transform:uppercase}.conf.high{color:var(--ok)}.conf.medium{color:var(--warn)}
.conf.low{color:var(--muted)}#igv{border:1px solid var(--line);border-radius:8px;min-height:120px;background:#fff}
select{font:inherit;padding:4px;margin:6px 0;max-width:100%}code{background:var(--code);padding:0 4px;border-radius:4px}
.flag{background:var(--warn-bg);color:var(--warn);padding:8px 12px;border-radius:6px;font-size:.9em}
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
               color: '#000000'}}]
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
         f"<style>{CSS}</style></head><body><main>"]
    p.append(f"<p class='muted'>PeakAnnoAgent report · {datetime.now():%Y-%m-%d %H:%M}</p>"
             f"<h1>{esc(narrative.get('title') or 'Peak annotation report')}</h1>")
    meta = [("Peaks", f"{Path(config.peaks).name} ({n_peaks})"), ("Genome", config.genome),
            ("Assay", bio.assay or "-"), ("Target", bio.target or "-"), ("Cell type", bio.cell_type or "-"),
            ("Condition", bio.condition or "-"),
            ("Mode", f"agent ({config.model})" if config.use_llm else "deterministic (no LLM)")]
    p.append("<div class='meta'>" + "".join(f"<div><b>{esc(k)}</b>{esc(v)}</div>" for k, v in meta) + "</div>")
    nav = ["summary", "hypotheses", "interpretation", "results", "browser", "references", "methods", "log"]
    p.append("<nav>" + "".join(f"<a href='#{n}'>{n.capitalize()}</a>" for n in nav) + "</nav>")
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
    p.append("</main></body></html>")
    out.write_text("\n".join(p))
    return {"path": str(out), "unverified_citations": unverified}
