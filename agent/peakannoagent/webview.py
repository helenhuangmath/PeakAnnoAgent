"""Server-side rendering of a run for the web app: the run is split into parts
(overview, summary, hypotheses, annotation, pathways, ...) that the page shows one at a
time. Reuses the report's markdown, table and list renderers; figures are served as
files instead of being embedded."""

from __future__ import annotations

import html
import os
from pathlib import Path
from typing import Callable, Dict, List, Optional

from .explore import explore_data, has_explore
from .evidence import Evidence
from .report import _dataset_list, _fmt, _paper_list, md, md_inline, methods_text
from .search import Run, _f

def esc(x) -> str:
    return html.escape("" if x is None else str(x))

PARTS = [
    # id, label, tools shown in it (None = narrative part)
    ("overview", "Overview", None),
    ("summary", "Summary", None),
    ("interpretation", "Interpretation", None),
    ("annotation", "Peaks & annotation", ("peak_qc", "annotate_peaks")),
    ("pathways", "Pathways", ("enrich_pathways", "network_viz", "pathway_consensus", "great_enrichment")),
    ("motifs", "Motifs", ("find_motifs",)),
    ("regulators", "Regulators", ("chipatlas_submit", "chipatlas_collect", "chipatlas_target_genes")),
    ("suggestions", "Suggestions", None),
    ("expression", "Expression & controls", ("integrate_rnaseq", "check_expectations")),
    ("literature", "Literature & data", ("search_literature", "search_geo")),
    ("methods", "Methods", None),
    ("log", "Analysis log", None),
]
HIDDEN = ("inspect_table", "chipatlas_status")
TOOL_LABEL = {"peak_qc": "Peak QC", "annotate_peaks": "Annotation", "enrich_pathways": "Pathway enrichment",
              "great_enrichment": "GREAT", "network_viz": "Pathway network",
              "pathway_consensus": "Pathway consensus", "find_motifs": "Motif enrichment", "chipatlas_submit": "ChIP-Atlas job",
              "chipatlas_collect": "ChIP-Atlas co-binding", "chipatlas_target_genes": "Known targets",
              "integrate_rnaseq": "RNA-seq integration", "check_expectations": "Positive controls",
              "search_literature": "Literature", "search_geo": "GEO datasets"}


def _table(rows: List[dict], show: int = 12, cap: int = 500) -> str:
    """Interactive table: all rows (up to `cap`) are rendered; rows beyond `show` are hidden until
    'show all'. Column headers sort (click), the search box filters. The page script does the rest."""
    if not rows:
        return ""
    cols = list(rows[0].keys())
    head = "".join(f'<th data-col="{i}" title="click to sort">{esc(c)}</th>' for i, c in enumerate(cols))
    body = []
    for i, r in enumerate(rows[:cap]):
        cls = ' class="xtra"' if i >= show else ""
        body.append(f"<tr{cls}>" + "".join(f"<td>{_fmt(r.get(c))}</td>" for c in cols) + "</tr>")
    n = len(rows)
    bar = ""
    if n > 5:
        more = (f'<button type="button" class="btn ghost small showall">show all {min(n, cap)}</button>' if n > show else "")
        capnote = f'<span class="muted small">first {cap} of {n} rows; download the file for all</span>' if n > cap else ""
        bar = (f'<div class="tbar"><input type="search" class="tfilter" placeholder="filter rows…" aria-label="filter rows">'
               f'<span class="tcount muted small">{n} rows</span>{more}{capnote}</div>')
    return (f'<div class="itable" data-show="{show}">{bar}<div class="scroll"><table class="data">'
            f'<thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div></div>')


def _summary_html(summary: dict) -> str:
    """Like report._summary_html, with interactive tables."""
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
            continue
        elif value not in (None, "", [], {}):
            scalars.append((label, _fmt(value)))
    kv = "".join(f"<tr><th>{esc(k)}</th><td>{v}</td></tr>" for k, v in scalars)
    return (f'<table class="kv">{kv}</table>' if kv else "") + "".join(blocks)


class View:
    def __init__(self, run: Run, file_url: Callable[[str], str]):
        self.run = run
        self.store = run.store
        self.file_url = file_url            # relative path inside the run dir -> URL

    # -- helpers -----------------------------------------------------------
    def _url(self, path: str) -> Optional[str]:
        p = self.run.resolve(path)
        if not p:
            return None
        try:
            return self.file_url(str(p.relative_to(self.run.dir)))
        except ValueError:
            return None

    def _figure(self, name: str, path: str) -> str:
        """Inline figure (PNG preview, or the PDF itself when no preview exists) with download links."""
        png = self._url(path)
        pdf_name = Path(path).with_suffix(".pdf").name
        pdf = self._url(str(Path(path).parent.parent / pdf_name)) or self._url(str(Path(path).parent / pdf_name))
        label = esc(name.replace("_", " "))
        links = " · ".join(f'<a href="{u}" download>{t}</a>' for u, t in ((png, "PNG"), (pdf, "PDF")) if u)
        cap = f'<figcaption><span>{label}</span><span class="dl">{links}</span></figcaption>'
        if png:
            return (f'<figure><a href="{png}" target="_blank" rel="noopener" title="open full size">'
                    f'<img loading="lazy" alt="{label}" src="{png}"></a>{cap}</figure>')
        if pdf:
            return (f'<figure class="pdf"><object data="{pdf}#toolbar=0&navpanes=0" type="application/pdf">'
                    f'<a href="{pdf}" target="_blank" rel="noopener">{label} (PDF)</a></object>{cap}</figure>')
        return ""

    def _tables(self, ev: Evidence) -> str:
        out = []
        for name, path in ev.tables.items():
            rows = self.run.table(ev.id, name)
            url = self._url(path)
            link = f' · <a href="{url}" download>download {esc(Path(path).name)}</a>' if url else ""
            out.append(f'<details class="tbl"><summary>{esc(name.replace("_", " "))} '
                       f'<span class="muted">({len(rows)} rows){link}</span></summary>{_table(rows, 15)}</details>')
        return "".join(out)

    def card(self, ev: Evidence) -> str:
        status = {"ok": "ok", "error": "err", "skipped": "skip"}.get(ev.status, "skip")
        h = [f'<article class="card ev" id="{ev.id}"><header><span class="eid">{ev.id}</span>'
             f'<h3>{esc(ev.title)}</h3><span class="pill {status}">{esc(ev.status)}</span></header>']
        if ev.rationale:
            h.append(f'<p class="why">{esc(ev.rationale)}</p>')
        if ev.error:
            h.append(f'<pre class="error">{esc(ev.error)}</pre>')
        summary = dict(ev.summary)
        papers, datasets = summary.pop("papers", None), summary.pop("datasets", None)
        for k in ("note", "clean_peaks", "gene_sets", "universe", "request_id", "cytoscape"):
            summary.pop(k, None)
        interactive = summary.pop("interactive", None) or {}
        if ev.tool == "search_literature" or ev.tool == "search_geo":
            summary.pop("query", None)
            h.append(f'<p class="muted">Query: <code>{esc(ev.params.get("query", ""))}</code></p>')
        if ev.tool == "enrich_pathways" and isinstance(summary.get("n_significant"), list):
            summary["n_significant"] = "; ".join(f"{x.get('method')} {x.get('group')}: {x.get('N')}"
                                                 for x in summary["n_significant"] if isinstance(x, dict))
            # lead with what matters: the top terms, then modules and hubs, then the counts
            order = ["top_terms", "pathway_modules", "hub_genes"]
            summary = {**{k: summary[k] for k in order if k in summary},
                       **{k: v for k, v in summary.items() if k not in order}}
        if has_explore(ev):
            h.append(f'<section class="explore-box" data-explore="{ev.id}"><h4>Explore the pathways</h4>'
                     f'<div class="explore" id="explore-{ev.id}"><p class="muted">Loading the interactive view…</p></div></section>')
        if summary:
            h.append(_summary_html(summary))
        if papers:
            h.append(_paper_list(papers))
        if datasets:
            h.append(_dataset_list(datasets))
        if interactive:
            links = [f'<a class="btn small" href="{self._url(v)}" target="_blank" rel="noopener">{esc(k.replace("_", " "))} ↗</a>'
                     for k, v in interactive.items() if self._url(v)]
            if links:
                h.append('<p class="links">Interactive: ' + " ".join(links) + "</p>")
        figs = "".join(self._figure(n, p) for n, p in ev.figures.items())
        if figs:
            h.append(f'<div class="figs">{figs}</div>')
        if ev.warnings:
            h.append('<ul class="warn">' + "".join(f"<li>{esc(w)}</li>" for w in ev.warnings) + "</ul>")
        h.append(self._tables(ev))
        h.append("</article>")
        return "".join(h)

    # -- parts -------------------------------------------------------------
    def stats(self) -> List[dict]:
        s, store = [], self.store
        qc, ann = store.latest("peak_qc"), store.latest("annotate_peaks")
        if qc:
            s.append({"label": "Peaks after QC", "value": qc.summary.get("n_clean_peaks"),
                      "sub": f"of {qc.summary.get('n_peaks')} · median {qc.summary.get('width', {}).get('median')} bp"})
        if ann:
            rc = ann.summary.get("region_class_counts", {})
            tot = sum(rc.values()) or 1
            s.append({"label": "Promoter peaks", "value": f"{100 * rc.get('promoter', 0) / tot:.0f}%",
                      "sub": f"{100 * (rc.get('distal', 0) + rc.get('distal_far', 0)) / tot:.0f}% distal"})
            s.append({"label": "Genes linked", "value": ann.summary.get("gene_set_sizes", {}).get("all"),
                      "sub": f"{ann.summary.get('gene_set_sizes', {}).get('promoter')} at promoters"})
        def count_sig(v) -> int:
            if isinstance(v, list):
                return sum(int(x.get("N") or 0) for x in v if isinstance(x, dict))
            try:
                return int(v or 0)
            except (TypeError, ValueError):
                return 0
        n_sig = sum(count_sig(e.summary.get("n_significant")) for e in store.items.values()
                    if e.tool == "enrich_pathways" and e.status == "ok")
        if n_sig:
            s.append({"label": "Enriched gene sets", "value": n_sig, "sub": "padj < 0.05, all collections"})
        mot = [e for e in store.items.values() if e.tool == "find_motifs" and e.status == "ok"]
        if mot:
            s.append({"label": "Significant motifs", "value": max(count_sig(e.summary.get("n_significant")) for e in mot),
                      "sub": f"{len(mot)} scan(s)"})
        pc = store.latest("pathway_consensus")
        if pc and pc.summary.get("lists"):
            lst = pc.summary["lists"].get("all") or next(iter(pc.summary["lists"].values()))
            tiers = lst.get("n_terms_by_tier") or {}
            s.append({"label": "Consensus themes", "value": len(lst.get("top_themes") or []),
                      "sub": f"{tiers.get('ROBUST', 0)} robust · {tiers.get('SUPPORTED', 0)} supported terms"})
        ca = store.latest("chipatlas_collect")
        if ca:
            s.append({"label": "Co-binding factors", "value": len(ca.summary.get("top_antigens") or []),
                      "sub": f"{ca.summary.get('n_experiments_tested')} public experiments"})
        exp = store.latest("check_expectations")
        if exp:
            s.append({"label": "Expected genes bound", "value": exp.summary.get("genes_recovered"),
                      "sub": f"TFs recovered {exp.summary.get('tfs_recovered')}"})
        if store.references:
            s.append({"label": "Papers", "value": len(store.references), "sub": f"{len(store.datasets)} datasets"})
        return s

    def overview(self) -> str:
        b = self.run.biology
        tiles = "".join(f'<div class="tile"><b>{esc(t["label"])}</b><span class="v">{esc(t["value"])}</span>'
                        f'<span class="s">{esc(t["sub"])}</span></div>' for t in self.stats())
        rows = [("Assay", b.get("assay")), ("Target", b.get("target")), ("Cell type", b.get("cell_type")),
                ("Condition", b.get("condition")), ("Expected genes", ", ".join(b.get("expected_genes") or [])),
                ("Expected TFs", ", ".join(b.get("expected_tfs") or []))]
        bio = "".join(f"<tr><th>{esc(k)}</th><td>{esc(v)}</td></tr>" for k, v in rows if v)
        notes = f'<p class="notes">{esc(b.get("description"))}</p>' if b.get("description") else ""
        steps = []
        for t in self.store.trace:
            ev = self.store.items.get(t["evidence"])
            if not ev or ev.tool in HIDDEN:
                continue
            steps.append(f'<li class="{esc(ev.status)}"><a href="#{ev.id}" class="eid">{ev.id}</a>'
                         f'<b>{esc(TOOL_LABEL.get(ev.tool, ev.tool))}</b>'
                         f'<span>{esc(ev.rationale or ev.title)}</span>'
                         f'<em>{esc(ev.status)}{f" · {ev.elapsed_sec:.0f}s" if ev.elapsed_sec else ""}</em></li>')
        exec_sum = self.run.narrative.get("executive_summary", "")
        lead = ""
        if self.run.config.get("use_llm") and exec_sum:
            first = exec_sum.split("</p>")[0] + "</p>" if self.run.narrative.get("html") else exec_sum.split("\n\n")[0]
            lead = f'<div class="lead">{self._md(first)}<p class="more"><a href="#summary" data-part="summary">Read the full summary →</a></p></div>'
        return (f'<div class="tiles">{tiles}</div>{lead}'
                f'<div class="two"><section><h3>What you told the agent</h3>'
                f'{"<table class=kv>" + bio + "</table>" if bio else ""}{notes}'
                f'{"" if bio or notes else "<p class=muted>No prior biology was given.</p>"}</section>'
                f'<section><h3>What the agent did</h3><ol class="steps">{"".join(steps)}</ol></section></div>')

    def _md(self, text: str) -> str:
        return text if self.run.narrative.get("html") else md(text, self.store)

    def _mdi(self, text: str) -> str:
        return text if self.run.narrative.get("html") else md_inline(text, self.store)

    def key_results(self) -> List[tuple]:
        """Deterministic one-line facts from the evidence, grouped by topic, each citing its card."""
        st = self.store
        cite = lambda ev: f' <a class="cite" href="#{ev.id}">{ev.id}</a>'
        groups: List[tuple] = []

        def add(topic, part, lines):
            lines = [ln for ln in lines if ln]
            if lines:
                groups.append((topic, part, lines))

        qc, ann = st.latest("peak_qc"), st.latest("annotate_peaks")
        lines = []
        if qc:
            w = qc.summary.get("width", {})
            lines.append(f"<b>{qc.summary.get('n_clean_peaks')}</b> peaks after QC (of {qc.summary.get('n_peaks')}; "
                         f"median width {w.get('median')} bp, {qc.summary.get('n_duplicated', 0)} duplicates, "
                         f"{qc.summary.get('n_blacklisted', 0)} blacklisted){cite(qc)}")
        if ann:
            rc = ann.summary.get("region_class_counts", {})
            tot = sum(rc.values()) or 1
            gs = ann.summary.get("gene_set_sizes", {})
            win = ann.summary.get("promoter_window", {})
            lines.append(f"<b>{100 * rc.get('promoter', 0) / tot:.0f}%</b> promoter-proximal (−{win.get('upstream')}/+{win.get('downstream')} bp), "
                         f"<b>{100 * (rc.get('distal', 0) + rc.get('distal_far', 0)) / tot:.0f}%</b> distal; "
                         f"median distance to a TSS {ann.summary.get('median_abs_distance_to_tss')} bp{cite(ann)}")
            lines.append(f"<b>{gs.get('all')}</b> linked genes ({gs.get('promoter')} promoter, {gs.get('distal')} distal), "
                         f"e.g. {', '.join((ann.summary.get('example_promoter_genes') or [])[:8])}{cite(ann)}")
        add("Peaks and genes", "annotation", lines)

        lines = []
        for ev in st.items.values():
            if ev.tool == "enrich_pathways" and ev.status == "ok":
                by_method: Dict[str, List[str]] = {}
                for t in ev.summary.get("top_terms") or []:
                    by_method.setdefault(t["method"], [])
                    if t["term"] not in by_method[t["method"]] and len(by_method[t["method"]]) < 3:
                        by_method[t["method"]].append(t["term"])
                for m, terms in by_method.items():
                    lines.append(f"<b>{esc(m)}</b>: {esc('; '.join(terms))}{cite(ev)}")
            if ev.tool == "network_viz" and ev.status == "ok":
                cl = ev.summary.get("clusters") or []
                lines.append(f"Network themes ({esc(ev.summary.get('method'))}, {esc(ev.summary.get('group'))} peaks): "
                             + "; ".join(f"<b>{esc(c['label'])}</b> ({c['n_pathways']})" for c in cl[:5]) + cite(ev))
            if ev.tool == "pathway_consensus" and ev.status == "ok":
                lists = ev.summary.get("lists") or {}
                lname = "all" if "all" in lists else next(iter(lists), None)
                if lname:
                    th = lists[lname].get("top_themes") or []
                    kg = lists[lname].get("key_genes") or []
                    if th:
                        lines.append(f"Consensus themes ({esc(lname)} peaks): " + "; ".join(esc(t["theme"]) for t in th[:4])
                                     + (f" · key genes <b>{esc(', '.join(k['gene'] for k in kg[:6]))}</b>" if kg else "") + cite(ev))
        add("Pathways", "pathways", lines)

        lines = []
        for ev in st.items.values():
            if ev.tool == "find_motifs" and ev.status == "ok":
                top = ev.summary.get("top_motifs") or []
                reg = ev.params.get("regions", "all")
                if top:
                    lines.append(f"{esc(reg)} peaks: {ev.summary.get('n_significant')} significant motifs; top "
                                 + ", ".join(f"<b>{esc(m['tf'])}</b>" for m in top[:6]) + cite(ev))
                else:
                    lines.append(f"{esc(reg)} peaks: no motif passed FDR{cite(ev)}")
        add("Motifs", "motifs", lines)

        lines = []
        ca = st.latest("chipatlas_collect")
        if ca and ca.summary.get("top_antigens"):
            lines.append("Public ChIP-seq overlapping these peaks most: "
                         + ", ".join(f"<b>{esc(a['antigen'])}</b>" for a in ca.summary["top_antigens"][:8])
                         + f" ({ca.summary.get('cell_class')}){cite(ca)}")
        for ev in st.items.values():
            if ev.tool == "chipatlas_target_genes" and ev.status == "ok":
                a = ev.summary.get("any", {})
                lines.append(f"Known <b>{esc(ev.summary.get('tf'))}</b> targets: {a.get('pct')}% of peak genes vs "
                             f"{a.get('background_pct')}% background (p {_f(a.get('pvalue'))}){cite(ev)}")
        add("Regulators", "regulators", lines)

        lines = []
        rn = st.latest("integrate_rnaseq")
        if rn:
            up, dn = rn.summary.get("up_any_peak", {}), rn.summary.get("down_any_peak", {})
            lines.append(f"RNA-seq: {up.get('pct_with_peak')}% of up-regulated and {dn.get('pct_with_peak')}% of "
                         f"down-regulated genes carry a peak (background {up.get('background_pct')}%){cite(rn)}")
        exp = st.latest("check_expectations")
        if exp:
            missing = [g["gene"] for g in exp.summary.get("genes") or [] if not g.get("bound")]
            lines.append(f"Expected genes bound <b>{exp.summary.get('genes_recovered')}</b>"
                         + (f" (not bound: {esc(', '.join(missing))})" if missing else "")
                         + f"; expected TFs recovered <b>{exp.summary.get('tfs_recovered')}</b>{cite(exp)}")
        add("Expression and controls", "expression", lines)

        lines = []
        if st.references:
            lines.append(f"{len(st.references)} papers and {len([d for d in st.datasets.values() if d.get('accession')])} "
                         f"GEO series retrieved")
        add("Literature", "literature", lines)
        return groups

    def summary(self) -> str:
        n = self.run.narrative
        out = []
        if self.run.config.get("use_llm") and n.get("executive_summary"):
            out.append(f'<section class="block"><h3>The agent\'s reading</h3>'
                       f'<div class="prose lead-prose">{self._md(n.get("executive_summary", ""))}</div></section>')
        groups = self.key_results()
        if groups:
            out.append('<section class="block"><h3>Key results</h3><div class="kr">')
            for topic, part, lines in groups:
                out.append(f'<div class="kr-row"><div class="kr-topic"><a href="#{part}" data-part="{part}">{esc(topic)}</a></div>'
                           f'<ul>' + "".join(f"<li>{ln}</li>" for ln in lines) + "</ul></div>")
            out.append("</div></section>")
        cav = n.get("caveats") or []
        if cav:
            out.append('<section class="block"><h3>Caveats</h3><ul class="prose">'
                       + "".join(f"<li>{self._mdi(c)}</li>" for c in cav) + "</ul></section>")
        return "".join(out)

    def suggestions(self) -> str:
        hyps = self.run.narrative.get("hypotheses") or []
        if not hyps:
            return ('<p class="muted">No suggestions: this run used the deterministic plan. Run in agent mode '
                    'to get testable hypotheses, each with a proposed next experiment.</p>')
        out = ['<p class="muted">Testable hypotheses the agent drew from the results, each with its evidence and '
               'the experiment that would test it.</p>']
        for i, h in enumerate(hyps, 1):
            cites = h.get("evidence_html") or " ".join(md_inline(f"[{e}]", self.store) for e in h.get("evidence", []))
            out.append(f'<article class="card hyp"><header><span class="hn">{i}</span>'
                       f'<h3>{self._mdi(h.get("statement", ""))}</h3>'
                       f'<span class="pill conf {esc(h.get("confidence"))}">{esc(h.get("confidence"))}</span></header>'
                       f'<p class="muted">Evidence {cites}</p>'
                       f'<p><b>Next experiment.</b> {self._mdi(h.get("next_experiment", ""))}</p></article>')
        return "".join(out)

    def interpretation(self) -> str:
        secs = self.run.narrative.get("sections") or []
        if not secs:
            return '<p class="muted">No interpretation: deterministic run.</p>'
        return "".join(f'<section class="prose"><h3>{self._mdi(s.get("heading", ""))}</h3>{self._md(s.get("body", ""))}</section>'
                       for s in secs)

    def tool_part(self, tools) -> str:
        cards = [self.card(e) for e in self.store.items.values() if e.tool in tools and e.tool not in HIDDEN]
        return "".join(cards) if cards else '<p class="muted">This step was not run.</p>'

    def literature(self) -> str:
        parts = [self.tool_part(("search_literature", "search_geo"))]
        exps = [(k, d) for k, d in self.store.datasets.items() if not d.get("accession")]
        if exps:
            parts.append("<h3>Public experiments cited (ChIP-Atlas)</h3><ul class='refs'>" + "".join(
                f"<li><a href='{esc(d.get('url'))}' target='_blank' rel='noopener'>{esc(k)}</a> · "
                f"{esc(d.get('antigen'))}, {esc(d.get('cell'))}</li>" for k, d in exps) + "</ul>")
        return "".join(parts)

    def methods(self) -> str:
        class Cfg:  # methods_text needs config attributes
            blacklist = self.run.config.get("blacklist")
            use_llm = self.run.config.get("use_llm")
            model = self.run.config.get("model")
        text, versions = methods_text(self.store, Cfg)
        out = f'<div class="prose"><p>{esc(text)}</p></div>'
        if versions:
            out += "<h3>Software versions</h3>" + _table([{"package": k, "version": v} for k, v in sorted(versions.items())], 100)
        return out

    def log(self) -> str:
        return _table([{"step": i + 1, "time": t["time"], "id": t["evidence"], "tool": t["tool"], "status": t["status"],
                        "why": t["rationale"]} for i, t in enumerate(self.store.trace)], 500)

    def parts(self) -> List[Dict]:
        done_tools = {e.tool for e in self.store.items.values() if e.status == "ok"}
        out = []
        for pid, label, tools in PARTS:
            if tools is None:
                html_ = getattr(self, pid)()
                count = {"suggestions": len(self.run.narrative.get("hypotheses") or []),
                         "log": len(self.store.trace)}.get(pid)
                empty = (pid == "suggestions" and not self.run.narrative.get("hypotheses")) or \
                        (pid == "interpretation" and not self.run.narrative.get("sections"))
            else:
                html_ = self.tool_part(tools) if pid != "literature" else self.literature()
                count = sum(1 for e in self.store.items.values() if e.tool in tools and e.tool not in HIDDEN)
                empty = not (set(tools) & done_tools)
            out.append({"id": pid, "label": label, "html": html_, "count": count, "empty": empty})
        return out
