"""Search over a finished run: genes, peaks, pathways, motifs, regulators, papers and
datasets, all pointing back to the evidence card and table they came from.

A `Run` wraps an output directory (<prefix>_manifest.json, <prefix>_evidence.json,
<prefix>_report.html). `RunIndex` flattens the run into searchable records."""

from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from .analysis import read_table
from .evidence import EvidenceStore
from .report import auto_narrative

KINDS = ("gene", "peak", "pathway", "module", "motif", "regulator", "expression", "control", "paper",
         "dataset", "hypothesis", "section", "result")
MAX_ROWS = 50000


@dataclass
class Record:
    kind: str
    title: str
    text: str                      # what is matched against the query (lower-cased on search)
    evidence: str = ""             # E-id the record comes from
    detail: str = ""               # one line shown under the title
    data: dict = field(default_factory=dict)
    score: float = 0.0


class Run:
    """A run output directory, finished or in progress."""

    def __init__(self, run_dir: str, prefix: Optional[str] = None):
        self.dir = Path(run_dir).expanduser().resolve()
        if prefix is None:
            manifests = sorted(self.dir.glob("*_manifest.json")) or sorted(self.dir.glob("*_evidence.json"))
            if not manifests:
                raise FileNotFoundError(f"No PeakAnnoAgent run found in {self.dir}")
            prefix = re.sub(r"_(manifest|evidence)\.json$", "", manifests[0].name)
        self.prefix = prefix
        self.manifest: dict = {}
        mpath = self.dir / f"{prefix}_manifest.json"
        if mpath.is_file():
            self.manifest = json.loads(mpath.read_text())
        epath = self.dir / f"{prefix}_evidence.json"
        self.store = EvidenceStore.load(str(epath)) if epath.is_file() else EvidenceStore()
        if not self.manifest.get("config") and self.report:
            self.manifest["config"] = self._config_from_report()
        self.narrative = self.manifest.get("narrative") or (self._narrative_from_report() if self.report else None) \
            or auto_narrative(self.store)

    def _narrative_from_report(self) -> Optional[dict]:
        """Older runs saved the agent's narrative only in the report: pull its HTML back out.
        Fields are already HTML (flag 'html': True); citation anchors (#E4) are kept."""
        text = self.report.read_text(errors="replace")
        if "automated summary" in text.lower() and "no llm interpretation" in text.lower():
            return None
        m = re.search(r"<div class='summary'>(.*?)</div>\s*<h2 id='hypotheses'>", text, re.S)
        if not m:
            return None
        summary = m.group(1)
        hyps = [{"statement": h.group(1), "confidence": h.group(2), "evidence_html": h.group(3),
                 "next_experiment": h.group(4), "evidence": re.findall(r">(E\d+)<", h.group(3))}
                for h in re.finditer(r"<div class='card'><header><h3>H\d+\. (.*?)</h3><span class='conf (\w+)'>"
                                     r".*?Evidence: (.*?)</p><p><b>Test:</b> (.*?)</p></div>", text, re.S)]
        inter = re.search(r"<h2 id='interpretation'>Interpretation</h2>(.*?)<h2 id='results'>", text, re.S)
        sections, caveats = [], []
        if inter:
            body = inter.group(1)
            cav = re.search(r"<h3>Caveats</h3><ul>(.*?)</ul>", body, re.S)
            if cav:
                caveats = re.findall(r"<li>(.*?)</li>", cav.group(1), re.S)
                body = body[: cav.start()]
            for sec in re.finditer(r"<h3>(.*?)</h3>(.*?)(?=<h3>|$)", body, re.S):
                sections.append({"heading": sec.group(1), "body": sec.group(2)})
        title = re.search(r"<title>(.*?)</title>", text, re.S)
        return {"html": True, "title": title.group(1) if title else "", "executive_summary": summary,
                "hypotheses": hyps, "sections": sections, "caveats": caveats}

    def _config_from_report(self) -> dict:
        """Runs without a manifest (older outputs): recover genome, biology and mode from the report header."""
        import html as _html
        text = self.report.read_text(errors="replace")[:20000]
        meta = {m.group(1): _html.unescape(m.group(2)).strip()
                for m in re.finditer(r"<div><b>([^<]+)</b>([^<]*)</div>", text)}
        blank = lambda v: "" if v in ("-", None) else v
        title = re.search(r"<title>(.*?)</title>", text, re.S)
        return {"genome": meta.get("Genome", ""), "peaks": re.sub(r"\s*\(\d+\)$", "", meta.get("Peaks", "")),
                "use_llm": meta.get("Mode", "").startswith("agent"),
                "model": re.sub(r"^agent \(|\)$", "", meta.get("Mode", "")) if meta.get("Mode", "").startswith("agent") else None,
                "biology": {"assay": blank(meta.get("Assay")), "target": blank(meta.get("Target")),
                            "cell_type": blank(meta.get("Cell type")), "condition": blank(meta.get("Condition")),
                            "description": "", "expected_genes": [], "expected_tfs": []},
                "title": _html.unescape(title.group(1)) if title else ""}

    @property
    def config(self) -> dict:
        return self.manifest.get("config") or {}

    @property
    def biology(self) -> dict:
        return self.config.get("biology") or {}

    @property
    def genome(self) -> str:
        return self.config.get("genome", "")

    @property
    def report(self) -> Optional[Path]:
        p = self.dir / f"{self.prefix}_report.html"
        return p if p.is_file() else None

    def resolve(self, path: str) -> Optional[Path]:
        """Tables are stored with the absolute path of the machine that ran the analysis;
        fall back to the same file name under this run's results/ when it moved."""
        p = Path(path)
        if p.is_file():
            return p
        for cand in (self.dir / "results" / p.name, self.dir / "results" / ".preview" / p.name, self.dir / p.name):
            if cand.is_file():
                return cand
        return None

    def table(self, evidence_id: str, name: Optional[str] = None) -> List[dict]:
        ev = self.store.items.get(evidence_id)
        if not ev or not ev.tables:
            return []
        key = name or next(iter(ev.tables))
        path = self.resolve(ev.tables.get(key, ""))
        return read_table(str(path)) if path else []

    def summary(self) -> dict:
        """Short JSON used by the web UI and as chat context."""
        qc = self.store.latest("peak_qc")
        return {
            "dir": str(self.dir), "prefix": self.prefix, "genome": self.genome,
            "peaks": Path(self.config.get("peaks", "")).name,
            "n_peaks": qc.summary.get("n_peaks") if qc else None,
            "mode": "agent" if self.config.get("use_llm") else "deterministic",
            "model": self.config.get("model") if self.config.get("use_llm") else None,
            "biology": self.biology, "started_at": self.config.get("started_at"),
            "finished_at": self.manifest.get("finished_at"), "title": self.narrative.get("title"),
            "n_steps": len(self.store.trace), "report": self.report.name if self.report else None,
            "n_references": len(self.store.references), "n_datasets": len(self.store.datasets),
        }


# ------------------------------------------------------------------ index

def _tokens(q: str) -> List[str]:
    return [t for t in re.split(r"[\s,;]+", q.lower().strip()) if t]


def _f(v, digits=3) -> str:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return str(v)
    return f"{x:.2e}" if x and abs(x) < 1e-3 else f"{x:.{digits}g}"


class RunIndex:
    def __init__(self, run: Run):
        self.run = run
        self.records: List[Record] = []
        self.genes: Dict[str, dict] = defaultdict(lambda: {"peaks": [], "pathways": [], "expression": None,
                                                           "control": None, "hub_in": []})
        self._build()

    # -- building ----------------------------------------------------------
    def _add(self, rec: Record) -> None:
        self.records.append(rec)

    def _build(self) -> None:
        run, store = self.run, self.run.store
        for ev in store.items.values():
            if ev.tool in ("inspect_table", "chipatlas_status"):
                continue
            self._add(Record("result", f"{ev.id} {ev.title}", f"{ev.id} {ev.title} {ev.tool} {ev.rationale}",
                             ev.id, f"{ev.tool} · {ev.status}" + (f" · {ev.error[:120]}" if ev.error else ""),
                             {"tool": ev.tool, "status": ev.status}))
        self._peaks()
        self._pathways()
        self._motifs()
        self._regulators()
        self._expression()
        self._controls()
        for key, p in store.references.items():
            self._add(Record("paper", p.get("title", key),
                             " ".join(str(p.get(k, "")) for k in ("title", "authors", "journal", "abstract", "id")),
                             "", f"{p.get('authors', '')} · {p.get('journal', '')} {p.get('year', '')} · {key}",
                             {"pmid": p.get("pmid"), "doi": p.get("doi"), "year": p.get("year"),
                              "url": f"https://pubmed.ncbi.nlm.nih.gov/{p['pmid']}/" if p.get("pmid") else
                              (f"https://doi.org/{p['doi']}" if p.get("doi") else "")}))
        for key, d in store.datasets.items():
            if d.get("accession"):
                self._add(Record("dataset", f"{key} {d.get('title', '')}",
                                 " ".join(str(d.get(k, "")) for k in ("accession", "title", "summary", "type")),
                                 "", f"{d.get('type', '')} · {d.get('n_samples', '')} samples · {d.get('date', '')}",
                                 {"url": d.get("url"), "accession": key}))
            else:
                self._add(Record("dataset", f"{key} {d.get('antigen', '')} ({d.get('cell', '')})",
                                 f"{key} {d.get('antigen', '')} {d.get('cell', '')} chip-atlas",
                                 "", "ChIP-Atlas experiment", {"url": d.get("url"), "accession": key}))
        n = run.narrative or {}
        if n.get("html"):                       # narrative recovered from a report: index plain text
            import html as _html
            plain = lambda t: _html.unescape(re.sub(r"<[^>]+>", "", t or ""))
            n = {**n, "executive_summary": plain(n.get("executive_summary")),
                 "hypotheses": [{**h, "statement": plain(h.get("statement")), "next_experiment": plain(h.get("next_experiment"))}
                                for h in n.get("hypotheses") or []],
                 "sections": [{"heading": plain(s_.get("heading")), "body": plain(s_.get("body"))} for s_ in n.get("sections") or []]}
        for i, h in enumerate(n.get("hypotheses") or [], 1):
            self._add(Record("hypothesis", f"Suggestion {i}. {h.get('statement', '')}",
                             f"{h.get('statement', '')} {h.get('next_experiment', '')}",
                             (h.get("evidence") or [""])[0],
                             f"{h.get('confidence', '')} confidence · test: {h.get('next_experiment', '')}",
                             {"evidence": h.get("evidence", []), "anchor": "suggestions"}))
        for s in n.get("sections") or []:
            self._add(Record("section", s.get("heading", ""), f"{s.get('heading', '')} {s.get('body', '')}", "",
                             (s.get("body") or "")[:160], {"anchor": "interpretation"}))
        if n.get("executive_summary"):
            self._add(Record("section", "Summary", n["executive_summary"], "", n["executive_summary"][:160],
                             {"anchor": "summary"}))
        for gene, info in self.genes.items():
            bits = []
            if info["peaks"]:
                bits.append(f"{len(info['peaks'])} peak(s): " + ", ".join(
                    f"{p['PeakID']} ({p['region_class']}, {p['distance_to_tss']} bp)" for p in info["peaks"][:4]))
            if info["pathways"]:
                bits.append(f"{len(info['pathways'])} pathway(s)")
            if info["expression"]:
                bits.append(f"RNA-seq {info['expression']['direction']} log2FC {info['expression']['log2fc']}")
            if info["control"] is not None:
                bits.append("expected gene: " + ("bound" if info["control"] else "NOT bound"))
            self._add(Record("gene", gene, f"{gene} " + " ".join(
                p.get("PeakID", "") for p in info["peaks"]), info.get("evidence", ""), " · ".join(bits),
                {"n_peaks": len(info["peaks"]), "pathways": info["pathways"][:12],
                 "expression": info["expression"], "control": info["control"], "hub_in": info["hub_in"]}))

    def _peaks(self) -> None:
        ann = self.run.store.latest("annotate_peaks")
        if not ann:
            return
        rows = self.run.table(ann.id, "peak_genes")[:MAX_ROWS]
        for r in rows:
            genes = [g for g in (r.get("promoter_symbols") or "").split(";") if g and g != "NA"]
            if r.get("nearest_symbol") and r["nearest_symbol"] != "NA":
                genes.append(r["nearest_symbol"])
            genes = list(dict.fromkeys(genes))
            locus = f"{r.get('Chr')}:{r.get('Start')}-{r.get('End')}"
            self._add(Record("peak", f"{r.get('PeakID')} {locus}",
                             f"{r.get('PeakID')} {locus} {r.get('Chr')} " + " ".join(genes),
                             ann.id, f"{r.get('region_class')} · nearest {r.get('nearest_symbol')} "
                                     f"({r.get('distance_to_tss')} bp to TSS)" +
                             (f" · promoter of {', '.join(genes[:-1])}" if len(genes) > 1 else ""),
                             {"locus": locus, "genes": genes, "region_class": r.get("region_class")}))
            for g in genes:
                self.genes[g]["peaks"].append(r)
                self.genes[g]["evidence"] = ann.id

    def _pathways(self) -> None:
        for ev in self.run.store.items.values():
            if ev.tool != "enrich_pathways" or ev.status != "ok":
                continue
            for r in self.run.table(ev.id, "enrichment")[:MAX_ROWS]:
                genes = [g for g in (r.get("genes") or "").split(";") if g]
                self._add(Record("pathway", r.get("term", ""),
                                 f"{r.get('term', '')} {r.get('term_id', '')} {r.get('method', '')} " + " ".join(genes),
                                 ev.id, f"{r.get('method')} · {r.get('group')} peaks · padj {_f(r.get('padj'))} · "
                                        f"{r.get('count')} genes", {"method": r.get("method"), "group": r.get("group"),
                                                                    "padj": r.get("padj"), "genes": genes[:40],
                                                                    "term_id": r.get("term_id")}))
                for g in genes:
                    self.genes[g]["pathways"].append(f"{r.get('term')} ({r.get('method')}, {r.get('group')})")
            for method, mods in (ev.summary.get("pathway_modules") or {}).items():
                if not isinstance(mods, list):
                    continue
                for m in mods:
                    self._add(Record("module", f"{method} module {m.get('module')}: {m.get('lead_term')}",
                                     f"{m.get('lead_term', '')} {m.get('other_terms', '')} {method}", ev.id,
                                     f"{m.get('n_terms')} terms · best padj {_f(m.get('best_padj'))} · "
                                     f"{m.get('other_terms', '')}"[:200], {"method": method}))
            for method, hubs in (ev.summary.get("hub_genes") or {}).items():
                if isinstance(hubs, list):
                    for h in hubs:
                        self.genes[h["gene"]]["hub_in"].append(f"{method} ({h.get('n_pathways')} pathways)")
        for ev in self.run.store.items.values():
            if ev.tool == "network_viz" and ev.status == "ok":
                for c in ev.summary.get("clusters") or []:
                    self._add(Record("module", f"Network cluster {c.get('cluster')}: {c.get('label')}",
                                     f"{c.get('label', '')} {c.get('members', '')} network", ev.id,
                                     f"{c.get('n_pathways')} pathways · best FDR {_f(c.get('best_fdr'))} · {c.get('members', '')}"[:220]))
            if ev.tool == "pathway_consensus" and ev.status == "ok":
                for lname, lst in (ev.summary.get("lists") or {}).items():
                    for t in lst.get("top_themes") or []:
                        self._add(Record("module", f"Consensus theme ({lname}): {t.get('theme')}",
                                         f"{t.get('theme', '')} {t.get('core_genes', '')} consensus pathwayagent", ev.id,
                                         f"{t.get('tier')} · score {_f(t.get('score'))} · {t.get('n_terms')} terms · "
                                         f"{t.get('databases')} · {t.get('core_genes')}"))
                    for k in lst.get("key_genes") or []:
                        self.genes[k["gene"]]["hub_in"].append(f"PathwayAgent key gene ({lname}, score {_f(k.get('score'))})")
        for ev in self.run.store.items.values():
            if ev.tool == "great_enrichment" and ev.status == "ok":
                for name in ev.tables:
                    for r in self.run.table(ev.id, name)[:2000]:
                        term = r.get("term") or r.get("description") or r.get("id") or ""
                        self._add(Record("pathway", term, f"{term} GREAT " + " ".join(r.values()), ev.id,
                                         "GREAT · " + ", ".join(f"{k} {_f(v)}" for k, v in r.items()
                                                                if k in ("padj", "p_adjust", "fold_enrichment"))))

    def _motifs(self) -> None:
        for ev in self.run.store.items.values():
            if ev.tool != "find_motifs" or ev.status != "ok":
                continue
            regions = ev.params.get("regions", "all")
            for r in self.run.table(ev.id, "motifs")[:MAX_ROWS]:
                self._add(Record("motif", f"{r.get('tf')} ({regions} peaks)",
                                 f"{r.get('motif', '')} {r.get('tf', '')} {r.get('family', '')} motif",
                                 ev.id, f"{r.get('family', '')} · log2 enrichment {_f(r.get('log2_enrichment'))} · "
                                        f"padj {_f(r.get('padj'))} · {r.get('pct_target')}% of peaks",
                                 {"regions": regions, "padj": r.get("padj"), "log2_enrichment": r.get("log2_enrichment"),
                                  "family": r.get("family"), "motif": r.get("motif")}))

    def _regulators(self) -> None:
        for ev in self.run.store.items.values():
            if ev.tool == "chipatlas_collect" and ev.status == "ok":
                for a in ev.summary.get("top_antigens") or []:
                    self._add(Record("regulator", a.get("antigen", ""),
                                     f"{a.get('antigen', '')} chip-atlas co-binding " + " ".join(a.get("cells") or []),
                                     ev.id, f"{a.get('n_significant_experiments')}/{a.get('n_experiments')} public "
                                            f"experiments · best log10 q {a.get('best_log10_q')} · fold "
                                            f"{a.get('best_fold_enrichment')} · {a.get('best_cell')}",
                                     {"best_experiment": a.get("best_experiment"), "cells": a.get("cells")}))
            if ev.tool == "chipatlas_target_genes" and ev.status == "ok":
                tf = ev.summary.get("tf", ev.params.get("tf"))
                self._add(Record("regulator", f"{tf} known targets", f"{tf} target genes chip-atlas", ev.id,
                                 f"{ev.summary.get('any', {}).get('n_known_targets')} peak genes are known {tf} "
                                 f"targets (p {_f(ev.summary.get('any', {}).get('pvalue'))})"))
                for g in ev.summary.get("strongest_shared_targets") or []:
                    self.genes[g]["hub_in"].append(f"known {tf} target")

    def _expression(self) -> None:
        ev = self.run.store.latest("integrate_rnaseq")
        if not ev:
            return
        for r in self.run.table(ev.id, "rnaseq_integration")[:MAX_ROWS]:
            g = r.get("gene", "")
            self.genes[g]["expression"] = {"direction": r.get("direction"), "log2fc": r.get("log2fc"),
                                           "padj": r.get("padj")}
            self._add(Record("expression", f"{g} {r.get('direction')}", f"{g} rna-seq {r.get('direction')} regulated",
                             ev.id, f"log2FC {r.get('log2fc')} · padj {r.get('padj')} · {r.get('n_peaks')} peak(s) "
                                    f"({r.get('classes')})"))

    def _controls(self) -> None:
        ev = self.run.store.latest("check_expectations")
        if not ev:
            return
        for g in ev.summary.get("genes") or []:
            self.genes[g["gene"]]["control"] = bool(g.get("bound"))
            self._add(Record("control", f"Expected gene {g['gene']}", f"{g['gene']} positive control expected",
                             ev.id, ("bound: %s peak(s), %s" % (g.get("n_peaks"), ", ".join(g.get("classes") or [])))
                             if g.get("bound") else "not bound"))
        for t in ev.summary.get("tfs") or []:
            self._add(Record("control", f"Expected TF {t['tf']}", f"{t['tf']} positive control expected tf motif",
                             ev.id, f"motif rank {t.get('motif_rank') or '-'} · ChIP-Atlas rank "
                                    f"{t.get('chipatlas_rank') or '-'}"))

    # -- querying ----------------------------------------------------------
    def search(self, query: str, kinds: Optional[Iterable[str]] = None, limit: int = 60) -> dict:
        toks = _tokens(query)
        wanted = set(kinds) if kinds else None
        hits: List[Record] = []
        if not toks:
            return {"query": query, "total": 0, "groups": {}}
        for rec in self.records:
            if wanted and rec.kind not in wanted:
                continue
            title, text = rec.title.lower(), rec.text.lower()
            score = 0.0
            for t in toks:
                if t == title or t in re.split(r"[\s(),:;/]+", title):
                    score += 3
                elif t in title:
                    score += 2
                elif t in text:
                    score += 1
                else:
                    score = 0
                    break
            if score:
                hits.append(Record(**{**asdict(rec), "score": score}))
        hits.sort(key=lambda r: (-r.score, KINDS.index(r.kind) if r.kind in KINDS else 99, r.title))
        groups: Dict[str, List[dict]] = defaultdict(list)
        for h in hits[: max(limit, 1) * 4]:
            if len(groups[h.kind]) < limit:
                groups[h.kind].append(asdict(h))
        return {"query": query, "total": len(hits),
                "groups": {k: groups[k] for k in KINDS if k in groups}}

    def gene(self, symbol: str) -> Optional[dict]:
        key = next((g for g in self.genes if g.lower() == symbol.lower()), None)
        if key is None:
            return None
        info = self.genes[key]
        return {"gene": key, "peaks": [{k: p.get(k) for k in ("PeakID", "Chr", "Start", "End", "region_class",
                                                               "distance_to_tss", "nearest_symbol")}
                                       for p in info["peaks"]],
                "pathways": info["pathways"], "expression": info["expression"], "control": info["control"],
                "hub_in": info["hub_in"], "evidence": info.get("evidence")}

    def counts(self) -> Dict[str, int]:
        out: Dict[str, int] = defaultdict(int)
        for r in self.records:
            out[r.kind] += 1
        return dict(out)
