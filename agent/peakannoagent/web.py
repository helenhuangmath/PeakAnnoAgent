"""Web app: submit a peak set with its biology, watch the agent work, then read the report,
search every result, and ask Claude about the run.

    peakanno-agent web --runs-dir runs --port 8766

Runs are directories under --runs-dir; existing result folders can be opened in place.
Requires `pip install 'peakannoagent[web]'` (Flask)."""

from __future__ import annotations

import io
import json
import os
import re
import shutil
import threading
import zipfile
import time
import traceback
import hashlib
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from flask import Flask, abort, jsonify, redirect, render_template, request, send_file, send_from_directory, url_for

from . import __version__
from .artifacts import write_json
from .config import DEFAULT_MODEL, SUPPORTED_GENOMES, RunConfig, build_biology
from .rbridge import RBridge, repo_root
from .report import md
from .search import KINDS, Run, RunIndex
from .explore import explore_data
from .webview import View

REPO = repo_root()
EXAMPLE_DIR = REPO / "examples" / "data"
# Example peak sets offered on the home page. Each has a BED file and a biology JSON; the first also has
# gene signatures and an RNA-seq table. Files are resolved against EXAMPLE_DIR.
EXAMPLES: Dict[str, dict] = {
    "tex": {
        "title": "Exhausted CD8 T cells, ATAC-seq",
        "genome": "mm10", "prefix": "tex", "n_peaks": 351, "kind": "Synthetic",
        "blurb": "351 ATAC peaks at exhaustion and effector genes, with three signatures and an RNA-seq table. "
                 "The full demo: pathways, motifs, regulators, positive controls, literature.",
        "source": "Built by examples/scripts/make_example_tex_peaks.R",
        "files": {"peaks": "example_tex_peaks_mm10.bed", "biology": "example_tex_biology.json",
                  "gmt": "example_tex_signatures.gmt", "rnaseq": "example_tex_rnaseq_de.tsv"},
    },
    "ctcf": {
        "title": "CTCF ChIP-seq, K562",
        "genome": "hg38", "prefix": "ctcf", "n_peaks": 1000, "kind": "ENCODE",
        "blurb": "The 1000 strongest CTCF peaks in K562. A TF with a known answer: the CTCF motif should "
                 "dominate and cohesin should co-bind.",
        "source": "ENCODE ENCSR000BPJ (Myers lab), file ENCFF582SNT, IDR thresholded peaks, GRCh38",
        "files": {"peaks": "example_ctcf_k562_hg38.bed", "biology": "example_ctcf_k562_biology.json"},
    },
    "k27ac": {
        "title": "H3K27ac ChIP-seq, human CD8 T cells",
        "genome": "hg38", "prefix": "k27ac", "n_peaks": 1000, "kind": "ENCODE",
        "blurb": "The 1000 strongest active-chromatin peaks in primary CD8 T cells. An enhancer mark: "
                 "distal motif scan for lineage TFs, T-cell gene sets.",
        "source": "ENCODE ENCSR835OJV (Bernstein lab), file ENCFF283CAQ, replicated peaks, GRCh38",
        "files": {"peaks": "example_h3k27ac_cd8_hg38.bed", "biology": "example_h3k27ac_cd8_biology.json"},
    },
    "k4me3": {
        "title": "H3K4me3 ChIP-seq, mouse macrophages",
        "genome": "mm10", "prefix": "k4me3", "n_peaks": 800, "kind": "ENCODE",
        "blurb": "The 800 strongest promoter-mark peaks in unstimulated bone-marrow macrophages. Nearly all "
                 "promoter: tests the TSS window and gene-set choice; motifs are skipped for a broad mark.",
        "source": "ENCODE ENCSR916IRY (Hoffmann lab), file ENCFF813FPW, replicated peaks, mm10",
        "files": {"peaks": "example_h3k4me3_bmdm_mm10.bed", "biology": "example_h3k4me3_bmdm_biology.json"},
    },
}
DEFAULT_EXAMPLE = "tex"
# Backwards compatibility for code that imported the single example.
EXAMPLE = {k: EXAMPLE_DIR / v for k, v in EXAMPLES[DEFAULT_EXAMPLE]["files"].items()}
# Test data shown on the home page: what each file is, so a new user knows what to upload.
TEST_DATA = [
    ("peaks", "example_tex_peaks_mm10.bed", "351 ATAC-seq peaks (mm10) from exhausted-like CD8 T cells; BED with peak ids"),
    ("biology", "example_tex_biology.json", "Known biology: assay, cell type, condition, expected genes and TFs, notes"),
    ("gmt", "example_tex_signatures.gmt", "Three exhaustion signatures (progenitor, terminal, effector-like) to test"),
    ("rnaseq", "example_tex_rnaseq_de.tsv", "Differential expression (gene, log2FoldChange, padj) for 47 genes"),
    ("peaks", "example_ctcf_k562_hg38.bed", "1000 strongest CTCF ChIP-seq peaks in K562 (ENCODE ENCFF582SNT, GRCh38)"),
    ("biology", "example_ctcf_k562_biology.json", "Known biology for the CTCF example"),
    ("peaks", "example_h3k27ac_cd8_hg38.bed", "1000 strongest H3K27ac peaks in human CD8 T cells (ENCODE ENCFF283CAQ, GRCh38)"),
    ("biology", "example_h3k27ac_cd8_biology.json", "Known biology for the H3K27ac example"),
    ("peaks", "example_h3k4me3_bmdm_mm10.bed", "800 strongest H3K4me3 peaks in mouse macrophages (ENCODE ENCFF813FPW, mm10)"),
    ("biology", "example_h3k4me3_bmdm_biology.json", "Known biology for the H3K4me3 example"),
]


def example_paths(key: str) -> Optional[Dict[str, Path]]:
    """Resolved file paths of an example, or None if the key is unknown or a file is missing."""
    if key == "1":  # old links: /?example=1
        key = DEFAULT_EXAMPLE
    ex = EXAMPLES.get(key)
    if not ex:
        return None
    paths = {k: EXAMPLE_DIR / v for k, v in ex["files"].items()}
    return paths if all(p.is_file() for p in paths.values()) else None


def example_form(key: str) -> Optional[dict]:
    """Form values that fill the new-run form with an example, so it can be edited before running."""
    paths = example_paths(key)
    if not paths:
        return None
    if key == "1":
        key = DEFAULT_EXAMPLE
    ex = EXAMPLES[key]
    bio = json.loads(paths["biology"].read_text())
    return {**{k: (", ".join(v) if isinstance(v, list) else v) for k, v in bio.items()},
            "name": f"example-{key}", "prefix": ex["prefix"], "genome": ex["genome"], "example": key,
            "example_title": ex["title"], "has_gmt": "gmt" in paths, "has_rnaseq": "rnaseq" in paths}
DEMO_RESULTS = [("demo-with-agent", REPO / "docs" / "demo" / "with_agent", "Demo: example analysed by the agent"),
                ("demo-without-agent", REPO / "docs" / "demo" / "without_agent", "Demo: example, deterministic plan")]
MAX_LOG = 2000


def slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", text).strip("-")[:40] or "run"


# ------------------------------------------------------------------ jobs

@dataclass
class Job:
    id: str
    dir: Path
    prefix: str
    status: str = "queued"               # queued | running | done | partial | failed | error | interrupted
    log: List[str] = field(default_factory=list)
    started: str = ""
    finished: str = ""
    error: str = ""
    name: str = ""
    run: Optional[Run] = None
    index: Optional[RunIndex] = None
    chat: object = None
    lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    @property
    def status_file(self) -> Path:
        return self.dir / "status.json"

    def add_log(self, msg: str) -> None:
        with self.lock:
            self.log.append(f"[{time.strftime('%H:%M:%S')}] {msg}")
            del self.log[:-MAX_LOG]
        self.persist()

    def persist(self) -> None:
        try:
            write_json(self.status_file, {"status": self.status, "started": self.started,
                                                    "finished": self.finished, "error": self.error,
                                                    "name": self.name, "prefix": self.prefix, "log": self.log[-400:]})
        except OSError:
            pass

    def load_run(self) -> Optional[Run]:
        if self.run is None:
            try:
                self.run = Run(str(self.dir), self.prefix)
            except FileNotFoundError:
                return None
            self.index = None
        return self.run

    def load_index(self) -> Optional[RunIndex]:
        if self.index is None and self.load_run() is not None:
            self.index = RunIndex(self.run)
        return self.index

    def to_dict(self) -> dict:
        d = {"id": self.id, "name": self.name, "status": self.status, "started": self.started,
             "finished": self.finished, "error": self.error, "prefix": self.prefix, "dir": str(self.dir),
             "n_log": len(self.log)}
        run = self.load_run() if self.status in ("done", "partial", "failed", "interrupted", "error") else None
        if run is not None:
            d["summary"] = run.summary()
        return d


class JobManager:
    def __init__(self, runs_dir: Path, rscript: Optional[str] = None):
        self.runs_dir = runs_dir
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.rscript = rscript or os.environ.get("PEAKANNOAGENT_RSCRIPT", "Rscript")
        self.jobs: Dict[str, Job] = {}
        self.registry = runs_dir / "_runs.json"       # id -> external directory (opened in place)
        self.lock = threading.Lock()
        self.discover()
        self.add_demos()

    def add_demos(self) -> None:
        """Committed demo outputs (docs/demo) appear in the run list so search and chat can be tried at once."""
        for job_id, path, name in DEMO_RESULTS:
            if path.is_dir() and job_id not in self.jobs and not any(j.dir == path.resolve() for j in self.jobs.values()):
                job = self._from_dir(job_id, path.resolve())
                if job:
                    job.name = name
                    with self.lock:
                        self.jobs[job_id] = job

    # -- discovery ---------------------------------------------------------
    def _register(self, job_id: str, path: Path) -> None:
        with self.lock:
            data = {}
            if self.registry.is_file():
                try:
                    data = json.loads(self.registry.read_text())
                except ValueError:
                    data = {}
            data[job_id] = str(path)
            write_json(self.registry, data)

    def _from_dir(self, job_id: str, path: Path) -> Optional[Job]:
        manifests = sorted(path.glob("*_manifest.json")) or sorted(path.glob("*_evidence.json"))
        status = {}
        if (path / "status.json").is_file():
            try:
                status = json.loads((path / "status.json").read_text())
            except ValueError:
                status = {}
        if not manifests and not status:
            return None
        prefix = status.get("prefix") or (re.sub(r"_(manifest|evidence)\.json$", "", manifests[0].name)
                                         if manifests else "run")
        job = Job(id=job_id, dir=path, prefix=prefix, name=status.get("name") or path.name,
                  started=status.get("started", ""), finished=status.get("finished", ""),
                  error=status.get("error", ""), log=list(status.get("log") or []))
        completion = {}
        manifest_path = path / f"{prefix}_manifest.json"
        if manifest_path.is_file():
            try:
                completion = json.loads(manifest_path.read_text())
            except ValueError:
                pass
        # A persisted analysis outcome takes precedence over stale worker status.
        st = {"success": "done", "partial": "partial", "failed": "failed"}.get(
            completion.get("status"), status.get("status", "interrupted"))
        if completion.get("reasons"):
            job.error = "; ".join(completion["reasons"])
        job.status = "interrupted" if st in ("queued", "running") else st
        return job

    def discover(self) -> None:
        with self.lock:
            for path in sorted(self.runs_dir.iterdir()):
                if path.is_dir() and path.name not in self.jobs:
                    job = self._from_dir(path.name, path)
                    if job:
                        self.jobs[path.name] = job
            if self.registry.is_file():
                try:
                    for job_id, p in json.loads(self.registry.read_text()).items():
                        if job_id not in self.jobs and Path(p).is_dir():
                            job = self._from_dir(job_id, Path(p))
                            if job:
                                self.jobs[job_id] = job
                except ValueError:
                    pass

    def open_existing(self, path: str) -> Job:
        p = Path(path).expanduser().resolve()
        if not p.is_dir():
            raise FileNotFoundError(f"Not a directory: {p}")
        for job in self.jobs.values():
            if job.dir == p:
                return job
        job_id = f"ext-{slug(p.name)}-{hashlib.sha256(str(p).encode()).hexdigest()[:16]}"
        job = self._from_dir(job_id, p)
        if job is None:
            raise FileNotFoundError(f"No PeakAnnoAgent run (<prefix>_manifest.json / _evidence.json) in {p}")
        with self.lock:
            self.jobs[job_id] = job
        self._register(job_id, p)
        return job

    def ordered(self) -> List[Job]:
        return sorted(self.jobs.values(), key=lambda j: (j.started or j.finished or j.id), reverse=True)

    # -- running -----------------------------------------------------------
    def new_job(self, name: str, prefix: str) -> Job:
        job_id = f"{datetime.now():%Y%m%d-%H%M%S}-{slug(name)}-{uuid.uuid4().hex}"
        d = self.runs_dir / job_id
        d.mkdir(exist_ok=False)
        (d / "input").mkdir()
        job = Job(id=job_id, dir=d, prefix=prefix, name=name)
        with self.lock:
            self.jobs[job_id] = job
        job.persist()
        return job

    def start(self, job: Job, config: RunConfig) -> None:
        from .agent import execute

        def work():
            job.status, job.started = "running", datetime.now().isoformat(timespec="seconds")
            job.persist()
            try:
                manifest = execute(config, log=job.add_log)
                job.status = {"success": "done", "partial": "partial", "failed": "failed"}[manifest["status"]]
                job.error = "; ".join(manifest.get("reasons", []))
            except Exception as err:
                job.status, job.error = "error", f"{type(err).__name__}: {err}"
                job.add_log("ERROR " + job.error)
                job.add_log(traceback.format_exc()[-2000:])
            job.finished = datetime.now().isoformat(timespec="seconds")
            job.run = job.index = job.chat = None
            job.persist()

        threading.Thread(target=work, name=f"peakanno-{job.id}", daemon=True).start()


# ------------------------------------------------------------------ app

def environment(rscript: str) -> dict:
    bridge = RBridge(rscript)
    try:
        import anthropic
        sdk = anthropic.__version__
    except ImportError:
        sdk = None
    from .external import available
    return {"rscript": rscript, "r_found": bool(bridge.available()), "anthropic_sdk": sdk, **available(),
            "api_key": bool(os.environ.get("ANTHROPIC_API_KEY") or os.environ.get("ANTHROPIC_AUTH_TOKEN")),
            "homer": bool(shutil.which("findMotifsGenome.pl"))}


def create_app(runs_dir: str = "runs", rscript: Optional[str] = None, model: str = DEFAULT_MODEL,
               chat_client=None) -> Flask:
    app = Flask(__name__)
    app.config["MAX_CONTENT_LENGTH"] = 512 * 1024 * 1024
    manager = JobManager(Path(runs_dir).expanduser().resolve(), rscript)
    app.manager = manager          # type: ignore[attr-defined]

    def get_job(job_id: str) -> Job:
        job = manager.jobs.get(job_id)
        if job is None:
            manager.discover()
            job = manager.jobs.get(job_id)
        if job is None:
            abort(404)
        return job

    def save_upload(job: Job, key: str) -> Optional[str]:
        f = request.files.get(key)
        if not f or not f.filename:
            return None
        dest = job.dir / "input" / f"{key}-{slug(Path(f.filename).name)}"
        f.save(dest)
        return str(dest)

    # -- pages -------------------------------------------------------------
    def page_args(**extra) -> dict:
        manager.discover()
        args = dict(jobs=[j.to_dict() for j in manager.ordered()], genomes=SUPPORTED_GENOMES, model=model,
                    env=environment(manager.rscript), version=__version__,
                    examples=[(k, ex) for k, ex in EXAMPLES.items() if example_paths(k)],
                    test_data=[(k, n, d, (EXAMPLE_DIR / n).stat().st_size) for k, n, d in TEST_DATA
                               if (EXAMPLE_DIR / n).is_file()])
        args.update(extra)
        return args

    @app.get("/")
    def index():
        # /?example=<key> fills the form with that example so it can be edited before running
        form = example_form(request.args.get("example", "")) if request.args.get("example") else None
        return render_template("index.html", **page_args(form=form))

    @app.get("/examples/<name>")
    def example_file(name):
        if name not in {n for _, n, _ in TEST_DATA}:
            abort(404)
        return send_from_directory(EXAMPLE_DIR, name, as_attachment=True)

    @app.post("/runs")
    def submit():
        form = request.form
        # use_example=<key>: run that example as is, in one click. example=<key>: the form was filled from
        # the example; its files stand in for anything the user did not upload.
        one_click = example_paths(form.get("use_example", "")) if form.get("use_example") else None
        prefill = example_paths(form.get("example", "")) if form.get("example") else None
        ex_key = form.get("use_example") or form.get("example") or ""
        ex_key = DEFAULT_EXAMPLE if ex_key == "1" else ex_key
        ex = EXAMPLES.get(ex_key) if (one_click or prefill) else None
        name = form.get("name") or (f"example-{ex_key}" if one_click else "run")
        prefix = slug(form.get("prefix") or (ex["prefix"] if ex else "run"))
        job = manager.new_job(name, prefix)
        genome = form.get("genome") or (ex["genome"] if ex else "mm10")
        try:
            if one_click:
                peaks = str(shutil.copy(one_click["peaks"], job.dir / "input"))
                file_fields = json.loads(one_click["biology"].read_text())
                gmt = str(shutil.copy(one_click["gmt"], job.dir / "input")) if "gmt" in one_click else None
                rnaseq = str(shutil.copy(one_click["rnaseq"], job.dir / "input")) if "rnaseq" in one_click else None
                blacklist = None
            else:
                peaks = save_upload(job, "peaks")
                if not peaks and form.get("peaks_text", "").strip():
                    peaks = str(job.dir / "input" / "peaks.bed")
                    Path(peaks).write_text(form["peaks_text"].strip() + "\n")
                if not peaks and prefill:
                    peaks = str(shutil.copy(prefill["peaks"], job.dir / "input"))
                if not peaks:
                    raise ValueError("A peak file (BED) is required.")
                file_fields = {}
                gmt, rnaseq, blacklist = save_upload(job, "gmt"), save_upload(job, "rnaseq"), save_upload(job, "blacklist")
                if prefill:
                    if "gmt" in prefill:
                        gmt = gmt or str(shutil.copy(prefill["gmt"], job.dir / "input"))
                    if "rnaseq" in prefill:
                        rnaseq = rnaseq or str(shutil.copy(prefill["rnaseq"], job.dir / "input"))
            biology = build_biology(file_fields, description=form.get("description"), assay=form.get("assay"),
                                    target=form.get("target"), cell_type=form.get("cell_type"),
                                    condition=form.get("condition"), expected_genes=form.get("expected_genes"),
                                    expected_tfs=form.get("expected_tfs"), custom_gmt=gmt, rnaseq=rnaseq)
            config = RunConfig(peaks=peaks, genome=genome, out_dir=str(job.dir), prefix=prefix,
                               biology=biology, blacklist=blacklist, offline=form.get("offline") == "1",
                               use_llm=form.get("mode", "agent") == "agent", model=form.get("model") or model,
                               effort=form.get("effort") or "high", max_turns=int(form.get("max_turns") or 40),
                               rscript=manager.rscript, chipatlas_wait_min=float(form.get("chipatlas_wait") or 25))
            if config.use_llm:
                try:
                    import anthropic  # noqa: F401
                except ImportError:
                    raise ValueError("Agent mode needs the Anthropic SDK (pip install 'peakannoagent[agent]'); "
                                     "choose the deterministic mode instead.")
        except (ValueError, FileNotFoundError, json.JSONDecodeError) as err:
            shutil.rmtree(job.dir, ignore_errors=True)
            manager.jobs.pop(job.id, None)
            return render_template("index.html", **page_args(error=str(err), form=form)), 400
        (job.dir / "input" / "biology.json").write_text(json.dumps(config.biology.__dict__, indent=1))
        manager.start(job, config)
        return redirect(url_for("run_page", job_id=job.id))

    @app.post("/runs/open")
    def open_existing():
        try:
            job = manager.open_existing(request.form.get("path", ""))
        except FileNotFoundError as err:
            return render_template("index.html", **page_args(error=str(err))), 400
        return redirect(url_for("run_page", job_id=job.id))

    @app.get("/runs/<job_id>")
    def run_page(job_id):
        job = get_job(job_id)
        run = job.load_run() if job.status != "running" else None
        report = run.report.name if run and run.report else None
        parts = []
        if run is not None:
            parts = View(run, lambda rel: url_for("run_file", job_id=job.id, filename=rel)).parts()
        return render_template("run.html", job=job.to_dict(), job_log=list(job.log), report=report, kinds=KINDS,
                               summary=run.summary() if run else None, env=environment(manager.rscript),
                               chat_model=model, parts=parts, title=(run.narrative.get("title") if run else None))

    @app.get("/runs/<job_id>/files/<path:filename>")
    def run_file(job_id, filename):
        job = get_job(job_id)
        return send_from_directory(job.dir, filename)

    @app.get("/runs/<job_id>/download.zip")
    def run_zip(job_id):
        """Everything the run produced (report, evidence, results tables and figures) as one zip."""
        job = get_job(job_id)
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            for p in sorted(job.dir.rglob("*")):
                if p.is_file() and not any(part.startswith(".") and part != ".preview" for part in p.relative_to(job.dir).parts):
                    z.write(p, str(Path(slug(job.name)) / p.relative_to(job.dir)))
        buf.seek(0)
        return send_file(buf, mimetype="application/zip", as_attachment=True, download_name=f"{slug(job.name)}.zip")

    @app.get("/api/runs/<job_id>/table/<eid>/<name>")
    def api_table(job_id, eid, name):
        """A result table as JSON rows (first 5000)."""
        job = get_job(job_id)
        run = job.load_run()
        if run is None or eid not in run.store.items or name not in run.store.items[eid].tables:
            abort(404)
        return jsonify({"evidence": eid, "table": name, "rows": run.table(eid, name)[:5000]})

    @app.get("/api/runs/<job_id>/explore/<eid>")
    def api_explore(job_id, eid):
        """Data for the interactive pathway explorer of one enrichment result."""
        job = get_job(job_id)
        run = job.load_run()
        if run is None or eid not in run.store.items:
            abort(404)
        ev = run.store.items[eid]
        return jsonify(explore_data(ev, lambda name: run.table(eid, name)))

    # -- api ---------------------------------------------------------------
    @app.get("/api/runs")
    def api_runs():
        manager.discover()
        return jsonify([j.to_dict() for j in manager.ordered()])

    @app.get("/api/runs/<job_id>")
    def api_run(job_id):
        job = get_job(job_id)
        since = int(request.args.get("since") or 0)
        d = job.to_dict()
        with job.lock:
            d["log"] = job.log[since:]
        return jsonify(d)

    @app.get("/api/runs/<job_id>/search")
    def api_search(job_id):
        job = get_job(job_id)
        index = job.load_index()
        if index is None:
            return jsonify({"error": "This run has no results yet."}), 409
        kinds = [k for k in request.args.getlist("kind") if k in KINDS] or None
        res = index.search(request.args.get("q", ""), kinds, int(request.args.get("limit") or 40))
        for hits in res["groups"].values():
            for h in hits:
                h.pop("text", None)
        res["counts"] = index.counts()
        return jsonify(res)

    @app.get("/api/runs/<job_id>/gene/<symbol>")
    def api_gene(job_id, symbol):
        job = get_job(job_id)
        index = job.load_index()
        if index is None:
            return jsonify({"error": "This run has no results yet."}), 409
        return jsonify(index.gene(symbol) or {"gene": symbol, "found": False})

    def render_entry(job: Job, entry: dict) -> dict:
        run = job.load_run()
        report = run.report.name if run and run.report else None
        html = md(entry.get("text", ""), run.store) if run else ""
        return dict(entry, html=html)

    def chat_session(job: Job):
        if job.chat is None:
            from .chat import ChatSession
            run = job.load_run()
            if run is None:
                return None
            job.chat = ChatSession(run, job.load_index(), client=chat_client, model=model,
                                   offline=bool(run.config.get("offline")))
        return job.chat

    @app.get("/api/runs/<job_id>/chat")
    def api_chat_history(job_id):
        job = get_job(job_id)
        session = chat_session(job)
        if session is None:
            return jsonify({"error": "This run has no results yet."}), 409
        return jsonify({"model": session.model, "messages": [render_entry(job, e) for e in session.transcript]})

    @app.post("/api/runs/<job_id>/chat")
    def api_chat(job_id):
        job = get_job(job_id)
        session = chat_session(job)
        if session is None:
            return jsonify({"error": "This run has no results yet."}), 409
        if chat_client is None:
            try:
                import anthropic  # noqa: F401
            except ImportError:
                return jsonify({"error": "The chat needs the Anthropic SDK: pip install 'peakannoagent[agent]'"}), 503
        text = (request.get_json(silent=True) or {}).get("message", "")
        try:
            entry = session.ask(text)
        except ValueError as err:
            return jsonify({"error": str(err)}), 400
        return jsonify(render_entry(job, entry))

    @app.delete("/api/runs/<job_id>/chat")
    def api_chat_clear(job_id):
        job = get_job(job_id)
        session = chat_session(job)
        if session is not None:
            session.clear()
        return jsonify({"ok": True})

    @app.get("/api/environment")
    def api_env():
        return jsonify(environment(manager.rscript))

    return app


def serve(runs_dir: str = "runs", host: str = "127.0.0.1", port: int = 8766, rscript: Optional[str] = None,
          model: str = DEFAULT_MODEL, open_browser: bool = False, debug: bool = False) -> None:
    import socket
    with socket.socket(socket.AF_INET6 if ":" in host else socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            probe.bind((host, port))
        except OSError:
            raise SystemExit(f"Port {port} on {host} is already used by another program, so PeakAnnoAgent "
                             f"was not started. Stop that program or choose another port: "
                             f"peakanno-agent web --port {port + 1}")
    app = create_app(runs_dir, rscript, model)
    url = f"http://{host}:{port}/"
    print(f"PeakAnnoAgent web: {url}  (runs in {Path(runs_dir).resolve()})")
    if open_browser:
        import webbrowser
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    app.run(host=host, port=port, debug=debug, threaded=True, use_reloader=False)
