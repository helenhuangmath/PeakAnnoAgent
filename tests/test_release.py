"""Regression tests for release-blocking failures, independent of live backends."""

import json
import tempfile
import unittest
from datetime import datetime
from pathlib import Path
from unittest.mock import patch

from peakannoagent.config import BiologyContext, RunConfig
from peakannoagent.evidence import EvidenceStore
from peakannoagent.outcome import assess_run
from peakannoagent.agent import execute
from peakannoagent.cli import build_parser, cmd_run
from peakannoagent.tools import RunContext, TOOLS, run_tool
from peakannoagent.engine import Engine
from types import SimpleNamespace
from helpers import EXAMPLE_PEAKS, FakeBridge


class CompletionTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.config = RunConfig(str(EXAMPLE_PEAKS), "mm10", self.tmp.name,
                                use_llm=False, offline=True)
        self.store = EvidenceStore()

    def add(self, name, status="ok", **result):
        return self.store.add(name, name, {}, {"status": status, "result": result,
                                              "error": "backend unavailable" if status != "ok" else None})

    def complete(self):
        for name in ("peak_qc", "annotate_peaks", "enrich_pathways", "find_motifs"):
            self.add(name)

    def test_no_annotation_is_failed_even_with_qc(self):
        self.add("peak_qc")
        self.add("annotate_peaks", "error")
        result = assess_run(self.store, self.config)
        self.assertEqual(result["status"], "failed")
        self.assertIn("annotate_peaks", result["unsuccessful_steps"])

    def test_failed_motifs_is_partial_but_absent_homer_is_not(self):
        self.complete()
        self.add("find_motifs", "skipped")        # HOMER is optional: a clean skip is not a failure
        self.assertEqual(assess_run(self.store, self.config)["status"], "success")
        self.add("find_motifs", "error")
        self.assertEqual(assess_run(self.store, self.config)["status"], "partial")

    def test_optional_skip_and_empty_results_can_succeed(self):
        self.complete()
        self.add("network_viz", "skipped")
        self.assertEqual(assess_run(self.store, self.config)["status"], "success")
        self.add("find_motifs", "error")
        self.add("find_motifs", n_significant=0)
        self.assertEqual(assess_run(self.store, self.config)["status"], "success")

    def test_method_failures_are_partial(self):
        self.complete()
        self.add("enrich_pathways", failed_methods={"kegg": "network down"})
        self.assertEqual(assess_run(self.store, self.config)["status"], "partial")

    def test_broad_marks_do_not_require_motifs(self):
        self.config.biology = BiologyContext(target="H3K27me3")
        for name in ("peak_qc", "annotate_peaks", "enrich_pathways"):
            self.add(name)
        self.assertEqual(assess_run(self.store, self.config)["status"], "success")

    def test_execute_writes_failed_report_manifest_and_provenance(self):
        class BrokenAnnotation(FakeBridge):
            def run(self, tool, params):
                if tool == "annotate_peaks":
                    return {"status": "error", "error": "missing annotation package"}
                return super().run(tool, params)

        factory = lambda **kw: RunContext(bridge=BrokenAnnotation(self.tmp.name), **kw)
        with patch("peakannoagent.agent.RunContext", side_effect=factory):
            manifest = execute(self.config, log=lambda _: None)
        self.assertEqual(manifest["status"], "failed")
        self.assertEqual(manifest["version"], "1.1.0")
        self.assertEqual(len(manifest["inputs"]["peaks"]["sha256"]), 64)
        self.assertIn("Analysis status: <strong>failed</strong>", Path(manifest["report"]).read_text())
        saved = json.loads(Path(self.tmp.name, "run_manifest.json").read_text())
        self.assertEqual(saved["status"], "failed")

    def test_previews_are_embedded_in_the_report_and_then_removed(self):
        png = bytes.fromhex("89504e470d0a1a0a") + b"\0" * 8 + (300).to_bytes(4, "big") + (200).to_bytes(4, "big") + b"\0" * 20

        class WithFigure(FakeBridge):
            def run(self, tool, params):
                outcome = super().run(tool, params)
                if tool == "annotate_peaks":
                    results = Path(params["out_dir"])
                    (results / ".preview").mkdir(parents=True, exist_ok=True)
                    (results / "r_feature_pie.pdf").write_bytes(b"%PDF-1.4")
                    (results / ".preview" / "r_feature_pie.png").write_bytes(png)
                    outcome["result"]["figures"] = {"feature_pie": str(results / "r_feature_pie.pdf")}
                return outcome

        for keep in (False, True):
            with self.subTest(keep_previews=keep):
                out = Path(self.tmp.name) / f"keep{keep}"
                config = RunConfig(str(EXAMPLE_PEAKS), "mm10", str(out), use_llm=False, offline=True,
                                   keep_previews=keep)
                factory = lambda **kw: RunContext(bridge=WithFigure(str(out)), **kw)
                with patch("peakannoagent.agent.RunContext", side_effect=factory):
                    manifest = execute(config, log=lambda _: None)
                html = Path(manifest["report"]).read_text()
                self.assertEqual(html.count("data:image/png;base64,"), 1)
                self.assertIn('href="results/r_feature_pie.pdf" download', html)
                self.assertNotIn(".preview", html)
                self.assertEqual((out / "results" / ".preview").exists(), keep)
                self.assertTrue((out / "results" / "r_feature_pie.pdf").is_file())

    def test_palette_is_validated_and_passed_to_the_engine(self):
        self.assertEqual(self.config.palette, "chippeakanno")
        with self.assertRaisesRegex(ValueError, "Unknown palette"):
            RunConfig(str(EXAMPLE_PEAKS), "mm10", self.tmp.name, palette="rainbows")
        config = RunConfig(str(EXAMPLE_PEAKS), "mm10", self.tmp.name, palette="Spectral")
        params = RunContext(config, bridge=FakeBridge(self.tmp.name), log=lambda _: None).base_params()
        self.assertEqual((params["palette"], params["previews"]), ("spectral", True))
        args = build_parser().parse_args(["run", "--peaks", "p.bed", "--genome", "mm10", "--out", "o",
                                          "--palette", "viridis"])
        self.assertEqual(args.palette, "viridis")

    def test_cli_exit_codes(self):
        args = build_parser().parse_args(["run", "--peaks", str(EXAMPLE_PEAKS), "--genome", "mm10",
                                          "--out", self.tmp.name, "--no-llm"])
        for status, code in (("success", 0), ("partial", 3), ("failed", 1)):
            with self.subTest(status=status), patch("peakannoagent.agent.execute",
                    return_value={"report": "report.html", "status": status}):
                self.assertEqual(cmd_run(args), code)

    def test_engine_turns_exceptions_into_tool_errors(self):
        result = Engine().run("annotate_peaks", {"peaks": str(Path(self.tmp.name) / "missing.bed"),
                                                 "out_dir": self.tmp.name})
        self.assertEqual(result["status"], "error")
        self.assertIn("missing.bed", result["error"])
        self.assertEqual(result["tool"], "annotate_peaks")

    def test_engine_reports_an_absent_optional_backend_as_skipped(self):
        with patch("peakannoagent.engine.motifs.shutil.which", return_value=None):
            result = Engine().run("find_motifs", {"peaks": str(EXAMPLE_PEAKS), "out_dir": self.tmp.name})
        self.assertEqual(result["status"], "skipped")
        self.assertIn("HOMER", result["error"])

    def test_unexpected_tool_exceptions_are_retained_as_failed_evidence(self):
        context = RunContext(self.config, log=lambda _: None)
        with patch.object(TOOLS["peak_qc"], "handler", side_effect=RuntimeError("unexpected backend crash")):
            evidence = run_tool(context, "peak_qc", {})
        self.assertEqual(evidence.status, "error")
        self.assertIn("unexpected backend crash", evidence.error)
        self.assertEqual(assess_run(context.evidence, self.config)["status"], "failed")


try:
    from peakannoagent.web import JobManager
except ImportError:
    JobManager = None


@unittest.skipIf(JobManager is None, "Flask not installed")
class JobRegressionTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.manager = JobManager(Path(self.tmp.name))

    def test_same_name_submissions_have_isolated_directories(self):
        with patch("peakannoagent.web.datetime") as clock:
            clock.now.return_value = datetime(2026, 11, 12, 12, 0, 0)
            first = self.manager.new_job("same", "run")
            second = self.manager.new_job("same", "run")
        self.assertNotEqual(first.id, second.id)
        self.assertNotEqual(first.dir, second.dir)
        self.assertTrue((first.dir / "input").is_dir())
        self.assertIn(first.id, self.manager.jobs)
        self.assertIn(second.id, self.manager.jobs)

    def test_manifest_outcome_survives_restart(self):
        job = self.manager.new_job("partial", "run")
        (job.dir / "run_manifest.json").write_text(json.dumps({"status": "partial", "reasons": ["motifs failed"]}))
        restored = self.manager._from_dir(job.id, job.dir)
        self.assertEqual(restored.status, "partial")
        self.assertIn("motifs failed", restored.error)

    def test_manifestless_job_does_not_crash_discovery_or_become_done(self):
        job = self.manager.new_job("queued", "run")
        restored = self.manager._from_dir(job.id, job.dir)
        self.assertEqual(restored.status, "interrupted")

    def test_worker_uses_analysis_outcome(self):
        job = self.manager.new_job("failed", "run")
        with patch("peakannoagent.agent.execute", return_value={"status": "failed", "reasons": ["annotation failed"]}), \
             patch("peakannoagent.web.threading.Thread") as thread:
            self.manager.start(job, None)
            thread.call_args.kwargs["target"]()
        self.assertEqual(job.status, "failed")
        self.assertIn("annotation failed", job.error)
