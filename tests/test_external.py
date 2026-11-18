"""Companion tools (SimpNetworkViz, PathwayAgent): input conversion, discovery, and the tool
wrappers' skip behaviour. The real programs run only when they are installed."""

import os
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from peakannoagent import external
from peakannoagent.config import BiologyContext, RunConfig
from peakannoagent.tools import RunContext, run_tool

from helpers import EXAMPLE_PEAKS, REPO, FakeBridge

DEMO_RESULTS = REPO / "docs" / "demo" / "without_agent" / "results"
ENRICHMENT = DEMO_RESULTS / "tex_enrich1_enrichment.tsv"
GENE_SETS = DEMO_RESULTS / "tex_peak_gene_sets.tsv"
UNIVERSE = DEMO_RESULTS / "tex_gene_universe.tsv"


@unittest.skipUnless(ENRICHMENT.is_file(), "demo tables not present")
class ConversionTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_enrichment_table_to_clusterprofiler_layout(self):
        out = self.tmp / "e.csv"
        n = external.enrichment_to_clusterprofiler(str(ENRICHMENT), "go_bp", "all", str(out))
        self.assertGreater(n, 100)
        lines = out.read_text().splitlines()
        self.assertEqual(lines[0], "ID,Description,pvalue,p.adjust,geneID,Count")
        genes = [line.split(",")[4] for line in lines[1:]]
        self.assertTrue(any("/" in g for g in genes))          # clusterProfiler separates genes with '/'
        self.assertFalse(any(";" in g for g in genes))
        methods = external.enrichment_methods(str(ENRICHMENT))
        self.assertIn("go_bp", methods)
        self.assertIn("all", methods["go_bp"])

    def test_gene_lists_one_file_per_group(self):
        files = external.write_gene_lists(str(GENE_SETS), str(self.tmp), ["promoter", "distal", "all", "nope"])
        self.assertEqual(set(files), {"promoter", "distal", "all"})
        genes = Path(files["promoter"]).read_text().split()
        self.assertGreater(len(genes), 50)
        self.assertEqual(len(genes), len(set(genes)))

    def test_discovery_prefers_environment(self):
        with mock.patch.dict(os.environ, {"PEAKANNOAGENT_SNV": "/opt/x/snv --quiet",
                                          "PEAKANNOAGENT_PATHWAYAGENT": "/opt/y/python -m pathwayagent"}):
            self.assertEqual(external.find_snv(), ["/opt/x/snv", "--quiet"])
            self.assertEqual(external.find_pathwayagent(), ["/opt/y/python", "-m", "pathwayagent"])


class ToolSkipTest(unittest.TestCase):
    """Without the programs, the tools are recorded as skipped (never an error) and point at the install."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        cfg = RunConfig(peaks=str(EXAMPLE_PEAKS), genome="mm10", out_dir=self.tmp, prefix="t",
                        biology=BiologyContext(), offline=True)
        self.ctx = RunContext(config=cfg, bridge=FakeBridge(self.tmp), log=lambda m: None)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_skipped_before_prerequisites(self):
        ev = run_tool(self.ctx, "network_viz", {})
        self.assertEqual(ev.status, "skipped")
        self.assertIn("enrich_pathways", ev.error)
        ev = run_tool(self.ctx, "pathway_consensus", {})
        self.assertEqual(ev.status, "skipped")
        self.assertIn("annotate_peaks", ev.error)

    def test_skipped_when_not_installed(self):
        run_tool(self.ctx, "annotate_peaks", {})
        with mock.patch.object(external, "find_pathwayagent", return_value=None):
            ev = run_tool(self.ctx, "pathway_consensus", {})
        self.assertEqual(ev.status, "skipped")
        self.assertIn("PathwayAgent", ev.error)

    def test_empty_enrichment_is_a_valid_network_skip(self):
        table = Path(self.tmp) / "empty.tsv"
        table.write_text("method\tgroup\tterm\tpadj\tgenes\n")
        self.ctx.evidence.add("enrich_pathways", "Pathways", {}, {"status": "ok", "result": {
            "tables": {"enrichment": str(table)}}})
        with mock.patch.object(external, "find_snv", return_value=["snv"]):
            ev = run_tool(self.ctx, "network_viz", {})
        self.assertEqual(ev.status, "skipped")
        self.assertIn("No enriched terms", ev.error)

    def test_too_few_significant_terms_is_a_network_skip_not_an_error(self):
        # Found on real data: small peak sets can have every term tested and none (or one) significant.
        table = Path(self.tmp) / "weak.tsv"
        table.write_text("method\tgroup\tterm\tpadj\tgenes\n"
                         "go_bp\tall\tGOBP_A\t0.04\tTox;Tcf7\n"
                         "go_bp\tall\tGOBP_B\t0.21\tTox\n"
                         "go_bp\tall\tGOBP_C\t0.90\tTcf7\n")
        self.ctx.evidence.add("enrich_pathways", "Pathways", {}, {"status": "ok", "result": {
            "tables": {"enrichment": str(table)}}})
        with mock.patch.object(external, "find_snv", return_value=["snv"]), \
                mock.patch.object(external, "run_snv", side_effect=AssertionError("must not be called")):
            ev = run_tool(self.ctx, "network_viz", {})
        self.assertEqual(ev.status, "skipped")
        self.assertIn("1 significant go_bp term ", ev.error)
        self.assertIn("at least two", ev.error)


@unittest.skipUnless(os.environ.get("PEAKANNOAGENT_LIVE_TESTS") == "1" and ENRICHMENT.is_file()
                     and external.find_snv(), "set PEAKANNOAGENT_LIVE_TESTS=1 for SimpNetworkViz integration")
class SnvLiveTest(unittest.TestCase):
    def test_networks_from_enrichment_table(self):
        tmp = Path(tempfile.mkdtemp())
        try:
            res = external.run_snv(str(ENRICHMENT), str(tmp / "snv"), "go_bp", "all", max_pathways=30,
                                   gene_network="cnet", title="test")
            self.assertEqual(res["n_pathways"], 30)
            self.assertGreaterEqual(res["n_clusters"], 1)
            self.assertTrue(res["clusters"][0]["label"])
            self.assertIn("pathway_network", res["figures"])
            self.assertTrue(Path(res["figures"]["pathway_network"]).is_file())
            self.assertIn("pathway_network", res["interactive"])
            self.assertIn("clusters", res["tables"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
