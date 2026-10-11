import json
import math
import os
import shutil
import tempfile
import unittest
from pathlib import Path

from peakannoagent import analysis, public_data
from peakannoagent.config import BiologyContext, RunConfig, build_biology, normalize_genome
from peakannoagent.evidence import EvidenceStore
from peakannoagent.rbridge import RBridge
from peakannoagent.report import auto_narrative, build_report, md

from helpers import EXAMPLE_PEAKS, MOTIFS, PEAK_TABLE

HERE = Path(__file__).parent


def write(tmp, name, text):
    path = Path(tmp) / name
    path.write_text(text)
    return str(path)


class ConfigTest(unittest.TestCase):
    def test_genome_aliases(self):
        self.assertEqual(normalize_genome("mouse"), "mm10")
        self.assertEqual(normalize_genome("GRCh38"), "hg38")
        with self.assertRaises(ValueError):
            normalize_genome("dm6")

    def test_biology_parsing(self):
        bio = build_biology({"assay": "ATAC-seq", "expected_genes": "Tox, Pdcd1"}, expected_tfs="Batf;Irf4")
        self.assertEqual(bio.expected_genes, ["Tox", "Pdcd1"])
        self.assertEqual(bio.expected_tfs, ["Batf", "Irf4"])
        self.assertFalse(bio.is_tf_assay)
        with self.assertRaises(ValueError):
            build_biology({"colour": "red"})

    def test_assay_classification(self):
        self.assertTrue(BiologyContext(assay="ChIP-seq", target="Tox").is_tf_assay)
        self.assertTrue(BiologyContext(assay="CUT&Tag", target="H3K27me3").is_broad_mark)
        self.assertFalse(BiologyContext(assay="CUT&Tag", target="H3K27ac").is_tf_assay)


class EvidenceTest(unittest.TestCase):
    def test_citation_checking(self):
        store = EvidenceStore()
        store.add("peak_qc", "QC", {}, {"status": "ok", "result": {"n_peaks": 3}})
        store.add_reference("123", {"title": "x"})
        store.add_dataset("GSE5", {"accession": "GSE5"})
        text = "Fine [E1], [PMID:123], [GSE5]; bad [E2] and [PMID: 999, E1] and [SRX1]."
        self.assertEqual(store.check_citations(text), ["E2", "PMID:999", "SRX1"])

    def test_compact_stays_valid_json_when_shrunk(self):
        store = EvidenceStore()
        rows = [{"term": f"term {i}", "genes": "A;B;C" * 20} for i in range(200)]
        ev = store.add("x", "X", {}, {"status": "ok", "result": {"top_terms": rows}})
        text = store.compact(ev, limit=3000)
        self.assertLess(len(text), 3000)
        data = json.loads(text)
        self.assertIn("more (use inspect_table)", data["top_terms"][-1])

    def test_single_warning_string_is_not_split(self):
        store = EvidenceStore()
        ev = store.add("x", "X", {}, {"status": "ok", "result": {"warnings": "Only one."},
                                      "r_warnings": "R said this."})
        self.assertEqual(ev.warnings, ["Only one.", "R said this."])


class AnalysisTest(unittest.TestCase):
    def test_hypergeom_matches_closed_form(self):
        # P(X>=2) drawing 3 from 10 with 4 successes = (C(4,2)C(6,1) + C(4,3)) / C(10,3)
        expected = (math.comb(4, 2) * 6 + math.comb(4, 3)) / math.comb(10, 3)
        self.assertAlmostEqual(analysis.hypergeom_sf(2, 10, 4, 3), expected)
        self.assertEqual(analysis.hypergeom_sf(0, 10, 4, 3), 1.0)
        self.assertEqual(analysis.hypergeom_sf(5, 10, 4, 3), 0.0)

    def test_peak_gene_links_skip_far_distal(self):
        tmp = tempfile.mkdtemp()
        links = analysis.peak_gene_links(write(tmp, "p.tsv", PEAK_TABLE))
        self.assertEqual(links["Tox"]["n_peaks"], 2)
        self.assertEqual(links["Tox"]["classes"], {"promoter", "distal"})
        self.assertNotIn("Gm123", links)

    def test_rnaseq_integration(self):
        tmp = tempfile.mkdtemp()
        de = write(tmp, "de.csv", "gene,log2FoldChange,padj\nTox,2,0.001\nPdcd1,1.5,0.01\nTcf7,-2,0.001\n"
                                  "Actb,0.1,0.9\nGapdh,0,0.8\n")
        res = analysis.integrate_rnaseq(de, write(tmp, "p.tsv", PEAK_TABLE), str(Path(tmp) / "o.tsv"))
        self.assertEqual(res["n_up"], 2)
        self.assertEqual(res["up_promoter_peak"]["n_with_peak"], 2)
        self.assertEqual(res["down_any_peak"]["n_with_peak"], 1)
        rows = analysis.read_table(res["tables"]["rnaseq_integration"])
        self.assertEqual({r["gene"] for r in rows}, {"Tox", "Pdcd1", "Tcf7"})

    def test_rnaseq_rejects_unknown_columns(self):
        tmp = tempfile.mkdtemp()
        with self.assertRaises(ValueError):
            analysis.integrate_rnaseq(write(tmp, "d.tsv", "a\tb\n1\t2\n"), write(tmp, "p.tsv", PEAK_TABLE), "x")

    def test_expectations(self):
        tmp = tempfile.mkdtemp()
        res = analysis.check_expectations(["tox", "Ifng"], ["Batf", "Ctcf"], write(tmp, "p.tsv", PEAK_TABLE),
                                          write(tmp, "m.tsv", MOTIFS), [{"antigen": "CTCF"}])
        self.assertEqual(res["genes_recovered"], "1/2")
        tfs = {t["tf"]: t for t in res["tfs"]}
        self.assertTrue(tfs["Batf"]["motif_significant"])
        self.assertFalse(tfs["Ctcf"]["motif_significant"])     # not enriched (log2 = 0)
        self.assertEqual(tfs["Ctcf"]["chipatlas_rank"], 1)


class ChipAtlasTest(unittest.TestCase):
    def test_parse_and_summarize_real_result(self):
        rows = public_data.parse_enrichment_tsv((HERE / "fixtures_chipatlas_result.tsv").read_text())
        self.assertGreater(len(rows), 30)
        first = rows[0]
        self.assertEqual(first["experiment"], "SRX2399084")
        self.assertEqual(first["antigen"], "Bcl11b")
        self.assertAlmostEqual(first["log10_q"], -8.12321)
        summary = public_data.summarize_by_antigen(rows)
        antigens = [s["antigen"] for s in summary]
        self.assertIn("Ep300", antigens)
        self.assertEqual(len(antigens), len(set(antigens)))

    def test_bed_payload_sanitises(self):
        tmp = tempfile.mkdtemp()
        bed = write(tmp, "x.bed", "track name=x\nchr1\t10\t20\tpeak-1\nchrUn_GL456.1\t5\t9\nbad\n")
        self.assertEqual(public_data.bed_payload(bed), "chr1\t10\t20\nchrUn_GL456_1\t5\t9")

    def test_submit_validates_inputs(self):
        with self.assertRaises(ValueError):
            public_data.submit_enrichment("x.bed", "mm39")


class ReportTest(unittest.TestCase):
    def test_markdown_and_citation_links(self):
        store = EvidenceStore()
        store.add("peak_qc", "QC", {}, {"status": "ok", "result": {}})
        html = md("**Bold** claim [E1].\n\n- one [PMID:42]\n- two <script>", store)
        self.assertIn("<strong>Bold</strong>", html)
        self.assertIn('href="#E1"', html)
        self.assertIn('class="cite unverified"', html)
        self.assertIn("&lt;script&gt;", html)
        self.assertIn("<ul><li>", html)

    def test_report_without_llm(self):
        tmp = tempfile.mkdtemp()
        cfg = RunConfig(peaks=str(EXAMPLE_PEAKS), genome="mm10", out_dir=tmp, use_llm=False)
        store = EvidenceStore()
        store.add("peak_qc", "Peak QC", {}, {"status": "ok", "result": {"n_peaks": 351, "n_clean_peaks": 350,
                                                                        "width": {"median": 350}}})
        store.add("great_enrichment", "GREAT", {}, {"status": "error", "error": "rGREAT missing"})
        out = build_report(store, cfg, None, str(Path(tmp) / "r.html"))
        html = Path(out["path"]).read_text()
        self.assertIn("Automated summary", html)
        self.assertIn("rGREAT missing", html)
        self.assertIn("igv.min.js", html)
        self.assertEqual(out["unverified_citations"], [])
        self.assertIn("351 peaks", auto_narrative(store)["executive_summary"])

    def test_auto_summary_lists_pathway_modules_and_hubs(self):
        store = EvidenceStore()
        store.add("enrich_pathways", "Enrichment", {}, {"status": "ok", "result": {
            "top_terms": [],
            "pathway_modules": {"go_bp": [{"module": 1, "lead_term": "T cell activation", "n_terms": 6},
                                          {"module": 2, "lead_term": "cytokine production", "n_terms": 3}]},
            "hub_genes": {"go_bp": [{"gene": "Batf", "n_pathways": 4}]}}})
        text = auto_narrative(store)["executive_summary"]
        self.assertIn("M1 T cell activation (6); M2 cytokine production (3)", text)
        self.assertIn("hub genes shared across top pathways: Batf", text)


@unittest.skipUnless(shutil.which(os.environ.get("PEAKANNOAGENT_RSCRIPT", "Rscript")), "Rscript not found")
class RBridgeIntegrationTest(unittest.TestCase):
    def test_peak_qc_through_r(self):
        tmp = tempfile.mkdtemp()
        bridge = RBridge(os.environ.get("PEAKANNOAGENT_RSCRIPT", "Rscript"))
        out = bridge.run("peak_qc", {"peaks": str(EXAMPLE_PEAKS), "out_dir": tmp, "prefix": "t"})
        if out["status"] != "ok" and "data.table" in out.get("error", ""):
            self.skipTest("data.table not installed")
        self.assertEqual(out["status"], "ok", out.get("error"))
        self.assertEqual(out["result"]["n_peaks"], 351)
        self.assertTrue(Path(out["result"]["figures"]["width"]).is_file())

    def test_unknown_tool_is_reported(self):
        out = RBridge(os.environ.get("PEAKANNOAGENT_RSCRIPT", "Rscript")).run("nope", {})
        self.assertEqual(out["status"], "error")
        self.assertIn("Unknown tool", out["error"])


if __name__ == "__main__":
    unittest.main()
