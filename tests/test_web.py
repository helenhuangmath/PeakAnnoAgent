"""Web app, run index and chat, against the committed demo run (docs/demo/without_agent).
The chat uses a scripted client, so no API key or network is needed."""

import json
import shutil
import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace

from peakannoagent.chat import ChatSession, run_context
from peakannoagent.evidence import EvidenceStore
from peakannoagent.search import Run, RunIndex

from helpers import REPO, block, response

DEMO = REPO / "docs" / "demo" / "without_agent"

try:
    import flask  # noqa: F401
    HAVE_FLASK = True
except ImportError:
    HAVE_FLASK = False


class ChatClient:
    """Non-beta messages.create with scripted responses (the chat uses the stable endpoint)."""

    def __init__(self, responses):
        self.responses = list(responses)
        self.requests = []
        self.messages = SimpleNamespace(create=self._create)

    def _create(self, **kwargs):
        self.requests.append(json.loads(json.dumps(kwargs, default=lambda o: o.__dict__)))
        return self.responses.pop(0)


@unittest.skipUnless(DEMO.is_dir(), "demo run not present")
class RunIndexTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # work on a copy so chat files never land in the repository
        cls.tmp = Path(tempfile.mkdtemp())
        cls.demo = cls.tmp / "demo"
        shutil.copytree(DEMO, cls.demo)
        cls.demo_run = Run(str(cls.demo))
        cls.index = RunIndex(cls.demo_run)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_store_round_trip(self):
        store = EvidenceStore.load(str(self.demo / "tex_evidence.json"))
        self.assertEqual(store.latest("annotate_peaks").id, "E2")
        self.assertTrue(store.known_citation("E4"))
        self.assertTrue(any(k.startswith("PMID:") for k in store.references))

    def test_tables_resolve_after_moving_the_run(self):
        rows = self.demo_run.table("E2", "peak_genes")
        self.assertGreater(len(rows), 300)
        self.assertIn("nearest_symbol", rows[0])

    def test_search_groups_and_gene_lookup(self):
        res = self.index.search("Tox")
        self.assertIn("gene", res["groups"])
        self.assertIn("peak", res["groups"])
        gene = res["groups"]["gene"][0]
        self.assertEqual(gene["title"], "Tox")
        self.assertEqual(gene["evidence"], "E2")
        self.assertEqual(gene["data"]["control"], True)
        info = self.index.gene("tox")
        self.assertEqual(info["gene"], "Tox")
        self.assertEqual(len(info["peaks"]), 3)
        self.assertIsNone(self.index.gene("NotAGene"))
        self.assertEqual(self.index.search("")["total"], 0)
        self.assertEqual(set(self.index.search("Tox", kinds=["motif"])["groups"]) - {"motif"}, set())

    def test_regulators_papers_datasets_indexed(self):
        counts = self.index.counts()
        kinds = ["regulator", "paper", "dataset", "pathway", "module", "expression", "control"]
        if any(self.demo.glob("results/*_motifs_*.tsv")):          # the demo has motifs only where HOMER ran
            kinds.append("motif")
        for kind in kinds:
            self.assertGreater(counts.get(kind, 0), 0, kind)
        self.assertIn("regulator", self.index.search("Bcl11b")["groups"])

    def test_chat_loop_runs_tools_and_checks_citations(self):
        client = ChatClient([
            response([block(type="tool_use", id="t1", name="gene_lookup", input={"gene": "Tox"}),
                      block(type="tool_use", id="t2", name="inspect_table",
                            input={"evidence_id": "E4", "sort_by": "padj", "n": 2})]),     # E4: pathway enrichment
            response([block(type="text", text="**Tox** has 3 peaks [E2]; see also [E99].")], stop="end_turn"),
        ])
        session = ChatSession(self.demo_run, self.index, client=client, offline=True)
        session.clear()
        entry = session.ask("Is Tox bound?")
        self.assertEqual([c["tool"] for c in entry["tool_calls"]], ["gene_lookup", "inspect_table"])
        self.assertTrue(all(c["ok"] for c in entry["tool_calls"]))
        self.assertEqual(entry["unverified_citations"], ["E99"])
        self.assertIsNone(entry["error"])
        req = client.requests[0]
        self.assertEqual(req["thinking"], {"type": "adaptive"})
        self.assertNotIn("search_literature", [t["name"] for t in req["tools"]])   # offline run
        self.assertIn("## Results (evidence ids)", req["system"][0]["text"])
        results = client.requests[1]["messages"][-1]["content"]
        self.assertEqual([r["tool_use_id"] for r in results], ["t1", "t2"])
        self.assertIn('"Tox"', results[0]["content"])
        # persisted and reloadable
        self.assertTrue(session.path.is_file())
        again = ChatSession(self.demo_run, self.index, client=client, offline=True)
        self.assertEqual(len(again.transcript), 2)
        self.assertEqual(len(again.messages), 4)
        again.clear()
        self.assertFalse(session.path.is_file())

    def test_chat_api_error_leaves_history_consistent(self):
        class Boom:
            messages = SimpleNamespace(create=lambda **kw: (_ for _ in ()).throw(RuntimeError("down")))
        session = ChatSession(self.demo_run, self.index, client=Boom(), offline=True)
        session.clear()
        entry = session.ask("hello")
        self.assertIn("RuntimeError", entry["error"])
        self.assertEqual(session.messages, [])
        session.clear()

    def test_run_context_is_bounded(self):
        text = run_context(self.demo_run, max_chars=5000)
        self.assertLessEqual(len(text), 5100)


@unittest.skipUnless(DEMO.is_dir() and HAVE_FLASK, "demo run or Flask not present")
class WebAppTest(unittest.TestCase):
    def setUp(self):
        from peakannoagent.web import create_app
        self.tmp = Path(tempfile.mkdtemp())
        self.demo = self.tmp / "demo"
        shutil.copytree(DEMO, self.demo)
        demos = patch("peakannoagent.web.DEMO_RESULTS", [("demo-without-agent", self.demo, "Test demo")])
        demos.start()
        self.addCleanup(demos.stop)
        self.client_script = ChatClient([
            response([block(type="text", text="Promoter peaks dominate [E2].")], stop="end_turn")])
        self.app = create_app(str(self.tmp / "runs"), chat_client=self.client_script)
        self.c = self.app.test_client()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def open_demo(self) -> str:
        r = self.c.post("/runs/open", data={"path": str(self.demo)})
        self.assertEqual(r.status_code, 302)
        return r.headers["Location"].rstrip("/").split("/")[-1]

    def test_index_and_open_existing_run(self):
        self.assertEqual(self.c.get("/").status_code, 200)
        jid = self.open_demo()
        page = self.c.get(f"/runs/{jid}")
        self.assertEqual(page.status_code, 200)
        self.assertIn(b"tex_report.html", page.data)
        for part in (b'id="part-overview"', b'id="part-suggestions"', b'id="part-pathways"', b'id="part-chat"', b'id="E2"', b'class="tile"'):
            self.assertIn(part, page.data)
        with self.c.get(f"/runs/{jid}/files/tex_report.html") as response:
            self.assertEqual(response.status_code, 200)
        with self.c.get(f"/runs/{jid}/files/results/tex_peak_genes.tsv") as response:
            self.assertEqual(response.status_code, 200)
        self.assertEqual(self.c.get(f"/runs/{jid}/files/../../etc/passwd").status_code, 404)
        self.assertEqual(self.c.get("/runs/nope").status_code, 404)
        self.assertEqual(self.c.get(f"/api/runs/{jid}").get_json()["summary"]["n_peaks"], 351)
        # registered: a fresh app over the same runs dir finds it again
        from peakannoagent.web import create_app
        self.assertIn(jid, create_app(str(self.tmp / "runs")).manager.jobs)

    def test_search_and_gene_api(self):
        jid = self.open_demo()
        d = self.c.get(f"/api/runs/{jid}/search?q=Tox&kind=gene&kind=peak").get_json()
        self.assertEqual(set(d["groups"]), {"gene", "peak"})
        self.assertNotIn("text", d["groups"]["gene"][0])
        self.assertEqual(self.c.get(f"/api/runs/{jid}/gene/Pdcd1").get_json()["control"], True)

    def test_chat_api_renders_citations_as_report_links(self):
        jid = self.open_demo()
        r = self.c.post(f"/api/runs/{jid}/chat", json={"message": "Where are the peaks?"})
        self.assertEqual(r.status_code, 200)
        d = r.get_json()
        self.assertIn('href="#E2"', d["html"])      # opens the evidence card in the page
        self.assertEqual(len(self.c.get(f"/api/runs/{jid}/chat").get_json()["messages"]), 2)
        self.assertEqual(self.c.post(f"/api/runs/{jid}/chat", json={"message": "  "}).status_code, 400)
        self.c.delete(f"/api/runs/{jid}/chat")
        self.assertEqual(self.c.get(f"/api/runs/{jid}/chat").get_json()["messages"], [])

    def test_test_data_demos_and_example_prefill(self):
        home = self.c.get("/")
        self.assertIn(b"Test data", home.data)
        self.assertIn(b"example_tex_peaks_mm10.bed", home.data)
        with self.c.get("/examples/example_tex_peaks_mm10.bed") as response:
            self.assertEqual(response.status_code, 200)
        self.assertEqual(self.c.get("/examples/../pyproject.toml").status_code, 404)
        self.assertEqual(self.c.get("/examples/nope.bed").status_code, 404)
        self.assertIn(b"example_ctcf_k562_hg38.bed", home.data)
        self.assertIn(b'name="use_example" value="k4me3"', home.data)
        filled = self.c.get("/?example=1")                      # old link still fills the default example
        self.assertIn(b'name="example" value="tex"', filled.data)
        self.assertIn(b"CD8 T cells", filled.data)
        filled = self.c.get("/?example=ctcf")
        self.assertIn(b'name="example" value="ctcf"', filled.data)
        self.assertIn(b'value="hg38" selected', filled.data)
        self.assertIn(b"K562", filled.data)
        self.assertNotIn(b'name="example"', self.c.get("/?example=nope").data)
        ids = {j["id"] for j in self.c.get("/api/runs").get_json()}
        self.assertIn("demo-without-agent", ids)
        page = self.c.get("/runs/demo-without-agent")
        self.assertEqual(page.status_code, 200)
        self.assertIn(b"tex_report.html", page.data)

    def test_explore_api_table_api_and_zip_download(self):
        jid = self.open_demo()
        d = self.c.get(f"/api/runs/{jid}/explore/E4").get_json()
        self.assertEqual(d["evidence"], "E4")
        self.assertTrue(d["enrichment"])
        self.assertEqual(set(d["enrichment"][0]) & {"method", "group", "term", "padj", "genes"}, {"method", "group", "term", "padj", "genes"})
        self.assertIn("go_bp", d["networks"])
        self.assertTrue(d["networks"]["go_bp"]["nodes"] and d["networks"]["go_bp"]["edges"])
        self.assertEqual(self.c.get(f"/api/runs/{jid}/explore/E999").status_code, 404)
        t = self.c.get(f"/api/runs/{jid}/table/E4/enrichment").get_json()
        self.assertEqual(t["table"], "enrichment")
        self.assertIn("term", t["rows"][0])
        self.assertEqual(self.c.get(f"/api/runs/{jid}/table/E4/nope").status_code, 404)
        page = self.c.get(f"/runs/{jid}")
        self.assertIn(b'id="explore-E4"', page.data)
        self.assertIn(b"explore.js", page.data)
        self.assertIn(b"download.zip", page.data)
        z = self.c.get(f"/runs/{jid}/download.zip")
        self.assertEqual(z.status_code, 200)
        import io
        import zipfile
        names = zipfile.ZipFile(io.BytesIO(z.data)).namelist()
        self.assertTrue(any(n.endswith("tex_evidence.json") for n in names))
        self.assertTrue(any("results/" in n and n.endswith(".tsv") for n in names))

    def test_submit_without_peaks_is_rejected_and_leaves_no_job(self):
        r = self.c.post("/runs", data={"name": "x", "genome": "mm10", "mode": "deterministic"})
        self.assertEqual(r.status_code, 400)
        self.assertIn(b"required", r.data)
        self.assertEqual([p for p in (self.tmp / "runs").iterdir() if p.is_dir()], [])

    def test_one_click_example_copies_its_files_and_starts_a_job(self):
        started = []
        manager = self.app.manager
        real_start = manager.start
        manager.start = lambda job, config: started.append((job, config))
        try:
            r = self.c.post("/runs", data={"use_example": "ctcf", "mode": "deterministic"})
            self.assertEqual(r.status_code, 302)
            job, config = started[0]
            self.assertEqual(config.genome, "hg38")
            self.assertEqual(config.prefix, "ctcf")
            self.assertEqual(job.name, "example-ctcf")
            self.assertTrue((job.dir / "input" / "example_ctcf_k562_hg38.bed").is_file())
            self.assertEqual(config.biology.target, "CTCF")
            self.assertIn("RAD21", config.biology.expected_tfs)
            self.assertIsNone(config.biology.custom_gmt)
            # the default example carries its signatures and RNA-seq table; "1" is the old key
            self.c.post("/runs", data={"use_example": "1", "mode": "deterministic"})
            job, config = started[1]
            self.assertEqual((config.genome, config.prefix), ("mm10", "tex"))
            self.assertTrue(config.biology.custom_gmt and config.biology.rnaseq)
            # an edited form keeps the example peaks but takes the user's biology
            self.c.post("/runs", data={"example": "k4me3", "genome": "mm10", "mode": "deterministic",
                                       "target": "H3K4me3", "cell_type": "BMDM"})
            job, config = started[2]
            self.assertTrue((job.dir / "input" / "example_h3k4me3_bmdm_mm10.bed").is_file())
            self.assertEqual(config.biology.cell_type, "BMDM")
            self.assertEqual(self.c.post("/runs", data={"use_example": "nope", "mode": "deterministic"}).status_code, 400)
        finally:
            manager.start = real_start
            for j in [s[0] for s in started]:
                shutil.rmtree(j.dir, ignore_errors=True)
                manager.jobs.pop(j.id, None)


if __name__ == "__main__":
    unittest.main()
