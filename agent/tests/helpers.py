"""Fakes shared by the tests: an R bridge that writes canned tables, and a scripted Claude client."""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

REPO = Path(__file__).resolve().parents[2]
EXAMPLE_PEAKS = REPO / "examples" / "data" / "example_tex_peaks_mm10.bed"

PEAK_TABLE = """PeakID\tChr\tStart\tEnd\tnearest_entrez\tdistance_to_tss\tshortest_distance\tinside_feature\tnearest_symbol\tpromoter_symbols\tregion_class
p1\tchr1\t100\t400\t1\t0\t0\toverlapStart\tTox\tTox\tpromoter
p2\tchr1\t9000\t9300\t1\t8800\t8800\tupstream\tTox\tNA\tdistal
p3\tchr2\t500\t800\t2\t-100\t100\tupstream\tPdcd1\tPdcd1\tpromoter
p4\tchr3\t500\t800\t3\t40000\t40000\tdownstream\tTcf7\tNA\tdistal
p5\tchr4\t500\t800\t4\t500000\t500000\tdownstream\tGm123\tNA\tdistal_far
"""

MOTIFS = """motif\ttf\tfamily\tpct_target\tpct_background\tlog2_enrichment\tpvalue\tpadj
MA1\tNFATC1\tRel\t30\t10\t1.5\t1e-10\t1e-8
MA2\tBatf::Jun\tbZIP\t25\t12\t1.0\t1e-6\t1e-4
MA3\tCTCF\tZf\t5\t5\t0.0\t0.5\t0.9
"""


class FakeBridge:
    """Stands in for RBridge: writes small tables and returns engine-shaped outcomes."""

    def __init__(self, out_dir: str):
        self.out = Path(out_dir)
        self.calls = []

    def available(self):
        return True

    def run(self, tool: str, params: dict) -> dict:
        self.calls.append((tool, params))
        self.out.mkdir(parents=True, exist_ok=True)
        if tool == "peak_qc":
            return {"status": "ok", "result": {"n_peaks": 5, "n_clean_peaks": 5, "width": {"median": 300},
                                               "warnings": [], "clean_peaks": params["peaks"]}}
        if tool == "annotate_peaks":
            table = self.out / "peak_genes.tsv"
            table.write_text(PEAK_TABLE)
            uni = self.out / "universe.tsv"
            uni.write_text("SYMBOL\n" + "\n".join(["Tox", "Pdcd1", "Tcf7", "Gm123"] + [f"G{i}" for i in range(96)]))
            sets = self.out / "sets.tsv"
            sets.write_text("Cluster\tSYMBOL\npromoter\tTox\npromoter\tPdcd1\n")
            return {"status": "ok", "result": {
                "n_peaks": 5, "region_class_counts": {"promoter": 2, "distal": 2, "distal_far": 1},
                "example_promoter_genes": ["Tox", "Pdcd1"], "promoter_window": {"upstream": 2000, "downstream": 500},
                "tables": {"peak_genes": str(table)}, "gene_sets": str(sets), "universe": str(uni),
                "warnings": [], "versions": {"ChIPpeakAnno": "3.36.1"}}}
        if tool == "find_motifs":
            table = self.out / "motifs.tsv"
            table.write_text(MOTIFS)
            return {"status": "ok", "result": {"backend": "monaLisa", "database": "JASPAR2020", "window_bp": 200,
                                               "n_significant": 2, "top_motifs": [{"tf": "NFATC1", "padj": 1e-8}],
                                               "tables": {"motifs": str(table)}, "warnings": []}}
        if tool == "enrich_pathways":
            return {"status": "ok", "result": {"methods": params["methods"], "top_terms": [
                {"method": "go_bp", "group": "promoter", "term": "T cell activation", "padj": 1e-5, "count": 2,
                 "genes": "Tox;Pdcd1"}], "warnings": []}}
        return {"status": "error", "error": f"{tool} not available in FakeBridge"}


def block(**kw):
    return SimpleNamespace(**kw)


def response(blocks, stop="tool_use"):
    return SimpleNamespace(content=blocks, stop_reason=stop,
                           usage=SimpleNamespace(input_tokens=100, output_tokens=50, cache_read_input_tokens=0))


class ScriptedClient:
    """Returns pre-scripted responses and records every request."""

    def __init__(self, responses):
        self.responses = list(responses)
        self.requests = []
        self.beta = SimpleNamespace(messages=SimpleNamespace(create=self._create))

    def _create(self, **kwargs):
        self.requests.append(json.loads(json.dumps(kwargs, default=lambda o: o.__dict__)))
        return self.responses.pop(0)
