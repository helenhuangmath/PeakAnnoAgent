"""Validate real offline mm10/hg38 analyses, recording manifests and environments.

Needs network access the first time (gene models and gene sets are downloaded and cached) and
HOMER with its mm10 and hg38 genomes for the motif checks. A partial/failed run is a failed gate.
Use --live-llm explicitly to add a credentialed model check (incurs API usage).
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

from peakannoagent.config import RunConfig, build_biology, load_biology_file
from peakannoagent.agent import execute
from peakannoagent.evidence import EvidenceStore


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    parser.add_argument("--live-llm", action="store_true")
    args = parser.parse_args()
    import peakannoagent
    examples = Path(peakannoagent.__file__).parent / "examples"
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "python-packages.txt").write_text(subprocess.check_output(
        [sys.executable, "-m", "pip", "freeze"], text=True))
    results = []
    cases = [("tex", "mm10", "example_tex_peaks_mm10.bed", "example_tex_biology.json"),
             ("ctcf", "hg38", "example_ctcf_k562_hg38.bed", "example_ctcf_k562_biology.json")]
    for name, genome, peaks, bio in cases:
        biology = build_biology(load_biology_file(str(examples / bio)))
        if name == "tex":
            biology.custom_gmt = str(examples / "example_tex_signatures.gmt")
            biology.rnaseq = str(examples / "example_tex_rnaseq_de.tsv")
        config = RunConfig(str(examples / peaks), genome, str(out / name), prefix=name,
                           biology=biology,
                           offline=True, use_llm=args.live_llm)
        try:
            manifest = execute(config)
            store = EvidenceStore.load(str(out / name / f"{name}_evidence.json"))
            qc = store.latest("peak_qc")
            checks = {"expected_peak_count": bool(qc and qc.summary.get("n_peaks") == (351 if name == "tex" else 1000))}
            if name == "ctcf":
                motif = store.latest("find_motifs")
                top = motif.summary.get("top_motifs") or [] if motif else []
                if isinstance(top, dict):
                    top = [top]
                if motif and motif.status == "ok":       # HOMER is optional; checked only when it ran
                    checks["ctcf_top_motif"] = bool(top and top[0].get("tf", "").upper() == "CTCF")
            status = manifest["status"] if all(checks.values()) else "failed"
            results.append({"case": name, "status": status, "analysis_status": manifest["status"],
                            "checks": checks, "reasons": manifest["reasons"],
                            "report": manifest["report"]})
        except Exception as err:
            results.append({"case": name, "status": "failed", "error": f"{type(err).__name__}: {err}"})
    (out / "validation.json").write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))
    return 0 if all(result["status"] == "success" for result in results) else 1


if __name__ == "__main__":
    sys.exit(main())
