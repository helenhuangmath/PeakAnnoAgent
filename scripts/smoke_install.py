"""Check an installed package with no source checkout on the import path: the bundled examples,
the web app, and a real annotation with the engine."""

import os
import subprocess
import sys
import tempfile

CODE = r'''
import json, sys, tempfile
from pathlib import Path
import peakannoagent
from peakannoagent import __version__, build_reference, read_peak_file, annotate_peaks
from peakannoagent.cli import main
assert __version__ == "1.0.2", __version__
package = Path(peakannoagent.__file__).parent
assert "site-packages" in str(package) or "dist-packages" in str(package), package
examples = package / "examples"
assert (examples / "example_tex_peaks_mm10.bed").is_file(), examples
assert (package / "templates" / "index.html").is_file() and (package / "static" / "app.js").is_file()

gtf = Path(tempfile.mkdtemp()) / "toy.gtf"
attrs = 'gene_id "G1"; transcript_id "T1"; gene_name "Toy"; gene_type "protein_coding";'
gtf.write_text("".join("\t".join(["chr1", "toy", kind, str(s), str(e), ".", "+", ".", attrs]) + "\n"
                       for kind, s, e in (("exon", 51980000, 51981000), ("exon", 51990000, 51992000),
                                          ("CDS", 51980500, 51981000), ("CDS", 51990000, 51991000))))
out = Path(tempfile.mkdtemp())
code = main(["annotate", "--peaks", str(examples / "example_tex_peaks_mm10.bed"), "--genome", "mm10",
             "--gtf", str(gtf), "--out", str(out), "--prefix", "smoke", "--no-download"])
assert code == 0, code
summary = json.loads((out / "smoke_annotation.json").read_text())
assert summary["n_peaks"] == 351 and len(summary["figures"]) == 8, summary
assert all(Path(p).is_file() for p in summary["figures"].values())

try:
    from peakannoagent.web import create_app, example_paths
except ImportError:
    print("Flask not installed: web app not checked")
else:
    assert example_paths("tex")
    with create_app("runs").test_client() as client:
        response = client.get("/")
        assert response.status_code == 200, response.data[:200]
print("Installed 1.0.2 package: examples, templates, annotation engine and figures available")
'''

with tempfile.TemporaryDirectory(prefix="peakanno-wheel-") as directory:
    env = dict(os.environ)
    env.pop("PYTHONPATH", None)
    env["PEAKANNOAGENT_CACHE"] = os.path.join(directory, "cache")
    env["MPLCONFIGDIR"] = os.path.join(directory, "mpl")
    subprocess.run([sys.executable, "-I", "-c", CODE], cwd=directory, env=env, check=True)
