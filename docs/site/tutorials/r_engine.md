# Tutorial 5: Python and R

Every analysis step is an ordinary Python function, so you can use the engine in a script or a
notebook without the agent, and call the same program from R.

## Annotate in Python

```python
from peakannoagent import build_reference, read_peak_file, annotate_peaks, run_annotation

ref = build_reference("mm10", "gencode")          # downloaded and built once, then cached
print(ref)
# <Reference mm10, GENCODE vM25 basic: 55,401 genes (21,859 protein-coding), 81,540 transcripts>

peaks = read_peak_file("peakannoagent/examples/example_tex_peaks_mm10.bed")
table = annotate_peaks(peaks, ref, promoter_upstream=2000, promoter_downstream=500)
```

`table` is a dictionary of equal-length columns, one entry per peak. It goes straight into pandas:

```python
import pandas as pd
df = pd.DataFrame(table)
df[["PeakID", "feature", "nearest_symbol", "distance_to_tss", "region_class"]].head()
df["feature"].value_counts(normalize=True)
```

`run_annotation()` does the same and also writes the tables, the gene sets and the figures:

```python
res = run_annotation(peaks, "results/py", "tex", genome="mm10", reference=ref)
res["feature_pct"]        # {'Promoter (<=1kb)': 20.23, ...}
res["tables"]             # paths of the TSV files
res["figures"]            # paths of the PDF figures
```

## Pathway enrichment in Python

The gene sets and the background written by `run_annotation()` feed the enrichment:

```python
from peakannoagent import run_enrichment

enr = run_enrichment(res["gene_sets"], "results/py", "tex_enrich", genome="mm10",
                     methods=["go_bp", "hallmark", "kegg"], universe=res["universe"])
enr["top_terms"][:3]
enr["pathway_modules"]["go_bp"]     # modules of the enrichment map
enr["hub_genes"]["go_bp"]           # genes shared by several top pathways
```

Methods are `go_bp`, `go_mf`, `go_cc`, `hallmark`, `reactome`, `wikipathways`, `kegg`,
`tf_targets`, `immunesigdb`, `cell_types` and `custom` (your own GMT, passed as `custom_gmt`).

## The same tools the agent calls

`Engine.run()` runs a tool by name with a dictionary of settings and returns a status and a result,
which is exactly what the agent sees:

```python
from peakannoagent import Engine

out = Engine().run("peak_qc", {"peaks": "peaks.bed", "out_dir": "results/py", "prefix": "tex"})
out["status"], out["result"]["n_clean_peaks"]
```

## Call it from R

The R package in `r/PeakAnnoAgent` runs the Python program and returns its tables as data frames.

```r
library(PeakAnnoAgent)

res <- annotate_peaks("peaks.bed", genome = "hg38", reference = "gencode")
head(res$annotation)
res$feature_distribution
res$figures                       # paths of the PDF figures

full <- run_analysis("peaks.bed", genome = "hg38", out_dir = "results/k562",
  assay = "TF ChIP-seq", target = "CTCF")
browseURL(full$report)
```

`peaks` can also be a data frame with `Chr`, `Start` and `End`. See
[R interface](../api/r_interface.md) for every argument.
