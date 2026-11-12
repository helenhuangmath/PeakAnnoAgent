# Tutorial 5: The R engine

Every analysis step of the agent is an ordinary R function. You can use the R engine on its own, in a
script or in R Markdown, without Python or an API key.

## Run the agent tools from R

`run_agent_tool()` runs the same tools the agent calls. It takes the same arguments and returns the same
results:

```r
library(PeakAnnoAgent)

anno <- run_agent_tool("annotate_peaks", list(
  peaks = "examples/data/example_tex_peaks_mm10.bed", genome = "mm10",
  out_dir = "results/r", prefix = "tex"))

enr <- run_agent_tool("enrich_pathways", list(
  gene_sets = anno$gene_sets, universe = anno$universe, genome = "mm10",
  methods = c("go_bp", "msigdb_H"), out_dir = "results/r", prefix = "tex"))

mot <- run_agent_tool("find_motifs", list(
  peaks = "examples/data/example_tex_peaks_mm10.bed", genome = "mm10",
  out_dir = "results/r", prefix = "tex"))

enr$pathway_modules$go_bp     # modules of the enrichment map
enr$hub_genes$go_bp           # genes shared by several top pathways
```

Each call writes its tables and figures to `out_dir` and returns a list with the summary, the file
paths and any warnings.

## Pathway networks from any enrichment table

`run_pathway_networks()` draws the enrichment map and the gene–pathway network for each method in a
standard enrichment table (`method`, `group`, `term`, `padj`, `count`, `genes`):

```r
tab <- data.table::fread("results/r/tex_enrich1_enrichment.tsv")
net <- run_pathway_networks(tab, out_dir = "results/r", prefix = "tex")
net$figures     # PNG previews (PDFs are written next to them)
net$modules     # module tables per method
```

The table can also come from another tool, such as clusterProfiler, as long as it has these columns.
See [Tutorial 4](pathway_networks.md) for the settings.

## Run single steps from the command line

`scripts/peak_agent.R` runs each step as a standalone command:

```bash
Rscript scripts/peak_agent.R --help
Rscript scripts/peak_agent.R annotate --peaks peaks.bed --species mm10 --methods chipseeker --out results/anno
Rscript scripts/peak_agent.R enrich --genes clusters.tsv --universe universe.tsv \
  --methods clusterprofiler --species mm10 --out results/pathways
```

The [R command-line reference](../reference/r_cli.md) lists every command.

## Vignette

The package vignette (`vignettes/peak-annotation-workflow.Rmd`) walks through the complete R workflow,
from peaks to figures.
