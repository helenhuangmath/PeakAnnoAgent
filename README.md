# Peak Annotation Agent

A clean R command-line framework for genomic peak annotation, gene/pathway enrichment, and manuscript-ready figures.

The goal is to keep each method easy to swap:

- peak annotation: `ChIPseeker`, nearest `TxDb`, `ChIPpeakAnno`, `annotatr`, `rGREAT`, `HOMER`
- pathway analysis: `clusterProfiler` GO/KEGG, `ReactomePA`, MSigDB, `Enrichr`, g:Profiler, WebGestaltR, DAVID, Metascape-ready exports
- figures: annotation bar plots, pie charts, TSS distance plots, and enrichment dot plots

The reusable code lives in `R/`, and the main runner is `scripts/peak_agent.R`. Project-specific one-off scripts have been replaced by configurable method wrappers.

## Best Way to Present This as an Agent

Use an R package plus command-line agent as the foundation. This is the most reproducible form for HPC/server runs and manuscript methods because it preserves exact inputs, methods, package versions, tables, and PDFs.

Add a Shiny web interface later for easy file upload and method selection. Build a REST API only if another app needs to submit jobs programmatically. See [docs/agent_design.md](docs/agent_design.md) for the suggested layered design.

## Repository Layout

```text
R/
  annotation.R      # ChIPseeker, rGREAT, HOMER wrappers
  enrichment.R      # clusterProfiler, Enrichr, DAVID, Metascape export
  figures.R         # reusable publication figure builders
  utils.R           # input parsing and shared helpers
scripts/
  peak_agent.R      # command-line runner
examples/data/
  example_peaks.bed
  example_cluster_genes.tsv
  example_annotation.tsv
  example_enrichment.tsv
docs/
  methods.md
  agent_design.md
tests/
  run_tests.R
  testthat/
```

## Install R Packages

Minimum packages:

```r
install.packages(c("data.table", "ggplot2"))
```

Or use the repository installer:

```bash
Rscript scripts/install_dependencies.R
```

Bioconductor methods:

```r
install.packages("BiocManager")
BiocManager::install(c(
  "ChIPseeker",
  "ChIPpeakAnno",
  "clusterProfiler",
  "enrichplot",
  "AnnotationDbi",
  "GenomicRanges",
  "GenomicFeatures",
  "IRanges",
  "annotatr",
  "ReactomePA",
  "WebGestaltR",
  "org.Mm.eg.db",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
))
```

API-backed methods:

```r
install.packages(c("enrichR", "gprofiler2", "msigdbr"))
BiocManager::install("rGREAT")
```

HOMER is not an R package. Install HOMER separately and make sure `annotatePeaks.pl` is available in `PATH`.

## Input Formats

Peak BED file:

```text
Chr   Start   End   PeakID
chr1  3671498 3671800 peak_001
```

Cluster/gene table:

```text
Cluster   SYMBOL
Cluster1  Sox2
Cluster1  Pax6
Cluster2  Il6
```

For enrichment, use the gene universe from all peaks that were tested whenever possible. That makes `clusterProfiler` results more defensible for publication than API methods that use broad default backgrounds.

## Hands-On Examples

Show available commands:

```bash
Rscript scripts/peak_agent.R --help
```

Export a Metascape-ready grouped gene table:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes examples/data/example_cluster_genes.tsv \
  --out results/example \
  --methods metascape \
  --cluster-col Cluster \
  --gene-col SYMBOL \
  --prefix example
```

Create example annotation and pathway figures without external services:

```bash
Rscript scripts/peak_agent.R figures \
  --annotation examples/data/example_annotation.tsv \
  --enrichment examples/data/example_enrichment.tsv \
  --out results/example_figures \
  --prefix example
```

Run Enrichr on clustered genes:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes examples/data/example_cluster_genes.tsv \
  --out results/example \
  --methods enrichr \
  --species mm10 \
  --cluster-col Cluster \
  --gene-col SYMBOL \
  --prefix example
```

Run local GO enrichment with `clusterProfiler`:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes path/to/KClusterGroups.tsv \
  --universe path/to/KCluster_universe.tsv \
  --out results/kmeans_pathway \
  --methods clusterprofiler \
  --species mm10 \
  --cluster-col Cluster \
  --gene-col SYMBOL \
  --ont BP \
  --top 10 \
  --prefix WT_gain_kmeans
```

Run several pathway methods from the same gene table:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes examples/data/example_cluster_genes.tsv \
  --out results/example_pathways \
  --methods clusterprofiler,kegg,reactome,msigdb,enrichr,gprofiler,metascape,export \
  --species mm10 \
  --cluster-col Cluster \
  --gene-col SYMBOL \
  --prefix example
```

Annotate peaks with ChIPseeker:

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks examples/data/example_peaks.bed \
  --out results/example \
  --methods chipseeker \
  --species mm10 \
  --prefix union_peaks
```

Run multiple annotation methods:

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks examples/data/example_peaks.bed \
  --out results/example_annotation \
  --methods chipseeker,nearest,chipeakanno,annotatr,rgreat,homer,webexport \
  --species mm10 \
  --prefix example
```

Annotate peaks with HOMER:

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks path/to/peaks.bed \
  --out results/homer \
  --methods homer \
  --species mm10 \
  --prefix union_peaks
```

Submit peak-to-gene association with rGREAT:

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks path/to/peaks.bed \
  --out results/great \
  --methods rgreat \
  --species mm10 \
  --great-rule oneClosest \
  --great-distance 2000 \
  --prefix union_peaks
```

## Testing

Run the offline unit tests (parsing, species mapping, method validation, and the
export paths that need no Bioconductor packages):

```bash
Rscript tests/run_tests.R
```

Unknown `--methods` values now fail fast with the list of valid methods instead of
running nothing, and `Rscript scripts/peak_agent.R --help` documents every flag.

## Recommended Workflow

1. Generate or collect the final peak set as BED.
2. Annotate peaks with `ChIPseeker` and optionally validate peak-gene links with `rGREAT` or HOMER.
3. Build cluster-specific gene lists from the peak annotation table or your k-means output.
4. Run `clusterProfiler` with a real tested-gene universe.
5. Run Enrichr, DAVID, or Metascape as complementary pathway discovery.
6. Use the generated TSV files as the source of truth for figures.

## Notes for Real Projects

- Use `--species mm10` or `--species hg38`; shortcuts `mouse` and `human` also work.
- Override annotation databases with `--txdb` and `--orgdb` when needed.
- Set `DAVID_EMAIL` before using DAVID:

```bash
export DAVID_EMAIL="your_email@example.com"
```

- Keep raw data outside the repository when files are large. The `.gitignore` excludes common BAM, bigWig, FASTQ, compressed BED, and result files.
