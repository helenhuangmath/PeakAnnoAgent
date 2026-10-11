# R command-line reference

`scripts/peak_agent.R` runs each analysis step of the R engine on its own, without the agent. Every
command writes TSV tables plus Nature-style PDF figures to `--out`.

```bash
Rscript scripts/peak_agent.R --help        # every command and flag
```

Unknown `--methods` values fail immediately and print the list of valid methods.

## Input formats

Peak BED file. Either use a header, or give the first three columns as chr/start/end:

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

For enrichment, use the genes from every peak tested as the gene universe whenever possible. This makes
`clusterProfiler` results more defensible than API methods that use broad default backgrounds.

## Peak QC and motifs

```bash
Rscript scripts/peak_agent.R qc --peaks peaks.bed --blacklist mm10-blacklist.v2.bed --out results/qc
Rscript scripts/peak_agent.R motifs --peaks peaks.bed --species mm10 --backend monalisa --out results/motifs
```

## Peak annotation

ChIPseeker:

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks examples/data/example_peaks.bed --out results/example \
  --methods chipseeker --species mm10 --prefix union_peaks
```

Several annotation methods at once:

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks examples/data/example_peaks.bed --out results/example_annotation \
  --methods chipseeker,nearest,chipeakanno,annotatr,rgreat,homer,webexport \
  --species mm10 --prefix example
```

HOMER (`annotatePeaks.pl` must be on `PATH`):

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks path/to/peaks.bed --out results/homer --methods homer --species mm10 --prefix union_peaks
```

rGREAT peak-to-gene association:

```bash
Rscript scripts/peak_agent.R annotate \
  --peaks path/to/peaks.bed --out results/great --methods rgreat --species mm10 \
  --great-rule oneClosest --great-distance 2000 --prefix union_peaks
```

## Pathway enrichment

Local GO enrichment with `clusterProfiler` and an explicit universe:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes path/to/KClusterGroups.tsv --universe path/to/KCluster_universe.tsv \
  --out results/kmeans_pathway --methods clusterprofiler --species mm10 \
  --cluster-col Cluster --gene-col SYMBOL --ont BP --top 10 --prefix WT_gain_kmeans
```

Several pathway methods from the same gene table:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes examples/data/example_cluster_genes.tsv --out results/example_pathways \
  --methods clusterprofiler,kegg,reactome,msigdb,enrichr,gprofiler,metascape,export \
  --species mm10 --cluster-col Cluster --gene-col SYMBOL --prefix example
```

Enrichr only:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes examples/data/example_cluster_genes.tsv --out results/example \
  --methods enrichr --species mm10 --cluster-col Cluster --gene-col SYMBOL --prefix example
```

Metascape-ready grouped gene table:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes examples/data/example_cluster_genes.tsv --out results/example \
  --methods metascape --cluster-col Cluster --gene-col SYMBOL --prefix example
```

DAVID needs a registered e-mail:

```bash
export DAVID_EMAIL="your_email@example.com"
```

## Figures from existing tables

These need no external services:

```bash
Rscript scripts/peak_agent.R figures \
  --annotation examples/data/example_annotation.tsv \
  --enrichment examples/data/example_enrichment.tsv \
  --out results/example_figures --prefix example
```

## Recommended workflow

1. Generate or collect the final peak set as BED.
2. Annotate the peaks. Optionally validate peak–gene links with `rGREAT` or HOMER.
3. Build gene lists for each cluster, from the annotation table or from k-means output.
4. Run `clusterProfiler` with a real tested-gene universe.
5. Run Enrichr, DAVID or Metascape as complementary pathway discovery.
6. Treat the generated TSV files as the source of truth for figures.

## Notes

- `--species` accepts `mm10`, `mm39`, `hg38`, `hg19`, `mouse` and `human`.
- Override the annotation databases with `--txdb` and `--orgdb`.
- Keep large raw data outside the repository. `.gitignore` excludes BAM, bigWig, FASTQ, compressed BED
  and result files.
