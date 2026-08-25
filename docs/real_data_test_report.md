# Real-Data Test Report

End-to-end run of `scripts/peak_agent.R` on published ChIP-seq peak sets, used to
check that the `annotate`, `enrich`, and `figures` commands behave on real input
rather than only on the small tables in `examples/data/`.

## Test data

Five MACS peak sets distributed with `ChIPseeker` (`inst/extdata/GEO_sample_data`),
all human hg19:

| Label      | GEO sample  | Assay                              | Peaks |
| ---------- | ----------- | ---------------------------------- | ----- |
| `CBX6`     | GSM1295076  | CBX6 ChIP (Polycomb CBX protein)   | 1331  |
| `CBX7`     | GSM1295077  | CBX7 ChIP (Polycomb CBX protein)   | 1663  |
| `AR_0M`    | GSM1174480  | Androgen receptor ChIP, vehicle    | 812   |
| `AR_1nM`   | GSM1174481  | Androgen receptor ChIP, 1 nM DHT   | 2296  |
| `AR_100nM` | GSM1174482  | Androgen receptor ChIP, 100 nM DHT | 1359  |

The two groups are a useful contrast: Polycomb CBX proteins are promoter/CpG-island
binders, androgen receptor is a distal enhancer binder. A correct annotation
pipeline has to separate them.

## Annotation databases

`TxDb.Hsapiens.UCSC.hg19.knownGene` could not be installed in the test environment,
so annotation ran against a TxDb built locally from the UCSC hg19 refGene GTF
(AWS iGenomes, `Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf`), with `gene_id`
remapped from symbol to Entrez via `org.Hs.eg.db` 3.18.0 so that ChIPseeker's
`annoDb` lookup resolves. The resulting TxDb has 24,643 genes and 52,749
transcripts, and was installed as `TxDb.Hsapiens.UCSC.hg19.refGene` and passed with
`--txdb`. Gene-level results are therefore refGene-based rather than knownGene-based.

## Commands

```bash
# annotate (per sample)
Rscript scripts/peak_agent.R annotate \
  --peaks GSM1295077_CBX7_BF_ChipSeq_mergedReps_peaks.bed \
  --out results/CBX7 \
  --methods chipseeker,nearest,chipeakanno,webexport \
  --species human \
  --txdb TxDb.Hsapiens.UCSC.hg19.refGene \
  --orgdb org.Hs.eg.db \
  --prefix CBX7

# enrich (genes with a peak within 10 kb of a TSS; universe = every gene
# assigned to any tested peak, as recommended in the README)
Rscript scripts/peak_agent.R enrich \
  --genes peak_gene_clusters.tsv \
  --universe peak_gene_universe.tsv \
  --out results/enrich \
  --methods clusterprofiler,msigdb,metascape,export \
  --species human --cluster-col Cluster --gene-col SYMBOL \
  --ont BP --top 8 --prefix CBX_AR

# figures
Rscript scripts/peak_agent.R figures \
  --annotation results/CBX7/CBX7_chipseeker_annotation.tsv \
  --enrichment results/enrich/CBX_AR_clusterprofiler_GO_BP.tsv \
  --out results/figures --prefix CBX7_real --padj 0.05 --top 8
```

## Results

Genomic distribution (`*_chipseeker_annotation_distribution.tsv`, % of peaks):

| Feature           | CBX6 | CBX7 | AR_0M | AR_1nM | AR_100nM |
| ----------------- | ---- | ---- | ----- | ------ | -------- |
| Promoter          | 47.3 | 43.1 | 1.6   | 2.8    | 1.8      |
| Distal Intergenic | 31.4 | 35.0 | 55.5  | 49.1   | 47.8     |
| Intron            | 11.7 | 12.2 | 40.4  | 43.2   | 46.0     |
| Exon              | 3.9  | 4.0  | 1.5   | 2.8    | 2.2      |
| 3' UTR            | 3.9  | 3.7  | 0.5   | 1.2    | 1.5      |
| Downstream        | 1.3  | 1.3  | 0.5   | 0.8    | 0.7      |
| 5' UTR            | 0.6  | 0.7  | -     | 0.1    | -        |

The promoter fraction separates the two assay types by roughly 25-fold, which is
what these factors are expected to do.

GO BP enrichment (`clusterProfiler`, BH-adjusted, peak-derived universe) returned
128 significant terms for CBX6 and 374 for CBX7, led by transcription regulation,
`pattern specification process` (p.adj 5.3e-16), `embryo development` (1.7e-14),
and `central nervous system development` (1.1e-13) — the canonical Polycomb target
program. The AR gene sets are small (122-363 genes) and returned no significant GO
BP terms against this universe.

MSigDB hallmark enrichment returned `HALLMARK_XENOBIOTIC_METABOLISM` (p.adj 0.011)
and `HALLMARK_ANDROGEN_RESPONSE` (p.adj 0.032) for `AR_1nM`, recovering the expected
positive control for the androgen-treated sample.

## Defects found and fixed

1. **`read_gene_table()` silently mis-assigned genes to clusters.** It called
   `split(unique(dt[[gene_col]]), dt[[cluster_col]])`; whenever a gene appears in
   more than one cluster the deduplicated gene vector is shorter than the cluster
   vector, so `split()` recycled it. On the real gene table (2281 rows, 1595 unique
   genes) three of five clusters came back empty or near-empty and a fourth kept
   only 307 of its 842 genes — with nothing but an R recycling warning. Fixed by
   splitting the full gene column and deduplicating per group, which the following
   `lapply(groups, unique)` already did. Covered by a regression test in
   `tests/testthat/test-io.R`.

2. **`run_msigdb_enrichment()` always failed.** `compareCluster(fun = "enricher")`
   resolves unrecognised function names with `get(fun, envir = .GlobalEnv)`;
   because the package only calls `requireNamespace("clusterProfiler")` and never
   attaches it, `enricher` is not on the search path and the run aborted with
   `object 'enricher' not found`. Passing the function itself
   (`clusterProfiler::enricher`) fixes it. `enrichGO` and `enrichKEGG` are
   special-cased inside `compareCluster` and were unaffected.

3. **Bar-plot percentage labels were clipped.** In `plot_annotation_summary_base()`,
   slices under 8% were labelled at `y_mid + 0.55`, outside the plot region: on the
   real CBX7 summary the four small labels were cut off at the top of the page and
   overlapped each other. Their percentages now go in the legend.

4. **ChIPseeker pie title was clipped.** `ChIPseeker::plotAnnoPie()` sets
   `par(mai = c(0, 0, 0, 0))`, so the `main =` title passed by
   `run_chipseeker_annotation()` had no margin and was cut in half. The title is now
   written into a reserved outer margin with `mtext()`.

## Not covered

`rgreat`, `homer`, `annotatr`, `enrichr`, `gprofiler`, `webgestalt`, and `david` all
need network access to external services (or, for HOMER, a separate installation)
and could not be exercised in this environment. `reactome` was not run because
`reactome.db` was unavailable. `kegg` needs `rest.kegg.jp`.
