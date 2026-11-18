# Agent tools

These are the tools the agent calls. In agent mode, Claude chooses which ones to run and with which
settings. With `--no-llm`, a fixed plan runs them in order. Tools marked *network* contact a
public service and are skipped with `--offline`. `peakanno-agent tools` prints the same list.

Every tool runs in the analysis engine or calls a public API. Its result is stored with an evidence ID
(`E1`, `E2`, …), which the report cites.

## `peak_qc`

Peak-set QC: counts, widths, chromosomes, duplicates, non-standard contigs and blacklist overlap. Writes a cleaned peak set that later tools use automatically. Run first.

## `annotate_peaks`

Peak annotation against the run's gene models (GENCODE by default; RefSeq, UCSC or a GTF if the user chose one). For every peak: the nearest TSS and its gene, all genes whose promoter window the peak overlaps (bidirectional and alternative promoters), and one genomic feature by priority (promoter, 5' UTR, 3' UTR, exon, intron, downstream, distal intergenic). Returns the feature-type proportions (feature_pct) and region distribution, classifies peaks as promoter / distal / distal_far and builds gene sets (promoter, distal, all) plus a genome background for enrichment. Draws the feature-type pie chart, feature bar, distance-to-TSS bars, feature-overlap upset plot, region distribution, TSS profile, TSS-distance histogram and peaks over chromosomes. Takes seconds. Choose the promoter window to fit the assay (tight for promoter-bound TFs, standard 2000/500 otherwise).

| Argument | Type | Description |
|---|---|---|
| `promoter_upstream` | integer | bp upstream of TSS (default 2000) |
| `promoter_downstream` | integer | bp downstream of TSS (default 500) |
| `distal_max_distance` | integer | Distal peaks further than this from any TSS are not linked to a gene (default 100000) |
| `gene_types` | string: `auto`, `protein_coding`, `all` | Genes to annotate against: protein-coding only (default) or all, including lncRNAs and pseudogenes |

## `enrich_pathways`

Over-representation analysis (hypergeometric test, Benjamini-Hochberg) of peak-linked gene sets against the genome background. Methods: go_bp, go_mf, go_cc, hallmark, reactome, wikipathways, kegg, tf_targets, immunesigdb, cell_types, custom. The public collections come from MSigDB (go_bp, go_mf, go_cc, hallmark, reactome, wikipathways, tf_targets = GTRD TF targets, immunesigdb, cell_types) and KEGG; 'custom' uses the user's GMT signatures. Pick collections that fit the biology (e.g. immunesigdb for immune cells). Also draws pathway networks per method: an enrichment map (terms linked by shared genes, grouped into modules - returned as pathway_modules) and a gene-pathway network (hub genes shared by several top pathways - returned as hub_genes). Interpret modules rather than redundant individual terms.

| Argument | Type | Description |
|---|---|---|
| `methods` (required) | array of string |  |
| `groups` | array of string: `promoter`, `distal`, `all` |  |
| `padj` | number |  |
| `top` | integer | Top terms per group to return |

## `network_viz`

Draw EnrichmentMap/AutoAnnotate-style networks from this run's enrichment results (SimpNetworkViz): pathways linked by shared genes, clustered into labelled themes, plus a gene-pathway network; interactive HTML and Cytoscape files. Run after enrich_pathways. Pick the collection and peak group to draw (default: go_bp, all peaks).

| Argument | Type | Description |
|---|---|---|
| `method` | string | A method present in the enrichment table, e.g. go_bp, hallmark, immunesigdb, custom |
| `group` | string: `promoter`, `distal`, `all` |  |
| `max_pathways` | integer | Most significant pathways to draw (default 40) |
| `cutoff` | number | Similarity needed for an edge (default 0.6 on hit genes; raise it if the network is one dense blob, lower it if nothing connects) |
| `gene_network` | string: `cnet`, `shared`, `none` |  |

## `pathway_consensus`

Run PathwayAgent on the peak-linked gene lists: several enrichment engines (local ORA on GO/KEGG/Reactome/WikiPathways/Hallmark/DisGeNET; with network access also Enrichr, g:Profiler, Reactome, STRING, clusterProfiler, fgsea) combined into cross-tool consensus themes with tiers (ROBUST / SUPPORTED / SINGLE-TOOL) and key driver genes. Complements enrich_pathways with independent engines. Run after annotate_peaks.

| Argument | Type | Description |
|---|---|---|
| `preset` | string: `quick`, `standard`, `everything` | quick = local engines only (seconds); standard adds web tools (~1-2 min) |
| `lists` | array of string: `promoter`, `distal`, `all` |  |
| `local_only` | boolean | Never send genes to web services |

## `find_motifs`

Known-motif enrichment at peak centres with HOMER. Skipped cleanly when HOMER is not installed. Use regions='distal' to find enhancer-driving TFs without promoter GC bias. Not meaningful for broad histone domains.

| Argument | Type | Description |
|---|---|---|
| `regions` | string: `all`, `promoter`, `distal` |  |
| `size` | integer | Window around peak centre, bp (default 200) |
| `max_peaks` | integer | Cap on peaks scanned (default 5000) |
| `padj` | number |  |

## `chipatlas_submit` *(network)*

Submit the peaks to ChIP-Atlas Enrichment Analysis: overlap with every public ChIP-seq/ATAC-seq experiment for this genome, vs random permutation. Identifies TFs/marks that co-bind and the most similar public experiments. Asynchronous: submit early, collect later.

| Argument | Type | Description |
|---|---|---|
| `antigen_class` | string: `TFs and others`, `Histone`, `ATAC-Seq`, `DNase-seq`, `RNA polymerase`, `Input control`, `Bisulfite-Seq` |  |
| `cell_class` | string | ChIP-Atlas cell-type class, e.g. 'All cell types', 'Blood', 'Neural', 'Liver', 'Pluripotent stem cell' |
| `regions` | string: `all`, `promoter`, `distal` |  |
| `threshold` | integer: `50`, `100`, `200`, `500` | MACS2 -10log10(q) threshold for public peaks |
| `distance` | integer |  |

## `chipatlas_collect` *(network)*

Fetch results of a ChIP-Atlas job (summarised per antigen). Returns ready=false if still running; set wait_minutes to block.

| Argument | Type | Description |
|---|---|---|
| `request_id` | string |  |
| `wait_minutes` | number |  |

## `chipatlas_target_genes` *(network)*

Known targets of a TF from public ChIP-seq (ChIP-Atlas, peaks within +/- 1/5/10 kb of TSS) and their overlap with this run's peak-linked genes (hypergeometric test).

| Argument | Type | Description |
|---|---|---|
| `tf` (required) | string | Official symbol as in ChIP-Atlas, e.g. Tox, CTCF |
| `distance_kb` | integer: `1`, `5`, `10` |  |

## `integrate_rnaseq`

Integrate the user's differential-expression table: are up/down-regulated genes enriched for peaks (promoter or any)? Lists bound DE genes as candidate direct targets.

| Argument | Type | Description |
|---|---|---|
| `padj` | number |  |
| `lfc` | number | abs log2FC cutoff |

## `check_expectations`

Positive controls: are the user's expected genes bound, and are expected TFs (and the ChIPped TF itself) recovered by motifs / ChIP-Atlas? Run after motifs and ChIP-Atlas.

| Argument | Type | Description |
|---|---|---|
| `extra_genes` | array of string |  |
| `extra_tfs` | array of string |  |

## `search_literature` *(network)*

Search Europe PMC for papers (with abstracts). Use specific queries combining TFs/genes/pathways found here with the cell type and condition. Returns PMIDs you may cite.

| Argument | Type | Description |
|---|---|---|
| `query` (required) | string |  |
| `max_results` | integer |  |
| `from_year` | integer |  |
| `sort` | string: `cited`, `recent` | 'cited' (default) finds landmark studies; 'recent' the newest work |

## `search_geo` *(network)*

Search NCBI GEO for related public datasets (series) for this organism.

| Argument | Type | Description |
|---|---|---|
| `query` (required) | string |  |
| `max_results` | integer |  |

## `inspect_table`

Read rows of a table produced earlier (by evidence id), optionally filtered and sorted - e.g. check a specific gene's peaks or a motif's statistics.

| Argument | Type | Description |
|---|---|---|
| `evidence_id` (required) | string |  |
| `table` | string |  |
| `filter_column` | string |  |
| `filter_value` | string |  |
| `sort_by` | string |  |
| `descending` | boolean |  |
| `n` | integer |  |
