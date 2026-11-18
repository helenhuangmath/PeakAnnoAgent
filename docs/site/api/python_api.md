# Python API

```python
from peakannoagent import (build_reference, read_peak_file, annotate_peaks, run_annotation,
                           run_peak_qc, run_enrichment, Engine)
```

All coordinates given to and returned by these functions are zero-based and half-open, as in BED.

## Gene models

### `build_reference(genome="mm10", reference=None, gtf=None, cache=None, rebuild=False)`

Build the annotation reference for a genome, or read it from the cache. Returns a `Reference`.

| Argument | Description |
|---|---|
| `genome` | `"mm10"`, `"mm39"`, `"hg38"` or `"hg19"` (`"mouse"` and `"human"` mean mm10 and hg38) |
| `reference` | `"gencode"` (default), `"gencode:<release>"`, `"refseq"` or `"ucsc"`; or the path or URL of a GTF |
| `gtf` | Path or URL of a GTF file; overrides `reference` |
| `cache` | Cache folder; default `PEAKANNOAGENT_CACHE` or `~/.cache/peakannoagent` |
| `rebuild` | Rebuild even if a cached copy exists |

A `Reference` holds `meta` (genome, source, counts, cache file), `genes` (ID, symbol, type,
position) and sorted arrays for transcription start sites, exons, introns, UTRs and transcript ends.
`Reference.symbols(gene_types)` lists the gene symbols used for annotation.

## Peaks

### `read_peak_file(path)`

Read a BED, narrowPeak or broadPeak file, or a table with `Chr`, `Start`, `End` columns (plain or
gzipped). Returns `Peaks`, with `ids`, `chrom`, `start` and `end`.

### `Peaks(ids, chrom, start, end)`

Build peaks from your own arrays, for example from a data frame:

```python
import numpy as np
from peakannoagent import Peaks
peaks = Peaks(list(df["name"]), df["chrom"].to_numpy(dtype=object),
              df["start"].to_numpy(), df["end"].to_numpy())
```

### `run_peak_qc(peaks, out_dir, prefix, blacklist=None)`

Count duplicates, non-standard contigs and blacklisted peaks, write the cleaned BED file and two
figures. Returns the summary.

## Annotation

### `annotate_peaks(peaks, reference, promoter_upstream=2000, promoter_downstream=500, downstream_bp=300, distal_max_distance=100000, gene_types="auto")`

Annotate peaks. Returns a dictionary of equal-length columns:

| Column | Content |
|---|---|
| `PeakID`, `Chr`, `Start`, `End` | The peak, as given |
| `feature` | One feature per peak: promoter (with a distance bin), 5' UTR, 3' UTR, 1st or other exon, 1st or other intron, downstream, distal intergenic |
| `feature_overlaps` | Every feature type the peak overlaps, separated by `;` |
| `nearest_symbol`, `nearest_gene_id`, `nearest_gene_type`, `nearest_transcript` | The gene of the nearest TSS |
| `distance_to_tss` | From the peak centre to that TSS; negative upstream of the gene |
| `promoter_symbols`, `promoter_gene_ids` | Every gene with a promoter window overlapping the peak |
| `region_class` | `promoter`, `distal` or `distal_far` |

`gene_types` is `"auto"` or `"protein_coding"` (protein-coding genes) or `"all"`.

### `run_annotation(peaks, out_dir, prefix, genome="mm10", reference=None, gtf=None, ..., draw_figures=True, cache=None)`

Annotate and write the tables, the gene sets and background for enrichment, and the eight figures.
`reference` may be a `Reference` or a name. Returns a JSON-serialisable summary with
`feature_pct`, `region_distribution_pct`, `region_class_counts`, `gene_set_sizes`, `timing_sec`,
`warnings`, and the paths in `tables`, `gene_sets`, `universe` and `figures`.

## Enrichment and networks

### `run_enrichment(gene_sets, out_dir, prefix, genome="mm10", methods=("go_bp",), universe=None, groups=None, padj=0.05, top=8, custom_gmt=None, cache=None)`

Over-representation analysis of each gene list in `gene_sets` (a table with `Cluster` and
`SYMBOL`), with dot plots, enrichment maps and gene-pathway networks. Returns `n_significant`,
`top_terms`, `pathway_modules`, `hub_genes`, `gene_set_sources`, and the table and figure paths.

| Method | Gene sets |
|---|---|
| `go_bp`, `go_mf`, `go_cc` | Gene Ontology (MSigDB C5 / M5) |
| `hallmark` | MSigDB Hallmark |
| `reactome`, `wikipathways` | MSigDB canonical pathways |
| `kegg` | KEGG (MSigDB KEGG Medicus for human, Enrichr KEGG 2019 for mouse) |
| `tf_targets` | GTRD transcription-factor targets |
| `immunesigdb` | ImmuneSigDB (human sets; matched to mouse genes by symbol) |
| `cell_types` | Cell-type signatures (MSigDB C8 / M8) |
| `custom` | Your own GMT file, `custom_gmt` |

## The engine

### `Engine().run(tool, params)`

Run one of `peak_qc`, `annotate_peaks`, `enrich_pathways`, `find_motifs` or `session_info` with a
dictionary of settings. Returns `{"status": "ok" | "error" | "skipped", "result": ..., "error": ...}`
and never raises for an analysis failure. This is the interface the agent uses.
