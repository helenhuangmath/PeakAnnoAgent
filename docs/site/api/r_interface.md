# R interface

The R package `PeakAnnoAgent` (in `r/PeakAnnoAgent`) lets you call the program from R. It builds a
command line, runs `peakanno-agent`, and reads the files it writes into data frames. It does no
analysis itself, and it needs only `jsonlite`.

## Install

Install the Python package first (see [Install](../installation.md)), then:

```r
# install.packages("remotes")
remotes::install_github("helenhuangmath/PeakAnnoAgent", subdir = "r/PeakAnnoAgent")
library(PeakAnnoAgent)
paa_version()
paa_check()
```

## Finding the program

`paa_command()` looks for the program in this order:

1. the R option `peakannoagent.command`, a character vector such as
   `c("/path/to/venv/bin/python", "-m", "peakannoagent")`;
2. the environment variable `PEAKANNOAGENT_COMMAND`;
3. `peakanno-agent` on the `PATH`;
4. `python3 -m peakannoagent`.

In RStudio the `PATH` often differs from your terminal's. If the program is not found, set the
option to the full path:

```r
options(peakannoagent.command = "/Users/you/venvs/peakanno/bin/peakanno-agent")
```

## Functions

### `annotate_peaks(peaks, genome, out_dir, prefix, reference, gtf, promoter_upstream, promoter_downstream, distal_max_distance, gene_types, figures)`

Runs [`annotate`](../subcommands/annotate.md). `peaks` is a BED file path or a data frame with
`Chr`, `Start`, `End` (zero-based, half-open) and optionally `PeakID`.

Returns a list:

| Element | Content |
|---|---|
| `annotation` | Data frame, one row per peak: feature, nearest gene, distance to TSS, promoter genes, region class |
| `feature_distribution`, `region_distribution`, `tss_distance_bins` | Summary data frames |
| `gene_sets`, `universe` | Gene sets for enrichment and the background |
| `figures` | Named paths of the PDF figures |
| `summary` | The run summary (reference used, percentages, timings, warnings) |
| `out_dir` | Where everything was written |

```r
res <- annotate_peaks("peaks.bed", genome = "mm10", reference = "gencode:vM25")
table(res$annotation$feature)
```

### `run_analysis(peaks, genome, out_dir, prefix, llm, offline, reference, gtf, blacklist, assay, target, cell_type, condition, biology, expected_genes, expected_tfs, gmt, rnaseq, echo)`

Runs [`run`](../peakannoagent/tool_usage.md). `llm = FALSE` (the default) is the deterministic
pipeline; `llm = TRUE` lets Claude plan and interpret and needs `ANTHROPIC_API_KEY`.

Returns a list with `status` (`"success"` or `"partial"`), `reasons`, `report` (the HTML file),
`annotation` and `enrichment` data frames, the paths of all `tables` and `figures`, and the
`manifest`. A failed run raises an error with the end of the program's log.

### `paa_command()`, `paa_check(offline = TRUE)`, `paa_version()`

Locate the program, print its installation check, and return its version.
