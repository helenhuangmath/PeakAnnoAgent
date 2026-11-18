# Changelog

## 1.1.0 — add-features

- **Interpretation without the agent**: runs with `--no-llm` now include an automated, rule-based
  reading of the results (where the peaks are, pathway themes, candidate regulators, positive
  controls, expression), each statement linked to its result.
- **Interactive pathway views** in the web app and report show readable pathway names wrapped onto
  several rows, keep labels apart and use the same colours and module names as the PDF figures.
- **Sign in with the browser**: agent mode and the chat accept an `ant auth login` profile as well
  as `ANTHROPIC_API_KEY`; messages say how to sign in when no credentials are found.
- **Figures**: module names are written out ("Module 1: ...") and explained in the figure; long
  names wrap and are never cut; gene-pathway networks grow to fit pathways with many genes.
- **Lookups work on Python builds without root certificates** (python.org on macOS): the HTTP
  helper falls back to the system certificate bundle; verification stays on.
- **Documentation**: Tutorial 6 walks through a real ENCODE run with every figure; the bundled
  without-agent demo is regenerated with the current engine.

## 1.1.0 — add-features

PeakAnnoAgent is now one Python package. The R package and every Bioconductor dependency are gone.

- **Python analysis engine** (`peakannoagent.engine`), replacing the R engine:
  - fast peak annotation: gene models are parsed once and cached, then peaks are annotated with
    array lookups; about a second for the whole step on the bundled example, well under a second
    of annotation for 100,000 peaks, 6-10 s once to build a reference
  - one genomic feature per peak (promoter, 5' UTR, 3' UTR, exon, intron, downstream, distal
    intergenic), the nearest TSS and gene, and all promoter-bound genes
  - peak QC, pathway over-representation (hypergeometric test, Benjamini-Hochberg), enrichment
    maps and gene-pathway networks, all figures drawn with Matplotlib as vector PDF
- **Gene models**: `--reference gencode` (default; pinned release per build), `gencode:<release>`,
  `refseq`, `ucsc`, or `--gtf <file>`. Peaks are annotated against protein-coding genes by default
  (`--gene-types all` for every gene). The web form has the same choice.
- **Figure colours**: `--palette chippeakanno` (ChIPpeakAnno's default colours, the default), `npg`
  (Nature Publishing Group), `spectral`, or any Matplotlib colormap name. Also in the web form, the MCP `start_run` tool and the R functions.
- **Figures**: larger type throughout; enrichment-map modules are laid out by their links inside
  rounded outlines instead of on a grid. Figures are written as PDF only: the PNG copies the report
  embeds (300 dpi) are removed when the report is written, so there is no `.preview` folder in the output.
- **`peakanno-agent annotate`**: annotation only, with its tables, eight figures and a JSON summary.
- **Python API**: `build_reference`, `read_peak_file`, `annotate_peaks`, `run_annotation`,
  `run_peak_qc`, `run_enrichment`, `Engine`.
- **R interface** (`r/PeakAnnoAgent`): `annotate_peaks()` and `run_analysis()` call the program and
  return data frames. It needs only `jsonlite`.
- **Gene sets** come from MSigDB 2024.1 (human and mouse) and Enrichr (mouse KEGG), downloaded once
  and cached. Methods: `go_bp`, `go_mf`, `go_cc`, `hallmark`, `reactome`, `wikipathways`, `kegg`,
  `tf_targets`, `immunesigdb`, `cell_types`, `custom`. The old `msigdb_*` names still work.
- Reference data is cached in `~/.cache/peakannoagent` (`PEAKANNOAGENT_CACHE` to change).
  `--offline` still skips ChIP-Atlas, literature and GEO; reference data is fetched once if missing.
- Reports embed each figure once instead of twice, which roughly halves their size.
- Documentation site in the MACS3 layout, published to GitHub Pages by a workflow.

Removed or changed behaviour:

- ChIPpeakAnno, ChIPseeker, clusterProfiler, ReactomePA, msigdbr, monaLisa and rGREAT are no
  longer used. The `great_enrichment` tool is removed.
- Motif enrichment needs HOMER. Without it the step is skipped, and a skipped motif step no longer
  makes a run "partial".
- The default gene models are GENCODE instead of the UCSC knownGene package; the annotation output
  files are named `<prefix>_peak_genes.tsv`, `<prefix>_feature_*`, `<prefix>_tss_*` and
  `<prefix>_region_*` instead of `<prefix>_chippeakanno_*`.
- `--rscript`, `PEAKANNOAGENT_RSCRIPT` and the `chipseeker_crosscheck` tool argument are removed.
- The repository layout is a standard Python project: `peakannoagent/`, `tests/`, `pyproject.toml`
  at the root; examples are bundled in `peakannoagent/examples/`.

## 1.0.2 — update-network-plots

Release preparation, 2026-11-12. See `docs/release_1.0.2.md` for validation and external checks.

- Correct BED-to-GRanges conversion and half-open blacklist overlaps; reject fractional coordinates and empty QC results.
- Import data.table in the installed R package, repair fixture discovery, move the source test runner out of installed test discovery, and document previously undocumented exports.
- Record successful, partial and failed analyses consistently in manifests, reports, CLI exit codes, web jobs and MCP reports. Preserve partial outputs and expose pathway-method failures.
- Isolate web jobs with unique, exclusively created directories; preserve outcomes across restarts, separate uploaded inputs, and write status/evidence/manifests atomically.
- Bundle the R engine and examples in Python wheels and source distributions. Align R/Python versions at 1.0.2 and require Python 3.10 or later.
- Add input SHA-256 fingerprints, runtime versions, unique tool logs, local dependency checks and explicit integration-test switches.
- Add GitHub Actions checks, resolved Python constraints, a Bioconductor container recipe and mouse/human biological validation with a CTCF motif positive control.
- Treat empty pathway results as a valid network skip and preserve genomic boundary/sequence warnings in reports.

- Pathway networks redrawn in the SimpNetworkViz / AutoAnnotate style, so they are no longer crowded:
  - each module is laid out on its own, and the modules are packed in rows with room for their labels
  - every module has a tinted outline and its name (`M1: ...`) above it, with one colour per module
  - node labels sit under the nodes; links between modules are faint
  - terms that join no module are grey, with no outline
  - the legend is drawn at the bottom, and the figure size follows the layout
- The gene–pathway network uses the same module colours. Genes found in a single pathway take that
  pathway's tint, and genes shared by several pathways are dark.
- The interactive pathway network uses the same layout, outlines and colours as the figure.
- Network node tables gain `x`/`y` columns (the figure layout, in inches).
- Pathway modules are now reproducible: Louvain clustering is seeded.
- `ggrepel` is no longer needed.

## 1.0.1

Public release.

- Four example peak sets on the home page, each runnable in one click or editable first:
  - Exhausted CD8 T cells, ATAC-seq (synthetic, mm10, 351 peaks, with signatures and RNA-seq)
  - CTCF ChIP-seq, K562 (ENCODE ENCFF582SNT, hg38, 1000 strongest peaks)
  - H3K27ac ChIP-seq, human CD8 T cells (ENCODE ENCFF283CAQ, hg38, 1000 strongest peaks)
  - H3K4me3 ChIP-seq, mouse bone-marrow macrophages (ENCODE ENCFF813FPW, mm10, 800 strongest peaks)
- Interactive pathway explorer in the Pathways card of the run page and in the HTML report: a clickable
  top-terms chart, a draggable pathway network (enrichment map, coloured by module) and a gene–pathway
  map. Click a pathway to list its genes, click a gene to see its pathways, peaks and expression; every
  view downloads as SVG, PNG or TSV. No external libraries.
- Every figure has PNG and PDF download links; the run page has a "Download all" zip of the report,
  evidence and every table and figure. New JSON routes: `/api/runs/<id>/explore/<E>` and
  `/api/runs/<id>/table/<E>/<name>`.
- MIT licence, author and contact recorded in the Python package, R package and README.

## 0.2.0

- Analysis agent, web app, search, chat, external tools (ChIP-Atlas, Europe PMC, GEO), HTML report.

## 0.1.0

- R package scaffold.
