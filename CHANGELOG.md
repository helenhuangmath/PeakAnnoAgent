# Changelog

## 1.0.2

Released 2026-11-12. See `docs/release_1.0.2.md` for validation and external checks.

- Correct BED-to-GRanges conversion and half-open blacklist overlaps; reject fractional coordinates and empty QC results.
- Import data.table in the installed R package and document previously undocumented exports. Testing files are excluded from release distributions.
- Record successful, partial and failed analyses consistently in manifests, reports, CLI exit codes, web jobs and MCP reports. Preserve partial outputs and expose pathway-method failures.
- Isolate web jobs with unique, exclusively created directories; preserve outcomes across restarts, separate uploaded inputs, and write status/evidence/manifests atomically.
- Bundle the R engine and examples in Python wheels and source distributions. Align R/Python versions at 1.0.2 and require Python 3.10 or later.
- Add input SHA-256 fingerprints, runtime versions, unique tool logs, local dependency checks and explicit integration-test switches.
- Add resolved Python constraints and a Bioconductor container recipe.
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
