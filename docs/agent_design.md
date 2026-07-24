# Agent Presentation

The best presentation for this project is a layered agent:

1. R package core
2. command-line agent for HPC and reproducible manuscript runs
3. optional Shiny web app for guided method selection
4. optional REST API after the CLI/package contract is stable

## Recommendation

Start with the `PeakAnnoAgent` R package plus CLI agent. This is the best fit for peak annotation and pathway analysis because the work often runs on HPC, depends on Bioconductor/HOMER installations, and needs exact input/output provenance for publication.

Use a web/Shiny layer as the interface, not as the core engine. Shiny can collect inputs, species, methods, cutoffs, and figure options, then call `scripts/peak_agent.R` behind the scenes.

Build an API only if another application needs to submit jobs programmatically. API deployment is heavier because annotation databases, HOMER, R/Bioconductor versions, and API credentials have to be managed on a server.

## Agent Contract

The agent should expose the same concepts in every interface:

- input peaks: BED or TSV with `Chr`, `Start`, `End`
- input genes: TSV with gene and optional cluster columns
- species: `mm10` or `hg38`
- annotation methods: `chipseeker`, `nearest`, `chipeakanno`, `annotatr`, `rgreat`, `homer`, `webexport`
- pathway methods: `clusterprofiler`, `kegg`, `reactome`, `msigdb`, `enrichr`, `gprofiler`, `webgestalt`, `metascape`, `david`, `export`
- output directory
- prefix
- p-adjust cutoff
- top terms to plot

## Interfaces

CLI:

```bash
Rscript scripts/peak_agent.R enrich \
  --genes clusters.tsv \
  --universe universe.tsv \
  --out results/pathway \
  --methods clusterprofiler,kegg,reactome,msigdb,enrichr,metascape \
  --species mm10 \
  --cluster-col Cluster \
  --gene-col SYMBOL \
  --prefix experiment1
```

Shiny:

- file upload controls for peaks, genes, universe
- species selector
- method checkboxes
- figure style controls
- run button
- results table and PDF download links

REST API:

- `POST /jobs`: submit files and method choices
- `GET /jobs/{id}`: check status
- `GET /jobs/{id}/results`: download archive

The API should call the same `PeakAnnoAgent` R functions as the CLI, not reimplement analysis logic.
