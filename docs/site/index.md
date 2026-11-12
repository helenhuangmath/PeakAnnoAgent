# PeakAnnoAgent

**From a peak file to an evidence-linked biological report: annotation, pathways, motifs, regulators
and literature in one run.**

PeakAnnoAgent annotates and interprets genomic peak sets from ChIP-seq, ATAC-seq, CUT&RUN or
CUT&Tag. You give it three things:

- a peak file (BED);
- a reference genome (`mm10`, `mm39`, `hg38` or `hg19`);
- what you already know about the biology.

It runs a full Bioconductor workflow, searches public data and the literature, and writes an HTML
report with testable hypotheses. **Every statement in the report links to the result it came from.**

```{image} ../figures/report.png
:alt: PeakAnnoAgent HTML report
:width: 760px
```

## Run it in your browser

PeakAnnoAgent runs on your own computer; the web app is a local page, not a hosted service. After
[installation](installation.md), start it with one command:

```bash
peakanno-agent web --port 8766 --open
```

**Run page: [http://127.0.0.1:8766](http://127.0.0.1:8766)** (works while the command above is running)

On the page you upload a peak file (or pick an example), choose the genome, and press **Run**.
[Tutorial 3](tutorials/web_app.md) walks through the page.

## What it does

| Step | What it does | Built on |
|---|---|---|
| **Peak QC** | Duplicates, non-standard contigs, blacklist removal, width and chromosome checks | data.table |
| **Annotation** | Nearest gene, bidirectional promoters, genomic-region distribution, TSS profiles | ChIPpeakAnno, ChIPseeker, rGREAT, HOMER |
| **Pathways** | GO, KEGG, Reactome, MSigDB and your own GMT signatures, against a genome background | clusterProfiler, ReactomePA, msigdbr |
| **Pathway networks** | Enrichment maps with labelled pathway modules, gene–pathway networks with hub genes, Cytoscape tables | igraph, ggplot2, SimpNetworkViz |
| **Pathway consensus** | Several independent enrichment engines combined into consensus themes | PathwayAgent |
| **Motifs** | Known-motif enrichment against a GC-matched background, genome-wide or distal-only | monaLisa + JASPAR, HOMER |
| **Regulators** | Co-binding with every public ChIP-seq/ATAC-seq experiment; known TF targets | ChIP-Atlas |
| **Your biology** | Positive-control genes and TFs, RNA-seq integration, assay-aware settings | — |
| **Context** | Related papers and public datasets | Europe PMC, NCBI GEO |
| **Report** | Self-contained HTML with summary, hypotheses, figures, tables, genome browser and Methods | — |

## Four ways to run it

| Mode | Command | Needs |
|---|---|---|
| **Web app**: form, live log, report, search, chat | `peakanno-agent web` | Flask; `ANTHROPIC_API_KEY` for agent mode and chat |
| **Agent**: Claude plans, runs tools and interprets | `peakanno-agent run …` | `ANTHROPIC_API_KEY` |
| **Deterministic**: fixed plan with a factual summary | `peakanno-agent run … --no-llm` | nothing |
| **R package**: every step as an R function | `run_agent_tool("annotate_peaks", …)` | R + Bioconductor |

## Where to start

- New here? Read the [Introduction](peakannoagent/introduction.md), then follow
  [Tutorial 1](tutorials/exhausted_cd8.md) with the bundled example. It takes about five minutes.
- Want every option? See [Tool usage](peakannoagent/tool_usage.md).
- Want to check your installation on real data? Run one of the
  [ENCODE examples](tutorials/encode_examples.md).

```{toctree}
:caption: Getting started
:maxdepth: 2

installation
```

```{toctree}
:caption: PeakAnnoAgent
:maxdepth: 2

peakannoagent/introduction
peakannoagent/tool_usage
peakannoagent/outputs
```

```{toctree}
:caption: Tutorials
:maxdepth: 2

tutorials/exhausted_cd8
tutorials/encode_examples
tutorials/web_app
tutorials/pathway_networks
tutorials/r_engine
```

```{toctree}
:caption: Reference
:maxdepth: 2

reference/agent_tools
reference/r_cli
reference/methods
reference/agent_design
reference/changelog
```

## Citation and contact

If you use PeakAnnoAgent, please cite the software:

> Huang H. *PeakAnnoAgent: evidence-linked annotation and interpretation of genomic peak sets.*
> GitHub, 2026. <https://github.com/helenhuangmath/PeakAnnoAgent>

Please also cite the tools whose results you report (ChIPpeakAnno, clusterProfiler, monaLisa, JASPAR,
ChIP-Atlas, MSigDB and others). The **Methods** section of every report lists them with their package
versions.

**Author and maintainer:** Hua Huang, University of Pennsylvania. Bug reports and feature requests:
[GitHub issues](https://github.com/helenhuangmath/PeakAnnoAgent/issues). Released under the MIT
License.
