# Installation

PeakAnnoAgent has two parts:

- an **R engine**, which does all the statistics and draws the figures;
- a **Python agent**, which plans the analysis, calls public resources and writes the report.

## Requirements

| Component | Version | Needed for |
|---|---|---|
| R | ≥ 4.1 with Bioconductor | everything (annotation, enrichment, motifs, figures) |
| Python | ≥ 3.10 | the command line, report and web app |
| Claude API key (`ANTHROPIC_API_KEY`) | — | agent mode and the chat; not needed with `--no-llm` |

## 1. Get the code

```bash
git clone https://github.com/helenhuangmath/PeakAnnoAgent.git
cd PeakAnnoAgent
```

## 2. Install the R engine

```bash
Rscript scripts/install_dependencies.R
```

This installs ChIPpeakAnno, ChIPseeker, clusterProfiler, ReactomePA, msigdbr, monaLisa, JASPAR2024
and igraph, plus the TxDb, OrgDb and BSgenome packages for mm10 and hg38.

To install only the R package:

```r
# install.packages("remotes")
remotes::install_github("helenhuangmath/PeakAnnoAgent")
```

## 3. Install the Python agent

```bash
pip install -e "agent[agent]"    # Claude agent mode (Anthropic SDK)
pip install -e "agent[web]"      # web app (Flask)
pip install -e "agent[mcp]"      # optional: use from Claude Code / Claude Desktop
```

The deterministic mode (`--no-llm`) needs only the Python standard library.

## 4. Check the installation

```bash
peakanno-agent check --genome mm10
```

`check` reports which R packages, optional tools, API credentials and network services were found.
Unavailable optional tools are reported as skipped. Missing required analysis dependencies produce a partial or failed outcome.

:::{tip}
If your Bioconductor packages live in a different R installation, point the agent at it:

```bash
export PEAKANNOAGENT_RSCRIPT=/path/to/Rscript
```
:::

## Optional tools

| Tool | Adds | How to install |
|---|---|---|
| [HOMER](http://homer.ucsd.edu/homer/) | HOMER known-motif enrichment and annotation | put `findMotifsGenome.pl` on `PATH` |
| rGREAT | region-based GREAT enrichment | `BiocManager::install("rGREAT")` |
| [SimpNetworkViz](https://github.com/helenhuangmath/SimpNetworkViz) | clustered pathway networks with interactive HTML and Cytoscape files | `pip install git+https://github.com/helenhuangmath/SimpNetworkViz` |
| [PathwayAgent](https://github.com/helenhuangmath/PathwayAgent) | cross-tool pathway consensus | clone next to this repository with its `.venv`, or set `PEAKANNOAGENT_PATHWAYAGENT` |

## Building this documentation

```bash
pip install -r docs/site/requirements.txt
sphinx-build -b html docs/site docs/site/_build/html
open docs/site/_build/html/index.html
```

The commands above build the documentation locally. Hosted documentation configuration is managed separately.
