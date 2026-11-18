# Install PeakAnnoAgent

PeakAnnoAgent is one Python package. There is no R or Bioconductor dependency.

## Requirements

| | |
|---|---|
| Python | 3.10 or later |
| Packages | NumPy and Matplotlib, installed automatically |
| Operating system | macOS and Linux are tested; Windows is not |
| Disk | About 0.2 MB for the package; 20-40 MB per genome build for the downloaded gene models, plus up to 10 MB per gene-set library |

## Install from source

```bash
git clone https://github.com/helenhuangmath/PeakAnnoAgent.git
cd PeakAnnoAgent
pip install -e .                 # annotation, enrichment and the deterministic pipeline
```

Optional extras:

```bash
pip install -e ".[agent]"        # Claude agent mode (Anthropic SDK)
pip install -e ".[web]"          # web app (Flask)
pip install -e ".[mcp]"          # use from Claude Code or Claude Desktop
```

A virtual environment keeps it separate from your other Python packages:

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[web,agent]"
```

## Check the installation

```bash
peakanno-agent check
```

It lists the Python packages, whether HOMER and the optional companion tools were found, the cache
location, your API credentials, and whether the public services can be reached. See
[check](subcommands/check.md).

Then run the bundled example, which needs no API key:

```bash
peakanno-agent annotate --peaks peakannoagent/examples/example_tex_peaks_mm10.bed \
  --genome mm10 --out results/anno --prefix tex
```

The first run downloads the mouse GENCODE gene models (about 20 MB) and builds the reference in a
few seconds. Later runs take about a second.

## Reference data and the cache

Gene models and gene-set libraries are public files that are downloaded the first time they are
needed, then kept in a cache:

- `~/.cache/peakannoagent` by default (or `$XDG_CACHE_HOME/peakannoagent`);
- the folder named by the `PEAKANNOAGENT_CACHE` environment variable, if set.

On a machine without internet access, copy the cache folder from another machine, or pass your own
GTF with `--gtf` and your own gene sets with `--gmt`. `peakanno-agent annotate --no-download`
refuses to download and uses the cache only.

## Optional tools

| Tool | Used for | Without it |
|---|---|---|
| `ANTHROPIC_API_KEY` and the Anthropic SDK | Agent mode and the chat in the web app | Use `--no-llm` |
| [HOMER](http://homer.ucsd.edu/homer/) with its genome package, `findMotifsGenome.pl` on `PATH` | Motif enrichment | The motif step is skipped |
| [SimpNetworkViz](https://github.com/helenhuangmath/SimpNetworkViz) | Clustered, interactive pathway networks | The step is skipped; the built-in networks are still drawn |
| [PathwayAgent](https://github.com/helenhuangmath/PathwayAgent) | Cross-tool pathway consensus | The step is skipped |

## Use it from R

A small R package wraps the command line for people who work in R. It runs the Python program and
returns its tables as data frames, so the Python package must be installed first.

```r
# install.packages("remotes")
remotes::install_github("helenhuangmath/PeakAnnoAgent", subdir = "r/PeakAnnoAgent")
library(PeakAnnoAgent)
paa_check()
```

See [R interface](api/r_interface.md).

## Docker

```bash
docker build -t peakannoagent .
docker run --rm -v "$PWD:/work" -v peakanno-cache:/cache peakannoagent \
  annotate --peaks peaks.bed --genome hg38 --out results
```

## Building this documentation

```bash
pip install -r docs/site/requirements.txt
sphinx-build -b html docs/site docs/site/_build/html
```
