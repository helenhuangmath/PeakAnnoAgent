# PeakAnnoAgent: fast peak annotation and evidence-linked interpretation

[![Release checks](https://github.com/helenhuangmath/PeakAnnoAgent/actions/workflows/ci.yml/badge.svg)](https://github.com/helenhuangmath/PeakAnnoAgent/actions/workflows/ci.yml)
[![Documentation](https://github.com/helenhuangmath/PeakAnnoAgent/actions/workflows/docs.yml/badge.svg)](https://github.com/helenhuangmath/PeakAnnoAgent/actions/workflows/docs.yml)
![License](https://img.shields.io/github/license/helenhuangmath/PeakAnnoAgent)
![Top language](https://img.shields.io/github/languages/top/helenhuangmath/PeakAnnoAgent)
![Last commit](https://img.shields.io/github/last-commit/helenhuangmath/PeakAnnoAgent)

**Latest version:** {{ release }} · source on [GitHub](https://github.com/helenhuangmath/PeakAnnoAgent) · one Python package, released
under the MIT License.

## Introduction

Peak callers such as MACS turn ChIP-seq, ATAC-seq, CUT&RUN or CUT&Tag reads into a list of genomic
intervals. The next questions are biological: which genes and genomic features the peaks sit on,
which pathways those genes belong to, which transcription factors bind there, and what all of it
suggests about the experiment. PeakAnnoAgent answers these in one run.

You give it a peak file (BED), a genome build (`mm10`, `mm39`, `hg38` or `hg19`) and, optionally,
what you already know about the biology. It then:

- annotates every peak against the gene models you choose (UCSC knownGene, GENCODE, RefSeq,
  Ensembl or any GTF) in about a second for 100,000 peaks;
- tests the linked genes for pathway enrichment and draws pathway networks;
- finds enriched transcription-factor motifs and co-bound public ChIP-seq experiments;
- searches the literature and public datasets;
- writes a self-contained HTML report in which every statement links to the result behind it.

An LLM agent (Claude) can plan the analysis and interpret the results, or a fixed, deterministic
plan can run with no API key. Every step is also an R function.

```{image} ../figures/report.png
:alt: PeakAnnoAgent HTML report
:width: 760px
```

## Changes in the development version

### Features and performance

- **One Python package.** The analysis engine is now Python; R and Bioconductor are no longer
  needed. `pip install` is the whole installation.
- **Fast peak annotation.** Gene models are parsed once and cached, and peaks are annotated with
  array lookups: about a second for the whole step on the bundled example and well under a second
  of annotation for 100,000 peaks.
- **`peakanno-agent annotate`**: annotation only, from the command line.
- **Choice of gene models**: `--reference gencode` (default), `gencode:<release>`, `refseq`,
  `ucsc`, or `--gtf <file>`. See [Gene models and the annotation method](execution.md).
- **Annotation figures**, each as PDF and PNG: feature-type pie chart, feature bar,
  distance-to-TSS bars, feature-overlap upset plot, peaks over chromosomes, genomic-region
  distribution, TSS profile and TSS-distance histogram.
- **Pathway enrichment** with a built-in hypergeometric test on MSigDB and KEGG gene sets, which
  are downloaded once and cached.

### Removed

- The R package and its Bioconductor backends (ChIPpeakAnno, ChIPseeker, clusterProfiler,
  monaLisa, rGREAT). Motif enrichment now needs HOMER and is skipped without it; GREAT region
  enrichment is gone.

### Documentation

- This site, a Python API reference, and a small R interface for people who work in R.

The full history is in the [changelog](reference/changelog.md).

## Install

PeakAnnoAgent is a Python package (Python 3.10 or later). From a source checkout:

```bash
git clone https://github.com/helenhuangmath/PeakAnnoAgent.git
cd PeakAnnoAgent
pip install -e ".[web,agent]"        # engine, web app and Claude SDK
peakanno-agent check                 # verify the installation
```

See [Install PeakAnnoAgent](installation.md) for the options, the reference-data cache and the
optional tools. macOS and Linux are tested; Windows is not.

## Usage

Example for an ATAC-seq peak set, with the agent interpreting the results:

```bash
peakanno-agent run --peaks peaks.bed --genome mm10 --out results/tex \
  --assay "ATAC-seq" --cell-type "CD8 T cells"
```

Example for the same analysis as a fixed pipeline, offline, with no API key:

```bash
peakanno-agent run --peaks peaks.bed --genome mm10 --out results/tex --no-llm --offline
```

Example with RefSeq gene models, or your own GTF (the default is GENCODE):

```bash
peakanno-agent run --peaks peaks.bed --genome hg38 --out results/k562 --reference refseq
peakanno-agent run --peaks peaks.bed --genome hg38 --out results/k562 --gtf my_annotation.gtf.gz
```

Example for the web app, where you upload peaks and watch the run in your browser:

```bash
peakanno-agent web --open
```

Example for annotation only, in about a second:

```bash
peakanno-agent annotate --peaks peaks.bed --genome hg38 --out results/anno
```

Example in Python:

```python
from peakannoagent import build_reference, read_peak_file, annotate_peaks
ref = build_reference("hg38", "gencode")
table = annotate_peaks(read_peak_file("peaks.bed"), ref)
```

Example in R, through the [R interface](api/r_interface.md):

```r
library(PeakAnnoAgent)
res <- annotate_peaks("peaks.bed", genome = "hg38")
```

The subcommands of `peakanno-agent`:

| Subcommand | Description |
|---|---|
| [`run`](peakannoagent/tool_usage.md) | Main function: annotate peaks, find pathways, motifs and regulators, write the HTML report |
| [`annotate`](subcommands/annotate.md) | Annotate peaks only: per-peak table, summaries and figures |
| [`web`](subcommands/web.md) | Start the local web app: submit runs, watch progress, search results, ask questions |
| [`check`](subcommands/check.md) | Check packages, optional tools, the cache, credentials and network access |
| [`tools`](subcommands/tools.md) | List the analysis tools the agent can call |

For questions that come up often, see [Common Q & A](qa.md).

## Contribute

Please read the [Code of Conduct](CODE_OF_CONDUCT.md) and the
[contributing guide](CONTRIBUTING.md). Bug reports and feature requests go to
[GitHub issues](https://github.com/helenhuangmath/PeakAnnoAgent/issues).

## Citation

If you use PeakAnnoAgent, please cite the software:

> Huang H. *PeakAnnoAgent: evidence-linked annotation and interpretation of genomic peak sets.*
> GitHub, 2026. <https://github.com/helenhuangmath/PeakAnnoAgent>

Please also cite the resources whose results you report (GENCODE, MSigDB, KEGG, HOMER, ChIP-Atlas
and others). The **Methods** section of every report lists
them with their versions.

## Other useful links

- [MACS](https://macs3-project.github.io/MACS/): peak calling
- [ChIPseeker](https://bioconductor.org/packages/ChIPseeker) and
  [ChIPpeakAnno](https://bioconductor.org/packages/ChIPpeakAnno): peak annotation in R, whose conventions this tool follows
- [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb): gene-set collections
- [GENCODE](https://www.gencodegenes.org/): reference gene annotation
- [ChIP-Atlas](https://chip-atlas.org/): public ChIP-seq and ATAC-seq experiments
- [HOMER](http://homer.ucsd.edu/homer/) and [JASPAR](https://jaspar.elixir.no/): motif analysis
- [deepTools](https://deeptools.readthedocs.io/) and [bedtools](https://bedtools.readthedocs.io/):
  signal tracks and interval arithmetic

```{toctree}
:hidden:
:maxdepth: 2

installation
execution
subcommands_index
fileformats_index
tutorial
qa
CODE_OF_CONDUCT
CONTRIBUTING
api/index
```
