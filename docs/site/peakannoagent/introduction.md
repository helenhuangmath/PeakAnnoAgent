# Introduction

## Method

PeakAnnoAgent turns a set of genomic peaks into a biological interpretation that can be checked line by
line. An **R engine** runs the statistics: peak QC, annotation with ChIPpeakAnno, over-representation
analysis with clusterProfiler, pathway networks, known-motif enrichment with monaLisa or HOMER. A
**Python agent** decides which analyses to run and with which settings, given what you already know
about the experiment. It queries ChIP-Atlas, Europe PMC and NCBI GEO, and then writes the report.

In agent mode, Claude plans the analysis and interprets the results. In deterministic mode
(`--no-llm`), a fixed plan runs and the report gives a factual summary. The statistics are the same in
both modes, because they always come from the R engine.

```text
           peaks + genome + known biology
                          │
                ┌─────────▼─────────┐        Claude plans the analysis, chooses
                │  Python agent     │◄──────  settings, and interprets results
                │  (agent/)         │        (or a fixed plan with --no-llm)
                └──┬──────────┬─────┘
       JSON bridge │          │ HTTPS
     ┌─────────────▼───┐   ┌──▼──────────────────────────────┐
     │ R engine (R/)   │   │ ChIP-Atlas · Europe PMC · GEO    │
     │ QC, ChIPpeakAnno│   └─────────────────────────────────┘
     │ enrichment,     │
     │ networks, motifs│
     └────────┬────────┘
              ▼
   evidence store (E1, E2, …) ──► HTML report with checked citations
```

Every tool result is stored with an evidence ID (`E1`, `E2`, …). The report may cite only those IDs and
the PMIDs or GEO accessions actually retrieved in the run. Any other citation is sent back to the model
once; if it remains, it is struck through in the report.

## Basic usage

The repository includes a small synthetic example: **351 ATAC-seq peaks from exhausted-like CD8 T
cells** on mm10. It comes with the known biology, three exhaustion gene signatures and an RNA-seq
differential-expression table:

| File | Content |
|---|---|
| `examples/data/example_tex_peaks_mm10.bed` | 351 peaks (BED) |
| `examples/data/example_tex_biology.json` | assay, cell type, condition, expected genes and TFs |
| `examples/data/example_tex_signatures.gmt` | progenitor, terminal and effector-like exhaustion signatures |
| `examples/data/example_tex_rnaseq_de.tsv` | DESeq2-style table (`gene`, `log2FoldChange`, `padj`) |

### Run PeakAnnoAgent

```bash
peakanno-agent run \
  --peaks examples/data/example_tex_peaks_mm10.bed --genome mm10 \
  --biology-file examples/data/example_tex_biology.json \
  --gmt examples/data/example_tex_signatures.gmt \
  --rnaseq examples/data/example_tex_rnaseq_de.tsv \
  --out results/tex --prefix tex
```

Add `--no-llm` to run without an API key. Add `--offline` to skip the ChIP-Atlas, Europe PMC and GEO
calls.

Only `--peaks`, `--genome` and `--out` are required. The other inputs make the analysis fit your
experiment; see [Tool usage](tool_usage.md) for every option.

### Results

Open `results/tex/tex_report.html` in a browser. The run folder contains:

| File | Description |
|---|---|
| `tex_report.html` | Self-contained report: summary, hypotheses, interpretation, every result with figures and tables, genome browser, references, Methods, analysis log |
| `results/*.pdf` | Every figure as a publication-ready vector PDF (editable Helvetica text) |
| `results/*.tsv` | Every table: annotation, enrichment, network nodes and edges, motifs, ChIP-Atlas, RNA-seq overlap |
| `tex_evidence.json` | Every tool result with its evidence ID, papers and datasets |
| `tex_manifest.json` | Parameters, timing, model usage, unverified citations |

[Outputs](outputs.md) describes each file in detail. [Tutorial 1](../tutorials/exhausted_cd8.md)
walks through this example's results.

```{image} ../../figures/enrichment_map.png
:alt: Pathway enrichment map of the example
:width: 560px
```

*Pathway enrichment map of the example (GO biological process). Terms that share genes are grouped
into labelled modules.*

## Citation

> Huang H. *PeakAnnoAgent: evidence-linked annotation and interpretation of genomic peak sets.*
> GitHub, 2026. <https://github.com/helenhuangmath/PeakAnnoAgent>
