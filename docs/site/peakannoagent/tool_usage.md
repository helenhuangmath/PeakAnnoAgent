# run

The main function. The other subcommands are listed on the [Subcommands](../subcommands_index.md) page.

```text
peakanno-agent run --peaks PEAKS --genome GENOME --out OUT [options]
```

## Required input

| Option | Description |
|---|---|
| `--peaks` | Peak file: BED, or a TSV with `Chr`/`Start`/`End` columns |
| `--genome` | `mm10`, `mm39`, `hg38` or `hg19` (`mouse` and `human` also work) |
| `--out` | Output directory |

The peak file can have a header, or give chromosome, start and end as its first three columns:

```text
chr1    3671498    3671800    peak_001
chr1    4785526    4786010    peak_002
```

## Options

| Option | Default | Description |
|---|---|---|
| `--prefix` | `run` | Prefix of every output file |
| `--blacklist` | none | Blacklist BED to remove, e.g. the ENCODE mm10 blacklist |
| `--reference` | `gencode` | Gene models for annotation: `gencode`, `gencode:<release>`, `refseq` or `ucsc` (see [Gene models](../execution.md)) |
| `--gtf` | none | Annotate against this GTF file instead (GENCODE, Ensembl, RefSeq or custom) |
| `--palette` | `chippeakanno` | Figure colours: `chippeakanno` (ChIPpeakAnno's default colours, default), `npg` (Nature Publishing Group), `spectral`, or a Matplotlib colormap name such as `Set2`, `Dark2`, `tab10`, `viridis` or `RdYlBu_r` |
| `--gmt` | none | Your own gene signatures (GMT) to test next to the public collections |
| `--rnaseq` | none | Differential-expression table: a gene column, a log2 fold-change column (`log2FoldChange`, `logFC`, `log2FC`, `lfc` or `avg_log2FC`) and an adjusted p-value (`padj`, `adj.P.Val`, `FDR` or `qvalue`) |
| `--no-llm` | off | Fixed analysis plan and a factual summary; no API key needed |
| `--offline` | off | No ChIP-Atlas, Europe PMC or GEO calls. Gene models and gene sets are still downloaded once if they are not cached |

## Known biology

What you know about the experiment changes the analysis, not just the text of the report.

| Option | Example | Effect |
|---|---|---|
| `--assay` | `"ATAC-seq"`, `"TF ChIP-seq"`, `"H3K27ac CUT&Tag"` | Sets the promoter window; broad histone marks skip motif scanning; enhancer assays get a distal-only motif scan |
| `--target` | `Tox`, `H3K27ac` | The ChIPped TF becomes a positive control for the motif and ChIP-Atlas results |
| `--cell-type` | `"CD8 T cells"` | Chooses the gene-set collections (for example ImmuneSigDB for immune cells) and the literature queries |
| `--condition` | `"chronic LCMV, day 30"` | Used in the literature and GEO searches |
| `--expected-genes` | `Tox,Pdcd1,Havcr2` | Genes expected to be bound, reported as positive controls |
| `--expected-tfs` | `Nfatc1,Batf` | TFs expected to drive the peaks, checked against motifs and ChIP-Atlas |
| `--biology` | free text | Background for planning and interpretation |
| `--biology-file` | `biology.json` or `notes.md` | The same fields as a JSON file, or free-text notes |

### Biology file

A JSON file with the same fields as the options above:

```json
{
  "description": "Peaks from chronically stimulated (exhausted-like) CD8 T cells. Exhaustion is driven by TOX and NR4A family factors downstream of persistent TCR/NFAT signalling.",
  "assay": "ATAC-seq",
  "target": "",
  "cell_type": "CD8 T cells",
  "condition": "chronic LCMV clone 13 infection, day 30",
  "expected_genes": ["Tox", "Pdcd1", "Havcr2", "Entpd1", "Tcf7", "Cd160"],
  "expected_tfs": ["Nfatc1", "Nr4a1", "Batf", "Tcf7"]
}
```

A text or Markdown file is read as free-text notes. Notes are always passed to the model as data, never
as instructions.

### Signature file (GMT)

One gene set per line: name, description, then the genes, separated by tabs.

```text
TEX_PROGENITOR    synthetic    Tcf7    Slamf6    Cxcr5    Id3    Sell    Ccr7    Lef1
TEX_TERMINAL      synthetic    Havcr2  Entpd1    Cd160    Cd244a Tox     Pdcd1   Lag3
```

## Advanced options

| Option | Default | Description |
|---|---|---|
| `--model` | `claude-opus-5-5` | Claude model for planning and interpretation |
| `--effort` | `high` | Reasoning effort: `low`, `medium`, `high`, `xhigh` or `max` |
| `--max-turns` | 40 | Maximum agent turns |
| `--chipatlas-wait` | 25 | Minutes the deterministic plan waits for ChIP-Atlas |

ChIP-Atlas jobs run on the DDBJ supercomputer and usually take 5–30 minutes. The agent submits its
job early and collects the results later.

## Environment variables

| Variable | Description |
|---|---|
| `ANTHROPIC_API_KEY` | Claude credentials for agent mode and the chat |
| `PEAKANNOAGENT_CACHE` | Folder for downloaded gene models and built references |
| `PEAKANNOAGENT_SNV` | SimpNetworkViz command, if `snv` is not on `PATH` |
| `PEAKANNOAGENT_PATHWAYAGENT` | PathwayAgent command, if it is not found automatically |

## Advanced example

A CD8 T-cell TF ChIP-seq with a blacklist, expected targets and an RNA-seq table, using the agent at a
higher reasoning effort:

```bash
peakanno-agent run \
  --peaks tox_chip_peaks.bed --genome mm10 --prefix tox \
  --blacklist mm10-blacklist.v2.bed \
  --assay "TF ChIP-seq" --target Tox --cell-type "CD8 T cells" \
  --condition "chronic LCMV clone 13, day 30" \
  --expected-genes Pdcd1,Havcr2,Entpd1 --expected-tfs Nfatc1,Nr4a1 \
  --rnaseq tox_ko_vs_wt_deseq2.tsv \
  --effort xhigh --out results/tox
```

## Use Claude Code as the agent (MCP)

```bash
pip install -e ".[mcp]"
claude mcp add peakanno -- python -m peakannoagent.mcp_server
```

Claude Code then calls `start_run`, the analysis tools and `write_report` itself.
