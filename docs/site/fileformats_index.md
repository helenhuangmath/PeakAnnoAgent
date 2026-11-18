# File formats

What PeakAnnoAgent reads and what it writes.

## Input files

### Peaks (BED)

A BED file, or a tab-separated table with `Chr`, `Start` and `End` columns. Coordinates are
zero-based and half-open, as in BED. A fourth column is used as the peak name; without one, peaks
are named `peak_00001`, `peak_00002` and so on. narrowPeak and broadPeak files from MACS work as
they are.

```text
chr1    3671498    3671800    peak_001
chr1    4785526    4786010    peak_002
```

### Gene models (GTF)

Optional (`--gtf`). Any GTF with `gene_id` and `transcript_id` attributes and `exon` records, plain
or gzipped. `CDS` records are used to find UTRs, `gene_name` gives the gene symbols, and
`gene_type` or `gene_biotype` gives the gene types. GFF3 is not supported. See
[Gene models and the annotation method](execution.md).

### Blacklist (BED)

Optional (`--blacklist`). Peaks that overlap a blacklisted region are removed during QC, for
example with the ENCODE blacklist for the genome build.

### Gene signatures (GMT)

Optional (`--gmt`). One gene set per line: name, description, then the genes, separated by tabs.

```text
TEX_PROGENITOR    synthetic    Tcf7    Slamf6    Cxcr5    Id3    Sell    Ccr7    Lef1
TEX_TERMINAL      synthetic    Havcr2  Entpd1    Cd160    Cd244a Tox     Pdcd1   Lag3
```

### Differential expression table

Optional (`--rnaseq`). A table with a gene column, a log2 fold-change column (`log2FoldChange`,
`logFC`, `log2FC`, `lfc` or `avg_log2FC`) and an adjusted p-value column (`padj`, `adj.P.Val`,
`FDR` or `qvalue`). DESeq2, edgeR, limma and Seurat outputs work as they are.

### Biology file

Optional (`--biology-file`). A JSON file with the fields `description`, `assay`, `target`,
`cell_type`, `condition`, `expected_genes` and `expected_tfs`, or a text or Markdown file of
free-text notes. See [run](peakannoagent/tool_usage.md#biology-file).

## Output files

```{toctree}
:maxdepth: 2

peakannoagent/outputs
```
