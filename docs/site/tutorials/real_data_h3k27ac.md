# Tutorial 6: A real data set from start to finish

This tutorial runs PeakAnnoAgent on real ENCODE data and reads every part of the output. All
numbers and figures on this page come from one run of the command below, with no agent and no API
key, so you can reproduce them.

| Item | Details |
|---|---|
| **Data** | H3K27ac ChIP-seq in primary human CD8 T cells, ENCODE experiment ENCSR835OJV, file ENCFF283CAQ (replicated peaks, GRCh38) |
| **Peaks used** | The 1000 strongest peaks, bundled as `peakannoagent/examples/example_h3k27ac_cd8_hg38.bed` |
| **Question** | Which genes, pathways and regulators do the strongest active-chromatin regions in resting CD8 T cells point to? |
| **Run time** | About two minutes the first time (the gene models are downloaded and cached), about one minute afterwards |

## 1. Run the analysis

```bash
peakanno-agent run \
  --peaks peakannoagent/examples/example_h3k27ac_cd8_hg38.bed --genome hg38 \
  --biology-file peakannoagent/examples/example_h3k27ac_cd8_biology.json \
  --out results/h3k27ac --prefix h3k27ac --no-llm
```

The biology file tells the program what the experiment is and which genes and factors you expect
to find. They are used later as positive controls:

```json
{
  "assay": "H3K27ac ChIP-seq",
  "target": "H3K27ac",
  "cell_type": "CD8 T cells",
  "condition": "primary human blood CD8 T cells, resting",
  "expected_genes": ["CD8A", "CD8B", "CD3E", "CD3D", "LCK", "ZAP70", "IL7R", "TCF7", "RUNX3", "EOMES", "TBX21"],
  "expected_tfs": ["RUNX3", "ETS1", "TCF7", "LEF1", "TBX21", "EOMES", "FLI1"]
}
```

The log shows each step and its evidence ID. The same IDs appear in the report, so every statement
can be traced to the step that produced it:

```text
[E1]  peak_qc            -> ok
[E2]  annotate_peaks     -> ok
[E3]  chipatlas_submit   -> ok
[E4]  enrich_pathways    -> ok
[E5]  network_viz        -> ok
[E7]  find_motifs        -> skipped: HOMER is not installed
[E9]  chipatlas_collect  -> ok
[E10] check_expectations -> ok
[E11] search_literature  -> ok
[E12] search_geo         -> ok
Analysis status: success
```

Motif enrichment needs [HOMER](http://homer.ucsd.edu/homer/). It was not installed for this run, so
the step was skipped and the run still finished; install HOMER and its `hg38` package to get it.

## 2. Open the report

`results/h3k27ac/h3k27ac_report.html` is one self-contained file. It starts with the run settings
and a summary, and every result below it has its figure, its table and download links.

```{image} ../../figures/tutorial_h3k27ac/report.png
:alt: Top of the HTML report for the H3K27ac run
:width: 760px
```

The figures are also written as vector PDFs and the tables as TSV files in `results/h3k27ac/results/`.

## 3. Peak QC (E1)

All 1000 peaks passed: none duplicated, none on non-standard chromosomes. The median width is
2221 bp (range 497 to 18,349 bp), which is typical of a histone mark and much broader than a
transcription-factor peak.

## 4. Where the peaks are (E2)

Peaks were annotated against GENCODE v46 protein-coding genes with a promoter window of 2000 bp
upstream to 500 bp downstream of a TSS.

```{image} ../../figures/tutorial_h3k27ac/feature_pie.png
:alt: Feature types of the H3K27ac peaks
:width: 620px
```

| Feature | Peaks |
|---|---|
| Promoter (within 1 kb of a TSS) | 90.0 % |
| Promoter (1 to 2 kb) | 0.4 % |
| Distal intergenic | 6.8 % |
| Introns, exons, UTRs, downstream | 2.8 % |

The median distance to the nearest TSS is 224 bp. The distance plot shows the same thing: 84 % of
peaks lie within 1 kb of a TSS, with more signal downstream of the TSS (51 %) than upstream (33 %),
because H3K27ac spreads into the first nucleosomes of the gene body.

```{image} ../../figures/tutorial_h3k27ac/tss_distance_bins.png
:alt: Distance of peaks to the nearest TSS
:width: 700px
```

A peak this broad usually overlaps several feature types at once. The upset plot counts peaks by the
combination they overlap: 746 peaks cover a promoter, 5' UTR, exon and intron together.

```{image} ../../figures/tutorial_h3k27ac/feature_upset.png
:alt: Peaks by the combination of features they overlap
:width: 640px
```

**What this tells you.** H3K27ac marks both active promoters and enhancers, but the *strongest*
peaks are almost all promoters. Keep this in mind for the next sections: this subset describes
highly active promoters, not the enhancers that define the cell type.

## 5. Pathways (E4)

The genes linked to the peaks (1120 from promoters, 74 from distal peaks, 1178 in total) are tested
against GO Biological Process, KEGG and Hallmark gene sets with a hypergeometric test, using the
20,036 annotated protein-coding genes as background.

```{image} ../../figures/tutorial_h3k27ac/go_bp_dotplot.png
:alt: GO Biological Process dot plot
:width: 640px
```

| Collection | Significant terms (promoter genes) | Top terms |
|---|---|---|
| GO BP | 109 | RNA splicing (adjusted p = 1.9e-18, 85 genes), mRNA processing (3.1e-17, 87 genes) |
| Hallmark | 3 | TNFα signalling via NF-κB (1.0e-03), MYC targets V1 (3.5e-03), protein secretion (4.7e-02) |
| KEGG | 1 | Antigen processing and presentation by MHC class I (2.1e-02) |

```{image} ../../figures/tutorial_h3k27ac/hallmark_dotplot.png
:alt: Hallmark dot plot
:width: 520px
```

**What this tells you.** The top pathways are housekeeping functions: RNA splicing, mRNA
metabolism, translation, protein folding. That is the expected result for the strongest promoters
of any cell. The immune signal is smaller but present: NF-κB target genes (`CD69`, `IL7R`, `JUN`,
`FOSB`, `NFKB2`, `REL`, `RELB`, `KLF2`, `NR4A2`) and MHC class I antigen presentation.

## 6. Pathway networks (E4, E5)

Many of the 109 GO terms describe the same genes. The enrichment map links terms that share genes
and groups them into modules, so you read a handful of themes instead of a long list. A module is
numbered by significance and named after its most significant term.

```{image} ../../figures/tutorial_h3k27ac/enrichment_map.png
:alt: Enrichment map of GO terms
:width: 720px
```

Here the top GO terms form RNA splicing (Modules 1 and 2), mRNA decay and translation (Module 3),
apoptotic signalling (Module 4), protein ubiquitination (Module 5) and smaller modules.

The gene-pathway network shows which genes drive each pathway. For Hallmark, the NF-κB module holds
the immediate-early and T-cell activation genes named above:

```{image} ../../figures/tutorial_h3k27ac/gene_pathway_network.png
:alt: Gene-pathway network for the Hallmark terms
:width: 640px
```

If [SimpNetworkViz](https://github.com/helenhuangmath/SimpNetworkViz) is installed, the run also
writes a simplified pathway network with automatically named clusters, an interactive HTML version
and Cytoscape files, in `results/h3k27ac_snv_go_bp_all/`. It found nine clusters in the 40 most
significant GO terms:

```{image} ../../figures/tutorial_h3k27ac/snv_pathway_network.png
:alt: Simplified pathway network from SimpNetworkViz
:width: 620px
```

See [Tutorial 4](pathway_networks.md) for how to read these networks in detail.

## 7. Public ChIP-seq that binds the same regions (E3, E9)

The peaks are sent to [ChIP-Atlas](https://chip-atlas.org), which tests them against 29,254 public
transcription-factor and chromatin-factor experiments. The table lists, for each factor, its best
experiment:

| Factor | Best experiment in | Fold enrichment | Significant experiments |
|---|---|---|---|
| KMT2A | Cord blood mononuclear cells | 52.1 | 106 of 129 |
| MED12 | T cells | 30.2 | 44 of 60 |
| MED1 | T cells | 57.7 | 209 of 244 |
| HDAC1 | PBMC | 51.6 | 71 of 89 |
| BCL11B | TALL-1 (T-ALL) | 33.0 | 14 of 14 |
| BRD4 | PBMC | 64.9 | 657 of 778 |
| TBX21 (T-bet) | Th1 cells | 30.0 | 6 of 6 |
| RUNX1 | Acute myeloid leukaemia | 46.2 | 142 of 169 |
| ETS1 | CUTLL1 (T-ALL) | 42.2 | 62 of 67 |

ChIP-Atlas is a live database, so these numbers shift slightly from run to run as experiments are
added; everything else on this page is deterministic.

**What this tells you.** This is where the cell type shows. The general co-activators of active
chromatin come first (Mediator, BRD4, the H3K4 methyltransferase KMT2A), and their best-matching
experiments are in T cells and blood cells. T-cell lineage factors follow: T-bet, BCL11B, RUNX1 and
ETS1.

## 8. Positive controls (E10)

The run checks the genes and factors you listed in the biology file:

| Check | Result |
|---|---|
| Expected genes with a peak | 1 of 11 (`IL7R`, two promoter peaks) |
| Expected factors among the ChIP-Atlas results | 2 of 7 (TBX21 at rank 9, ETS1 at rank 16) |

A low score here is a result, not a failure. `CD8A`, `CD3E`, `LCK` and the other identity genes are
acetylated in CD8 T cells, but their peaks are not among the 1000 strongest in the genome, which
are dominated by housekeeping promoters. The check tells you that this subset is the wrong one for
a question about T-cell identity.

## 9. Literature and data sets (E11, E12)

The run searches Europe PMC and GEO with the cell type and assay. It returned eight papers on
H3K27ac in CD8 T cells and related GEO series, listed with links in the *References* section of the
report. Without the agent this is one generic query, so treat the list as a starting point.

## 10. What the run shows, and what to do next

1. The 1000 strongest H3K27ac peaks in resting CD8 T cells are highly active promoters (90 %).
2. Their genes are mostly housekeeping (RNA processing, translation), with an NF-κB and
   early-activation component.
3. The cell type is visible in the public ChIP-seq overlap (T-cell Mediator, T-bet, BCL11B, ETS1),
   not in the pathway results.
4. Most T-cell identity genes are missing from this subset.

To study T-cell identity with this mark:

- **Use the full peak set**, not the strongest 1000. Download ENCFF283CAQ from ENCODE and run it
  the same way; PeakAnnoAgent annotates 100,000 peaks in about a second.
- **Look at the distal peaks.** The run writes them to `results/h3k27ac_distal_peaks.bed`.
  Enhancers carry more cell-type information than promoters.
- **Add your own gene signatures** with `--gmt` and an expression table with `--rnaseq`, as in
  [Tutorial 1](exhausted_cd8.md).
- **Install HOMER** to add motif enrichment.

## Running the same data with the agent

Drop `--no-llm` and set `ANTHROPIC_API_KEY`:

```bash
export ANTHROPIC_API_KEY=...
peakanno-agent run \
  --peaks peakannoagent/examples/example_h3k27ac_cd8_hg38.bed --genome hg38 \
  --biology-file peakannoagent/examples/example_h3k27ac_cd8_biology.json \
  --out results/h3k27ac_agent --prefix h3k27ac
```

The statistics are computed by the same engine, so the tables above do not change for the same
settings. The agent changes three things:

- **The plan.** It reads the biology file and chooses steps and settings, for example adding the
  ImmuneSigDB collection because the cells are T cells, restricting ChIP-Atlas to blood cells, and
  testing promoter and distal peaks separately.
- **The literature search.** Several targeted queries instead of one generic query.
- **The report.** A written summary, an interpretation and hypotheses with a proposed experiment
  each. Every statement cites an evidence ID, and citations that cannot be verified are flagged.

This page documents the run without the agent. Agent output varies between runs, so it is not
reproduced here; [Tutorial 1](exhausted_cd8.md) walks through an agent run on the bundled example.
