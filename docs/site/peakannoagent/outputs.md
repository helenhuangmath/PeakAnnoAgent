# Outputs

A run writes one folder (`--out`). Its top level holds the report and two JSON records. All figures and
tables are under `results/`. Names start with `--prefix`. The examples below use `tex`, from the
bundled example.

```text
results/tex/
├── tex_report.html          # the report (self-contained, open in any browser)
├── tex_evidence.json        # every tool result with its evidence ID
├── tex_manifest.json        # parameters, timing, model usage, unverified citations
└── results/
    ├── tex_qc_*                         # peak QC
    ├── tex_chippeakanno_*               # annotation
    ├── tex_enrich1_enrichment*          # pathway enrichment and dot plots
    ├── tex_enrich1_network_*            # enrichment maps (+ Cytoscape node/edge tables)
    ├── tex_enrich1_cnet_*               # gene–pathway networks
    ├── tex_*_motifs_*                   # motif enrichment
    ├── tex_chipatlas_*                  # ChIP-Atlas co-binding
    └── tex_rnaseq_integration.tsv       # RNA-seq overlap
```

## The report

`<prefix>_report.html` is a single file with everything embedded, so it can be e-mailed or archived.

| Section | Content |
|---|---|
| Summary | What the peaks are, the main pathways and regulators, and caveats. In agent mode this is a titled narrative; with `--no-llm` it is a list of key statistics |
| Hypotheses | Agent mode only. Testable hypotheses, each with a confidence, evidence IDs and a proposed experiment |
| Results | One card per analysis step: key statistics, figures, table previews and download links. The Pathways card includes the interactive [pathway explorer](../tutorials/pathway_networks.md) |
| Genome browser | igv.js view of the peaks |
| References | Papers and datasets retrieved in the run |
| Methods | Auto-written methods with package versions |
| Analysis log | Every tool call, its settings and its status |

## Evidence and manifest

`<prefix>_evidence.json` holds one record per tool call: evidence ID, tool, settings, status, summary
statistics, and the paths of its figures and tables. It is the source of truth for the report and the
web app.

`<prefix>_manifest.json` records the configuration, timing, Claude token usage and any citation that
could not be verified.

## Tables

| File | Content |
|---|---|
| `*_qc_clean_peaks.bed` | Peaks after QC; later steps use this file |
| `*_qc_chromosomes.tsv` | Peaks per chromosome |
| `*_chippeakanno_peak_genes.tsv` | Every peak with its nearest gene, distance to TSS, region class (promoter / distal / distal_far) and all overlapping promoters |
| `*_chippeakanno_region_distribution.tsv` | Genomic-region distribution |
| `*_peak_gene_sets.tsv`, `*_gene_universe.tsv` | Gene sets tested for enrichment, and the genome background |
| `*_enrich1_enrichment.tsv` | All enrichment results: `method`, `group`, `term_id`, `term`, `gene_ratio`, `bg_ratio`, `pvalue`, `padj`, `count`, `genes` |
| `*_enrich1_network_<method>_nodes.tsv` | Enrichment-map nodes: `term`, `module`, `module_label`, `padj`, `count`, `groups`, `genes`, and the figure layout `x`, `y` (inches) |
| `*_enrich1_network_<method>_edges.tsv` | Enrichment-map edges: `from`, `to`, `jaccard` |
| `*_all_motifs_monalisa.tsv`, `*_distal_motifs_monalisa.tsv` | Known-motif enrichment (all peaks, distal peaks) |
| `*_chipatlas_<id>.tsv` | ChIP-Atlas co-binding, one row per public experiment |
| `*_rnaseq_integration.tsv` | Differentially expressed genes with and without peaks |

The node and edge tables import directly into Cytoscape; see
[Pathway networks](../tutorials/pathway_networks.md#open-the-network-in-cytoscape).

## Figures

Every figure is written as a vector PDF that follows Nature's artwork guidelines (Helvetica text that
stays editable in Illustrator, no Dingbats glyphs). A PNG preview is embedded in the report.

| File | Figure |
|---|---|
| `*_qc_peak_width.pdf`, `*_qc_chromosomes.pdf` | Peak widths and peaks per chromosome |
| `*_chippeakanno_region_distribution.pdf` | Genomic-region distribution |
| `*_chippeakanno_tss_distance.pdf`, `*_chippeakanno_tss_profile.pdf` | Distance to TSS and TSS binding profile |
| `*_enrich1_enrichment_<method>_dotplot.pdf` | Top terms per peak group |
| `*_enrich1_network_<method>.pdf` | Enrichment map with labelled modules |
| `*_enrich1_cnet_<method>.pdf` | Gene–pathway network |
