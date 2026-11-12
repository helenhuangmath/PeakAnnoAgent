# Tutorial 1: Exhausted CD8 T cells (ATAC-seq)

This tutorial runs PeakAnnoAgent on the bundled example and reads the report section by section. You
will:

1. run the analysis with and without the agent;
2. check the peaks and their annotation;
3. read the pathway results as modules instead of long lists;
4. compare motifs, ChIP-Atlas regulators and positive controls;
5. connect the peaks to gene expression.

**Time:** about 5 minutes with `--no-llm --offline`. Agent mode with ChIP-Atlas takes 15–35 minutes,
mostly waiting for the ChIP-Atlas job.

## The data

The example is synthetic. It contains 351 ATAC-seq peaks (mm10) placed around genes of chronically
stimulated, exhausted-like CD8 T cells:

```text
examples/data/
├── example_tex_peaks_mm10.bed        # 351 peaks
├── example_tex_biology.json          # what we know about the experiment
├── example_tex_signatures.gmt        # progenitor / terminal / effector-like exhaustion signatures
└── example_tex_rnaseq_de.tsv         # 47 genes with log2 fold change and padj
```

The biology file tells PeakAnnoAgent what to expect:

```json
{
  "assay": "ATAC-seq",
  "cell_type": "CD8 T cells",
  "condition": "chronic LCMV clone 13 infection, day 30",
  "expected_genes": ["Tox", "Pdcd1", "Havcr2", "Entpd1", "Tcf7", "Cd160"],
  "expected_tfs": ["Nfatc1", "Nr4a1", "Batf", "Tcf7"]
}
```

## Step 1: Run the analysis

**With the agent** (needs `ANTHROPIC_API_KEY`):

```bash
peakanno-agent run \
  --peaks examples/data/example_tex_peaks_mm10.bed --genome mm10 \
  --biology-file examples/data/example_tex_biology.json \
  --gmt examples/data/example_tex_signatures.gmt \
  --rnaseq examples/data/example_tex_rnaseq_de.tsv \
  --out results/tex --prefix tex
```

**Without the agent** (no key needed). Add `--no-llm`, and `--offline` to skip the web services:

```bash
peakanno-agent run \
  --peaks examples/data/example_tex_peaks_mm10.bed --genome mm10 \
  --biology-file examples/data/example_tex_biology.json \
  --gmt examples/data/example_tex_signatures.gmt \
  --rnaseq examples/data/example_tex_rnaseq_de.tsv \
  --out results/tex_fixed --prefix tex --no-llm --offline
```

Then open the report:

```bash
open results/tex/tex_report.html
```

:::{tip}
Generate either example locally with the commands above, then open its HTML report to follow along.
:::

### What the agent changes

The statistics are identical in both modes. The agent changes which analyses run and with which
settings:

| | With agent | Without agent (`--no-llm`) |
|---|---|---|
| Promoter window | −1000/+500 bp, chosen for sharp ATAC peaks | −2000/+500 bp default |
| Gene sets | GO BP, Hallmark, **ImmuneSigDB** and your signatures, promoter vs distal | GO BP, **KEGG**, Hallmark and your signatures |
| Motifs | Distal-only scan first, then all peaks | All peaks, then distal |
| Regulators | ChIP-Atlas restricted to blood experiments; Tcf7 and Batf target tests | ChIP-Atlas, all cell classes |
| Summary | Narrative with three testable hypotheses | List of key statistics |

## Step 2: Peaks and annotation

**Peak QC (E1).** All 351 peaks pass QC: no duplicates, no non-standard contigs, median width 350 bp.

**Annotation (E2).** ChIPpeakAnno assigns each peak to its nearest gene and classifies it:

| Class | Peaks |
|---|---|
| promoter | 71 |
| distal | 277 |
| distal_far | 3 |

Most peaks (79 %) are distal, which is typical of ATAC-seq in T cells. The agent therefore runs a
distal-only motif scan as well. The gene sets used for enrichment are **promoter** (72 genes),
**distal** (178) and **all** (192). The background is all 24,511 annotated genes.

```{image} ../../figures/region_distribution.png
:alt: Genomic distribution of peaks
:width: 520px
```

## Step 3: Pathways, read as modules

`enrich_pathways` (E4) tests each gene set against the genome background. Hundreds of GO terms are
significant, and many are near-duplicates ("T cell activation", "regulation of T cell activation",
"lymphocyte activation", …). Rather than reading a long list, read the **enrichment map**. It links terms
that share genes and groups them into modules:

```{image} ../../figures/enrichment_map.png
:alt: Pathway enrichment map
:width: 640px
```

Each tinted outline is a module, named after its most significant term (`M1`, `M2`, …). Here the
example's GO terms fall into five themes: lymphocyte differentiation, T-cell differentiation, adaptive
immune response, alpha-beta T-cell differentiation, and leukocyte cell–cell adhesion. The report and
the agent refer to these modules instead of to the individual terms. See
[Pathway networks](pathway_networks.md) for how to read the figure.

Your own signatures are tested the same way. The terminal-exhaustion signature is strongest at
promoters, while the progenitor and effector-like signatures appear in both groups:

```{image} ../../figures/demo/with_agent_signatures_dotplot.png
:alt: Exhaustion signature dot plot
:width: 520px
```

## Step 4: Motifs, regulators and positive controls

**Motifs (E5, E6).** With only 351 peaks, no known motif passes the multiple-testing correction. The
best nominal hits are STAT3, SP1/SP4 and Foxd3. The report says this plainly and does not build claims
on them.

**ChIP-Atlas co-binding (E16).** The peaks are compared with 4,119 public blood ChIP-seq experiments.
**Ep300**, a marker of active enhancers, is the top factor: 75 of 125 Ep300 experiments are
significant. Bcl11b comes second.

**Known TF targets (E7, E8).** Peak-linked promoter genes are strongly enriched among public Batf
targets (93 % vs 56 % background, p = 4 × 10⁻¹²) and Tcf7 targets (89 % vs 66 %, p = 1 × 10⁻⁵).

**Positive controls (E22).** Every expected gene (Tox, Pdcd1, Havcr2, Entpd1, Tcf7, Cd160) has at least
one peak, most at the promoter. This is a quick sanity check that the peaks reflect the expected cell
state.

## Step 5: Connect peaks to expression

`integrate_rnaseq` (E9) asks whether differentially expressed genes carry peaks. All 17 up-regulated
and all 10 down-regulated genes have a promoter peak (about 1.7-fold over background), so they are
candidate direct targets. The table `tex_rnaseq_integration.tsv` lists them.

## Step 6: Hypotheses

In agent mode, the report ends with testable hypotheses. Each one cites the evidence above by ID and
proposes the experiment that would test it. In the demo run the report is titled:

> *Accessible chromatin of chronically stimulated CD8 T cells: an Ep300/Runx3/Eomes enhancer landscape
> spanning progenitor and terminal exhaustion genes*

Click any evidence badge (`E16`, `E7`, …) in the report to jump to the result it refers to.

## Next steps

- Run the same analysis from the browser: [Tutorial 3: Web app](web_app.md).
- Try real data with a known answer: [Tutorial 2: ENCODE examples](encode_examples.md).
- Explore the pathway network interactively: [Tutorial 4: Pathway networks](pathway_networks.md).
