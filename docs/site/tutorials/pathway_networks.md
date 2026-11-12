# Tutorial 4: Pathway networks

Enrichment tests return long lists of overlapping terms. "T cell activation", "regulation of T cell
activation", "lymphocyte activation" and "leukocyte activation" mostly share the same genes.
PeakAnnoAgent therefore draws two networks for every gene-set collection. Together they show the few
themes behind such a list.

## The enrichment map

```{image} ../../figures/enrichment_map.png
:alt: Pathway enrichment map
:width: 680px
```

### How it is built

1. **Nodes.** The 30 most significant terms (`padj ≤ 0.05`), pooled over the peak groups (promoter,
   distal, all). Each term keeps its best `padj` and the union of its genes.
2. **Edges.** Two terms are linked when their gene sets overlap (Jaccard index ≥ 0.25). Each term keeps
   only its 3 strongest links, so redundant ontologies such as GO do not collapse into one hairball.
3. **Modules.** Louvain community detection on the weighted network groups the terms into modules. Each
   module is named after its most significant term, and modules are numbered by significance (`M1` is
   the best). The clustering is seeded, so the same input always gives the same modules.
4. **Layout.** As in Cytoscape's AutoAnnotate and
   [SimpNetworkViz](https://github.com/helenhuangmath/SimpNetworkViz), each module is laid out on its
   own. The modules are then packed in rows, with room for every label, so they never overlap.

### How to read it

| Element | Meaning |
|---|---|
| Tinted outline and coloured title | One module (theme), with its name and number |
| Node size | Number of genes in the term |
| Node colour | Module. Grey terms share too few genes to join a module |
| Dark edge, width | Shared genes within a module (Jaccard index) |
| Faint edge | A link between two different modules |

Read the modules first, then open a module to see which terms it contains. The report and the agent
describe the results in terms of modules (`M1: Lymphocyte differentiation`, 8 terms), not as 30
separate terms.

## The gene–pathway network

```{image} ../../figures/gene_pathway_network.png
:alt: Gene-pathway network
:width: 600px
```

The gene–pathway network shows the lead term of each of the top modules (up to five pathways, up to 25
genes each) and links each pathway to its genes:

- Pathway nodes have the same colour as their module in the enrichment map.
- A gene found in one pathway takes that pathway's tint.
- **Dark genes are shared by two or more pathways.** They are the hubs that connect the themes, and
  are listed as `hub_genes` in the results. In the example these include *Tox*, *Tcf7*, *Lef1*, *Batf*
  and *Prdm1*.

## The interactive explorer

The Pathways card in the web app and the HTML report has an interactive version:

- **Pathway network.** It uses the same layout, outlines and colours as the figure, and you can drag
  nodes to rearrange them. Click a pathway to list its genes. Hover to see its module, `padj` and gene
  count.
- **Gene–pathway map.** Click a gene to see its pathways, its peaks and its RNA-seq direction.
- **Download.** Each view downloads as SVG (editable), PNG or TSV.

## Open the network in Cytoscape

Each enrichment map comes with two tables in `results/`:

```text
<prefix>_enrich1_network_<method>_nodes.tsv   # term, module, module_label, padj, count, groups, genes, x, y
<prefix>_enrich1_network_<method>_edges.tsv   # from, to, jaccard
```

1. *File ▸ Import ▸ Network from File*: choose the edges table, with `from` as the source and `to` as the
   target.
2. *File ▸ Import ▸ Table from File*: choose the nodes table and match `term` to the node name.
3. Style the nodes by `module` (colour) and `count` (size), and the edges by `jaccard` (width).
4. Optional: run AutoAnnotate using the existing `module` column.

The `x` and `y` columns hold the figure layout in inches. Use them to reproduce the same arrangement.

## Clustered networks with SimpNetworkViz

If [SimpNetworkViz](https://github.com/helenhuangmath/SimpNetworkViz) is installed, the agent can also
run `network_viz` on the enrichment results. This adds MCL-clustered pathway networks with
automatic word labels, a gene network, an interactive HTML viewer, and `.cyjs` files for Cytoscape:

```bash
pip install git+https://github.com/helenhuangmath/SimpNetworkViz
peakanno-agent check --genome mm10      # confirms that snv is found
```

## Tune the network in R

`build_enrichment_map()` and `plot_enrichment_map()` are ordinary R functions:

```r
library(PeakAnnoAgent)
tab <- data.table::fread("results/tex/results/tex_enrich1_enrichment.tsv")[method == "go_bp"]

em <- build_enrichment_map(tab, padj_cutoff = 0.05, top_n = 40, min_jaccard = 0.3, k = 3)
em$nodes[, .(term, module, module_label, padj, count)]

p <- plot_enrichment_map(em, title = "GO biological process")
size <- attr(p, "fig_size")                      # inches; the layout is drawn at 1 unit = 1 inch
ggplot2::ggsave("enrichment_map.pdf", p, width = size[1], height = size[2])
```

| Argument | Default | Effect |
|---|---|---|
| `top_n` | 30 | Number of terms drawn |
| `min_jaccard` | 0.25 | Minimum gene overlap for an edge. Raise it if modules merge; lower it if nothing connects |
| `k` | 3 | Strongest links kept per term |
| `padj_cutoff` | 0.05 | Significance threshold for a term to be drawn |
