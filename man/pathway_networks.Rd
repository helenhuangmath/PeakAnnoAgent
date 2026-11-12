\name{build_enrichment_map}
\alias{build_enrichment_map}
\alias{plot_enrichment_map}
\alias{plot_gene_concept_network}
\alias{run_pathway_networks}
\title{Build and Draw Pathway Networks}
\description{Group enriched pathways by shared genes and draw pathway and gene-pathway networks.}
\usage{
build_enrichment_map(tab, padj_cutoff = 0.05, top_n = 30,
  min_jaccard = 0.25, k = 3)
plot_enrichment_map(em, title = "Pathway enrichment map")
plot_gene_concept_network(tab, padj_cutoff = 0.05, n_terms = 5,
  max_genes_per_term = 25, title = "Gene-pathway network",
  lead_terms = NULL, term_colours = NULL)
run_pathway_networks(tab, out_dir, prefix, padj_cutoff = 0.05,
  top_n = 30, min_jaccard = 0.25)
}
\arguments{
  \item{tab}{Enrichment table with method, group, term, padj, count and semicolon-separated genes.}
  \item{padj_cutoff}{Adjusted p-value threshold.}
  \item{top_n}{Maximum number of pathways for an enrichment map.}
  \item{min_jaccard}{Minimum shared-gene Jaccard similarity for edges.}
  \item{k}{Number of strongest edges retained per node before taking their union.}
  \item{em}{An enrichment map returned by build_enrichment_map.}
  \item{title}{Figure title.}
  \item{n_terms}{Number of top pathways to show in the gene-pathway network.}
  \item{max_genes_per_term}{Maximum displayed genes per pathway.}
  \item{lead_terms}{Optional vector of pathway terms selected as module representatives.}
  \item{term_colours}{Optional named vector of colours for pathway terms.}
  \item{out_dir}{Output directory.}
  \item{prefix}{Output file prefix.}
}
\value{The builder returns a graph with node and edge tables, or NULL for fewer than three
terms. Plot functions return plot objects or NULL when insufficient data exist.
The file-producing wrapper returns figure/table paths, module and hub-gene summaries, and warnings.}
\details{Modules use Louvain community detection with a fixed random seed. Networks summarize
gene-set overlap and are not molecular interaction or causal regulatory networks.}
