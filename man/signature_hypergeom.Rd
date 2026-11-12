\name{signature_hypergeom}
\alias{signature_hypergeom}
\title{Test Custom Gene Signatures Against an Explicit Universe}
\description{Compute one-sided hypergeometric enrichment for every gene-list/signature pair.}
\usage{signature_hypergeom(gene_lists, term2gene, universe, min_size = 3)}
\arguments{
  \item{gene_lists}{Named list of gene-symbol vectors.}
  \item{term2gene}{Table with term and gene columns.}
  \item{universe}{Unique eligible background genes; all tests are restricted to this universe.}
  \item{min_size}{Minimum signature size after intersection with the universe.}
}
\value{A data.table with cluster, signature, overlap counts, query/background ratios,
raw p-values and Benjamini-Hochberg adjusted p-values computed within each query group.}
\details{Zero-hit signatures are retained. The background is not reduced to the union
of genes present in the supplied signatures.}
