#!/usr/bin/env Rscript

cran <- c("data.table", "ggplot2", "enrichR", "gprofiler2", "msigdbr")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc <- c(
  "ChIPseeker",
  "ChIPpeakAnno",
  "clusterProfiler",
  "enrichplot",
  "AnnotationDbi",
  "GenomicRanges",
  "GenomicFeatures",
  "IRanges",
  "annotatr",
  "ReactomePA",
  "rGREAT",
  "RDAVIDWebService",
  "WebGestaltR",
  "org.Mm.eg.db",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene"
)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
