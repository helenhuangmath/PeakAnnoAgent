#!/usr/bin/env Rscript
# Builds examples/data/example_tex_peaks_mm10.bed: a small synthetic mm10 peak set
# placed at promoters and nearby distal sites of T-cell exhaustion / effector
# genes plus random background, so the full agent can be exercised end to end.
suppressPackageStartupMessages({
  library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db); library(GenomicRanges)
})
set.seed(7)
genes <- c("Tox", "Tox2", "Pdcd1", "Havcr2", "Entpd1", "Cd160", "Cd244a", "Lag3", "Tigit", "Ctla4",
  "Nr4a1", "Nr4a2", "Nr4a3", "Batf", "Irf4", "Prdm1", "Eomes", "Tbx21", "Tcf7", "Id3", "Id2",
  "Gzmb", "Gzmk", "Prf1", "Ifng", "Tnf", "Il10", "Il21", "Cxcr5", "Cx3cr1", "Klrg1", "Il7r",
  "Sell", "Ccr7", "Bcl6", "Myb", "Nfatc1", "Nfatc2", "Junb", "Fos", "Bach2", "Runx3", "Cd8a",
  "Cd8b1", "Ikzf2", "Lef1", "Zeb2", "Klf2", "S1pr1", "Ptger4", "Rgs16", "Ccl3", "Ccl4", "Ccl5",
  "Cish", "Socs1", "Socs3", "Stat3", "Stat4", "Il2ra", "Il2rb", "Cd69", "Egr2", "Egr1", "Dusp6")
ids <- mapIds(org.Mm.eg.db, genes, "ENTREZID", "SYMBOL")
g <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
g <- g[names(g) %in% ids]
tss <- resize(g, 1, fix = "start")
promoter <- resize(shift(tss, sample(-300:100, length(tss), TRUE) * ifelse(strand(tss) == "-", -1, 1)), 400, fix = "center")
distal <- unlist(GRangesList(lapply(seq_along(tss), function(i) {
  k <- sample(2:5, 1)
  resize(shift(tss[rep(i, k)], sample(c(-1, 1), k, TRUE) * sample(3000:60000, k)), 350, fix = "center")
})))
bg_genes <- sample(g <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene), 60)
background <- resize(shift(resize(bg_genes, 1, fix = "start"), sample(-200000:200000, 60)), 300, fix = "center")
all <- c(granges(promoter), granges(distal), granges(background))
strand(all) <- "*"
all <- sort(all[seqnames(all) %in% paste0("chr", c(1:19, "X"))])
out <- data.frame(chr = as.character(seqnames(all)), start = start(all) - 1L, end = end(all),
  name = sprintf("tex_peak_%04d", seq_along(all)), score = sample(50:1000, length(all), TRUE))
write.table(out, "examples/data/example_tex_peaks_mm10.bed", sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = FALSE)
cat("Wrote", nrow(out), "peaks\n")
