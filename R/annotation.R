peaks_to_granges <- function(peaks) {
  require_pkg("GenomicRanges")
  require_pkg("IRanges")
  GenomicRanges::GRanges(
    seqnames = peaks$Chr,
    ranges = IRanges::IRanges(start = peaks$Start, end = peaks$End)
  )
}

run_chipseeker_annotation <- function(peaks, out_dir, prefix, species = "mm10",
                                      txdb_pkg = NULL, orgdb_pkg = NULL,
                                      tss_upstream = 1000, tss_downstream = 1000) {
  require_pkg("ChIPseeker")
  require_pkg("ggplot2")
  txdb <- load_txdb(species, txdb_pkg)
  orgdb_pkg <- if (is.null(orgdb_pkg)) infer_species_defaults(species)$orgdb else orgdb_pkg
  require_pkg(orgdb_pkg)
  gr <- peaks_to_granges(peaks)
  anno <- ChIPseeker::annotatePeak(
    gr,
    tssRegion = c(-as.integer(tss_upstream), as.integer(tss_downstream)),
    TxDb = txdb,
    annoDb = orgdb_pkg,
    assignGenomicAnnotation = TRUE,
    genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")
  )

  anno_df <- as.data.frame(anno)
  names(anno_df)[names(anno_df) == "seqnames"] <- "Chr"
  names(anno_df)[names(anno_df) == "start"] <- "Start"
  names(anno_df)[names(anno_df) == "end"] <- "End"
  out_table <- file.path(out_dir, paste0(prefix, "_chipseeker_annotation.tsv"))
  write_tsv(anno_df, out_table)

  anno_simple <- sub(" \\(.*", "", anno_df$annotation)
  dist_tab <- as.data.frame(table(Annotation = anno_simple))
  dist_tab$Pct <- round(100 * dist_tab$Freq / sum(dist_tab$Freq), 2)
  dist_tab <- dist_tab[order(-dist_tab$Freq), ]
  write_tsv(dist_tab, file.path(out_dir, paste0(prefix, "_chipseeker_annotation_distribution.tsv")))

  grDevices::pdf(file.path(out_dir, paste0(prefix, "_chipseeker_pie.pdf")), height = 4.8, width = 7)
  # ChIPseeker::plotAnnoPie() sets par(mai = c(0, 0, 0, 0)), so a title passed as
  # main= has no margin to draw into and gets clipped. Reserve an outer margin
  # and write the title there instead.
  graphics::par(oma = c(0, 0, 2, 0))
  ChIPseeker::plotAnnoPie(anno)
  graphics::mtext(prefix, side = 3, outer = TRUE, line = 0.2, font = 2, cex = 1.1)
  grDevices::dev.off()

  p_bar <- ChIPseeker::plotAnnoBar(anno, title = prefix) + ggplot2::theme_classic()
  ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_chipseeker_bar.pdf")), p_bar, height = 3.2, width = 9)

  p_tss <- ChIPseeker::plotDistToTSS(anno, title = paste(prefix, "distance to TSS")) + ggplot2::theme_classic()
  ggplot2::ggsave(file.path(out_dir, paste0(prefix, "_chipseeker_tss_distance.pdf")), p_tss, height = 3.2, width = 9)

  invisible(list(annotation = out_table, object = anno))
}

run_rgreat_annotation <- function(peaks, out_dir, prefix, species = "mm10",
                                  rule = "oneClosest", one_distance = 2000) {
  require_pkg("rGREAT")
  gr <- peaks_to_granges(peaks)
  job <- rGREAT::submitGreatJob(gr, rule = rule, adv_oneDistance = as.integer(one_distance), species = species)
  association <- as.data.frame(rGREAT::getRegionGeneAssociations(job))
  out_table <- file.path(out_dir, paste0(prefix, "_rgreat_region_gene_association.tsv"))
  write_tsv(association, out_table)
  invisible(list(annotation = out_table, object = job))
}

run_homer_annotation <- function(peaks_path, out_dir, prefix, genome = "mm10") {
  exe <- Sys.which("annotatePeaks.pl")
  if (exe == "") {
    stop("HOMER annotatePeaks.pl was not found in PATH.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_table <- file.path(out_dir, paste0(prefix, "_homer_annotation.tsv"))
  err <- system2(exe, args = c(peaks_path, genome), stdout = out_table, stderr = TRUE)
  status <- attr(err, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("HOMER annotatePeaks.pl failed:\n%s", paste(err, collapse = "\n")), call. = FALSE)
  }
  invisible(list(annotation = out_table))
}

run_nearest_txdb_annotation <- function(peaks, out_dir, prefix, species = "mm10",
                                        txdb_pkg = NULL, orgdb_pkg = NULL) {
  require_pkg("GenomicRanges")
  require_pkg("GenomicFeatures")
  require_pkg("AnnotationDbi")
  require_pkg("data.table")
  txdb <- load_txdb(species, txdb_pkg)
  orgdb <- load_orgdb(species, orgdb_pkg)

  gr <- peaks_to_granges(peaks)
  gene_gr <- GenomicFeatures::genes(txdb)
  hits <- GenomicRanges::distanceToNearest(gr, gene_gr, ignore.strand = TRUE)
  out <- data.table::as.data.table(as.data.frame(hits))
  out[, PeakIndex := queryHits]
  out[, GeneIndex := subjectHits]
  gene_ids <- names(gene_gr)[out$GeneIndex]
  symbols <- AnnotationDbi::mapIds(orgdb, keys = gene_ids, column = "SYMBOL",
    keytype = "ENTREZID", multiVals = "first")
  out[, ENTREZID := gene_ids]
  out[, SYMBOL := unname(symbols)]
  out <- cbind(peaks[out$PeakIndex, ], out[, .(ENTREZID, SYMBOL, distance)])
  out_table <- file.path(out_dir, paste0(prefix, "_nearest_txdb_annotation.tsv"))
  write_tsv(out, out_table)
  invisible(list(annotation = out_table))
}

run_chipeakanno_annotation <- function(peaks, out_dir, prefix, species = "mm10",
                                       txdb_pkg = NULL) {
  require_pkg("ChIPpeakAnno")
  require_pkg("GenomicFeatures")
  txdb <- load_txdb(species, txdb_pkg)
  gr <- peaks_to_granges(peaks)
  anno_data <- GenomicFeatures::genes(txdb)
  anno <- ChIPpeakAnno::annotatePeakInBatch(gr, AnnotationData = anno_data)
  out_table <- file.path(out_dir, paste0(prefix, "_chIPpeakAnno_annotation.tsv"))
  write_tsv(as.data.frame(anno), out_table)
  invisible(list(annotation = out_table, object = anno))
}

run_chippeakananno_annotation <- run_chipeakanno_annotation

run_annotatr_annotation <- function(peaks, out_dir, prefix, species = "mm10",
                                    annotations = NULL) {
  require_pkg("annotatr")
  defaults <- infer_species_defaults(species)
  annot_species <- defaults$annotatr
  if (is.null(annotations) || length(annotations) == 0) {
    annotations <- paste0(annot_species, c("_genes_promoters", "_genes_5UTRs",
      "_genes_exons", "_genes_introns", "_genes_3UTRs"))
  }
  gr <- peaks_to_granges(peaks)
  annots <- annotatr::build_annotations(genome = annot_species, annotations = annotations)
  anno <- annotatr::annotate_regions(regions = gr, annotations = annots, ignore.strand = TRUE)
  out_table <- file.path(out_dir, paste0(prefix, "_annotatr_annotation.tsv"))
  write_tsv(as.data.frame(anno), out_table)
  invisible(list(annotation = out_table, object = anno))
}

export_annotation_web_inputs <- function(peaks, out_dir, prefix) {
  require_pkg("data.table")
  bed <- peaks[, intersect(c("Chr", "Start", "End"), names(peaks)), with = FALSE]
  bed_path <- file.path(out_dir, paste0(prefix, "_webtools_peaks.bed"))
  data.table::fwrite(bed, bed_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(list(bed = bed_path))
}
