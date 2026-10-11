# Full ChIPpeakAnno workflow used by the agent: nearest-gene and promoter
# (TSS-window) assignment, genomic-region distribution, TSS-centred binding
# profile, and peak-derived gene sets with a matching genome background.

chippeakanno_gene_annotation <- function(species, txdb_pkg = NULL) {
  require_pkg("ChIPpeakAnno")
  txdb <- load_txdb(species, txdb_pkg)
  suppressMessages(ChIPpeakAnno::toGRanges(txdb, feature = "gene"))
}

add_symbols <- function(entrez, species, orgdb_pkg = NULL) {
  require_pkg("AnnotationDbi")
  orgdb <- load_orgdb(species, orgdb_pkg)
  ids <- unique(stats::na.omit(as.character(entrez)))
  if (!length(ids)) return(rep(NA_character_, length(entrez)))
  map <- suppressMessages(AnnotationDbi::mapIds(orgdb, keys = ids, column = "SYMBOL",
    keytype = "ENTREZID", multiVals = "first"))
  unname(map[as.character(entrez)])
}

run_chippeakanno_pipeline <- function(peaks, out_dir, prefix, species = "mm10",
                                      txdb_pkg = NULL, orgdb_pkg = NULL,
                                      promoter_upstream = 2000, promoter_downstream = 500,
                                      distal_max_distance = 100000) {
  require_pkg("ChIPpeakAnno")
  require_pkg("GenomicRanges")
  require_pkg("data.table")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  promoter_upstream <- as.integer(promoter_upstream)
  promoter_downstream <- as.integer(promoter_downstream)
  txdb <- load_txdb(species, txdb_pkg)
  anno_data <- chippeakanno_gene_annotation(species, txdb_pkg)
  gr <- peaks_to_granges(peaks)
  warnings <- character()

  # 1. Nearest gene (by TSS) for every peak.
  nearest <- ChIPpeakAnno::annotatePeakInBatch(gr, AnnotationData = anno_data,
    output = "nearestLocation", FeatureLocForDistance = "TSS", PeakLocForDistance = "middle",
    select = "first")
  nearest_dt <- data.table::data.table(
    PeakID = nearest$peak,
    nearest_entrez = as.character(nearest$feature),
    distance_to_tss = nearest$distancetoFeature,
    shortest_distance = nearest$shortestDistance,
    inside_feature = as.character(nearest$insideFeature)
  )
  nearest_dt[, nearest_symbol := add_symbols(nearest_entrez, species, orgdb_pkg)]

  # 2. Every gene whose TSS window (-upstream/+downstream) overlaps the peak;
  #    keeps both genes of a bidirectional promoter.
  promoter <- suppressMessages(ChIPpeakAnno::annotatePeakInBatch(gr, AnnotationData = anno_data,
    output = "overlapping", FeatureLocForDistance = "TSS",
    bindingRegion = c(-promoter_upstream, promoter_downstream)))
  promoter_dt <- if (length(promoter)) {
    data.table::data.table(PeakID = promoter$peak, entrez = as.character(promoter$feature))
  } else {
    data.table::data.table(PeakID = character(), entrez = character())
  }
  promoter_dt[, symbol := add_symbols(entrez, species, orgdb_pkg)]
  promoter_by_peak <- promoter_dt[!is.na(symbol), .(promoter_symbols = paste(unique(symbol), collapse = ";")),
    by = PeakID]

  peak_table <- merge(peaks[, .(PeakID, Chr, Start, End)], nearest_dt, by = "PeakID", all.x = TRUE)
  peak_table <- merge(peak_table, promoter_by_peak, by = "PeakID", all.x = TRUE)
  peak_table[, region_class := data.table::fifelse(!is.na(promoter_symbols), "promoter",
    data.table::fifelse(abs(distance_to_tss) <= distal_max_distance, "distal", "distal_far"))]
  data.table::setorder(peak_table, Chr, Start)
  peak_table_path <- file.path(out_dir, paste0(prefix, "_chippeakanno_peak_genes.tsv"))
  write_tsv(peak_table, peak_table_path)

  # 3. Genomic-region distribution (ChIPpeakAnno::assignChromosomeRegion).
  region <- ChIPpeakAnno::assignChromosomeRegion(gr, nucleotideLevel = FALSE,
    precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),
    TxDb = txdb, proximal.promoter.cutoff = c(upstream = promoter_upstream, downstream = promoter_downstream))
  region_dt <- data.table::data.table(Region = names(region$percentage),
    Pct = round(as.numeric(region$percentage), 2))
  region_path <- file.path(out_dir, paste0(prefix, "_chippeakanno_region_distribution.tsv"))
  write_tsv(region_dt, region_path)
  region_fig <- save_figure(function() {
    graphics::par(mar = c(6, 17, 3.5, 1))
    ord <- rev(seq_len(nrow(region_dt)))
    graphics::barplot(region_dt$Pct[ord], names.arg = region_dt$Region[ord], horiz = TRUE, las = 1,
      xlim = range(0, pretty(region_dt$Pct)), col = publication_colors()$bar, border = NA, xlab = "% of peaks", main = "Genomic region distribution")
  }, out_dir, paste0(prefix, "_chippeakanno_region_distribution"), width = 9, height = 5.5)

  # 4. Binding profile around TSSs.
  tss_fig <- tryCatch(save_figure(function() {
    graphics::par(col = publication_colors()$primary, mar = c(6, 7, 3.5, 1))   # binOverFeature has no colour argument
    ChIPpeakAnno::binOverFeature(GenomicRanges::granges(gr), annotationData = anno_data,
      radius = 5000, nbins = 50, FUN = length, errFun = 0, xlab = "Distance from TSS (bp)",
      ylab = "Peak count", main = "Peak distribution around TSS")
  }, out_dir, paste0(prefix, "_chippeakanno_tss_profile"), width = 8, height = 5.5),
  error = function(e) {
    warnings <<- c(warnings, paste("binOverFeature failed:", conditionMessage(e)))
    NULL
  })

  dist_fig <- save_figure(function() {
    d <- peak_table$distance_to_tss
    d <- d[!is.na(d)]
    clipped <- pmax(pmin(d, 100000), -100000) / 1000
    graphics::par(mar = c(6, 7, 3.5, 1))
    graphics::hist(clipped, breaks = 80, col = publication_colors()$bar, border = "white",
      main = "Distance to nearest TSS", xlab = "kb (clipped to +/-100 kb)")
    graphics::abline(v = 0, lty = 2, col = publication_colors()$accent)
  }, out_dir, paste0(prefix, "_chippeakanno_tss_distance"), width = 8, height = 5.5)

  # 5. Gene sets for enrichment plus a genome-wide background.
  gene_sets <- data.table::rbindlist(list(
    data.table::data.table(Cluster = "promoter", SYMBOL = unique(stats::na.omit(promoter_dt$symbol))),
    data.table::data.table(Cluster = "distal",
      SYMBOL = unique(stats::na.omit(peak_table[region_class == "distal", nearest_symbol]))),
    data.table::data.table(Cluster = "all", SYMBOL = unique(stats::na.omit(c(promoter_dt$symbol,
      peak_table[region_class != "distal_far", nearest_symbol]))))
  ))
  gene_sets_path <- file.path(out_dir, paste0(prefix, "_peak_gene_sets.tsv"))
  write_tsv(gene_sets, gene_sets_path)
  universe <- data.table::data.table(SYMBOL = sort(unique(stats::na.omit(
    add_symbols(names(anno_data), species, orgdb_pkg)))))
  universe_path <- file.path(out_dir, paste0(prefix, "_gene_universe.tsv"))
  write_tsv(universe, universe_path)

  class_counts <- peak_table[, .N, by = region_class]
  if (peak_table[, mean(region_class == "promoter")] > 0.6) {
    warnings <- c(warnings, "Most peaks are promoter-proximal: typical of promoter-bound factors or H3K4me3.")
  }
  if (peak_table[, mean(region_class == "distal_far")] > 0.3) {
    warnings <- c(warnings, sprintf("Over 30%% of peaks are > %d kb from any TSS; nearest-gene links are weak - prefer GREAT.",
      as.integer(distal_max_distance / 1000)))
  }
  top_promoter <- utils::head(peak_table[region_class == "promoter", unique(unlist(strsplit(promoter_symbols, ";")))], 40)

  list(
    n_peaks = nrow(peak_table),
    region_class_counts = stats::setNames(as.list(class_counts$N), class_counts$region_class),
    region_distribution_pct = stats::setNames(as.list(region_dt$Pct), region_dt$Region),
    median_abs_distance_to_tss = stats::median(abs(peak_table$distance_to_tss), na.rm = TRUE),
    gene_set_sizes = as.list(gene_sets[, .N, by = Cluster][, stats::setNames(N, Cluster)]),
    example_promoter_genes = top_promoter,
    promoter_window = list(upstream = promoter_upstream, downstream = promoter_downstream),
    warnings = warnings,
    tables = list(peak_genes = peak_table_path, region_distribution = region_path),
    gene_sets = gene_sets_path,
    universe = universe_path,
    figures = Filter(Negate(is.null), list(region_distribution = region_fig[["png"]],
      tss_profile = if (!is.null(tss_fig)) tss_fig[["png"]], tss_distance = dist_fig[["png"]])),
    versions = package_versions(c("ChIPpeakAnno", "GenomicRanges", txdb_pkg %||% infer_species_defaults(species)$txdb,
      orgdb_pkg %||% infer_species_defaults(species)$orgdb))
  )
}
