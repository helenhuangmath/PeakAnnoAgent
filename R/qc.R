# Peak-set QC that needs only data.table and base graphics, so it runs before
# (and independently of) any Bioconductor annotation backend.

standard_chromosomes <- function(chr) {
  grepl("^chr([0-9]+|X|Y|M|MT)$", chr)
}

# Overlap of peaks with a blacklist BED, using data.table::foverlaps.
blacklist_overlaps <- function(peaks, blacklist_path) {
  require_pkg("data.table")
  bl <- data.table::fread(blacklist_path, header = FALSE, select = 1:3,
    col.names = c("Chr", "Start", "End"))
  bl[, `:=`(Chr = as.character(Chr), Start = as.integer(Start), End = as.integer(End))]
  q <- data.table::copy(peaks[, .(PeakID, Chr, Start, End)])
  data.table::setkey(bl, Chr, Start, End)
  hits <- data.table::foverlaps(q, bl, by.x = c("Chr", "Start", "End"), nomatch = NULL)
  unique(hits$PeakID)
}

run_peak_qc <- function(peaks, out_dir, prefix, blacklist = NULL) {
  require_pkg("data.table")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  peaks <- data.table::copy(peaks)
  peaks[, Width := End - Start]

  warnings <- character()
  n <- nrow(peaks)
  if (n == 0) stop("The peak file contains no peaks.", call. = FALSE)
  bad_coords <- peaks[is.na(Start) | is.na(End) | End <= Start, .N]
  if (bad_coords > 0) warnings <- c(warnings, sprintf("%d peaks have End <= Start or missing coordinates.", bad_coords))
  dups <- sum(duplicated(peaks[, .(Chr, Start, End)]))
  if (dups > 0) warnings <- c(warnings, sprintf("%d duplicated peak intervals.", dups))
  nonstd <- peaks[!standard_chromosomes(Chr), .N]
  if (nonstd > 0) warnings <- c(warnings, sprintf("%d peaks on non-standard contigs (random/Un/alt).", nonstd))
  if (!any(startsWith(peaks$Chr, "chr"))) {
    warnings <- c(warnings, "Chromosome names lack the 'chr' prefix; UCSC-based annotation expects chr1, chr2, ...")
  }
  width_q <- stats::quantile(peaks$Width, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  if (width_q[[3]] > 5000) {
    warnings <- c(warnings, "Median peak width > 5 kb: looks like broad domains; use region-based (GREAT) enrichment and avoid summit-centred motif scans.")
  }
  if (n < 200) {
    warnings <- c(warnings, sprintf("Only %d peaks: enrichment and motif statistics will have low power.", n))
  }

  blacklisted <- character()
  if (!is.null(blacklist) && nzchar(blacklist)) {
    blacklisted <- blacklist_overlaps(peaks, blacklist)
    if (length(blacklisted)) {
      warnings <- c(warnings, sprintf("%d peaks (%.1f%%) overlap the blacklist.",
        length(blacklisted), 100 * length(blacklisted) / n))
    }
  }

  chrom <- peaks[, .N, by = Chr][order(-N)]
  chrom[, Pct := round(100 * N / n, 2)]
  write_tsv(chrom, file.path(out_dir, paste0(prefix, "_qc_chromosomes.tsv")))

  clean <- peaks[standard_chromosomes(Chr) & !is.na(Start) & End > Start & !PeakID %in% blacklisted]
  clean <- unique(clean, by = c("Chr", "Start", "End"))
  clean_path <- file.path(out_dir, paste0(prefix, "_qc_clean_peaks.bed"))
  data.table::fwrite(clean[, .(Chr, Start, End, PeakID)], clean_path, sep = "\t",
    col.names = FALSE, quote = FALSE)

  width_fig <- save_figure(function() {
    graphics::par(mar = c(6, 7, 3.5, 1))
    graphics::hist(log10(pmax(peaks$Width, 1)), breaks = 40, col = publication_colors()$bar,
      border = "white", main = "Peak width distribution", xlab = "log10(peak width, bp)")
    graphics::abline(v = log10(width_q[[3]]), lty = 2, col = publication_colors()$accent)
  }, out_dir, paste0(prefix, "_qc_peak_width"), width = 8, height = 5.5)

  chrom_fig <- save_figure(function() {
    top <- chrom[standard_chromosomes(Chr)]
    ord <- order(suppressWarnings(as.numeric(sub("^chr", "", top$Chr))), top$Chr)
    top <- top[ord]
    graphics::par(mar = c(5.5, 7.5, 3.5, 1))
    graphics::barplot(top$N, names.arg = sub("^chr", "", top$Chr), las = 2, col = publication_colors()$bar,
      border = NA, main = "Peaks per chromosome", ylab = "Peaks")
  }, out_dir, paste0(prefix, "_qc_chromosomes"), width = 10, height = 5.5)

  list(
    n_peaks = n,
    n_clean_peaks = nrow(clean),
    width = as.list(stats::setNames(as.numeric(width_q), c("min", "q25", "median", "q75", "max"))),
    n_duplicated = dups,
    n_nonstandard_chrom = nonstd,
    n_blacklisted = length(blacklisted),
    top_chromosomes = utils::head(chrom, 5),
    warnings = warnings,
    tables = list(chromosomes = file.path(out_dir, paste0(prefix, "_qc_chromosomes.tsv"))),
    clean_peaks = clean_path,
    figures = list(width = width_fig[["png"]], chromosomes = chrom_fig[["png"]]),
    pdfs = list(width_fig[["pdf"]], chrom_fig[["pdf"]])
  )
}
