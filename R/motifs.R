# Known-motif enrichment. HOMER is used when installed; otherwise monaLisa with
# JASPAR CORE vertebrate PWMs against a GC-matched random genomic background.
# Both backends write one standardised table:
#   motif, tf, family, pct_target, pct_background, log2_enrichment, pvalue, padj

motif_tf_name <- function(x) {
  x <- sub("/.*$", "", x)            # HOMER: "CTCF(Zf)/CD4+-CTCF-ChIP-Seq/Homer"
  x <- sub("\\(.*$", "", x)
  trimws(x)
}

# Keep the strongest peaks when a score column exists, else a reproducible sample.
cap_peaks <- function(peaks, max_peaks = 5000) {
  if (nrow(peaks) <= max_peaks) return(peaks)
  score_col <- intersect(c("Score", "score", "signalValue", "V5", "V7"), names(peaks))[1]
  if (!is.na(score_col) && is.numeric(peaks[[score_col]])) {
    return(peaks[order(-get(score_col))][seq_len(max_peaks)])
  }
  set.seed(1)
  peaks[sort(sample.int(nrow(peaks), max_peaks))]
}

get_jaspar_pwms <- function(species = "mm10", tax_group = "vertebrates") {
  require_pkg("TFBSTools")
  opts <- list(matrixtype = "PWM", tax_group = tax_group, collection = "CORE")
  if (requireNamespace("JASPAR2024", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {
    con <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPAR2024::JASPAR2024()))
    on.exit(RSQLite::dbDisconnect(con), add = TRUE)
    return(list(pwms = TFBSTools::getMatrixSet(con, opts), source = "JASPAR2024"))
  }
  for (pkg in c("JASPAR2022", "JASPAR2020")) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      db <- get(pkg, envir = asNamespace(pkg))
      return(list(pwms = TFBSTools::getMatrixSet(db, opts), source = pkg))
    }
  }
  stop("Install JASPAR2024 (or JASPAR2022/JASPAR2020) for monaLisa motif enrichment.", call. = FALSE)
}

run_monalisa_motifs <- function(peaks, out_dir, prefix, species = "mm10", size = 200,
                                max_peaks = 5000, bsgenome_pkg = NULL) {
  require_pkg("monaLisa")
  require_pkg("Biostrings")
  require_pkg("SummarizedExperiment")
  require_pkg("BiocParallel")
  bsg <- load_bsgenome(species, bsgenome_pkg)
  peaks <- cap_peaks(peaks, max_peaks)
  gr <- GenomicRanges::resize(peaks_to_granges(peaks), width = as.integer(size), fix = "center")
  seqs <- Biostrings::getSeq(bsg, gr)
  jaspar <- get_jaspar_pwms(species)
  se <- monaLisa::calcBinnedMotifEnrR(seqs = seqs, pwmL = jaspar$pwms, background = "genome",
    genome = bsg, genome.regions = NULL, genome.oversample = 2,
    BPPARAM = BiocParallel::SerialParam(RNGseed = 42), verbose = FALSE)
  a <- function(name) as.numeric(SummarizedExperiment::assay(se, name)[, 1])
  cd <- SummarizedExperiment::colData(se)
  rd <- SummarizedExperiment::rowData(se)
  res <- data.table::data.table(
    motif = paste(rd$motif.id, rd$motif.name, sep = " "),
    tf = rd$motif.name,
    family = vapply(jaspar$pwms[rd$motif.id], function(m) paste(TFBSTools::tags(m)$family, collapse = ";"),
      character(1)),
    pct_target = round(100 * a("sumForegroundWgtWithHits") / cd$totalWgtForeground[1], 2),
    pct_background = round(100 * a("sumBackgroundWgtWithHits") / cd$totalWgtBackground[1], 2),
    log2_enrichment = round(a("log2enr"), 3),
    pvalue = 10^(-a("negLog10P")),
    padj = 10^(-a("negLog10Padj"))
  )
  list(table = res[order(padj, -log2_enrichment)], backend = "monaLisa", database = jaspar$source,
    n_sequences = length(seqs), versions = package_versions(c("monaLisa", "TFBSTools", jaspar$source,
      bsgenome_pkg %||% infer_species_defaults(species)$bsgenome)))
}

read_homer_known_results <- function(path) {
  dt <- data.table::fread(path, sep = "\t", check.names = FALSE)
  pct <- function(x) as.numeric(sub("%", "", x, fixed = TRUE))
  data.table::data.table(
    motif = dt[["Motif Name"]],
    tf = motif_tf_name(dt[["Motif Name"]]),
    family = ifelse(grepl("\\(", dt[["Motif Name"]]), sub("^[^(]*\\(([^)]*)\\).*$", "\\1", dt[["Motif Name"]]), ""),
    pct_target = pct(dt[[7]]),
    pct_background = pct(dt[[9]]),
    log2_enrichment = round(log2((pct(dt[[7]]) + 0.1) / (pct(dt[[9]]) + 0.1)), 3),
    pvalue = as.numeric(dt[["P-value"]]),
    padj = as.numeric(dt[["q-value (Benjamini)"]])
  )
}

run_homer_motifs <- function(peaks, out_dir, prefix, species = "mm10", size = 200,
                             max_peaks = 5000, genome = NULL, threads = 4, denovo = FALSE) {
  exe <- Sys.which("findMotifsGenome.pl")
  if (exe == "") stop("HOMER findMotifsGenome.pl was not found in PATH.", call. = FALSE)
  genome <- genome %||% infer_species_defaults(species)$homer
  peaks <- cap_peaks(peaks, max_peaks)
  homer_dir <- file.path(out_dir, paste0(prefix, "_homer_motifs"))
  dir.create(homer_dir, recursive = TRUE, showWarnings = FALSE)
  bed <- file.path(homer_dir, "input_peaks.bed")
  data.table::fwrite(peaks[, .(Chr, Start, End, PeakID, ".", "+")], bed, sep = "\t",
    col.names = FALSE, quote = FALSE)
  args <- c(bed, genome, homer_dir, "-size", as.character(size), "-p", as.character(threads))
  if (!isTRUE(denovo)) args <- c(args, "-nomotif")
  log <- system2(exe, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status")
  known <- file.path(homer_dir, "knownResults.txt")
  if ((!is.null(status) && status != 0) || !file.exists(known)) {
    stop(sprintf("HOMER findMotifsGenome.pl failed:\n%s", paste(utils::tail(log, 20), collapse = "\n")),
      call. = FALSE)
  }
  list(table = read_homer_known_results(known), backend = "HOMER", database = "HOMER known motifs",
    n_sequences = nrow(peaks), homer_dir = homer_dir, versions = list())
}

run_motif_enrichment <- function(peaks, out_dir, prefix, species = "mm10", backend = "auto",
                                 size = 200, max_peaks = 5000, padj_cutoff = 0.05, top = 25, ...) {
  require_pkg("data.table")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (backend == "auto") {
    backend <- if (nzchar(Sys.which("findMotifsGenome.pl"))) "homer" else "monalisa"
  }
  res <- switch(backend,
    homer = run_homer_motifs(peaks, out_dir, prefix, species = species, size = size, max_peaks = max_peaks, ...),
    monalisa = run_monalisa_motifs(peaks, out_dir, prefix, species = species, size = size, max_peaks = max_peaks),
    stop(sprintf("Unknown motif backend '%s' (use auto, homer, monalisa).", backend), call. = FALSE)
  )
  tab <- res$table
  table_path <- file.path(out_dir, paste0(prefix, "_motifs_", backend, ".tsv"))
  write_tsv(tab, table_path)
  sig <- tab[padj <= padj_cutoff & log2_enrichment > 0]

  fig <- NULL
  if (nrow(sig)) {
    show <- utils::head(sig, min(top, 20))
    fig <- save_figure(function() {
      graphics::par(mar = c(6, 16, 3.5, 1))
      vals <- rev(pmin(-log10(pmax(show$padj, 1e-300)), 300))
      graphics::barplot(vals, names.arg = rev(show$tf), horiz = TRUE, las = 1, xlim = range(0, pretty(vals)), col = publication_colors()$bar,
        border = NA, xlab = "-log10 adjusted p", main = sprintf("Known motif enrichment (%s)", res$backend))
    }, out_dir, paste0(prefix, "_motifs_", backend), width = 9, height = max(5.5, nrow(show) * 0.45 + 2.5))
  }

  warnings <- character()
  nominal <- NULL
  if (!nrow(sig)) {
    warnings <- c(warnings, "No motif passed the adjusted p-value cutoff; top_nominal lists the best motifs by raw p-value (not significant).")
    nominal <- utils::head(tab[log2_enrichment > 0][order(pvalue), .(tf, log2_enrichment, pvalue, padj)], 10)
  }
  list(
    backend = res$backend,
    database = res$database,
    n_sequences = res$n_sequences,
    window_bp = size,
    n_significant = nrow(sig),
    top_motifs = utils::head(sig[, .(tf, family, pct_target, pct_background, log2_enrichment, padj)], top),
    top_nominal = nominal,
    warnings = warnings,
    tables = list(motifs = table_path),
    figures = if (!is.null(fig)) list(motifs = fig[["png"]]) else list(),
    versions = res$versions
  )
}
