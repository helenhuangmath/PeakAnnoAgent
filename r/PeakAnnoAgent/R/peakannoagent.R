# R interface to the peakannoagent Python program. Every function builds a command line for
# `peakanno-agent`, runs it, and reads the files it writes. No analysis happens in R.

#' Locate the peakanno-agent command
#'
#' Looked up in this order: the R option `peakannoagent.command`, the environment variable
#' `PEAKANNOAGENT_COMMAND`, `peakanno-agent` on the PATH, and `python3 -m peakannoagent`.
paa_command <- function() {
  cmd <- getOption("peakannoagent.command")
  if (is.null(cmd) && nzchar(Sys.getenv("PEAKANNOAGENT_COMMAND"))) {
    cmd <- strsplit(Sys.getenv("PEAKANNOAGENT_COMMAND"), " +")[[1]]
  }
  if (is.null(cmd)) {
    exe <- Sys.which("peakanno-agent")
    cmd <- if (nzchar(exe)) unname(exe) else c(unname(Sys.which("python3")), "-m", "peakannoagent")
  }
  if (!nzchar(cmd[[1]])) {
    stop("peakanno-agent was not found. Install it with `pip install peakannoagent`, or set ",
      "options(peakannoagent.command = ...).", call. = FALSE)
  }
  cmd
}

paa_call <- function(args, echo = FALSE) {
  cmd <- paa_command()
  err <- tempfile("paa_stderr_")
  on.exit(unlink(err), add = TRUE)
  out <- suppressWarnings(system2(cmd[[1]], c(cmd[-1], shQuote(args)), stdout = TRUE,
    stderr = if (echo) "" else err))
  status <- attr(out, "status")
  status <- if (is.null(status)) 0L else status
  messages <- if (!echo && file.exists(err)) readLines(err, warn = FALSE) else character()
  list(status = status, stdout = as.character(out), stderr = messages)
}

paa_fail <- function(res, what) {
  tail <- utils::tail(c(res$stderr, res$stdout), 15)
  stop(sprintf("%s failed (exit status %d):\n%s", what, res$status, paste(tail, collapse = "\n")), call. = FALSE)
}

read_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  utils::read.delim(path, check.names = FALSE, na.strings = "NA", stringsAsFactors = FALSE, quote = "")
}

# A peak file path, or a data frame with Chr/Start/End (or chromosome, start, end as the first
# three columns) written to a temporary BED file. Coordinates are zero-based and half-open.
peaks_path <- function(peaks) {
  if (is.character(peaks) && length(peaks) == 1) {
    if (!file.exists(peaks)) stop("Peak file not found: ", peaks, call. = FALSE)
    return(normalizePath(peaks))
  }
  if (!is.data.frame(peaks) || ncol(peaks) < 3) {
    stop("`peaks` must be a file path or a data frame with chromosome, start and end columns.", call. = FALSE)
  }
  cols <- match(c("chr", "start", "end"), tolower(names(peaks)))
  if (anyNA(cols)) cols <- 1:3
  bed <- data.frame(chr = peaks[[cols[1]]], start = format(as.numeric(peaks[[cols[2]]]), scientific = FALSE, trim = TRUE),
    end = format(as.numeric(peaks[[cols[3]]]), scientific = FALSE, trim = TRUE))
  id <- match("peakid", tolower(names(peaks)))
  if (!is.na(id)) bed$id <- peaks[[id]]
  path <- tempfile(fileext = ".bed")
  utils::write.table(bed, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  path
}

#' Annotate peaks
#'
#' Runs `peakanno-agent annotate` and returns its tables.
annotate_peaks <- function(peaks, genome, out_dir = tempfile("peakanno_"), prefix = "run", reference = "gencode",
                           gtf = NULL, promoter_upstream = 2000, promoter_downstream = 500,
                           distal_max_distance = 100000, gene_types = c("auto", "protein_coding", "all"),
                           figures = TRUE, palette = "chippeakanno") {
  gene_types <- match.arg(gene_types)
  args <- c("annotate", "--peaks", peaks_path(peaks), "--genome", genome, "--out", out_dir, "--prefix", prefix,
    "--reference", reference, "--promoter-upstream", as.integer(promoter_upstream),
    "--promoter-downstream", as.integer(promoter_downstream),
    "--distal-max-distance", format(distal_max_distance, scientific = FALSE), "--gene-types", gene_types, "--palette", palette)
  if (!is.null(gtf)) args <- c(args, "--gtf", normalizePath(gtf))
  if (!isTRUE(figures)) args <- c(args, "--no-figures")
  res <- paa_call(args)
  if (res$status != 0) paa_fail(res, "peakanno-agent annotate")
  table <- function(name) read_table(file.path(out_dir, paste0(prefix, "_", name, ".tsv")))
  summary <- jsonlite::fromJSON(file.path(out_dir, paste0(prefix, "_annotation.json")), simplifyVector = TRUE)
  pdfs <- if (length(summary$figures)) {
    stats <- names(summary$figures)
    structure(file.path(out_dir, paste0(prefix, "_", stats, ".pdf")), names = stats)
  } else character()
  list(annotation = table("peak_genes"), feature_distribution = table("feature_distribution"),
    region_distribution = table("region_distribution"), tss_distance_bins = table("tss_distance_bins"),
    gene_sets = table("peak_gene_sets"), universe = table("gene_universe")$SYMBOL,
    figures = pdfs, summary = summary, out_dir = normalizePath(out_dir))
}

#' Run the full analysis
#'
#' Runs `peakanno-agent run`: QC, annotation, pathway enrichment, motifs, public-data lookups
#' and the HTML report.
run_analysis <- function(peaks, genome, out_dir, prefix = "run", llm = FALSE, offline = FALSE, reference = "gencode",
                         gtf = NULL, blacklist = NULL, assay = NULL, target = NULL, cell_type = NULL,
                         condition = NULL, biology = NULL, expected_genes = NULL, expected_tfs = NULL, gmt = NULL,
                         rnaseq = NULL, palette = "chippeakanno", echo = interactive()) {
  args <- c("run", "--peaks", peaks_path(peaks), "--genome", genome, "--out", out_dir, "--prefix", prefix,
    "--reference", reference, "--palette", palette)
  if (!isTRUE(llm)) args <- c(args, "--no-llm")
  if (isTRUE(offline)) args <- c(args, "--offline")
  add <- function(flag, value, file = FALSE) {
    if (is.null(value)) return(invisible())
    args <<- c(args, flag, if (file) normalizePath(value) else paste(value, collapse = ","))
  }
  add("--gtf", gtf, TRUE); add("--blacklist", blacklist, TRUE); add("--gmt", gmt, TRUE); add("--rnaseq", rnaseq, TRUE)
  add("--assay", assay); add("--target", target); add("--cell-type", cell_type); add("--condition", condition)
  add("--biology", biology); add("--expected-genes", expected_genes); add("--expected-tfs", expected_tfs)
  res <- paa_call(args, echo = echo)
  # 0 = success, 3 = partial (some steps failed or were unavailable); anything else is an error.
  if (!res$status %in% c(0, 3)) paa_fail(res, "peakanno-agent run")
  manifest <- jsonlite::fromJSON(file.path(out_dir, paste0(prefix, "_manifest.json")), simplifyVector = FALSE)
  results <- file.path(out_dir, "results")
  list(status = manifest$status, reasons = unlist(manifest$reasons),
    report = normalizePath(file.path(out_dir, paste0(prefix, "_report.html"))),
    annotation = read_table(file.path(results, paste0(prefix, "_peak_genes.tsv"))),
    enrichment = read_table(file.path(results, paste0(prefix, "_enrich1_enrichment.tsv"))),
    tables = list.files(results, pattern = "\\.tsv$", full.names = TRUE),
    figures = list.files(results, pattern = "\\.pdf$", full.names = TRUE),
    manifest = manifest, out_dir = normalizePath(out_dir))
}

#' Check the installation
paa_check <- function(offline = TRUE) {
  res <- paa_call(c("check", if (isTRUE(offline)) "--offline"))
  cat(res$stdout, sep = "\n")
  invisible(res$status == 0)
}

#' Version of the peakannoagent program
paa_version <- function() {
  res <- paa_call("--version")
  if (res$status != 0) paa_fail(res, "peakanno-agent --version")
  sub("^.* ", "", res$stdout[[length(res$stdout)]])
}
