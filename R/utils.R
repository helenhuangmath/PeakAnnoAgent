`%||%` <- function(a, b) if (is.null(a)) b else a

msg <- function(...) {
  cat(sprintf(...), "\n")
}

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required. Install it before running this method.", pkg), call. = FALSE)
  }
}

split_csv <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "") {
    character()
  } else {
    trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  }
}

parse_cli_args <- function(args) {
  out <- list()
  pos <- character()
  i <- 1
  while (i <= length(args)) {
    token <- args[[i]]
    if (startsWith(token, "--")) {
      key <- sub("^--", "", token)
      if (i == length(args) || startsWith(args[[i + 1]], "--")) {
        out[[key]] <- TRUE
        i <- i + 1
      } else {
        out[[key]] <- args[[i + 1]]
        i <- i + 2
      }
    } else {
      pos <- c(pos, token)
      i <- i + 1
    }
  }
  out$command <- if (length(pos)) pos[[1]] else NULL
  out
}

arg_value <- function(args, key, default = NULL) {
  value <- args[[key]]
  if (is.null(value) || identical(value, "")) default else value
}

# Canonical method names accepted by each command, including aliases.
known_methods <- function(command) {
  switch(command,
    annotate = c("chipseeker", "nearest", "txdb", "chipeakanno", "chippeakanno",
      "annotatr", "rgreat", "homer", "webexport", "export"),
    enrich = c("clusterprofiler", "kegg", "clusterprofiler_kegg", "reactome",
      "reactomepa", "msigdb", "enrichr", "gprofiler", "gprofiler2", "webgestalt",
      "webgestaltr", "metascape", "david", "export", "genelist"),
    character()
  )
}

# Fail fast on typos so a mistyped --methods value never becomes a silent no-op.
validate_methods <- function(methods, command) {
  valid <- known_methods(command)
  unknown <- setdiff(methods, valid)
  if (length(unknown)) {
    stop(sprintf(
      "Unknown %s method(s): %s.\nValid methods: %s.",
      command, paste(unknown, collapse = ", "), paste(valid, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(methods)
}

read_peak_file <- function(path) {
  require_pkg("data.table")
  dt <- data.table::fread(path)
  lower <- tolower(names(dt))
  if (all(c("chr", "start", "end") %in% lower)) {
    names(dt)[match(c("chr", "start", "end"), lower)] <- c("Chr", "Start", "End")
  } else if (ncol(dt) >= 3) {
    names(dt)[1:3] <- c("Chr", "Start", "End")
  } else {
    stop("Peak file must contain at least three columns: Chr, Start, End.", call. = FALSE)
  }
  dt[, Start := as.integer(Start)]
  dt[, End := as.integer(End)]
  dt[, Chr := as.character(Chr)]
  if (!"PeakID" %in% names(dt)) {
    id_col <- intersect(c("name", "Name", "peak", "Peak", "V4"), names(dt))[1]
    ids <- if (!is.na(id_col)) as.character(dt[[id_col]]) else character()
    if (!length(ids) || anyNA(ids) || anyDuplicated(ids)) {
      ids <- sprintf("peak_%05d", seq_len(nrow(dt)))
    }
    dt[, PeakID := ids]
  }
  dt
}

read_gene_table <- function(path, gene_col = "SYMBOL", cluster_col = NULL) {
  require_pkg("data.table")
  dt <- data.table::fread(path)
  if (!gene_col %in% names(dt)) {
    stop(sprintf("Gene column '%s' was not found in %s.", gene_col, path), call. = FALSE)
  }
  dt <- dt[!is.na(get(gene_col)) & get(gene_col) != ""]
  if (!is.null(cluster_col) && cluster_col != "" && cluster_col %in% names(dt)) {
    groups <- split(dt[[gene_col]], as.character(dt[[cluster_col]]))
  } else {
    groups <- list(All = unique(dt[[gene_col]]))
  }
  lapply(groups, unique)
}

write_tsv <- function(x, path) {
  require_pkg("data.table")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(path)
}

# Genome builds the agent knows how to annotate. Each entry names the Bioconductor
# annotation packages plus the identifiers the external services expect.
genome_table <- function() {
  mouse <- list(orgdb = "org.Mm.eg.db", homer = "mm10", enrichr_suffix = "Mouse",
    kegg = "mmu", reactome = "mouse", gprofiler = "mmusculus", msigdb = "Mus musculus",
    organism = "mouse")
  human <- list(orgdb = "org.Hs.eg.db", homer = "hg38", enrichr_suffix = "Human",
    kegg = "hsa", reactome = "human", gprofiler = "hsapiens", msigdb = "Homo sapiens",
    organism = "human")
  list(
    mm10 = utils::modifyList(mouse, list(genome = "mm10",
      txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene", bsgenome = "BSgenome.Mmusculus.UCSC.mm10",
      annotatr = "mm10", chipatlas = "mm10")),
    mm39 = utils::modifyList(mouse, list(genome = "mm39", homer = "mm39",
      txdb = "TxDb.Mmusculus.UCSC.mm39.knownGene", bsgenome = "BSgenome.Mmusculus.UCSC.mm39",
      annotatr = "mm39", chipatlas = NA_character_)),
    hg38 = utils::modifyList(human, list(genome = "hg38",
      txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene", bsgenome = "BSgenome.Hsapiens.UCSC.hg38",
      annotatr = "hg38", chipatlas = "hg38")),
    hg19 = utils::modifyList(human, list(genome = "hg19", homer = "hg19",
      txdb = "TxDb.Hsapiens.UCSC.hg19.knownGene", bsgenome = "BSgenome.Hsapiens.UCSC.hg19",
      annotatr = "hg19", chipatlas = "hg19"))
  )
}

supported_genomes <- function() names(genome_table())

infer_species_defaults <- function(species) {
  key <- switch(tolower(species), mouse = "mm10", human = "hg38", tolower(species))
  defaults <- genome_table()[[key]]
  if (is.null(defaults)) {
    stop(sprintf("Supported species shortcuts: %s, mouse, human.",
      paste(supported_genomes(), collapse = ", ")), call. = FALSE)
  }
  defaults
}

load_orgdb <- function(species, orgdb_pkg = NULL) {
  defaults <- infer_species_defaults(species)
  orgdb_pkg <- if (is.null(orgdb_pkg)) defaults$orgdb else orgdb_pkg
  require_pkg(orgdb_pkg)
  get(orgdb_pkg, envir = asNamespace(orgdb_pkg))
}

load_txdb <- function(species, txdb_pkg = NULL) {
  defaults <- infer_species_defaults(species)
  txdb_pkg <- if (is.null(txdb_pkg)) defaults$txdb else txdb_pkg
  require_pkg(txdb_pkg)
  get(txdb_pkg, envir = asNamespace(txdb_pkg))
}

symbols_to_entrez <- function(genes, species = "mm10", orgdb_pkg = NULL) {
  require_pkg("AnnotationDbi")
  orgdb <- load_orgdb(species, orgdb_pkg)
  mapped <- AnnotationDbi::select(
    orgdb,
    keys = unique(genes),
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL"
  )
  unique(mapped$ENTREZID[!is.na(mapped$ENTREZID)])
}

load_bsgenome <- function(species, bsgenome_pkg = NULL) {
  bsgenome_pkg <- if (is.null(bsgenome_pkg)) infer_species_defaults(species)$bsgenome else bsgenome_pkg
  require_pkg(bsgenome_pkg)
  get(bsgenome_pkg, envir = asNamespace(bsgenome_pkg))
}

# PNG previews for the HTML report live in a hidden .preview/ folder so the results folder
# itself holds only the PDF figures.
preview_png_path <- function(out_dir, name) {
  dir <- file.path(out_dir, ".preview")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(name, ".png"))
}

# Write a ggplot or base-graphics figure as a Nature-style vector PDF (manuscript) plus a PNG
# preview for the HTML report.
# `draw` is either a ggplot object or a zero-argument function that draws on the
# current device. Returns the two paths.
save_figure <- function(draw, out_dir, name, width = 7, height = 4.5) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- preview_png_path(out_dir, name)
  render <- function() {
    if (inherits(draw, "ggplot")) {
      print(draw)
    } else {
      publication_par()
      draw()
    }
  }
  publication_pdf(pdf_path, width = width, height = height)
  tryCatch(render(), finally = grDevices::dev.off())
  grDevices::png(png_path, width = width, height = height, units = "in", res = 120, bg = "white")
  tryCatch(render(), finally = grDevices::dev.off())
  c(pdf = pdf_path, png = png_path)
}

# Method names such as "msigdb_C7:IMMUNESIGDB" -> safe file-name parts ("msigdb_C7-IMMUNESIGDB").
file_tag <- function(x) gsub("[^A-Za-z0-9_.-]+", "-", x)

# Versions of the packages a tool actually used, for the run manifest.
package_versions <- function(pkgs) {
  pkgs <- unique(pkgs[nzchar(pkgs)])
  installed <- pkgs[vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  stats::setNames(as.list(vapply(installed, function(p) as.character(utils::packageVersion(p)),
    character(1))), installed)
}
