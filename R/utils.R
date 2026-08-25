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
    groups <- split(dt[[gene_col]], dt[[cluster_col]])
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

infer_species_defaults <- function(species) {
  switch(tolower(species),
    mm10 = list(orgdb = "org.Mm.eg.db", txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene", homer = "mm10", enrichr_suffix = "Mouse", kegg = "mmu", reactome = "mouse", gprofiler = "mmusculus", msigdb = "Mus musculus", annotatr = "mm10"),
    mouse = list(orgdb = "org.Mm.eg.db", txdb = "TxDb.Mmusculus.UCSC.mm10.knownGene", homer = "mm10", enrichr_suffix = "Mouse", kegg = "mmu", reactome = "mouse", gprofiler = "mmusculus", msigdb = "Mus musculus", annotatr = "mm10"),
    hg38 = list(orgdb = "org.Hs.eg.db", txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene", homer = "hg38", enrichr_suffix = "Human", kegg = "hsa", reactome = "human", gprofiler = "hsapiens", msigdb = "Homo sapiens", annotatr = "hg38"),
    human = list(orgdb = "org.Hs.eg.db", txdb = "TxDb.Hsapiens.UCSC.hg38.knownGene", homer = "hg38", enrichr_suffix = "Human", kegg = "hsa", reactome = "human", gprofiler = "hsapiens", msigdb = "Homo sapiens", annotatr = "hg38"),
    stop("Supported species shortcuts: mm10, mouse, hg38, human.", call. = FALSE)
  )
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
