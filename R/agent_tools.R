# Tool layer used by the agent (scripts/peak_tools.R). Each tool_* function takes
# a named parameter list, writes its tables/figures under params$out_dir, and
# returns a compact, JSON-serialisable summary that the agent reasons over. The
# full tables on disk stay the source of truth for figures and the report.

agent_tool_names <- function() {
  c("peak_qc", "annotate_peaks", "enrich_pathways", "great_enrichment", "find_motifs", "session_info")
}

param <- function(params, key, default = NULL) {
  value <- params[[key]]
  if (is.null(value) || (length(value) == 1 && (is.na(value) || identical(value, "")))) default else value
}

tool_peak_qc <- function(params) {
  peaks <- read_peak_file(param(params, "peaks"))
  run_peak_qc(peaks, param(params, "out_dir"), param(params, "prefix", "run"),
    blacklist = param(params, "blacklist"))
}

tool_annotate_peaks <- function(params) {
  peaks <- read_peak_file(param(params, "peaks"))
  out <- run_chippeakanno_pipeline(peaks, param(params, "out_dir"), param(params, "prefix", "run"),
    species = param(params, "genome", "mm10"),
    txdb_pkg = param(params, "txdb"), orgdb_pkg = param(params, "orgdb"),
    promoter_upstream = param(params, "promoter_upstream", 2000),
    promoter_downstream = param(params, "promoter_downstream", 500),
    distal_max_distance = param(params, "distal_max_distance", 100000))
  if (isTRUE(param(params, "chipseeker", FALSE)) && requireNamespace("ChIPseeker", quietly = TRUE)) {
    cs <- tryCatch(run_chipseeker_annotation(peaks, param(params, "out_dir"), param(params, "prefix", "run"),
      species = param(params, "genome", "mm10")), error = function(e) {
        out$warnings <<- c(out$warnings, paste("ChIPseeker cross-check failed:", conditionMessage(e)))
        NULL
      })
    if (!is.null(cs)) out$tables$chipseeker <- cs$annotation
  }
  out
}

msigdb_term2gene <- function(species, collection, subcollection = NULL) {
  require_pkg("msigdbr")
  species_name <- infer_species_defaults(species)$msigdb
  args <- list(species = species_name)
  if ("collection" %in% names(formals(msigdbr::msigdbr))) {
    args$collection <- collection
    if (!is.null(subcollection)) args$subcollection <- subcollection
  } else {
    args$category <- collection
    if (!is.null(subcollection)) args$subcategory <- subcollection
  }
  msig <- do.call(msigdbr::msigdbr, args)
  unique(as.data.frame(msig)[, c("gs_name", "gene_symbol")])
}

read_gmt <- function(path) {
  lines <- strsplit(readLines(path, warn = FALSE), "\t", fixed = TRUE)
  data.table::rbindlist(lapply(lines, function(x) {
    if (length(x) < 3) return(NULL)
    data.table::data.table(term = x[[1]], gene = unique(x[-(1:2)]))
  }))
}

# Standardise any compareCluster-style result into the report's columns.
standardise_enrichment <- function(res, method) {
  if (is.null(res) || !nrow(res)) return(data.table::data.table())
  res <- data.table::as.data.table(res)
  data.table::data.table(
    method = method,
    group = as.character(res$Cluster),
    term_id = as.character(res$ID),
    term = sub(" - [A-Z][a-z]+ [a-z]+ \\([^)]*\\)$", "", as.character(res$Description)),
    gene_ratio = as.character(res$GeneRatio),
    bg_ratio = as.character(res$BgRatio),
    pvalue = as.numeric(res$pvalue),
    padj = as.numeric(res$p.adjust),
    count = as.integer(res$Count),
    genes = gsub("/", ";", as.character(res$geneID), fixed = TRUE)
  )
}

# Hypergeometric test of small user signatures against a genome-wide background.
# clusterProfiler::enricher would shrink the background to genes present in the
# GMT and drop sets under 10 genes, which is wrong for a handful of signatures.
signature_hypergeom <- function(gene_lists, term2gene, universe, min_size = 3) {
  universe <- unique(universe)
  sets <- lapply(split(term2gene$gene, term2gene$term), function(x) intersect(unique(x), universe))
  sets <- sets[lengths(sets) >= min_size]
  rows <- list()
  for (g in names(gene_lists)) {
    genes <- intersect(unique(gene_lists[[g]]), universe)
    for (term in names(sets)) {
      hit <- intersect(genes, sets[[term]])
      p <- stats::phyper(length(hit) - 1, length(sets[[term]]), length(universe) - length(sets[[term]]),
        length(genes), lower.tail = FALSE)
      rows[[length(rows) + 1]] <- data.table::data.table(Cluster = g, ID = term, Description = term,
        GeneRatio = sprintf("%d/%d", length(hit), length(genes)),
        BgRatio = sprintf("%d/%d", length(sets[[term]]), length(universe)),
        pvalue = p, geneID = paste(hit, collapse = "/"), Count = length(hit))
    }
  }
  res <- data.table::rbindlist(rows)
  if (nrow(res)) res[, p.adjust := stats::p.adjust(pvalue, "BH"), by = Cluster]
  res
}

compare_enricher <- function(gene_lists, term2gene, universe, padj_cutoff) {
  cc <- clusterProfiler::compareCluster(geneClusters = gene_lists, fun = clusterProfiler::enricher,
    TERM2GENE = term2gene, universe = universe, pAdjustMethod = "BH",
    pvalueCutoff = padj_cutoff, qvalueCutoff = 1)
  if (is.null(cc)) NULL else cc@compareClusterResult
}

tool_enrich_pathways <- function(params) {
  require_pkg("clusterProfiler")
  require_pkg("data.table")
  out_dir <- param(params, "out_dir")
  prefix <- param(params, "prefix", "run")
  species <- param(params, "genome", "mm10")
  padj_cutoff <- as.numeric(param(params, "padj", 0.05))
  top <- as.integer(param(params, "top", 8))
  methods <- unlist(param(params, "methods", c("go_bp")))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  gene_lists <- read_gene_table(param(params, "gene_sets"), gene_col = "SYMBOL", cluster_col = "Cluster")
  groups <- unlist(param(params, "groups", names(gene_lists)))
  gene_lists <- gene_lists[intersect(groups, names(gene_lists))]
  gene_lists <- gene_lists[lengths(gene_lists) >= 5]
  if (!length(gene_lists)) stop("No gene group has at least 5 genes to test.", call. = FALSE)
  universe <- NULL
  if (!is.null(param(params, "universe"))) {
    universe <- unique(unlist(read_gene_table(param(params, "universe"), gene_col = "SYMBOL")))
  }
  warnings <- character()
  if (is.null(universe)) warnings <- c(warnings, "No gene universe given: ORA uses the full annotation background.")
  orgdb <- load_orgdb(species)

  run_one <- function(method) {
    if (method %in% c("go_bp", "go_mf", "go_cc")) {
      ont <- toupper(sub("go_", "", method))
      cc <- clusterProfiler::compareCluster(geneClusters = gene_lists, fun = clusterProfiler::enrichGO, OrgDb = orgdb,
        keyType = "SYMBOL", ont = ont, universe = universe, pAdjustMethod = "BH",
        pvalueCutoff = padj_cutoff, qvalueCutoff = 1)
      return(if (is.null(cc)) NULL else cc@compareClusterResult)
    }
    if (method %in% c("kegg", "reactome")) {
      entrez <- lapply(gene_lists, symbols_to_entrez, species = species)
      uni_entrez <- if (!is.null(universe)) symbols_to_entrez(universe, species) else NULL
      defaults <- infer_species_defaults(species)
      cc <- if (method == "kegg") {
        clusterProfiler::compareCluster(geneClusters = entrez, fun = clusterProfiler::enrichKEGG,
          organism = defaults$kegg, universe = uni_entrez, pAdjustMethod = "BH", pvalueCutoff = padj_cutoff)
      } else {
        require_pkg("ReactomePA")
        clusterProfiler::compareCluster(geneClusters = entrez, fun = ReactomePA::enrichPathway,
          organism = defaults$reactome, universe = uni_entrez, pAdjustMethod = "BH",
          pvalueCutoff = padj_cutoff, readable = TRUE)
      }
      if (is.null(cc)) return(NULL)
      res <- cc@compareClusterResult
      if (method == "kegg") {
        sym <- AnnotationDbi::mapIds(orgdb, keys = unique(unlist(strsplit(res$geneID, "/"))),
          column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
        res$geneID <- vapply(strsplit(res$geneID, "/"), function(x) paste(stats::na.omit(sym[x]), collapse = "/"),
          character(1))
      }
      return(res)
    }
    if (startsWith(method, "msigdb_")) {
      parts <- strsplit(sub("^msigdb_", "", method), ":", fixed = TRUE)[[1]]
      t2g <- msigdb_term2gene(species, toupper(parts[[1]]), if (length(parts) > 1) parts[[2]] else NULL)
      return(compare_enricher(gene_lists, t2g, universe, padj_cutoff))
    }
    if (method == "custom") {
      gmt <- param(params, "custom_gmt")
      if (is.null(gmt)) stop("method 'custom' needs custom_gmt.", call. = FALSE)
      bg <- universe %||% AnnotationDbi::keys(orgdb, keytype = "SYMBOL")
      return(signature_hypergeom(gene_lists, read_gmt(gmt), bg))
    }
    stop(sprintf("Unknown enrichment method '%s'.", method), call. = FALSE)
  }

  all_res <- list()
  for (m in methods) {
    res <- tryCatch(run_one(m), error = function(e) {
      warnings <<- c(warnings, sprintf("%s failed: %s", m, conditionMessage(e)))
      NULL
    })
    all_res[[m]] <- standardise_enrichment(res, m)
  }
  combined <- data.table::rbindlist(all_res, fill = TRUE)
  table_path <- file.path(out_dir, paste0(prefix, "_enrichment.tsv"))
  write_tsv(combined, table_path)

  figures <- list()
  if (nrow(combined)) {
    for (m in unique(combined$method)) {
      sub <- combined[method == m & padj <= padj_cutoff]
      if (!nrow(sub)) next
      sub <- sub[, .(Cluster = group, Term = term, Adjusted.P.value = padj, Count = count)]
      name <- paste0(prefix, "_enrichment_", file_tag(m), "_dotplot")
      pdf <- file.path(out_dir, paste0(name, ".pdf"))
      png <- preview_png_path(out_dir, name)
      ok <- tryCatch({
        plot_enrichment_dotplot(sub, pdf, padj_cutoff = padj_cutoff, top_show = top,
          title = sprintf("Pathway enrichment: %s", m))
        if (requireNamespace("ggplot2", quietly = TRUE)) {
          plot_enrichment_dotplot(sub, png, padj_cutoff = padj_cutoff, top_show = top,
            title = sprintf("Pathway enrichment: %s", m))
        }
        TRUE
      }, error = function(e) FALSE)
      if (ok && file.exists(png)) figures[[m]] <- png
    }
  }

  networks <- if (nrow(combined)) {
    run_pathway_networks(combined, out_dir, prefix, padj_cutoff = padj_cutoff)
  } else {
    list(figures = list(), tables = list(), modules = list(), hub_genes = list(), warnings = character())
  }
  figures <- c(figures, networks$figures)
  warnings <- c(warnings, networks$warnings)

  top_terms <- if (nrow(combined)) {
    sig <- combined[padj <= padj_cutoff]
    sig[order(method, group, padj), utils::head(.SD, top), by = .(method, group)][,
      .(method, group, term, padj = signif(padj, 3), count,
        genes = vapply(strsplit(genes, ";"), function(g) paste(utils::head(g, 12), collapse = ";"), character(1)))]
  } else data.table::data.table()

  list(
    methods = methods,
    groups = as.list(lengths(gene_lists)),
    universe_size = length(universe),
    n_significant = if (nrow(combined)) combined[padj <= padj_cutoff, .N, by = .(method, group)] else list(),
    top_terms = top_terms,
    pathway_modules = networks$modules,
    hub_genes = networks$hub_genes,
    warnings = warnings,
    tables = c(list(enrichment = table_path), networks$tables),
    figures = figures,
    versions = package_versions(c("clusterProfiler", "msigdbr", "ReactomePA", "igraph", "ggrepel",
      infer_species_defaults(species)$orgdb))
  )
}

tool_great_enrichment <- function(params) {
  require_pkg("rGREAT")
  out_dir <- param(params, "out_dir")
  prefix <- param(params, "prefix", "run")
  species <- param(params, "genome", "mm10")
  gene_sets <- param(params, "gene_sets", "GO:BP")
  peaks <- read_peak_file(param(params, "peaks"))
  gr <- peaks_to_granges(peaks)
  res <- rGREAT::great(gr, gene_sets = gene_sets, tss_source = infer_species_defaults(species)$txdb,
    verbose = FALSE)
  tab <- data.table::as.data.table(rGREAT::getEnrichmentTable(res))
  table_path <- file.path(out_dir, paste0(prefix, "_great_", gsub("[^A-Za-z0-9]", "_", gene_sets), ".tsv"))
  write_tsv(tab, table_path)
  padj_cutoff <- as.numeric(param(params, "padj", 0.05))
  sig <- tab[p_adjust <= padj_cutoff][order(p_adjust)]
  list(
    gene_sets = gene_sets,
    n_significant = nrow(sig),
    top_terms = utils::head(sig[, .(id, description, fold_enrichment = round(fold_enrichment, 2),
      p_adjust = signif(p_adjust, 3), observed_region_hits)], as.integer(param(params, "top", 15))),
    warnings = if (!nrow(sig)) "No GREAT term passed the cutoff." else character(),
    tables = list(great = table_path),
    figures = list(),
    versions = package_versions(c("rGREAT"))
  )
}

tool_find_motifs <- function(params) {
  peaks <- read_peak_file(param(params, "peaks"))
  run_motif_enrichment(peaks, param(params, "out_dir"), param(params, "prefix", "run"),
    species = param(params, "genome", "mm10"), backend = param(params, "backend", "auto"),
    size = as.integer(param(params, "size", 200)), max_peaks = as.integer(param(params, "max_peaks", 5000)),
    padj_cutoff = as.numeric(param(params, "padj", 0.05)), top = as.integer(param(params, "top", 25)))
}

tool_session_info <- function(params) {
  genome <- param(params, "genome", "mm10")
  d <- infer_species_defaults(genome)
  pkgs <- c("ChIPpeakAnno", "ChIPseeker", "clusterProfiler", "ReactomePA", "msigdbr", "rGREAT", "monaLisa",
    "TFBSTools", "JASPAR2024", "JASPAR2022", "JASPAR2020", d$txdb, d$orgdb, d$bsgenome, "data.table", "jsonlite")
  list(
    r_version = R.version.string,
    genome = d$genome,
    packages = package_versions(pkgs),
    missing = setdiff(pkgs[!grepl("^JASPAR20(22|20)$", pkgs)], names(package_versions(pkgs))),
    homer = nzchar(Sys.which("findMotifsGenome.pl"))
  )
}

run_agent_tool <- function(name, params) {
  if (!name %in% agent_tool_names()) {
    stop(sprintf("Unknown tool '%s'. Valid: %s", name, paste(agent_tool_names(), collapse = ", ")), call. = FALSE)
  }
  get(paste0("tool_", name), mode = "function")(params)
}
