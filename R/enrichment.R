default_enrichr_libraries <- function(species = "mouse") {
  suffix <- infer_species_defaults(species)$enrichr_suffix
  if (suffix == "Mouse") {
    c("GO_Biological_Process_2021", "KEGG_2019_Mouse", "WikiPathways_2019_Mouse",
      "Reactome_2022", "MSigDB_Hallmark_2020", "TRRUST_Transcription_Factors_2019")
  } else {
    c("GO_Biological_Process_2021", "KEGG_2021_Human", "WikiPathways_2021_Human",
      "Reactome_2022", "MSigDB_Hallmark_2020", "TRRUST_Transcription_Factors_2019")
  }
}

run_clusterprofiler_go <- function(gene_lists, out_dir, prefix, species = "mm10",
                                   universe = NULL, ont = "BP", padj_cutoff = 0.05,
                                   top_show = 10) {
  require_pkg("clusterProfiler")
  require_pkg("enrichplot")
  require_pkg("ggplot2")
  require_pkg("data.table")
  orgdb_pkg <- infer_species_defaults(species)$orgdb
  require_pkg(orgdb_pkg)
  orgdb <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))

  cc <- clusterProfiler::compareCluster(
    geneClusters = gene_lists,
    fun = "enrichGO",
    OrgDb = orgdb,
    keyType = "SYMBOL",
    ont = ont,
    universe = universe,
    pAdjustMethod = "BH",
    pvalueCutoff = as.numeric(padj_cutoff),
    qvalueCutoff = 0.1,
    readable = FALSE
  )
  cc <- enrichplot::pairwise_termsim(cc)
  res <- data.table::as.data.table(cc@compareClusterResult)
  write_tsv(res, file.path(out_dir, paste0(prefix, "_clusterprofiler_GO_", ont, ".tsv")))

  if (nrow(res) > 0) {
    p <- enrichplot::dotplot(cc, showCategory = as.integer(top_show), label_format = 55) +
      scale_color_padj_publication() +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggsave_publication(file.path(out_dir, paste0(prefix, "_clusterprofiler_GO_", ont, "_dotplot.pdf")),
      p, width = max(10, length(gene_lists) * 1.1), height = 8.5)
  }
  invisible(res)
}

run_clusterprofiler_kegg <- function(gene_lists, out_dir, prefix, species = "mm10",
                                     padj_cutoff = 0.05, top_show = 10) {
  require_pkg("clusterProfiler")
  require_pkg("enrichplot")
  require_pkg("ggplot2")
  require_pkg("data.table")
  organism <- infer_species_defaults(species)$kegg
  entrez_lists <- lapply(gene_lists, symbols_to_entrez, species = species)
  cc <- clusterProfiler::compareCluster(
    geneClusters = entrez_lists,
    fun = "enrichKEGG",
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = as.numeric(padj_cutoff)
  )
  res <- data.table::as.data.table(cc@compareClusterResult)
  write_tsv(res, file.path(out_dir, paste0(prefix, "_clusterprofiler_KEGG.tsv")))
  if (nrow(res) > 0) {
    p <- enrichplot::dotplot(cc, showCategory = as.integer(top_show), label_format = 55) +
      scale_color_padj_publication() +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggsave_publication(file.path(out_dir, paste0(prefix, "_clusterprofiler_KEGG_dotplot.pdf")),
      p, width = max(10, length(gene_lists) * 1.1), height = 8.5)
  }
  invisible(res)
}

run_reactomepa <- function(gene_lists, out_dir, prefix, species = "mm10",
                           padj_cutoff = 0.05, top_show = 10) {
  require_pkg("ReactomePA")
  require_pkg("enrichplot")
  require_pkg("ggplot2")
  require_pkg("data.table")
  organism <- infer_species_defaults(species)$reactome
  entrez_lists <- lapply(gene_lists, symbols_to_entrez, species = species)
  cc <- clusterProfiler::compareCluster(
    geneClusters = entrez_lists,
    fun = ReactomePA::enrichPathway,
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = as.numeric(padj_cutoff),
    readable = TRUE
  )
  res <- data.table::as.data.table(cc@compareClusterResult)
  write_tsv(res, file.path(out_dir, paste0(prefix, "_reactomepa.tsv")))
  if (nrow(res) > 0) {
    p <- enrichplot::dotplot(cc, showCategory = as.integer(top_show), label_format = 55) +
      scale_color_padj_publication() +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggsave_publication(file.path(out_dir, paste0(prefix, "_reactomepa_dotplot.pdf")),
      p, width = max(10, length(gene_lists) * 1.1), height = 8.5)
  }
  invisible(res)
}

run_gprofiler2 <- function(gene_lists, out_dir, prefix, species = "mm10",
                           padj_cutoff = 0.05, top_show = 10) {
  require_pkg("gprofiler2")
  require_pkg("data.table")
  organism <- infer_species_defaults(species)$gprofiler
  all_res <- data.table::rbindlist(lapply(names(gene_lists), function(group) {
    gost <- gprofiler2::gost(query = unique(gene_lists[[group]]), organism = organism,
      correction_method = "fdr", significant = FALSE)
    if (is.null(gost$result) || !nrow(gost$result)) return(NULL)
    d <- data.table::as.data.table(gost$result)
    d[, Cluster := group]
    rename_old <- intersect(c("term_name", "p_value", "intersection"), names(d))
    data.table::setnames(d, old = rename_old,
      new = c("Term", "Adjusted.P.value", "Genes")[seq_along(rename_old)],
      skip_absent = TRUE)
    if ("Genes" %in% names(d) && is.list(d$Genes)) {
      d[, Genes := vapply(Genes, paste, collapse = ";", character(1))]
    }
    d
  }), fill = TRUE)
  write_tsv(all_res, file.path(out_dir, paste0(prefix, "_gprofiler2.tsv")))
  plot_enrichment_dotplot(all_res, file.path(out_dir, paste0(prefix, "_gprofiler2_dotplot.pdf")),
    padj_cutoff = as.numeric(padj_cutoff), top_show = as.integer(top_show))
  invisible(all_res)
}

run_msigdb_enrichment <- function(gene_lists, out_dir, prefix, species = "mm10",
                                  collection = "H", universe = NULL,
                                  padj_cutoff = 0.05, top_show = 10) {
  require_pkg("msigdbr")
  require_pkg("clusterProfiler")
  require_pkg("enrichplot")
  require_pkg("ggplot2")
  require_pkg("data.table")
  species_name <- infer_species_defaults(species)$msigdb
  msig <- msigdbr::msigdbr(species = species_name, category = collection)
  term2gene <- unique(msig[, c("gs_name", "gene_symbol")])
  cc <- clusterProfiler::compareCluster(
    geneClusters = gene_lists,
    fun = "enricher",
    TERM2GENE = term2gene,
    universe = universe,
    pAdjustMethod = "BH",
    pvalueCutoff = as.numeric(padj_cutoff)
  )
  res <- data.table::as.data.table(cc@compareClusterResult)
  write_tsv(res, file.path(out_dir, paste0(prefix, "_msigdb_", collection, ".tsv")))
  if (nrow(res) > 0) {
    p <- enrichplot::dotplot(cc, showCategory = as.integer(top_show), label_format = 55) +
      scale_color_padj_publication() +
      theme_publication() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggsave_publication(file.path(out_dir, paste0(prefix, "_msigdb_", collection, "_dotplot.pdf")),
      p, width = max(10, length(gene_lists) * 1.1), height = 8.5)
  }
  invisible(res)
}

run_webgestaltr <- function(gene_lists, out_dir, prefix, species = "mm10",
                            database = "geneontology_Biological_Process",
                            padj_cutoff = 0.05, top_show = 10) {
  require_pkg("WebGestaltR")
  require_pkg("data.table")
  organism <- if (infer_species_defaults(species)$kegg == "mmu") "mmusculus" else "hsapiens"
  all_res <- data.table::rbindlist(lapply(names(gene_lists), function(group) {
    out <- WebGestaltR::WebGestaltR(
      enrichMethod = "ORA",
      organism = organism,
      enrichDatabase = database,
      interestGene = unique(gene_lists[[group]]),
      interestGeneType = "genesymbol",
      isOutput = FALSE
    )
    if (is.null(out) || !nrow(out)) return(NULL)
    d <- data.table::as.data.table(out)
    d[, Cluster := group]
    rename_old <- intersect(c("description", "FDR", "userId"), names(d))
    data.table::setnames(d, old = rename_old,
      new = c("Term", "Adjusted.P.value", "Genes")[seq_along(rename_old)],
      skip_absent = TRUE)
    if ("Genes" %in% names(d) && is.list(d$Genes)) {
      d[, Genes := vapply(Genes, paste, collapse = ";", character(1))]
    }
    d
  }), fill = TRUE)
  write_tsv(all_res, file.path(out_dir, paste0(prefix, "_webgestaltr.tsv")))
  plot_enrichment_dotplot(all_res, file.path(out_dir, paste0(prefix, "_webgestaltr_dotplot.pdf")),
    padj_cutoff = as.numeric(padj_cutoff), top_show = as.integer(top_show))
  invisible(all_res)
}

run_enrichr <- function(gene_lists, out_dir, prefix, species = "mouse",
                        libraries = NULL, padj_cutoff = 0.05, top_show = 10) {
  require_pkg("enrichR")
  require_pkg("data.table")
  require_pkg("ggplot2")
  if (is.null(libraries) || length(libraries) == 0) {
    libraries <- default_enrichr_libraries(species)
  }
  enrichR::setEnrichrSite("Enrichr")
  all_res <- data.table::rbindlist(lapply(names(gene_lists), function(group) {
    genes <- unique(toupper(gene_lists[[group]]))
    e <- tryCatch(enrichR::enrichr(genes, libraries), error = function(err) {
      warning(sprintf("Enrichr failed for %s: %s", group, conditionMessage(err)))
      NULL
    })
    if (is.null(e)) return(NULL)
    data.table::rbindlist(lapply(names(e), function(db) {
      d <- data.table::as.data.table(e[[db]])
      if (!nrow(d)) return(NULL)
      d[, .(Cluster = group, Library = db, Term, Overlap, P.value, Adjusted.P.value, Combined.Score, Genes)]
    }), fill = TRUE)
  }), fill = TRUE)

  out_table <- file.path(out_dir, paste0(prefix, "_enrichr.tsv"))
  write_tsv(all_res, out_table)
  plot_enrichment_dotplot(all_res, file.path(out_dir, paste0(prefix, "_enrichr_dotplot.pdf")),
    padj_cutoff = as.numeric(padj_cutoff), top_show = as.integer(top_show))
  invisible(all_res)
}

export_metascape_input <- function(gene_lists, out_dir, prefix) {
  require_pkg("data.table")
  dt <- data.table::rbindlist(lapply(names(gene_lists), function(group) {
    data.table::data.table(Group = group, Gene = gene_lists[[group]])
  }))
  out_table <- file.path(out_dir, paste0(prefix, "_metascape_input.tsv"))
  write_tsv(dt, out_table)
  invisible(out_table)
}

export_gene_lists <- function(gene_lists, out_dir, prefix) {
  require_pkg("data.table")
  long <- data.table::rbindlist(lapply(names(gene_lists), function(group) {
    data.table::data.table(Group = group, Gene = gene_lists[[group]])
  }))
  wide_path <- file.path(out_dir, paste0(prefix, "_gene_lists_wide.tsv"))
  long_path <- file.path(out_dir, paste0(prefix, "_gene_lists_long.tsv"))
  max_len <- max(lengths(gene_lists))
  wide <- data.table::as.data.table(lapply(gene_lists, function(x) c(x, rep(NA_character_, max_len - length(x)))))
  write_tsv(long, long_path)
  write_tsv(wide, wide_path)
  invisible(list(long = long_path, wide = wide_path))
}

run_david_enrichment <- function(gene_lists, out_dir, prefix, species = "mouse",
                                 email = Sys.getenv("DAVID_EMAIL")) {
  require_pkg("RDAVIDWebService")
  require_pkg("data.table")
  if (email == "") {
    stop("Set DAVID_EMAIL or pass --david-email before running the DAVID backend.", call. = FALSE)
  }
  david <- RDAVIDWebService::DAVIDWebService$new(email = email)
  all_res <- data.table::rbindlist(lapply(names(gene_lists), function(group) {
    RDAVIDWebService::addList(david, gene_lists[[group]], idType = "OFFICIAL_GENE_SYMBOL",
      listName = group, listType = "Gene")
    chart <- RDAVIDWebService::getFunctionalAnnotationChart(david)
    d <- data.table::as.data.table(chart)
    if (!nrow(d)) return(NULL)
    d[, Cluster := group]
    d
  }), fill = TRUE)
  write_tsv(all_res, file.path(out_dir, paste0(prefix, "_david.tsv")))
  invisible(all_res)
}
