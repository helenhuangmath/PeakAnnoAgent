# Pathway networks built from standardised enrichment tables
# (columns: method, group, term, padj, count, genes).
#
# Enrichment map: significant terms are nodes, edges join terms whose gene sets
#   overlap (Jaccard), Louvain communities group redundant terms into modules.
# Gene-pathway network: the top terms linked to their genes; genes shared by
#   several terms are the hubs.
# Node/edge tables are written for Cytoscape alongside the figures.

split_genes <- function(x) {
  unique(trimws(unlist(strsplit(x, "[;/,]"))))
}

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (u == 0) 0 else length(intersect(a, b)) / u
}

# MSigDB-style names (HALLMARK_TNFA_SIGNALING_VIA_NFKB) -> "TNFA signaling via NFKB".
pretty_term <- function(x) {
  msig <- grepl("^[A-Z0-9]+(_[A-Z0-9]+)+$", x)
  y <- x
  y[msig] <- sub("^(HALLMARK|KEGG|REACTOME|WP|BIOCARTA|GOBP|GOMF|GOCC)_", "", x[msig])
  y[msig] <- gsub("_", " ", y[msig])
  y[msig] <- vapply(y[msig], function(s) {
    w <- strsplit(tolower(s), " ")[[1]]
    lower_words <- c("via", "and", "of", "in", "to", "the", "by", "for", "up", "dn", "vs", "or", "on", "at")
    keep_upper <- (nchar(w) <= 4 | grepl("[0-9]", w)) & !w %in% lower_words   # acronyms, GSE12345, CD8
    w[keep_upper] <- toupper(w[keep_upper])
    w[1] <- paste0(toupper(substr(w[1], 1, 1)), substr(w[1], 2, nchar(w[1])))
    paste(w, collapse = " ")
  }, character(1))
  y
}

wrap_label <- function(x, width = 28) {
  x <- pretty_term(x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1), USE.NAMES = FALSE)
}

# One node per term: genes pooled over groups, best padj, and the groups where it is significant.
enrichment_nodes <- function(tab, padj_cutoff = 0.05, top_n = 30) {
  sig <- data.table::as.data.table(tab)[padj <= padj_cutoff]
  if (!nrow(sig)) return(data.table::data.table())
  nodes <- sig[, .(padj = min(padj), groups = paste(sort(unique(group)), collapse = ";"),
    genes = paste(unique(unlist(lapply(genes, split_genes))), collapse = ";")), by = term]
  nodes[, count := lengths(strsplit(genes, ";"))]
  utils::head(nodes[order(padj)], top_n)
}

# Edges: gene-set overlap (Jaccard) >= min_jaccard, keeping only each term's k strongest
# links. Redundant ontologies (GO) otherwise connect almost every pair into a hairball.
build_enrichment_map <- function(tab, padj_cutoff = 0.05, top_n = 30, min_jaccard = 0.25, k = 3) {
  require_pkg("igraph")
  nodes <- enrichment_nodes(tab, padj_cutoff, top_n)
  if (nrow(nodes) < 3) return(NULL)
  gene_sets <- lapply(nodes$genes, split_genes)
  pairs <- utils::combn(nrow(nodes), 2)
  w <- apply(pairs, 2, function(p) jaccard(gene_sets[[p[1]]], gene_sets[[p[2]]]))
  all_edges <- data.table::data.table(i = pairs[1, ], j = pairs[2, ], jaccard = w)[jaccard >= min_jaccard]
  ranked <- data.table::rbindlist(list(all_edges[, .(node = i, i, j, jaccard)], all_edges[, .(node = j, i, j, jaccard)]))
  strongest <- ranked[order(-jaccard), utils::head(.SD, k), by = node]
  kept <- unique(strongest[, .(i, j, jaccard)])
  edges <- data.table::data.table(from = nodes$term[kept$i], to = nodes$term[kept$j], jaccard = round(kept$jaccard, 3))
  g <- igraph::graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  modules <- if (igraph::ecount(g) > 0) {
    igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$jaccard))
  } else {
    stats::setNames(seq_len(nrow(nodes)), nodes$term)
  }
  nodes[, module := as.integer(modules[term])]
  # Name each module after its most significant term; order modules by best padj.
  mod_names <- nodes[order(padj), .(module_label = term[1], best = min(padj), size = .N), by = module]
  mod_names <- mod_names[order(best)]
  mod_names[, module_id := seq_len(.N)]
  nodes <- merge(nodes, mod_names[, .(module, module_id, module_label)], by = "module")
  nodes[, module := module_id][, module_id := NULL]
  list(graph = g, nodes = nodes[order(module, padj)], edges = edges)
}

plot_enrichment_map <- function(em, title = "Pathway enrichment map") {
  require_pkg("ggplot2")
  pc <- publication_colors()
  set.seed(1)
  g <- em$graph
  # Lay out each connected module separately and pack them together, so terms without links
  # do not drift to the far edges of the canvas.
  fr <- function(sub) {
    igraph::layout_with_fr(sub, weights = if (igraph::ecount(sub)) igraph::E(sub)$jaccard * 5 else NULL)
  }
  lay <- igraph::layout_components(g, layout = fr)
  pos <- data.table::data.table(term = igraph::V(g)$name, x = lay[, 1], y = lay[, 2])
  nodes <- merge(em$nodes, pos, by = "term")
  nodes[, neglog := -log10(pmax(padj, 1e-300))]
  edges <- merge(merge(em$edges, pos, by.x = "from", by.y = "term"), pos, by.x = "to", by.y = "term",
    suffixes = c("", "end"))

  # Shaded hull behind each multi-term module, labelled with the module's lead term.
  hulls <- nodes[, {
    if (.N >= 3) {
      h <- grDevices::chull(x, y)
      .(x = x[h], y = y[h])
    } else {
      .(x = x, y = y)
    }
  }, by = module]
  mod_lab <- nodes[, .(x = mean(x), y = max(y) + 0.05 * diff(range(nodes$y)), n = .N,
    label = module_label[1]), by = module][n >= 2]

  mods <- nodes[, .(module_label = module_label[1], n_terms = .N), by = module][order(module)]
  caption <- paste(sprintf("M%d: %s (%d term%s)", mods$module, pretty_term(mods$module_label), mods$n_terms,
    ifelse(mods$n_terms == 1, "", "s")), collapse = "\n")

  label_fun <- if (requireNamespace("ggrepel", quietly = TRUE)) {
    ggrepel::geom_text_repel(data = nodes, ggplot2::aes(x, y, label = wrap_label(term, 26)), size = 5,
      colour = "black", lineheight = 0.85, max.overlaps = 40, seed = 1, min.segment.length = 0.3,
      segment.colour = "grey60")
  } else {
    ggplot2::geom_text(data = nodes, ggplot2::aes(x, y, label = wrap_label(term, 26)), size = 4.6,
      colour = "black", vjust = -1)
  }
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = hulls, ggplot2::aes(x, y, group = module), fill = "grey85", alpha = 0.45,
      colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend,
      linewidth = jaccard), colour = "grey55", alpha = 0.7) +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y, size = count, fill = neglog), shape = 21,
      colour = "black", stroke = 0.3) +
    label_fun +
    ggplot2::geom_label(data = mod_lab, ggplot2::aes(x, y, label = paste0("M", module)), fontface = "bold",
      size = 6.5, colour = "black", fill = "white", label.size = 0.3) +
    ggplot2::scale_fill_gradient(low = pc$light, high = pc$primary, name = "-log10 padj") +
    ggplot2::scale_size_continuous(range = c(3.5, 12), name = "# genes", breaks = count_breaks) +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 1.6), name = "Gene overlap\n(Jaccard)") +
    ggplot2::labs(title = title, caption = caption) +
    theme_network_publication()
  p
}

# `lead_terms` (e.g. one per enrichment-map module) avoids showing near-duplicate pathways.
plot_gene_concept_network <- function(tab, padj_cutoff = 0.05, n_terms = 5, max_genes_per_term = 25,
                                      title = "Gene-pathway network", lead_terms = NULL) {
  require_pkg("igraph")
  require_pkg("ggplot2")
  pc <- publication_colors()
  all_nodes <- enrichment_nodes(tab, padj_cutoff, top_n = Inf)
  chosen <- unique(c(intersect(lead_terms, all_nodes$term), all_nodes$term))
  nodes <- all_nodes[match(utils::head(chosen, n_terms), term)]
  if (nrow(nodes) < 2) return(NULL)
  links <- data.table::rbindlist(lapply(seq_len(nrow(nodes)), function(i) {
    data.table::data.table(term = nodes$term[i], gene = utils::head(split_genes(nodes$genes[i]), max_genes_per_term))
  }))
  degree <- links[, .N, by = gene]
  g <- igraph::graph_from_data_frame(links, directed = FALSE)
  set.seed(1)
  lay <- igraph::layout_with_fr(g)
  pos <- data.table::data.table(name = igraph::V(g)$name, x = lay[, 1], y = lay[, 2])
  pos[, type := data.table::fifelse(name %in% nodes$term, "pathway", "gene")]
  pos <- merge(pos, degree, by.x = "name", by.y = "gene", all.x = TRUE)
  pos <- merge(pos, nodes[, .(name = term, padj, count)], by = "name", all.x = TRUE)
  edges <- merge(merge(links, pos[, .(name, x, y)], by.x = "term", by.y = "name"),
    pos[, .(name, xend = x, yend = y)], by.x = "gene", by.y = "name")
  terms <- pos[type == "pathway"]
  genes <- pos[type == "gene"]
  genes[, hub := N >= 2]

  label_layer <- function(data, size, face) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      ggrepel::geom_text_repel(data = data, ggplot2::aes(x, y, label = label), size = size, fontface = face,
        colour = "black", max.overlaps = 60, seed = 1, lineheight = 0.85, segment.colour = "grey70",
        box.padding = if (face == "bold") 0.9 else 0.3, point.padding = if (face == "bold") 0.6 else 0.1)
    } else {
      ggplot2::geom_text(data = data, ggplot2::aes(x, y, label = label), size = size, fontface = face,
        colour = "black", vjust = -0.9)
    }
  }
  terms[, label := wrap_label(name, 24)]
  genes[, label := name]
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), colour = "grey75",
      linewidth = 0.3) +
    ggplot2::geom_point(data = genes, ggplot2::aes(x, y, size = hub), shape = 21, fill = pc$bar, colour = "white",
      stroke = 0.3) +
    ggplot2::geom_point(data = terms, ggplot2::aes(x, y, fill = -log10(pmax(padj, 1e-300))), size = 14, shape = 22,
      colour = "black", stroke = 0.4) +
    label_layer(genes, 6, "italic") +
    label_layer(terms, 5.6, "bold") +
    ggplot2::scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 5.5), labels = c("1 pathway", "2+ pathways"),
      name = "Gene") +
    ggplot2::scale_fill_gradient(low = pc$light, high = pc$primary, name = "-log10 padj") +
    ggplot2::labs(title = title) +
    theme_network_publication()
  list(plot = p, hub_genes = genes[hub == TRUE][order(-N), .(gene = name, n_pathways = N)], n_nodes = nrow(pos))
}

# Networks for every method in a standardised enrichment table. Returns figures,
# Cytoscape tables and a module summary for the agent.
run_pathway_networks <- function(tab, out_dir, prefix, padj_cutoff = 0.05, top_n = 30, min_jaccard = 0.25) {
  if (!requireNamespace("igraph", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    return(list(warnings = "Install igraph and ggplot2 to draw pathway networks.", figures = list(),
      tables = list(), modules = list()))
  }
  tab <- data.table::as.data.table(tab)
  figures <- list(); tables <- list(); modules <- list(); hubs <- list(); warnings <- character()
  for (m in unique(tab$method)) {
    sub <- tab[method == m]
    em <- tryCatch(build_enrichment_map(sub, padj_cutoff, top_n, min_jaccard), error = function(e) {
      warnings <<- c(warnings, sprintf("Enrichment map for %s failed: %s", m, conditionMessage(e)))
      NULL
    })
    if (!is.null(em)) {
      n_terms <- nrow(em$nodes)
      fig <- save_figure(plot_enrichment_map(em, sprintf("Pathway enrichment map: %s", m)), out_dir,
        paste0(prefix, "_network_", file_tag(m)), width = if (n_terms <= 10) 10 else 12.5,
        height = if (n_terms <= 10) 7.5 else max(9.5, 0.2 * n_terms + 6.5))
      figures[[paste0("map_", m)]] <- fig[["png"]]
      node_path <- file.path(out_dir, paste0(prefix, "_network_", file_tag(m), "_nodes.tsv"))
      edge_path <- file.path(out_dir, paste0(prefix, "_network_", file_tag(m), "_edges.tsv"))
      write_tsv(em$nodes[, .(term, module, module_label, padj, count, groups, genes)], node_path)
      write_tsv(em$edges, edge_path)
      tables[[paste0("network_nodes_", m)]] <- node_path
      tables[[paste0("network_edges_", m)]] <- edge_path
      mods <- em$nodes[, .(lead_term = module_label[1], n_terms = .N, best_padj = signif(min(padj), 3),
        other_terms = paste(utils::head(setdiff(term, module_label[1]), 3), collapse = "; ")), by = module]
      modules[[m]] <- mods[order(module)]
    }
    leads <- if (!is.null(em)) unique(em$nodes[order(module), module_label]) else NULL
    cn <- tryCatch(plot_gene_concept_network(sub, padj_cutoff, title = sprintf("Gene-pathway network: %s", m),
      lead_terms = leads),
      error = function(e) {
        warnings <<- c(warnings, sprintf("Gene-pathway network for %s failed: %s", m, conditionMessage(e)))
        NULL
      })
    if (!is.null(cn)) {
      big <- cn$n_nodes > 40
      fig <- save_figure(cn$plot, out_dir, paste0(prefix, "_cnet_", file_tag(m)), width = if (big) 13 else 11,
        height = if (big) 11 else 9)
      figures[[paste0("cnet_", m)]] <- fig[["png"]]
      if (nrow(cn$hub_genes)) hubs[[m]] <- utils::head(cn$hub_genes, 15)
    }
  }
  list(figures = figures, tables = tables, modules = modules, hub_genes = hubs, warnings = warnings)
}
