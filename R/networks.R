# Pathway networks built from standardised enrichment tables
# (columns: method, group, term, padj, count, genes).
#
# Enrichment map: significant terms are nodes, edges join terms whose gene sets
#   overlap (Jaccard), Louvain communities group redundant terms into modules.
# Gene-pathway network: the top terms linked to their genes; genes shared by
#   several terms are the hubs.
# Node/edge tables are written for Cytoscape alongside the figures (x, y = figure layout, inches).

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
    set.seed(1)   # Louvain breaks ties at random; fix it so modules are reproducible
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

# ---- Network drawing (EnrichmentMap / AutoAnnotate look, after SimpNetworkViz) -------------
# Figures are laid out in physical units: 1 data unit = 1 inch on the saved figure, so text
# sizes (pt) and the room reserved for labels agree. Each module is laid out on its own and the
# modules are packed in rows, every one in a tinted outline with its name above it.

net_style <- list(node_pt = 9.5, gene_pt = 9, title_pt = 13.5, legend_pt = 10.5, hull_pad = 0.12,
  block_gap = 0.45, label_gap = 0.03, wrap = 18, max_lines = 4)

net_palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#EE4C97", "#0072B5",
  "#BC3C29", "#925E9F", "#FDAF91", "#AD002A")
net_neutral <- "#BDBDBD"

net_colour <- function(i) ifelse(i <= 0, net_neutral, net_palette[(i - 1) %% length(net_palette) + 1])

lighten <- function(col, amount) {
  m <- grDevices::col2rgb(col) / 255
  m <- m + (1 - m) * amount
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

darken <- function(col, amount) {
  m <- grDevices::col2rgb(col) / 255 * (1 - amount)
  grDevices::rgb(m[1, ], m[2, ], m[3, ])
}

text_w <- function(n_chars, pt) n_chars * 0.58 * pt / 72
text_h <- function(n_lines, pt) n_lines * 1.12 * pt / 72

# "lymphocyte differentiation" -> "Lymphocyte differentiation"; leaves "mRNA ..." alone.
cap_first <- function(x) {
  first <- sub("[ ].*", "", x)
  low <- grepl("^[a-z][a-z-]*$", first)
  x[low] <- paste0(toupper(substr(x[low], 1, 1)), substr(x[low], 2, nchar(x[low])))
  x
}

# Wrapped node label with at most `max_lines` lines.
net_label <- function(x, width = net_style$wrap, max_lines = net_style$max_lines) {
  x <- cap_first(pretty_term(x))
  vapply(x, function(s) {
    l <- strwrap(s, width = width)
    if (length(l) > max_lines) {
      l <- l[seq_len(max_lines)]
      last <- strsplit(l[max_lines], " ", fixed = TRUE)[[1]]
      while (length(last) > 1 && nchar(paste(last, collapse = " ")) > width - 3) last <- utils::head(last, -1)
      l[max_lines] <- paste0(sub("[ ,;:-]*$", "", paste(last, collapse = " ")), "...")
    }
    paste(l, collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

label_dims <- function(label) {
  lines <- strsplit(label, "\n", fixed = TRUE)
  list(lines = lengths(lines), chars = vapply(lines, function(l) max(nchar(l), 0L), integer(1)))
}

# Node boxes: circle of radius r with its label underneath. w = half width; up/down = extent
# above/below the node centre.
node_boxes <- function(r, label, pt) {
  d <- label_dims(label)
  has <- nzchar(label)
  list(w = pmax(r, ifelse(has, text_w(d$chars, pt) / 2, 0)) + 0.02, up = r,
    down = r + ifelse(has, net_style$label_gap + text_h(d$lines, pt), 0))
}

# Push overlapping boxes apart along the axis that needs the smaller move.
remove_box_overlaps <- function(P, b, gap = 0.08, iters = 400) {
  n <- nrow(P)
  if (n < 2) return(P)
  hh <- (b$up + b$down) / 2
  off <- (b$up - b$down) / 2
  for (it in seq_len(iters)) {
    cx <- P[, 1]
    cy <- P[, 2] + off
    dx <- outer(cx, cx, "-")
    dy <- outer(cy, cy, "-")
    ox <- outer(b$w, b$w, "+") + gap - abs(dx)
    oy <- outer(hh, hh, "+") + gap - abs(dy)
    over <- ox > 0 & oy > 0
    over[lower.tri(over, diag = TRUE)] <- FALSE
    if (!any(over)) break
    idx <- which(over, arr.ind = TRUE)
    for (k in seq_len(nrow(idx))) {
      i <- idx[k, 1]
      j <- idx[k, 2]
      if (ox[i, j] <= oy[i, j]) {
        s <- if (dx[i, j] == 0) 1 else sign(dx[i, j])
        P[i, 1] <- P[i, 1] + s * (ox[i, j] / 2 + 1e-3)
        P[j, 1] <- P[j, 1] - s * (ox[i, j] / 2 + 1e-3)
      } else {
        s <- if (dy[i, j] == 0) 1 else sign(dy[i, j])
        P[i, 2] <- P[i, 2] + s * (oy[i, j] / 2 + 1e-3)
        P[j, 2] <- P[j, 2] - s * (oy[i, j] / 2 + 1e-3)
      }
    }
  }
  P
}

# Lay out `nodes` (a group of vertices of g) around the origin without overlapping boxes.
local_layout <- function(g, nodes, b, seed = 1) {
  n <- length(nodes)
  if (n == 1) return(matrix(0, 1, 2))
  sub <- igraph::induced_subgraph(g, nodes)
  if (igraph::ecount(sub) == 0) {
    cols <- ceiling(sqrt(n * 1.6))
    sx <- 2 * max(b$w) + 0.1
    sy <- max(b$up + b$down) + 0.12
    i <- seq_len(n) - 1
    P <- cbind((i %% cols) * sx, -(i %/% cols) * sy)
    return(sweep(P, 2, colMeans(P)))
  }
  set.seed(seed)
  lay <- if (n == 2) {
    matrix(c(-1, 1, 0, 0), 2)
  } else if (igraph::is_connected(sub)) {
    igraph::layout_with_kk(sub)
  } else {
    igraph::layout_components(sub, layout = igraph::layout_with_kk)
  }
  rownames(lay) <- igraph::V(sub)$name
  P <- lay[nodes, , drop = FALSE]
  P <- sweep(P, 2, colMeans(P))
  if (n >= 3) {
    ev <- eigen(stats::cov(P), symmetric = TRUE)$vectors   # long axis horizontal
    P <- P %*% ev
  }
  d <- as.matrix(stats::dist(P))
  diag(d) <- Inf
  nn <- stats::median(apply(d, 1, min))
  spacing <- stats::median(2 * b$w) + 0.12
  if (nn > 1e-9) P <- P * spacing / nn
  P <- remove_box_overlaps(P, b)
  ext <- block_extent(P, b, 0)
  sweep(P, 2, c((ext[1] + ext[2]) / 2, (ext[3] + ext[4]) / 2))
}

block_extent <- function(P, b, pad) {
  c(min(P[, 1] - b$w) - pad, max(P[, 1] + b$w) + pad, min(P[, 2] - b$down) - pad, max(P[, 2] + b$up) + pad)
}

# Pack blocks (each list(P, w, h, top), P centred on its body) in rows of roughly `aspect`.
pack_blocks <- function(blocks, gap = net_style$block_gap, aspect = 1.4) {
  area <- sum(vapply(blocks, function(k) (k$w + gap) * (k$h + gap), numeric(1)))
  row_w <- max(max(vapply(blocks, `[[`, numeric(1), "w")), sqrt(area * aspect))
  rows <- list(); cur <- list(); cur_w <- 0
  for (k in blocks) {
    if (length(cur) && cur_w + k$w > row_w) {
      rows[[length(rows) + 1]] <- cur
      cur <- list(); cur_w <- 0
    }
    cur[[length(cur) + 1]] <- k
    cur_w <- cur_w + k$w + gap
  }
  if (length(cur)) rows[[length(rows) + 1]] <- cur
  out <- list(); y_top <- 0
  for (row in rows) {
    rw <- sum(vapply(row, `[[`, numeric(1), "w")) + gap * (length(row) - 1)
    x <- -rw / 2
    for (k in row) {
      body_h <- k$h - k$top
      k$P <- sweep(k$P, 2, c(x + k$w / 2, y_top - k$top - body_h / 2), "+")
      out[[length(out) + 1]] <- k
      x <- x + k$w + gap
    }
    y_top <- y_top - max(vapply(row, `[[`, numeric(1), "h")) - gap
  }
  out
}

circle_pts <- function(x, y, r, n = 40) {
  a <- seq(0, 2 * pi, length.out = n + 1)[-1]
  data.frame(x = x + r * cos(a), y = y + r * sin(a))
}

circles_df <- function(x, y, r, id) {
  do.call(rbind, lapply(seq_along(x), function(i) cbind(circle_pts(x[i], y[i], r[i]), id = id[i])))
}

# Rounded convex outline around the nodes and their labels.
hull_df <- function(P, r, b, pad = net_style$hull_pad) {
  pts <- do.call(rbind, lapply(seq_len(nrow(P)), function(i) {
    corners <- cbind(c(P[i, 1] - b$w[i], P[i, 1] + b$w[i]), P[i, 2] - b$down[i])
    rbind(circle_pts(P[i, 1], P[i, 2], r[i] + pad, 32),
      do.call(rbind, lapply(1:2, function(k) circle_pts(corners[k, 1] + c(0.05, -0.05)[k], corners[k, 2] + 0.05,
        pad + 0.05, 12))))
  }))
  pts[grDevices::chull(pts$x, pts$y), ]
}

radius_scale <- function(v, rmin, rmax, lim = range(v)) {
  s <- sqrt(pmax(v, 1)); sl <- sqrt(pmax(lim, 1))
  if (diff(sl) <= 0) return(rep((rmin + rmax) / 2, length(v)))
  rmin + (rmax - rmin) * (s - sl[1]) / diff(sl)
}

pt_size <- function(pt) pt / ggplot2::.pt

# Node labels under the nodes on a translucent patch of the background, so edges passing
# behind them do not cut through the text.
net_text <- function(d, y, pt, bg, face = "plain") {
  d <- data.frame(x = d$x, y = y, label = d$label)
  ggplot2::geom_label(data = d, ggplot2::aes(x, y, label = label), vjust = 1, size = pt_size(pt),
    fontface = face, lineheight = 0.9, colour = "#1A1A1A", fill = grDevices::adjustcolor(bg, alpha.f = 0.8),
    label.size = 0, label.padding = ggplot2::unit(0.6, "pt"), label.r = ggplot2::unit(0, "pt"))
}

# Shared canvas: fixed aspect, 1 unit = 1 inch, figure size attached as attr "fig_size".
net_canvas <- function(p, xr, yr, title) {
  margin <- 0.25
  title_h <- if (nzchar(title)) 0.45 else 0
  p <- p + ggplot2::coord_fixed(xlim = xr, ylim = yr, expand = FALSE, clip = "off") +
    ggplot2::scale_linewidth_identity() + ggplot2::labs(title = title) +
    ggplot2::theme_void(base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 15, colour = "black", hjust = 0.5,
      margin = ggplot2::margin(0, 0, 8, 0)), plot.margin = ggplot2::margin(18, 18, 18, 18),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA))
  attr(p, "fig_size") <- c(diff(xr) + 2 * margin, diff(yr) + 2 * margin + title_h)
  p
}

plot_enrichment_map <- function(em, title = "Pathway enrichment map") {
  require_pkg("ggplot2")
  st <- net_style
  g <- em$graph
  nodes <- data.table::copy(em$nodes)
  nodes[, r := radius_scale(count, 0.07, 0.22)]
  nodes[, label := net_label(term)]
  sizes <- nodes[, .N, by = module]
  nodes[, clustered := module %in% sizes[N >= 2, module]]
  # Modules keep their numbers; colours follow module order so the gene network can match.
  nodes[, colour := data.table::fifelse(clustered, net_colour(module), net_neutral)]

  sig_order <- function(d) d[order(padj, -count, term)]
  groups <- c(lapply(sort(unique(nodes[clustered == TRUE, module])), function(m) sig_order(nodes[module == m])),
    if (any(!nodes$clustered)) list(sig_order(nodes[clustered == FALSE])) else NULL)
  blocks <- lapply(groups, function(d) {
    b <- node_boxes(d$r, d$label, st$node_pt)
    P <- local_layout(g, d$term, b)
    hulled <- d$clustered[1]
    ext <- block_extent(P, b, if (hulled) st$hull_pad else 0.02)
    P <- sweep(P, 2, c((ext[1] + ext[2]) / 2, (ext[3] + ext[4]) / 2))
    w <- ext[2] - ext[1]; h <- ext[4] - ext[3]; top <- 0; ttl <- ""
    if (hulled) {
      wrap <- as.integer(min(40, max(18, w / text_w(1, st$title_pt))))
      ttl <- paste(strwrap(sprintf("M%d: %s", d$module[1], cap_first(pretty_term(d$module_label[1]))), wrap),
        collapse = "\n")
      dd <- label_dims(ttl)
      top <- text_h(dd$lines, st$title_pt) + 0.1
      w <- max(w, text_w(dd$chars, st$title_pt) + 0.1)
    }
    list(d = d, P = P, b = b, w = w, h = h + top, top = top, title = ttl, hulled = hulled)
  })
  blocks <- pack_blocks(blocks)

  pos <- data.table::rbindlist(lapply(blocks, function(k) {
    data.table::data.table(term = k$d$term, x = k$P[, 1], y = k$P[, 2], down = k$b$down)
  }))
  nodes <- merge(nodes, pos, by = "term")
  hulls <- data.table::rbindlist(lapply(blocks[vapply(blocks, `[[`, logical(1), "hulled")], function(k) {
    h <- hull_df(k$P, k$d$r, k$b)
    data.table::data.table(module = k$d$module[1], x = h$x, y = h$y, colour = k$d$colour[1])
  }))
  titles <- data.table::rbindlist(lapply(blocks[vapply(blocks, `[[`, logical(1), "hulled")], function(k) {
    m <- k$d$module[1]
    data.table::data.table(x = mean(range(hulls[module == m, x])), y = max(hulls[module == m, y]) + 0.06,
      label = k$title, colour = darken(k$d$colour[1], 0.45))
  }))

  edges <- data.table::copy(em$edges)
  if (nrow(edges)) {
    edges <- merge(edges, nodes[, .(from = term, x, y, m1 = module)], by = "from")
    edges <- merge(edges, nodes[, .(to = term, xend = x, yend = y, m2 = module)], by = "to")
    jr <- range(edges$jaccard)
    edges[, lw := if (diff(jr) > 0) 0.3 + 1.5 * (jaccard - jr[1]) / diff(jr) else 0.9]
    # Links between modules stay visible but faint, so the modules read as separate themes.
    edges[, within := m1 == m2]
    edges[within == FALSE, lw := 0.35]
    edges[, alpha := data.table::fifelse(within, 0.6, 0.18)]
  }
  circles <- circles_df(nodes$x, nodes$y, nodes$r, nodes$term)
  circles$fill <- nodes$colour[match(circles$id, nodes$term)]

  # Legend along the bottom: gene-set size circles and edge widths.
  all_y <- c(nodes$y - nodes$down, hulls$y)
  xr <- range(c(nodes$x - 0.7, nodes$x + 0.7, hulls$x))
  ly <- min(all_y) - 0.45
  lim <- range(nodes$count)
  brk <- unique(round(c(lim[1], stats::median(nodes$count), lim[2])))
  lr <- radius_scale(brk, 0.07, 0.22, lim)
  lx <- xr[1] + 0.25 + cumsum(c(0, utils::head(lr, -1) + lr[-1] + 0.3))
  leg_c <- circles_df(lx, rep(ly - 0.55, length(brk)), lr, paste0("leg", seq_along(brk)))
  leg_t <- data.table::data.table(x = c(xr[1], lx, xr[1] + 3.0),
    y = c(ly, rep(ly - 0.55 - max(lr) - 0.12, length(brk)), ly),
    label = c("Gene-set size", brk, "Edge width = shared genes (Jaccard)"),
    vjust = c(1, 1, rep(1, length(brk) - 1), 1), hjust = c(0, rep(0.5, length(brk)), 0))
  leg_e <- if (nrow(edges)) {
    data.table::data.table(x = xr[1] + 3.0 + c(0, 1.0, 2.0), xend = xr[1] + 3.0 + c(0.7, 1.7, 2.7),
      y = ly - 0.55, lw = c(0.3, 1.05, 1.8),
      lab = sprintf("%.2f", c(jr[1], mean(jr), jr[2])))
  } else NULL
  yr <- c(ly - 0.55 - max(lr) - 0.45, max(c(titles$y + 0.5, nodes$y + nodes$r)) + 0.05)
  xr <- range(c(xr, if (!is.null(leg_e)) max(leg_e$xend) + 0.2))

  p <- ggplot2::ggplot()
  if (nrow(hulls)) {
    p <- p + ggplot2::geom_polygon(data = hulls, ggplot2::aes(x, y, group = module),
      fill = lighten(hulls$colour, 0.86), colour = lighten(hulls$colour, 0.25), linewidth = 0.6)
  }
  if (nrow(edges)) {
    p <- p + ggplot2::geom_segment(data = edges, ggplot2::aes(x, y, xend = xend, yend = yend, linewidth = lw),
      colour = "#8F8F8F", alpha = edges$alpha, lineend = "round")
  }
  p <- p + ggplot2::geom_polygon(data = circles, ggplot2::aes(x, y, group = id), fill = circles$fill,
    colour = darken(circles$fill, 0.3), linewidth = 0.3) +
    net_text(nodes, nodes$y - nodes$r - st$label_gap, st$node_pt,
      data.table::fifelse(nodes$clustered, lighten(nodes$colour, 0.86), "#FFFFFF"))
  if (nrow(titles)) {
    p <- p + ggplot2::geom_text(data = titles, ggplot2::aes(x, y, label = label), colour = titles$colour,
      vjust = 0, size = pt_size(st$title_pt), lineheight = 0.9)
  }
  p <- p + ggplot2::geom_polygon(data = leg_c, ggplot2::aes(x, y, group = id), fill = "#E6E6E6", colour = "#6E6E6E",
    linewidth = 0.3) +
    ggplot2::geom_text(data = leg_t, ggplot2::aes(x, y, label = label), hjust = leg_t$hjust, vjust = 1,
      size = pt_size(st$legend_pt))
  if (!is.null(leg_e)) {
    p <- p + ggplot2::geom_segment(data = leg_e, ggplot2::aes(x, y, xend = xend, yend = y, linewidth = lw),
      colour = "#8F8F8F") +
      ggplot2::geom_text(data = leg_e, ggplot2::aes((x + xend) / 2, y - 0.12, label = lab), vjust = 1,
        size = pt_size(st$legend_pt - 1))
  }
  p <- net_canvas(p, xr, yr, title)
  attr(p, "positions") <- nodes[, .(term, x = round(x, 3), y = round(y, 3), colour)]
  p
}

# `lead_terms` (e.g. one per enrichment-map module) avoids showing near-duplicate pathways;
# `term_colours` (term -> colour) matches the pathway hubs to the enrichment-map modules.
plot_gene_concept_network <- function(tab, padj_cutoff = 0.05, n_terms = 5, max_genes_per_term = 25,
                                      title = "Gene-pathway network", lead_terms = NULL, term_colours = NULL) {
  require_pkg("igraph")
  require_pkg("ggplot2")
  st <- net_style
  all_nodes <- enrichment_nodes(tab, padj_cutoff, top_n = Inf)
  chosen <- unique(c(intersect(lead_terms, all_nodes$term), all_nodes$term))
  terms <- all_nodes[match(utils::head(chosen, n_terms), term)]
  if (nrow(terms) < 2) return(NULL)
  links <- data.table::rbindlist(lapply(seq_len(nrow(terms)), function(i) {
    data.table::data.table(term = terms$term[i], gene = utils::head(split_genes(terms$genes[i]), max_genes_per_term))
  }))
  # Gene names can collide with term names in principle; keep the two vertex sets apart.
  links[, gid := paste0("gene:", gene)]
  degree <- links[, .(N = .N, first = term[1]), by = .(gene, gid)]
  tcol <- stats::setNames(net_colour(seq_len(nrow(terms))), terms$term)
  if (!is.null(term_colours)) {
    hit <- intersect(names(term_colours), terms$term)
    tcol[hit] <- term_colours[hit]
  }
  g <- igraph::graph_from_data_frame(links[, .(term, gid)], directed = FALSE)

  vt <- data.table::data.table(id = terms$term, label = net_label(terms$term),
    r = radius_scale(terms$count, 0.12, 0.26), pt = st$node_pt, type = "pathway", fill = tcol[terms$term],
    face = "plain")
  vg <- degree[, .(id = gid, label = gene, r = 0.05 + 0.012 * pmin(N, 4), pt = st$gene_pt, type = "gene",
    fill = data.table::fifelse(N >= 2, "#525252", lighten(tcol[first], 0.35)), face = "italic")]
  v <- rbind(vt, vg)
  b <- node_boxes(v$r, v$label, v$pt)
  P <- local_layout(g, v$id, b)
  v[, `:=`(x = P[, 1], y = P[, 2], down = b$down)]

  edges <- merge(links[, .(term, gid)], v[, .(term = id, x, y)], by = "term")
  edges <- merge(edges, v[, .(gid = id, xend = x, yend = y)], by = "gid")
  edges[, colour := tcol[term]]
  circles <- circles_df(v$x, v$y, v$r, v$id)
  circles$fill <- v$fill[match(circles$id, v$id)]

  xr <- range(c(v$x - b$w, v$x + b$w)) + c(-0.1, 0.1)
  ly <- min(v$y - v$down) - 0.45
  leg <- data.table::data.table(r = c(0.19, 0.062, 0.074), fill = c("#9E9E9E", lighten("#9E9E9E", 0.35), "#525252"),
    label = c("Pathway (size = # genes)", "Gene in 1 pathway (pathway tint)", "Gene in 2+ pathways"))
  item_w <- 2 * leg$r + 0.08 + text_w(nchar(leg$label), st$legend_pt) + 0.35
  leg[, x := xr[1] + r + cumsum(c(0, utils::head(item_w, -1)))]
  leg_c <- circles_df(leg$x, rep(ly - 0.2, 3), leg$r, leg$label)
  leg_c$fill <- leg$fill[match(leg_c$id, leg$label)]
  yr <- c(ly - 0.55, max(v$y + v$r) + 0.1)
  xr <- range(c(xr, xr[1] + sum(item_w)))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = edges, ggplot2::aes(x, y, xend = xend, yend = yend), colour = edges$colour,
      alpha = 0.4, linewidth = 0.35) +
    ggplot2::geom_polygon(data = circles, ggplot2::aes(x, y, group = id), fill = circles$fill,
      colour = darken(circles$fill, 0.3), linewidth = 0.3) +
    net_text(v, v$y - v$r - st$label_gap, v$pt, "#FFFFFF", v$face) +
    ggplot2::geom_polygon(data = leg_c, ggplot2::aes(x, y, group = id), fill = leg_c$fill,
      colour = darken(leg_c$fill, 0.3), linewidth = 0.3) +
    ggplot2::geom_text(data = leg, ggplot2::aes(x + r + 0.08, ly - 0.2, label = label), hjust = 0,
      size = pt_size(st$legend_pt))
  p <- net_canvas(p, xr, yr, title)
  genes <- v[type == "gene"]
  genes <- merge(genes, degree[, .(id = gid, N)], by = "id")
  list(plot = p, hub_genes = genes[N >= 2][order(-N, label), .(gene = label, n_pathways = N)], n_nodes = nrow(v))
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
    term_colours <- NULL
    if (!is.null(em)) {
      map <- plot_enrichment_map(em, sprintf("Pathway enrichment map: %s", m))
      size <- attr(map, "fig_size")
      fig <- save_figure(map, out_dir, paste0(prefix, "_network_", file_tag(m)), width = size[1], height = size[2])
      figures[[paste0("map_", m)]] <- fig[["png"]]
      xy <- attr(map, "positions")
      em$nodes <- merge(em$nodes, xy[, .(term, x, y)], by = "term", sort = FALSE)[order(module, padj)]
      term_colours <- xy[colour != net_neutral, stats::setNames(colour, term)]
      node_path <- file.path(out_dir, paste0(prefix, "_network_", file_tag(m), "_nodes.tsv"))
      edge_path <- file.path(out_dir, paste0(prefix, "_network_", file_tag(m), "_edges.tsv"))
      write_tsv(em$nodes[, .(term, module, module_label, padj, count, groups, genes, x, y)], node_path)
      write_tsv(em$edges, edge_path)
      tables[[paste0("network_nodes_", m)]] <- node_path
      tables[[paste0("network_edges_", m)]] <- edge_path
      mods <- em$nodes[, .(lead_term = module_label[1], n_terms = .N, best_padj = signif(min(padj), 3),
        other_terms = paste(utils::head(setdiff(term, module_label[1]), 3), collapse = "; ")), by = module]
      modules[[m]] <- mods[order(module)]
    }
    leads <- if (!is.null(em)) unique(em$nodes[order(module), module_label]) else NULL
    cn <- tryCatch(plot_gene_concept_network(sub, padj_cutoff, title = sprintf("Gene-pathway network: %s", m),
      lead_terms = leads, term_colours = term_colours),
      error = function(e) {
        warnings <<- c(warnings, sprintf("Gene-pathway network for %s failed: %s", m, conditionMessage(e)))
        NULL
      })
    if (!is.null(cn)) {
      size <- attr(cn$plot, "fig_size")
      fig <- save_figure(cn$plot, out_dir, paste0(prefix, "_cnet_", file_tag(m)), width = size[1], height = size[2])
      figures[[paste0("cnet_", m)]] <- fig[["png"]]
      if (nrow(cn$hub_genes)) hubs[[m]] <- utils::head(cn$hub_genes, 15)
    }
  }
  list(figures = figures, tables = tables, modules = modules, hub_genes = hubs, warnings = warnings)
}
