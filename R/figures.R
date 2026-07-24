plot_annotation_summary <- function(annotation_table, out_dir, prefix) {
  require_pkg("data.table")
  dt <- data.table::fread(annotation_table)
  anno_col <- if ("annotation" %in% names(dt)) "annotation" else if ("Annotation" %in% names(dt)) "Annotation" else NULL
  if (is.null(anno_col)) stop("Annotation table needs an 'annotation' or 'Annotation' column.", call. = FALSE)
  dt[, Annotation := sub(" \\(.*", "", get(anno_col))]
  tab <- dt[, .N, by = Annotation][order(-N)]
  tab[, Pct := 100 * N / sum(N)]
  write_tsv(tab, file.path(out_dir, paste0(prefix, "_annotation_summary.tsv")))

  plot_annotation_summary_base(tab, out_dir, prefix)
  invisible(tab)
}

plot_enrichment_dotplot <- function(enrichment_table, out_file, padj_cutoff = 0.05, top_show = 10) {
  require_pkg("data.table")
  dt <- if (is.character(enrichment_table)) data.table::fread(enrichment_table) else data.table::as.data.table(enrichment_table)
  if (!nrow(dt)) return(invisible(NULL))

  padj_col <- intersect(c("Adjusted.P.value", "p.adjust", "adj.P.Val"), names(dt))[1]
  term_col <- intersect(c("Term", "Description", "Category"), names(dt))[1]
  group_col <- intersect(c("Cluster", "Group"), names(dt))[1]
  if (is.na(padj_col) || is.na(term_col) || is.na(group_col)) {
    stop("Enrichment table must contain adjusted p-value, term, and cluster/group columns.", call. = FALSE)
  }
  dt <- dt[get(padj_col) <= padj_cutoff]
  if (!nrow(dt)) return(invisible(NULL))

  if ("Genes" %in% names(dt)) {
    dt[, nGene := lengths(strsplit(Genes, ";|,", perl = TRUE))]
  } else if ("Count" %in% names(dt)) {
    dt[, nGene := as.numeric(Count)]
  } else {
    dt[, nGene := 1]
  }
  dt[, negLogPadj := pmin(-log10(get(padj_col)), 50)]
  top <- dt[order(get(group_col), get(padj_col)), head(.SD, top_show), by = group_col]
  dt <- dt[get(term_col) %in% unique(top[[term_col]])]
  dt[, TermPlot := factor(get(term_col), levels = rev(unique(top[[term_col]])))]

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot_enrichment_dotplot_base(dt, out_file, group_col = group_col, term_col = term_col)
    return(invisible(dt))
  }

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = get(group_col), y = TermPlot, size = nGene, color = negLogPadj)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient(low = "#fdae61", high = "#7f0000", name = "-log10 padj") +
    ggplot2::scale_size_continuous(name = "# genes") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Pathway enrichment")
  ggplot2::ggsave(out_file, p, width = 10, height = max(4, length(unique(dt[[term_col]])) * 0.28 + 2),
    limitsize = FALSE, bg = "white")
  invisible(p)
}

plot_annotation_summary_base <- function(tab, out_dir, prefix) {
  feature_levels <- c("Promoter", "5' UTR", "3' UTR", "1st Exon", "Other Exon",
    "1st Intron", "Other Intron", "Downstream", "Downstream (<=300)", "Distal Intergenic",
    "Intron", "Exon")
  feature_colors <- c(
    "Promoter" = "#a6cee3",
    "5' UTR" = "#1f78b4",
    "3' UTR" = "#b2df8a",
    "1st Exon" = "#33a02c",
    "Other Exon" = "#fb9a99",
    "1st Intron" = "#e31a1c",
    "Other Intron" = "#fdbf6f",
    "Downstream" = "#ff7f00",
    "Downstream (<=300)" = "#ff7f00",
    "Distal Intergenic" = "#cab2d6",
    "Intron" = "#e31a1c",
    "Exon" = "#fb9a99"
  )
  extras <- setdiff(tab$Annotation, names(feature_colors))
  if (length(extras)) {
    extra_cols <- grDevices::hcl.colors(length(extras), "Set 2")
    names(extra_cols) <- extras
    feature_colors <- c(feature_colors, extra_cols)
  }
  tab$Annotation <- factor(tab$Annotation, levels = c(intersect(feature_levels, tab$Annotation), extras))
  tab <- tab[order(tab$Annotation), ]
  colors <- unname(feature_colors[as.character(tab$Annotation)])
  bar_pdf <- file.path(out_dir, paste0(prefix, "_annotation_bar.pdf"))
  pie_pdf <- file.path(out_dir, paste0(prefix, "_annotation_pie.pdf"))

  draw_bar <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4.2, 1.3))
    graphics::par(mar = c(4, 2, 1, 1), xpd = FALSE)
    widths <- matrix(tab$Pct, nrow = nrow(tab))
    y_mid <- graphics::barplot(widths, horiz = TRUE, beside = FALSE, col = colors, border = NA,
      xlim = c(0, 100), axes = FALSE, names.arg = "", xlab = "")
    graphics::axis(1, at = seq(0, 100, 25))
    x_left <- c(0, head(cumsum(tab$Pct), -1))
    x_mid <- x_left + tab$Pct / 2
    labels <- ifelse(tab$Pct %% 1 == 0, sprintf("%.0f%%", tab$Pct), sprintf("%.1f%%", tab$Pct))
    inside <- tab$Pct >= 8
    graphics::text(x_mid[inside], y_mid, labels[inside], cex = 0.8, col = "black")
    if (any(!inside)) {
      graphics::text(x_mid[!inside], y_mid + 0.55, labels[!inside], cex = 0.7, col = "black")
    }
    graphics::par(mar = c(1, 0, 1, 1))
    graphics::plot.new()
    graphics::legend("topleft", legend = as.character(tab$Annotation),
      fill = colors, title = "Feature", bty = "n", cex = 0.82)
  }
  grDevices::pdf(bar_pdf, width = 9, height = 3.2, bg = "white"); draw_bar(); grDevices::dev.off()

  draw_pie <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(3.2, 1.2))
    graphics::par(mar = c(2, 2, 2, 1), xpd = FALSE)
    graphics::pie(tab$N, labels = rep("", nrow(tab)), col = colors, border = "black", main = "")
    graphics::par(mar = c(1, 0, 1, 1))
    graphics::plot.new()
    graphics::legend("left", legend = as.character(tab$Annotation),
      fill = colors, title = "Feature", bty = "n", cex = 0.82)
  }
  grDevices::pdf(pie_pdf, width = 7, height = 4.5, bg = "white"); draw_pie(); grDevices::dev.off()
}

plot_enrichment_dotplot_base <- function(dt, out_file, group_col, term_col) {
  pdf_file <- out_file
  groups <- unique(dt[[group_col]])
  terms <- rev(levels(dt$TermPlot))
  term_labels <- vapply(terms, function(x) paste(strwrap(x, width = 34), collapse = "\n"), character(1))
  dt[, x := match(get(group_col), groups)]
  dt[, y := match(as.character(TermPlot), terms)]
  size <- 1.3 + 4 * sqrt(dt$nGene / max(dt$nGene))
  pal <- grDevices::colorRampPalette(c("#fdae61", "#d95f0e", "#7f0000"))(100)
  color_idx <- pmax(1, pmin(100, round(1 + 99 * dt$negLogPadj / max(dt$negLogPadj))))

  draw_dot <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4.5, 1))
    graphics::par(mar = c(7, 18, 4, 1))
    graphics::plot(
      dt$x, dt$y,
      pch = 16,
      cex = size,
      col = pal[color_idx],
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = "",
      xlim = c(0.5, length(groups) + 0.5),
      ylim = c(0.5, length(terms) + 0.5),
      main = "Pathway enrichment"
    )
    graphics::axis(1, at = seq_along(groups), labels = groups, las = 2)
    graphics::axis(2, at = seq_along(terms), labels = term_labels, las = 1, cex.axis = 0.85)
    graphics::grid(nx = length(groups), ny = length(terms), col = "grey90")
    graphics::box()

    graphics::par(mar = c(4, 0, 4, 1))
    graphics::plot.new()
    size_values <- unique(round(pretty(dt$nGene, n = 3)))
    size_values <- size_values[size_values > 0]
    size_values <- size_values[seq_len(min(3, length(size_values)))]
    graphics::legend("topleft",
      legend = as.character(size_values),
      title = "# genes",
      pch = 16,
      pt.cex = 1.3 + 4 * sqrt(size_values / max(dt$nGene)),
      col = "black",
      bty = "n",
      cex = 0.85
    )
    graphics::legend("left",
      legend = c("2", "4", "6"),
      title = "-log10 padj",
      pch = 16,
      pt.cex = 1.3,
      col = pal[round(seq(15, 90, length.out = 3))],
      bty = "n",
      cex = 0.85
    )
  }
  grDevices::pdf(pdf_file, width = 10, height = max(4.5, length(terms) * 0.35 + 2), bg = "white"); draw_dot(); grDevices::dev.off()
}
