# Nature-style figures: Helvetica/Arial sans serif, black text, thin L-shaped axes with outward
# ticks, no gridlines or chart junk, and vector PDFs with editable text. Colours are dark (Penn)
# red with Penn Blue as the secondary and the Nature Publishing Group (NPG) palette for
# categories. Every figure in the package takes its colours, fonts and PDF device from here.
publication_colors <- function() {
  list(
    primary = "#990000",   # dark red (Penn Red)
    accent = "#011F5B",    # Penn Blue, for reference lines and contrast
    mid = "#C24A4A",
    light = "#F2CFCF",
    grey = "#7F7F7F",
    bar = "#404040",       # dark grey for bar plots and histograms
    npg = c("#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F", "#8491B4", "#91D1C2",
      "#DC0000", "#7E6148", "#B09C85")
  )
}

# Whole-number legend breaks for gene counts (no "17.5 genes").
count_breaks <- function(limits) {
  b <- unique(round(pretty(limits, n = 4)))
  b[b >= max(1, floor(limits[1])) & b <= ceiling(limits[2])]
}

# Default ggplot2 base font size (pt) for package figures.
publication_base_size <- 18

# Light-to-dark-red ramp for significance (darker = more significant).
publication_gradient <- function(n = 100) {
  pc <- publication_colors()
  grDevices::colorRampPalette(c(pc$light, pc$mid, pc$primary))(n)
}

# PDF device used for every figure: Helvetica (Arial-metric) text that stays editable in
# Illustrator, and real circles instead of Dingbats glyphs, as Nature's artwork guide asks.
publication_pdf <- function(filename, width = 7, height = 4.5, bg = "white", ...) {
  grDevices::pdf(file = filename, width = width, height = height, family = "Helvetica", useDingbats = FALSE,
    bg = bg, ...)
}

# ggsave() to a Nature-style PDF.
ggsave_publication <- function(file, plot, width, height, ...) {
  ggplot2::ggsave(file, plot, width = width, height = height, device = publication_pdf, limitsize = FALSE, ...)
}

# Base-graphics defaults applied to every figure written by save_figure().
publication_par <- function() {
  graphics::par(family = "sans", las = 1, tcl = -0.3, mgp = c(3.4, 0.8, 0), cex.axis = 1.5,
    cex.lab = 1.8, cex.main = 1.7, font.main = 2, bty = "l", lwd = 0.8, col.main = "black",
    col.lab = "black", col.axis = "black", col = "black", fg = "black")
}

theme_publication <- function(base_size = publication_base_size) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size),
      axis.title = ggplot2::element_text(colour = "black", size = base_size + 2),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
      legend.text = ggplot2::element_text(colour = "black", size = base_size - 1),
      legend.title = ggplot2::element_text(colour = "black", size = base_size),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(0.5, "cm"),
      plot.title = ggplot2::element_text(face = "bold", colour = "black", size = base_size + 2),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "black", size = base_size),
      panel.grid = ggplot2::element_blank()
    )
}

# Network figures: Nature-style text on a blank canvas.
theme_network_publication <- function(base_size = publication_base_size) {
  ggplot2::theme_void(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      plot.title = ggplot2::element_text(face = "bold", colour = "black", size = base_size + 3),
      plot.caption = ggplot2::element_text(colour = "black", hjust = 0, size = base_size, lineheight = 1.15),
      plot.caption.position = "plot",
      legend.text = ggplot2::element_text(colour = "black", size = base_size - 1),
      legend.title = ggplot2::element_text(colour = "black", size = base_size),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

# enrichplot dot plots colour by p.adjust (smaller = more significant = darker).
scale_color_padj_publication <- function() {
  pc <- publication_colors()
  ggplot2::scale_color_gradient(low = pc$primary, high = pc$light, name = "p.adjust")
}

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

plot_enrichment_dotplot <- function(enrichment_table, out_file, padj_cutoff = 0.05, top_show = 10,
                                    title = "Pathway enrichment") {
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
  top <- dt[order(get(group_col), get(padj_col)), utils::head(.SD, top_show), by = group_col]
  dt <- dt[get(term_col) %in% unique(top[[term_col]])]
  # Readable, wrapped labels: MSigDB names such as GSE7218_UNSTIM_VS_... would otherwise push the
  # plot panel off the page.
  dt[, TermPlot := factor(wrap_label(get(term_col), 42), levels = rev(unique(wrap_label(top[[term_col]], 42))))]

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    plot_enrichment_dotplot_base(dt, out_file, group_col = group_col, term_col = term_col)
    return(invisible(dt))
  }

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = get(group_col), y = TermPlot, size = nGene, color = negLogPadj)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient(low = publication_colors()$light, high = publication_colors()$primary,
      name = "-log10 padj") +
    ggplot2::scale_size_continuous(name = "# genes", breaks = count_breaks, range = c(3, 10)) +
    theme_publication() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.3)) +
    ggplot2::labs(x = NULL, y = NULL, title = title)
  height <- max(5.5, sum(lengths(strsplit(levels(dt$TermPlot), "\n"))) * 0.32 + 3)
  if (grepl("\\.pdf$", out_file)) {
    ggsave_publication(out_file, p, width = 13, height = height)
  } else {
    ggplot2::ggsave(out_file, p, width = 13, height = height, limitsize = FALSE, bg = "white")
  }
  invisible(p)
}

plot_annotation_summary_base <- function(tab, out_dir, prefix) {
  feature_levels <- c("Promoter", "5' UTR", "3' UTR", "1st Exon", "Other Exon",
    "1st Intron", "Other Intron", "Downstream", "Downstream (<=300)", "Distal Intergenic",
    "Intron", "Exon")
  npg <- publication_colors()$npg
  feature_colors <- c(
    "Promoter" = publication_colors()$primary,
    "5' UTR" = npg[4],
    "3' UTR" = npg[3],
    "1st Exon" = npg[7],
    "Other Exon" = npg[6],
    "1st Intron" = npg[2],
    "Other Intron" = npg[5],
    "Downstream" = npg[9],
    "Downstream (<=300)" = npg[9],
    "Distal Intergenic" = npg[10],
    "Intron" = npg[2],
    "Exon" = npg[6]
  )
  extras <- setdiff(tab$Annotation, names(feature_colors))
  if (length(extras)) {
    extra_cols <- rep(c(publication_colors()$accent, npg[c(1, 8)], publication_colors()$grey),
      length.out = length(extras))
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
    publication_par()
    graphics::par(mar = c(5, 2, 1, 1), xpd = FALSE)
    widths <- matrix(tab$Pct, nrow = nrow(tab))
    y_mid <- graphics::barplot(widths, horiz = TRUE, beside = FALSE, col = colors, border = NA,
      xlim = c(0, 100), axes = FALSE, names.arg = "", xlab = "")
    graphics::axis(1, at = seq(0, 100, 25))
    x_left <- c(0, utils::head(cumsum(tab$Pct), -1))
    x_mid <- x_left + tab$Pct / 2
    labels <- ifelse(tab$Pct %% 1 == 0, sprintf("%.0f%%", tab$Pct), sprintf("%.1f%%", tab$Pct))
    inside <- tab$Pct >= 8
    graphics::text(x_mid[inside], y_mid, labels[inside], cex = 1.4, col = "black")
    if (any(!inside)) {
      graphics::text(x_mid[!inside], y_mid + 0.55, labels[!inside], cex = 1.2, col = "black")
    }
    graphics::par(mar = c(1, 0, 1, 1))
    graphics::plot.new()
    graphics::legend("topleft", legend = as.character(tab$Annotation),
      fill = colors, title = "Feature", bty = "n", cex = 1.4)
  }
  publication_pdf(bar_pdf, width = 11, height = 4.2); draw_bar(); grDevices::dev.off()

  draw_pie <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(3.2, 1.2))
    publication_par()
    graphics::par(mar = c(2, 2, 2, 1), xpd = FALSE)
    graphics::pie(tab$N, labels = rep("", nrow(tab)), col = colors, border = "black", main = "")
    graphics::par(mar = c(1, 0, 1, 1))
    graphics::plot.new()
    graphics::legend("left", legend = as.character(tab$Annotation),
      fill = colors, title = "Feature", bty = "n", cex = 1.4)
  }
  publication_pdf(pie_pdf, width = 9, height = 5.5); draw_pie(); grDevices::dev.off()
}

plot_enrichment_dotplot_base <- function(dt, out_file, group_col, term_col) {
  pdf_file <- out_file
  groups <- unique(dt[[group_col]])
  terms <- rev(levels(dt$TermPlot))
  term_labels <- vapply(terms, function(x) paste(strwrap(x, width = 34), collapse = "\n"), character(1))
  dt[, x := match(get(group_col), groups)]
  dt[, y := match(as.character(TermPlot), terms)]
  size <- 1.3 + 4 * sqrt(dt$nGene / max(dt$nGene))
  pal <- publication_gradient(100)
  color_idx <- pmax(1, pmin(100, round(1 + 99 * dt$negLogPadj / max(dt$negLogPadj))))

  draw_dot <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    publication_par()
    graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(4.5, 1))
    graphics::par(mar = c(8, 26, 4, 1))
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
    graphics::axis(2, at = seq_along(terms), labels = term_labels, las = 1, cex.axis = 1.3)
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
      cex = 1.4
    )
    graphics::legend("left",
      legend = c("2", "4", "6"),
      title = "-log10 padj",
      pch = 16,
      pt.cex = 1.3,
      col = pal[round(seq(15, 90, length.out = 3))],
      bty = "n",
      cex = 1.4
    )
  }
  publication_pdf(pdf_file, width = 13, height = max(6, length(terms) * 0.5 + 3)); draw_dot(); grDevices::dev.off()
}
