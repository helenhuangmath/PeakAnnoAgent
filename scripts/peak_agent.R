#!/usr/bin/env Rscript

`%||%` <- function(a, b) if (is.null(a)) b else a

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/peak_agent.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "R"))) repo_root <- normalizePath(getwd(), mustWork = TRUE)

source(file.path(repo_root, "R", "utils.R"))
source(file.path(repo_root, "R", "figures.R"))
source(file.path(repo_root, "R", "annotation.R"))
source(file.path(repo_root, "R", "enrichment.R"))

usage <- function() {
  cat("
PeakAnnoAgent

Usage:
  Rscript scripts/peak_agent.R <command> [options]

Commands:
  annotate   Annotate genomic peaks with gene/feature context.
  enrich     Run gene-set / pathway enrichment on clustered gene lists.
  figures    Build figures from existing annotation or enrichment tables.

Common options (all commands):
  --out DIR              Output directory                       (default: results)
  --prefix NAME          Output file prefix                     (default: today's date)
  --species NAME         mm10 | hg38 | mouse | human            (default: mm10)
  --help                 Show this message.

annotate options:
  --peaks FILE           BED/TSV with Chr, Start, End           (required)
  --methods LIST         Comma-separated. Choices:
                           chipseeker, nearest, chipeakanno, annotatr,
                           rgreat, homer, webexport             (default: chipseeker)
  --txdb PKG             Override TxDb package name.
  --orgdb PKG            Override OrgDb package name.
  --tss-upstream N       Promoter window upstream of TSS        (default: 1000)
  --tss-downstream N     Promoter window downstream of TSS      (default: 1000)
  --annotatr-annotations LIST  Override annotatr annotation codes.
  --great-rule RULE      rGREAT rule, e.g. oneClosest           (default: oneClosest)
  --great-distance N     rGREAT single-gene distance            (default: 2000)
  --homer-genome NAME    HOMER genome id (defaults from species).

enrich options:
  --genes FILE           TSV with a gene column (and optional cluster column)  (required)
  --methods LIST         Comma-separated. Choices:
                           clusterprofiler, kegg, reactome, msigdb, enrichr,
                           gprofiler, webgestalt, metascape, david, export
                                                                (default: clusterprofiler)
  --gene-col NAME        Gene symbol column                     (default: SYMBOL)
  --cluster-col NAME     Cluster/group column                   (default: Cluster)
  --universe FILE        Background gene universe (recommended for clusterprofiler/msigdb).
  --ont ONT              GO ontology: BP | MF | CC              (default: BP)
  --padj N               Adjusted p-value cutoff                (default: 0.05)
  --top N                Top terms per cluster to plot          (default: 10)
  --msigdb-collection C  MSigDB collection, e.g. H, C2          (default: H)
  --libraries LIST       Override Enrichr libraries.
  --webgestalt-database DB  WebGestaltR database name.
  --david-email ADDR     DAVID account email (or set DAVID_EMAIL).

figures options:
  --annotation FILE      Annotation table -> bar/pie summary figures.
  --enrichment FILE      Enrichment table -> dot plot.
  --padj N               Adjusted p-value cutoff for the dot plot  (default: 0.05)
  --top N                Top terms per cluster                     (default: 10)

Examples:
  Rscript scripts/peak_agent.R annotate --peaks examples/data/example_peaks.bed --out results/example --methods webexport --species mm10 --prefix example
  Rscript scripts/peak_agent.R enrich --genes examples/data/example_cluster_genes.tsv --out results/example --methods metascape,export --cluster-col Cluster --gene-col SYMBOL --prefix example
")
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$command) || isTRUE(args$help)) {
  usage()
  quit(status = if (isTRUE(args$help)) 0 else 1)
}

command <- args$command
out_dir <- arg_value(args, "out", "results")
prefix <- arg_value(args, "prefix", format(Sys.Date(), "%Y%m%d"))
species <- arg_value(args, "species", "mm10")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (command == "annotate") {
  peaks_path <- arg_value(args, "peaks")
  if (is.null(peaks_path)) stop("--peaks is required.", call. = FALSE)
  methods <- split_csv(arg_value(args, "methods", "chipseeker"))
  validate_methods(methods, "annotate")
  peaks <- read_peak_file(peaks_path)
  defaults <- infer_species_defaults(species)

  if ("chipseeker" %in% methods) {
    msg("Running ChIPseeker annotation on %d peaks.", nrow(peaks))
    run_chipseeker_annotation(
      peaks = peaks,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      txdb_pkg = arg_value(args, "txdb", NULL),
      orgdb_pkg = arg_value(args, "orgdb", NULL),
      tss_upstream = arg_value(args, "tss-upstream", 1000),
      tss_downstream = arg_value(args, "tss-downstream", 1000)
    )
  }
  if ("rgreat" %in% methods) {
    msg("Running rGREAT annotation on %d peaks.", nrow(peaks))
    run_rgreat_annotation(
      peaks = peaks,
      out_dir = out_dir,
      prefix = prefix,
      species = arg_value(args, "great-species", defaults$homer),
      rule = arg_value(args, "great-rule", "oneClosest"),
      one_distance = arg_value(args, "great-distance", 2000)
    )
  }
  if ("homer" %in% methods) {
    msg("Running HOMER annotatePeaks.pl on %s.", peaks_path)
    run_homer_annotation(peaks_path, out_dir, prefix, genome = arg_value(args, "homer-genome", defaults$homer))
  }
  if ("nearest" %in% methods || "txdb" %in% methods) {
    msg("Running nearest TxDb gene annotation on %d peaks.", nrow(peaks))
    run_nearest_txdb_annotation(
      peaks = peaks,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      txdb_pkg = arg_value(args, "txdb", NULL),
      orgdb_pkg = arg_value(args, "orgdb", NULL)
    )
  }
  if ("chipeakanno" %in% methods || "chippeakanno" %in% methods) {
    msg("Running ChIPpeakAnno annotation on %d peaks.", nrow(peaks))
    run_chipeakanno_annotation(
      peaks = peaks,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      txdb_pkg = arg_value(args, "txdb", NULL)
    )
  }
  if ("annotatr" %in% methods) {
    msg("Running annotatr annotation on %d peaks.", nrow(peaks))
    run_annotatr_annotation(
      peaks = peaks,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      annotations = split_csv(arg_value(args, "annotatr-annotations", ""))
    )
  }
  if ("webexport" %in% methods || "export" %in% methods) {
    msg("Exporting annotation web-tool BED input.")
    export_annotation_web_inputs(peaks, out_dir, prefix)
  }
} else if (command == "enrich") {
  genes_path <- arg_value(args, "genes")
  if (is.null(genes_path)) stop("--genes is required.", call. = FALSE)
  methods <- split_csv(arg_value(args, "methods", "clusterprofiler"))
  validate_methods(methods, "enrich")
  gene_col <- arg_value(args, "gene-col", "SYMBOL")
  cluster_col <- arg_value(args, "cluster-col", "Cluster")
  gene_lists <- read_gene_table(genes_path, gene_col = gene_col, cluster_col = cluster_col)
  msg("Loaded %d gene group(s): %s.", length(gene_lists), paste(names(gene_lists), collapse = ", "))

  universe <- NULL
  universe_path <- arg_value(args, "universe", NULL)
  if (!is.null(universe_path)) {
    universe <- unique(unlist(read_gene_table(universe_path, gene_col = gene_col, cluster_col = NULL)))
  }
  if ("clusterprofiler" %in% methods) {
    msg("Running clusterProfiler GO enrichment.")
    run_clusterprofiler_go(
      gene_lists = gene_lists,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      universe = universe,
      ont = arg_value(args, "ont", "BP"),
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10))
    )
  }
  if ("kegg" %in% methods || "clusterprofiler_kegg" %in% methods) {
    msg("Running clusterProfiler KEGG enrichment.")
    run_clusterprofiler_kegg(
      gene_lists = gene_lists,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10))
    )
  }
  if ("reactome" %in% methods || "reactomepa" %in% methods) {
    msg("Running ReactomePA pathway enrichment.")
    run_reactomepa(
      gene_lists = gene_lists,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10))
    )
  }
  if ("msigdb" %in% methods) {
    msg("Running MSigDB enrichment.")
    run_msigdb_enrichment(
      gene_lists = gene_lists,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      collection = arg_value(args, "msigdb-collection", "H"),
      universe = universe,
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10))
    )
  }
  if ("enrichr" %in% methods) {
    msg("Running Enrichr enrichment.")
    libraries <- split_csv(arg_value(args, "libraries", ""))
    run_enrichr(
      gene_lists = gene_lists,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      libraries = libraries,
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10))
    )
  }
  if ("gprofiler" %in% methods || "gprofiler2" %in% methods) {
    msg("Running g:Profiler enrichment.")
    run_gprofiler2(
      gene_lists = gene_lists,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10))
    )
  }
  if ("webgestalt" %in% methods || "webgestaltr" %in% methods) {
    msg("Running WebGestaltR enrichment.")
    run_webgestaltr(
      gene_lists = gene_lists,
      out_dir = out_dir,
      prefix = prefix,
      species = species,
      database = arg_value(args, "webgestalt-database", "geneontology_Biological_Process"),
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10))
    )
  }
  if ("metascape" %in% methods) {
    msg("Exporting Metascape-ready input table.")
    export_metascape_input(gene_lists, out_dir, prefix)
  }
  if ("export" %in% methods || "genelist" %in% methods) {
    msg("Exporting generic long and wide gene-list tables.")
    export_gene_lists(gene_lists, out_dir, prefix)
  }
  if ("david" %in% methods) {
    msg("Running DAVID enrichment.")
    run_david_enrichment(gene_lists, out_dir, prefix, species = species, email = arg_value(args, "david-email", Sys.getenv("DAVID_EMAIL")))
  }
} else if (command == "figures") {
  annotation <- arg_value(args, "annotation", NULL)
  enrichment <- arg_value(args, "enrichment", NULL)
  if (!is.null(annotation)) {
    plot_annotation_summary(annotation, out_dir, prefix)
  }
  if (!is.null(enrichment)) {
    plot_enrichment_dotplot(enrichment, file.path(out_dir, paste0(prefix, "_enrichment_dotplot.pdf")),
      padj_cutoff = as.numeric(arg_value(args, "padj", 0.05)),
      top_show = as.integer(arg_value(args, "top", 10)))
  }
} else {
  usage()
  stop(sprintf("Unknown command: %s", command), call. = FALSE)
}

msg("Done. Outputs are in %s", normalizePath(out_dir, mustWork = FALSE))
