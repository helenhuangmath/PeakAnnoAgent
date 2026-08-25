# Offline tests for file reading and the export paths that need no Bioconductor.

repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
example_dir <- file.path(repo_root, "examples", "data")

test_that("read_peak_file normalises Chr/Start/End and coerces coordinates", {
  peaks <- read_peak_file(file.path(example_dir, "example_peaks.bed"))
  expect_true(all(c("Chr", "Start", "End") %in% names(peaks)))
  expect_type(peaks$Start, "integer")
  expect_type(peaks$End, "integer")
  expect_gt(nrow(peaks), 0)
})

test_that("read_peak_file rejects files with fewer than three columns", {
  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t100", "chr2\t200"), tmp)
  expect_error(read_peak_file(tmp), "at least three columns")
})

test_that("read_gene_table splits genes by cluster and dedupes", {
  groups <- read_gene_table(file.path(example_dir, "example_cluster_genes.tsv"),
    gene_col = "SYMBOL", cluster_col = "Cluster")
  expect_true(length(groups) >= 2)
  expect_true(all(vapply(groups, function(g) !any(duplicated(g)), logical(1))))
})

test_that("read_gene_table keeps genes in the cluster they came from", {
  # Regression test: a gene shared by two clusters used to make the gene vector
  # shorter than the cluster vector, so split() recycled it and silently moved
  # genes into the wrong clusters (or emptied clusters entirely).
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(
    "Cluster\tSYMBOL",
    "Cluster1\tSox2",
    "Cluster1\tPax6",
    "Cluster1\tSox2",
    "Cluster2\tSox2",
    "Cluster2\tIl6",
    "Cluster3\tTnf"
  ), tmp)

  groups <- read_gene_table(tmp, gene_col = "SYMBOL", cluster_col = "Cluster")

  expect_equal(sort(names(groups)), c("Cluster1", "Cluster2", "Cluster3"))
  expect_setequal(groups$Cluster1, c("Sox2", "Pax6"))
  expect_setequal(groups$Cluster2, c("Sox2", "Il6"))
  expect_setequal(groups$Cluster3, "Tnf")
})

test_that("read_gene_table collapses to a single 'All' group without a cluster column", {
  groups <- read_gene_table(file.path(example_dir, "example_cluster_genes.tsv"),
    gene_col = "SYMBOL", cluster_col = NULL)
  expect_equal(names(groups), "All")
})

test_that("read_gene_table errors on a missing gene column", {
  expect_error(
    read_gene_table(file.path(example_dir, "example_cluster_genes.tsv"), gene_col = "GENE"),
    "was not found"
  )
})

test_that("export_metascape_input writes a Group/Gene table", {
  out <- tempfile("mscape")
  dir.create(out)
  gene_lists <- list(Cluster1 = c("Sox2", "Pax6"), Cluster2 = c("Il6", "Tnf"))
  path <- export_metascape_input(gene_lists, out, "unit")
  expect_true(file.exists(path))
  dt <- data.table::fread(path)
  expect_equal(names(dt), c("Group", "Gene"))
  expect_equal(nrow(dt), 4)
})

test_that("export_gene_lists writes both long and wide tables", {
  out <- tempfile("glists")
  dir.create(out)
  gene_lists <- list(Cluster1 = c("Sox2", "Pax6", "Hes1"), Cluster2 = c("Il6"))
  res <- export_gene_lists(gene_lists, out, "unit")
  expect_true(file.exists(res$long))
  expect_true(file.exists(res$wide))
  wide <- data.table::fread(res$wide, na.strings = "")
  expect_equal(nrow(wide), 3)                 # padded to the longest list
  expect_true(is.na(wide$Cluster2[3]))        # shorter list padded with NA
})
