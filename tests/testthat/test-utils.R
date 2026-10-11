# Offline unit tests for the parsing / configuration helpers in R/utils.R.

test_that("split_csv handles empty, single, and multi-value input", {
  expect_equal(split_csv(""), character())
  expect_equal(split_csv(NULL), character())
  expect_equal(split_csv(NA), character())
  expect_equal(split_csv("a"), "a")
  expect_equal(split_csv("a,b,c"), c("a", "b", "c"))
  expect_equal(split_csv(" a , b "), c("a", "b"))
})

test_that("parse_cli_args separates command, valued flags, and boolean flags", {
  args <- parse_cli_args(c("enrich", "--genes", "g.tsv", "--methods", "metascape", "--help"))
  expect_equal(args$command, "enrich")
  expect_equal(args$genes, "g.tsv")
  expect_equal(args$methods, "metascape")
  expect_true(isTRUE(args$help))
})

test_that("parse_cli_args returns NULL command when only flags are given", {
  args <- parse_cli_args(c("--help"))
  expect_null(args$command)
  expect_true(isTRUE(args$help))
})

test_that("arg_value falls back to the default for missing or empty values", {
  args <- list(out = "results", empty = "")
  expect_equal(arg_value(args, "out", "default"), "results")
  expect_equal(arg_value(args, "empty", "default"), "default")
  expect_equal(arg_value(args, "missing", "default"), "default")
  expect_null(arg_value(args, "missing"))
})

test_that("infer_species_defaults maps shortcuts and rejects unknown species", {
  mm <- infer_species_defaults("mm10")
  expect_equal(mm$orgdb, "org.Mm.eg.db")
  expect_equal(mm$kegg, "mmu")
  expect_equal(infer_species_defaults("mouse")$homer, "mm10")

  hg <- infer_species_defaults("hg38")
  expect_equal(hg$orgdb, "org.Hs.eg.db")
  expect_equal(hg$gprofiler, "hsapiens")
  expect_equal(infer_species_defaults("human")$msigdb, "Homo sapiens")

  expect_equal(infer_species_defaults("mm39")$txdb, "TxDb.Mmusculus.UCSC.mm39.knownGene")
  expect_equal(infer_species_defaults("hg19")$bsgenome, "BSgenome.Hsapiens.UCSC.hg19")
  expect_equal(infer_species_defaults("mm10")$bsgenome, "BSgenome.Mmusculus.UCSC.mm10")
  expect_error(infer_species_defaults("dog"), "Supported species")
})

test_that("validate_methods accepts known methods and rejects typos", {
  expect_silent(validate_methods(c("chipseeker", "rgreat"), "annotate"))
  expect_silent(validate_methods(c("clusterprofiler", "metascape", "export"), "enrich"))
  expect_silent(validate_methods(character(), "annotate"))

  expect_error(validate_methods("clusterprofilr", "enrich"), "Unknown enrich method")
  expect_error(validate_methods(c("chipseeker", "bogus"), "annotate"), "bogus")
})

test_that("save_figure keeps only the PDF in the results folder", {
  out <- tempfile("fig")
  paths <- save_figure(function() graphics::plot(1:3), out, "demo")
  expect_true(file.exists(paths[["pdf"]]))
  expect_equal(dirname(paths[["pdf"]]), out)
  expect_true(file.exists(paths[["png"]]))
  expect_equal(basename(dirname(paths[["png"]])), ".preview")
  expect_equal(list.files(out, pattern = "\\.png$"), character())
})
