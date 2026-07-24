#!/usr/bin/env Rscript
# Standalone test runner. The project is used as sourced scripts (not yet an
# installed package), so we source R/ into the global env and run testthat
# against the offline code paths.
#
#   Rscript tests/run_tests.R

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install 'testthat' to run the test suite: install.packages('testthat')", call. = FALSE)
}

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)

for (f in c("utils.R", "figures.R", "annotation.R", "enrichment.R")) {
  source(file.path(repo_root, "R", f))
}

reporter <- testthat::SummaryReporter$new()
results <- testthat::test_dir(
  file.path(repo_root, "tests", "testthat"),
  reporter = reporter,
  stop_on_failure = TRUE
)
