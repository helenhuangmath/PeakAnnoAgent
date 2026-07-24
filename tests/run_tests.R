#!/usr/bin/env Rscript
# Standalone test runner. The project is used as sourced scripts (not yet an
# installed package), so we source R/ into the global env and run testthat
# against the offline code paths.
#
#   Rscript tests/run_tests.R

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Install 'testthat' to run the test suite: install.packages('testthat')", call. = FALSE)
}

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg)) {
  this_file <- sub("^--file=", "", file_arg[1])
  repo_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
} else if (dir.exists(file.path(getwd(), "R"))) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
} else {
  repo_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
}

if (dir.exists(file.path(repo_root, "R"))) {
  for (f in c("utils.R", "figures.R", "annotation.R", "enrichment.R")) {
    source(file.path(repo_root, "R", f))
  }
} else {
  library(PeakAnnoAgent)
}

test_dir <- file.path(repo_root, "tests", "testthat")
if (!dir.exists(test_dir)) {
  test_dir <- file.path(getwd(), "tests", "testthat")
}

reporter <- testthat::SummaryReporter$new()
results <- testthat::test_dir(
  test_dir,
  reporter = reporter,
  stop_on_failure = TRUE
)
