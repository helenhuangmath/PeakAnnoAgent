# Offline tests for the agent's R engine (QC, signature tests, motif parsing).

repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."), mustWork = TRUE)
example_peaks <- file.path(repo_root, "examples", "data", "example_tex_peaks_mm10.bed")

test_that("read_peak_file assigns stable peak IDs", {
  peaks <- read_peak_file(example_peaks)
  expect_equal(peaks$PeakID[1], "tex_peak_0001")
  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t10\t20", "chr1\t30\t40"), tmp)
  expect_equal(read_peak_file(tmp)$PeakID, c("peak_00001", "peak_00002"))
})

test_that("run_peak_qc cleans contigs, duplicates and blacklisted peaks", {
  tmp <- tempfile(fileext = ".bed")
  writeLines(c("chr1\t100\t400\ta", "chr1\t100\t400\tb", "chrUn_JH584304\t5\t50\tc",
    "chr2\t1000\t1300\td", "chr3\t500\t900\te"), tmp)
  bl <- tempfile(fileext = ".bed")
  writeLines("chr2\t1200\t5000", bl)
  out <- tempfile("qc")
  res <- run_peak_qc(read_peak_file(tmp), out, "t", blacklist = bl)
  expect_equal(res$n_peaks, 5)
  expect_equal(res$n_duplicated, 1)
  expect_equal(res$n_nonstandard_chrom, 1)
  expect_equal(res$n_blacklisted, 1)
  expect_equal(res$n_clean_peaks, 2)
  expect_true(file.exists(res$figures$width))
  expect_true(any(grepl("Only 5 peaks", res$warnings)))
})

test_that("signature_hypergeom uses the full background and keeps small sets", {
  universe <- paste0("G", 1:1000)
  t2g <- data.table::data.table(term = c(rep("A", 5), rep("B", 3)),
    gene = c("G1", "G2", "G3", "G4", "G5", "G900", "G901", "G902"))
  res <- signature_hypergeom(list(x = c("G1", "G2", "G3", "G10")), t2g, universe)
  a <- res[ID == "A"]
  expect_equal(a$Count, 3)
  expect_equal(a$BgRatio, "5/1000")
  expect_equal(a$pvalue, stats::phyper(2, 5, 995, 4, lower.tail = FALSE))
  expect_equal(res[ID == "B"]$Count, 0)
})

test_that("motif names are reduced to TF symbols", {
  expect_equal(motif_tf_name("CTCF(Zf)/CD4+-CTCF-ChIP-Seq(Barski_et_al.)/Homer"), "CTCF")
  expect_equal(motif_tf_name("BATF(bZIP)/Th17-BATF-ChIP-Seq/Homer"), "BATF")
  expect_equal(motif_tf_name("Nr4a1"), "Nr4a1")
})

test_that("HOMER knownResults.txt is standardised", {
  tmp <- tempfile(fileext = ".txt")
  writeLines(c(paste("Motif Name", "Consensus", "P-value", "Log P-value", "q-value (Benjamini)",
    "# of Target Sequences with Motif(of 300)", "% of Target Sequences with Motif",
    "# of Background Sequences with Motif(of 40000)", "% of Background Sequences with Motif", sep = "\t"),
    paste("BATF(bZIP)/Th17-BATF-ChIP-Seq/Homer", "DATGASTCAT", "1e-20", "-46.05", "0.0000", "90", "30.00%",
      "4000", "10.00%", sep = "\t")), tmp)
  res <- read_homer_known_results(tmp)
  expect_equal(res$tf, "BATF")
  expect_equal(res$family, "bZIP")
  expect_equal(res$pct_target, 30)
  expect_gt(res$log2_enrichment, 1.5)
})

test_that("agent tool dispatch rejects unknown tools", {
  expect_error(run_agent_tool("nope", list()), "Unknown tool")
  expect_true("annotate_peaks" %in% agent_tool_names())
})
