# Pathway network construction (needs igraph; plotting is not exercised here).

network_table <- function() {
  data.table::data.table(
    method = "go_bp",
    group = c("promoter", "promoter", "promoter", "distal", "distal"),
    term = c("T cell activation", "T cell differentiation", "lymphocyte activation", "T cell activation",
      "cytokine production"),
    padj = c(1e-10, 1e-8, 1e-7, 1e-6, 1e-4),
    count = c(4L, 4L, 4L, 2L, 3L),
    genes = c("Tox;Pdcd1;Batf;Tcf7", "Tox;Pdcd1;Batf;Lef1", "Tox;Pdcd1;Batf;Cd8a", "Ifng;Tox",
      "Ifng;Tnf;Il2")
  )
}

test_that("enrichment nodes pool genes and groups per term", {
  nodes <- enrichment_nodes(network_table())
  ta <- nodes[term == "T cell activation"]
  expect_equal(ta$padj, 1e-10)
  expect_equal(ta$groups, "distal;promoter")
  expect_setequal(strsplit(ta$genes, ";")[[1]], c("Tox", "Pdcd1", "Batf", "Tcf7", "Ifng"))
})

test_that("enrichment map links overlapping terms and separates unrelated ones", {
  skip_if_not_installed("igraph")
  em <- build_enrichment_map(network_table(), min_jaccard = 0.25, k = 3)
  expect_equal(nrow(em$nodes), 4)
  expect_true(all(em$edges$jaccard >= 0.25))
  mod <- stats::setNames(em$nodes$module, em$nodes$term)
  expect_equal(mod[["T cell activation"]], mod[["T cell differentiation"]])
  expect_false(mod[["cytokine production"]] == mod[["T cell differentiation"]])
  expect_equal(em$nodes[module == 1, module_label][1], "T cell activation")   # module named by best term
})

test_that("k-nearest sparsification caps edges per term", {
  skip_if_not_installed("igraph")
  genes <- "A;B;C;D;E"
  tab <- data.table::data.table(method = "m", group = "g", term = paste0("t", 1:8), padj = 1e-5 * (1:8),
    count = 5L, genes = genes)
  em <- build_enrichment_map(tab, k = 2)
  degree <- table(c(em$edges$from, em$edges$to))
  expect_lt(nrow(em$edges), choose(8, 2))
  expect_true(max(degree) <= 2 * 2 + 2)
})

test_that("MSigDB names are made readable", {
  expect_equal(pretty_term("HALLMARK_TNFA_SIGNALING_VIA_NFKB"), "TNFA signaling via NFKB")
  expect_equal(pretty_term("T cell activation"), "T cell activation")
  expect_equal(pretty_term("GSE7218_UNSTIM_VS_ANTIGEN_STIM_BCELL_DN"), "GSE7218 unstim vs antigen STIM bcell dn")
  expect_equal(pretty_term("KAECH_NAIVE_VS_DAY15_EFF_CD8_TCELL_DN"), "Kaech naive vs DAY15 EFF CD8 tcell dn")
})

test_that("method names become safe file names", {
  expect_equal(file_tag("msigdb_C7:IMMUNESIGDB"), "msigdb_C7-IMMUNESIGDB")
  expect_equal(file_tag("go_bp"), "go_bp")
})
