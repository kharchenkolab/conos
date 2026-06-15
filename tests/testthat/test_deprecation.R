## Phase 5 (§7.4) — run* verbs are preferred; the old names are deprecated-but-working aliases.

## TRUE iff `expr` emits a deprecation warning (robust to other warnings, e.g. uwot's)
.warned_deprecated <- function(expr) {
  hit <- FALSE
  withCallingHandlers(force(expr),
    warning = function(w) { if (grepl("deprecat|run[A-Z]", conditionMessage(w))) hit <<- TRUE; invokeRestart("muffleWarning") })
  hit
}

test_that("deprecated method aliases warn but still do the work", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)

  expect_true(.warned_deprecated(con$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)))
  expect_gt(igraph::vcount(con$graph), 0)                                  # the alias still built the graph

  expect_true(.warned_deprecated(con$findCommunities(method = "leiden", verbose = FALSE)))
  expect_true("leiden" %in% names(con$clusters))

  expect_true(.warned_deprecated(de <- con$getDifferentialGenes(verbose = FALSE)))
  expect_true(is.list(de))

  expect_true(.warned_deprecated(con$embedGraph(method = "UMAP", verbose = FALSE)))
  expect_equal(nrow(con$embedding), length(getCellNames(con)))
})

test_that("the preferred run* verbs do not emit a deprecation warning", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  expect_false(.warned_deprecated(con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)))
  expect_false(.warned_deprecated(con$runClustering(verbose = FALSE)))
})
