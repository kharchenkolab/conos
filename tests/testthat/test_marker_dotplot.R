## Phase 4 — Conos$plotMarkerDotPlot (§6): top per-cluster markers via sccore::dotPlot, on the real panel.

test_that("plotMarkerDotPlot returns a ggplot of top per-cluster markers", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  con$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)
  con$findCommunities(method = "leiden", verbose = FALSE)

  p <- con$plotMarkerDotPlot(n.genes.per.group = 3, z.threshold = 0.5)
  expect_s3_class(p, "ggplot")

  ## explicit grouping also works
  g <- con$clusters[[1]]$groups
  expect_s3_class(con$plotMarkerDotPlot(groups = g, n.genes.per.group = 2, z.threshold = 0), "ggplot")
})

test_that("plotMarkerDotPlot errors clearly when there is no clustering", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  expect_error(con$plotMarkerDotPlot(), "no clustering available")
})
