## buildGraph(pairs.storage=) controls retention of the O(n^2) per-pair rotations (proposal §1.2/§5).
## Uses the real bundled MantonBM panel rather than synthetic data.

test_that("buildGraph pairs.storage keep/drop/disk manage the per-pair rotations", {
  skip_if_not_installed("pagoda2")  # panel samples are Pagoda2 objects
  data("small_panel.preprocessed", package = "conos", envir = environment())
  p <- small_panel.preprocessed
  bg <- function(con, storage)
    con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000,
                   verbose = FALSE, pairs.storage = storage)

  ## keep (default): rotations retained, graph built
  con <- Conos$new(p, n.cores = 1)
  bg(con, "keep")
  expect_gt(length(con$pairs[["PCA"]]), 0)
  expect_gt(igraph::vcount(con$graph), 0)

  ## drop: rotations freed, graph still built
  con2 <- Conos$new(p, n.cores = 1)
  bg(con2, "drop")
  expect_null(con2$pairs[["PCA"]])
  expect_gt(igraph::vcount(con2$graph), 0)

  ## disk: rotations off RAM but persisted to a scratch file; graph built
  con3 <- Conos$new(p, n.cores = 1)
  bg(con3, "disk")
  expect_null(con3$pairs[["PCA"]])
  dp <- con3$misc[["pairs.disk"]][["PCA"]]
  expect_true(!is.null(dp) && file.exists(dp))
  expect_gt(igraph::vcount(con3$graph), 0)

  ## re-run restores from disk (reused, not recomputed) and is retained under "keep"
  bg(con3, "keep")
  expect_gt(length(con3$pairs[["PCA"]]), 0)
})
