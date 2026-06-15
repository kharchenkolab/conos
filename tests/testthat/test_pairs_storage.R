## buildGraph(pairs.storage=) controls retention of the O(n^2) per-pair rotations (proposal §1.2).

test_that("buildGraph pairs.storage keeps or drops the per-pair rotations", {
  skip_if_not_installed("pagoda2")
  suppressMessages(library(pagoda2))
  mk <- function(seed, ncell = 90L) {
    set.seed(seed)
    n.pop <- 3L; cpp <- ncell / n.pop; nm <- 30L; nbg <- 200L; ng <- n.pop * nm + nbg
    pop <- rep(seq_len(n.pop), each = cpp)
    lam <- matrix(0.5, ng, ncell)
    for (k in seq_len(n.pop)) { gi <- ((k - 1) * nm + 1):(k * nm); lam[gi, pop == k] <- 4 }
    lam <- sweep(lam, 2, runif(ncell, 0.8, 1.2), "*")
    cnt <- matrix(rpois(length(lam), lam), ng, ncell)
    rownames(cnt) <- c(paste0("Mk", seq_len(n.pop * nm)), paste0("Bg", seq_len(nbg)))
    colnames(cnt) <- paste0("s", seed, "c", seq_len(ncell))
    cnt <- as(Matrix::Matrix(cnt, sparse = TRUE), "dgCMatrix")
    p2 <- Pagoda2$new(cnt, log.scale = TRUE, min.cells.per.gene = 5, min.transcripts.per.cell = 10,
                      n.cores = 1, verbose = FALSE)
    p2$runVariance(use.raw.variance = TRUE, verbose = FALSE)
    p2$runReduction(nPcs = 10, n.odgenes = 100, verbose = FALSE)
    p2
  }

  ## default "keep": rotations retained, graph built
  con <- Conos$new(list(s1 = mk(1), s2 = mk(2)), n.cores = 1)
  con$buildGraph(k = 10, k.self = 5, space = "PCA", ncomps = 10, n.odgenes = 100, verbose = FALSE)
  expect_gt(length(con$pairs[["PCA"]]), 0)
  expect_gt(igraph::vcount(con$graph), 0)

  ## "drop": rotations freed, graph still built
  con2 <- Conos$new(list(s1 = mk(1), s2 = mk(2)), n.cores = 1)
  con2$buildGraph(k = 10, k.self = 5, space = "PCA", ncomps = 10, n.odgenes = 100, verbose = FALSE,
                  pairs.storage = "drop")
  expect_null(con2$pairs[["PCA"]])
  expect_gt(igraph::vcount(con2$graph), 0)
})
