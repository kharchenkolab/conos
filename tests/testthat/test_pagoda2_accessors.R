## Verify conos's Pagoda2 sample accessors against pagoda2.1's name-keyed storage (proposal §4.1).
## The accessors must use the object's canonical keys (defaults$reduction, getGrouping, the embedding
## name) rather than a hardcoded "PCA", and still work for a standard run() object.

test_that("conos Pagoda2 accessors use pagoda2.1 canonical keys", {
  skip_if_not_installed("pagoda2")
  suppressMessages(library(pagoda2))
  set.seed(42)
  n.pop <- 3L; cpp <- 80L; nc <- n.pop * cpp; nm <- 40L; nbg <- 250L; ng <- n.pop * nm + nbg
  pop <- rep(seq_len(n.pop), each = cpp)
  lam <- matrix(0.5, ng, nc)
  for (k in seq_len(n.pop)) { gi <- ((k - 1) * nm + 1):(k * nm); lam[gi, pop == k] <- 4 }
  lam <- sweep(lam, 2, runif(nc, 0.8, 1.2), "*")
  cnt <- matrix(rpois(length(lam), lam), ng, nc)
  rownames(cnt) <- c(paste0("Mk", seq_len(n.pop * nm)), paste0("Bg", seq_len(nbg)))
  colnames(cnt) <- paste0("c", seq_len(nc))
  cnt <- as(Matrix::Matrix(cnt, sparse = TRUE), "dgCMatrix")

  p2 <- Pagoda2$new(cnt, log.scale = TRUE, min.cells.per.gene = 5, min.transcripts.per.cell = 10,
                    n.cores = 1, verbose = FALSE)
  p2$run(plots = "none", verbose = FALSE, variance = list(use.raw.variance = TRUE),
         pca = list(nPcs = 12, n.odgenes = 250), graph = list(k = 15),
         embedding = list(n_neighbors = 15, min_dist = 0.3))

  ## getPca returns the default reduction (canonical key, not a hardcoded "PCA")
  pca <- getPca(p2)
  expect_false(is.null(pca))
  expect_identical(pca, p2$reductions[[p2$defaults$reduction]])

  ## getEmbedding finds the run() embedding by name across reduction namespaces
  emb <- getEmbedding(p2, "UMAP")
  expect_false(is.null(emb))
  expect_equal(ncol(emb), 2L)

  ## getClustering returns the grouping factor (via getGrouping), matching the default grouping
  cl <- getClustering(p2, "leiden")
  expect_false(is.null(cl))
  expect_equal(length(cl), nrow(p2$getRawCounts()))
})
