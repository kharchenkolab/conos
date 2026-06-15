## Phase 3 — memory-lean merge (§1.4) + sparse SNN weighting (§1.3), and the snn.quantile && bug fix.

test_that("mergeCountMatrices is sparse and matches dense zero-padding for differing gene sets (§1.4)", {
  set.seed(1)
  a <- as(Matrix::Matrix(rpois(40 * 30, 1), 30, 40,
    dimnames = list(paste0("c", 1:30), paste0("g", 1:40)), sparse = TRUE), "CsparseMatrix")
  b <- as(Matrix::Matrix(rpois(40 * 30, 1), 30, 40,
    dimnames = list(paste0("d", 1:30), paste0("g", 21:60)), sparse = TRUE), "CsparseMatrix")
  ext_dense <- function(mtx, cn) { nn <- setdiff(cn, colnames(mtx)); e <- matrix(0, nrow(mtx), length(nn)); colnames(e) <- nn; cbind(mtx, e)[, cn] }
  gu <- Reduce(union, lapply(list(a, b), colnames))
  ref <- as(Reduce(rbind, lapply(list(a, b), ext_dense, gu)), "CsparseMatrix") # old dense-padded reference

  m <- conos:::mergeCountMatrices(list(a, b), transposed = TRUE)
  expect_s4_class(m, "CsparseMatrix")                                          # stays sparse (no densification)
  expect_equal(as.matrix(m[rownames(ref), colnames(ref)]), as.matrix(ref))     # values identical to dense padding
})

test_that(".conos_snn_jaccard equals the dense numerator/outer(pmin) form, but sparse (§1.3)", {
  set.seed(7); n1 <- 50; n2 <- 60
  m1 <- as(Matrix::sparseMatrix(i = sample(n1, 200, TRUE), j = sample(40, 200, TRUE), x = 1, dims = c(n1, 40)), "CsparseMatrix")
  m2 <- as(Matrix::sparseMatrix(i = sample(40, 200, TRUE), j = sample(n2, 200, TRUE), x = 1, dims = c(40, n2)), "CsparseMatrix")
  mnn1 <- as(Matrix::sparseMatrix(i = sample(n1, 150, TRUE), j = sample(n2, 150, TRUE), x = 1, dims = c(n1, n2)), "CsparseMatrix")
  dense <- as.matrix(((m1 %*% m2) * mnn1) / pmax(outer(Matrix::rowSums(m1), Matrix::colSums(m2), FUN = pmin), 1))
  sparse <- conos:::.conos_snn_jaccard(m1, m2, mnn1)
  expect_s4_class(sparse, "CsparseMatrix")
  expect_equal(as.matrix(sparse), dense)
})

test_that("buildGraph(snn=TRUE) runs for scalar and length-2 snn.quantile (regression: && / is.na length>1)", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, snn = TRUE, verbose = FALSE)
  expect_gt(igraph::ecount(con$graph), 0)

  con2 <- Conos$new(small_panel.preprocessed, n.cores = 1)
  con2$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000,
                  snn = TRUE, snn.quantile = c(0.1, 0.9), min.snn.jaccard = 0.1, verbose = FALSE)
  expect_gt(igraph::ecount(con2$graph), 0)
})
