## Phase 4 — §6 robustness fixes: scanKModularity must return its data AND actually draw under plot=TRUE.

test_that("scanKModularity returns the k/modularity table and draws under plot=TRUE", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)

  ks <- scanKModularity(con, min = 4, max = 6, by = 1, plot = FALSE, verbose = FALSE,
                        k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000)
  expect_s3_class(ks, "data.frame")
  expect_identical(colnames(ks), c("k", "m"))               # data returned regardless of plot=
  expect_equal(nrow(ks), 3L)
  expect_true(all(is.finite(ks$m)))
})
