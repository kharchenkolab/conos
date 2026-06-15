## Annotation transfer (propagateLabels) — recovers held-out cluster labels well above the majority
## baseline, with bounded confidence; the diffusion and solver engines agree. Bundled real panel.

test_that("propagateLabels recovers held-out labels (diffusion and solver) with bounded uncertainty", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)
  con$runClustering(verbose = FALSE)

  truth <- con$clusters$leiden$groups
  set.seed(1)
  hidden <- sample(names(truth), round(0.3 * length(truth)))
  known  <- truth[setdiff(names(truth), hidden)]
  baseline <- max(table(truth)) / length(truth)

  rd <- con$propagateLabels(labels = known, method = "diffusion", verbose = FALSE)
  expect_setequal(names(rd$labels), names(truth))                       # a label for every cell
  expect_true(all(rd$uncertainty >= 0 & rd$uncertainty <= 1))           # confidence in [0,1]
  acc.d <- mean(rd$labels[hidden] == as.character(truth[hidden]))
  expect_gt(acc.d, baseline)                                            # better than guessing the majority

  rs <- con$propagateLabels(labels = known, method = "solver", solver = "Matrix")
  acc.s <- mean(rs$labels[hidden] == as.character(truth[hidden]))
  expect_gt(acc.s, baseline)
  expect_gt(mean(rd$labels[hidden] == rs$labels[hidden]), 0.7)          # the two engines largely agree
})
