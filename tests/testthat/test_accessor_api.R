## Phase 1 — findCommunities(method=) string support (§7.1) and symmetric Conos accessors (§7.2).
## Uses the bundled real MantonBM panel.

test_that("findCommunities accepts a string method as well as a function (§7.1)", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  con$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)

  con$findCommunities(method = "leiden", verbose = FALSE)          # string
  expect_true("leiden" %in% names(con$clusters))
  expect_gt(length(unique(con$clusters[["leiden"]]$groups)), 1L)

  con$findCommunities(method = igraph::cluster_walktrap, verbose = FALSE) # function still works
  expect_true(length(names(con$clusters)) >= 2L)

  expect_error(con$findCommunities(method = "not_a_method"), "unknown community detection method")
})

test_that(".conos_resolve_community_method maps names (case/punctuation-insensitive) and rejects unknowns", {
  expect_identical(conos:::.conos_resolve_community_method("Leiden"), leiden.community)
  expect_identical(conos:::.conos_resolve_community_method("label.prop"), igraph::cluster_label_prop)
  expect_identical(conos:::.conos_resolve_community_method("multilevel"), igraph::cluster_louvain)
  expect_error(conos:::.conos_resolve_community_method("xyz"), "unknown")
})

test_that("Conos accessors are symmetric: getPca/getCountMatrix/getEmbedding (§7.2)", {
  skip_if_not_installed("pagoda2")
  skip_if_not_installed("uwot")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  con$buildGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)

  expect_null(getPca(con))                                          # no joint PCA -> NULL (not a dispatch error)

  cm <- getCountMatrix(con)                                         # joint (normalized) count matrix
  expect_true(inherits(cm, "Matrix") || is.matrix(cm))
  expect_equal(nrow(cm), length(getCellNames(con)))

  suppressWarnings(con$embedGraph(method = "UMAP", verbose = FALSE))
  expect_equal(nrow(getEmbedding(con)), length(getCellNames(con))) # default joint embedding
})
