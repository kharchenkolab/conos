## Phase 2 — planIntegration() commonality assessment + getModalities/getDefaultModality (§3.3 / §4.3).
## Multi-facet cases skip-gate on pagoda2.1 (facets); the single-modality case uses the real bundled panel.

test_that("planIntegration assesses the single common modality on the real panel (§3.3)", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(small_panel.preprocessed, n.cores = 1)
  plan <- con$planIntegration(verbose = FALSE)
  expect_s3_class(plan, "conosIntegrationPlan")
  expect_true("RNA" %in% plan$modality)
  rna <- plan[plan$modality == "RNA", ]
  expect_gt(rna$common.features, 0)
  expect_match(rna$verdict, "usable")
  expect_identical(con$misc$integration.plan$default.modality, "RNA") # resolved default recorded on the object
})

test_that("getModalities / getDefaultModality report RNA for legacy single-modality samples", {
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  s <- small_panel.preprocessed[[1]]
  expect_identical(getModalities(s), "RNA")
  expect_identical(getDefaultModality(s), "RNA")
})

## build a pagoda2.1 sample with an RNA + ADT facet; returns NULL on legacy pagoda2 (no facets) -> skip
.mk_multifacet <- function(seed, adt.feats) {
  set.seed(seed); ng <- 60; nc <- 90
  base <- matrix(rpois(ng * nc, 2), ng, nc, dimnames = list(paste0("g", 1:ng), paste0(seed, "c", 1:nc)))
  p <- tryCatch(pagoda2::Pagoda2$new(as(Matrix::Matrix(base, sparse = TRUE), "dgCMatrix"), verbose = FALSE,
        n.cores = 1, min.cells.per.gene = 0, min.transcripts.per.cell = 0, trim = 0, log.scale = TRUE),
        error = function(e) NULL)
  if (is.null(p) || !is.function(tryCatch(p$addFacet, error = function(e) NULL))) return(NULL)
  adt <- matrix(rpois(nc * length(adt.feats), 4) + 1L, nc, length(adt.feats),
    dimnames = list(paste0(seed, "c", 1:nc), adt.feats))
  p$addFacet("ADT", as(Matrix::Matrix(adt, sparse = TRUE), "dgCMatrix"), modelType = "plain", featureType = "protein")
  p
}

test_that("planIntegration flags multi-facet commonality: shared -> usable, disjoint -> not-integrable (§3.3)", {
  skip_if_not_installed("pagoda2")
  a <- .mk_multifacet(1, paste0("P", 1:10)); b <- .mk_multifacet(2, paste0("P", 1:10))
  skip_if(is.null(a) || is.null(b), "pagoda2.1 facets not available")

  con <- Conos$new(list(a = a, b = b), n.cores = 1)
  plan <- con$planIntegration(verbose = FALSE)
  expect_setequal(plan$modality, c("RNA", "ADT"))
  expect_match(plan$verdict[plan$modality == "ADT"], "usable") # 10 fully-shared proteins -> usable (overlap-primary)

  ad <- .mk_multifacet(1, paste0("Pa", 1:10)); bd <- .mk_multifacet(2, paste0("Pb", 1:10)) # disjoint ADT features
  con2 <- Conos$new(list(a = ad, b = bd), n.cores = 1)
  plan2 <- con2$planIntegration(verbose = FALSE)
  expect_match(plan2$verdict[plan2$modality == "ADT"], "not-integrable")
})
