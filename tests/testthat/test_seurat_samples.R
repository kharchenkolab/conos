## Conos on Seurat samples: the S4 accessor layer and the standard joint workflow run on
## CreateSeuratObject-built samples, including under SeuratObject >= 5 (the GetAssayData `layer=` path that
## replaced the now-defunct `slot=`). Skipped when Seurat is unavailable.

make_seurat_sample <- function(p2, npcs = 20) {
  cm <- getRawCountMatrix(p2, transposed = FALSE)          # genes x cells
  colnames(cm) <- make.unique(colnames(cm)); rownames(cm) <- make.unique(rownames(cm))
  s <- Seurat::CreateSeuratObject(counts = cm)
  s <- Seurat::NormalizeData(s, verbose = FALSE)
  s <- Seurat::FindVariableFeatures(s, nfeatures = 1000, verbose = FALSE)
  s <- Seurat::ScaleData(s, verbose = FALSE)
  Seurat::RunPCA(s, npcs = npcs, verbose = FALSE)
}

test_that("the accessor layer works on a Seurat sample (incl. SeuratObject >= 5)", {
  skip_if_not_installed("Seurat")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  s <- make_seurat_sample(small_panel.preprocessed[[1]])

  expect_equal(nrow(getPca(s)), length(getCellNames(s)))            # cells x PCs
  expect_equal(ncol(getRawCountMatrix(s)), length(getCellNames(s))) # genes x cells -> defunct slot= would error
  expect_equal(ncol(getCountMatrix(s)), length(getCellNames(s)))    # normalized expression
  expect_true(length(getOverdispersedGenes(s, 500)) > 0)
  expect_identical(getModalities(s), "RNA")
})

test_that("Conos integrates a panel of Seurat objects", {
  skip_if_not_installed("Seurat")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  seurat.panel <- lapply(small_panel.preprocessed, make_seurat_sample)

  con <- Conos$new(seurat.panel, n.cores = 1)
  con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)
  con$runClustering(verbose = FALSE)
  truth <- con$clusters$leiden$groups
  expect_gt(nlevels(truth), 1)

  set.seed(1)
  known <- truth[sample(names(truth), round(0.5 * length(truth)))]
  r <- con$propagateLabels(labels = known, verbose = FALSE)
  expect_setequal(names(r$labels), names(truth))                   # a label for every cell
})

test_that("runMarkers works on a Seurat panel (shared sccore::matrixDE core)", {
  skip_if_not_installed("Seurat")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  con <- Conos$new(lapply(small_panel.preprocessed, make_seurat_sample), n.cores = 1)
  con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)
  con$runClustering(verbose = FALSE)

  de <- con$runMarkers(z.threshold = 1, verbose = FALSE)     # low threshold: tiny bundled panel
  expect_type(de, "list")
  nonempty <- Filter(function(d) is.data.frame(d) && nrow(d) > 0, de)
  expect_gt(length(nonempty), 0)                             # markers found on Seurat samples
  expect_true(all(c("Gene", "M", "Z", "AUC", "Specificity") %in% colnames(nonempty[[1]])))
  expect_s3_class(con$plotMarkerDotPlot(n.genes.per.group = 3, z.threshold = 0.5, min.auc = 0.55), "ggplot")
})

test_that("a mixed panel (pagoda2 + Seurat) integrates and finds markers", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("pagoda2")
  data("small_panel.preprocessed", package = "conos", envir = environment())
  ## a current pagoda2 object (its own disk-backed marker path) alongside a Seurat object (sccore::matrixDE)
  raw1 <- getRawCountMatrix(small_panel.preprocessed[[1]], transposed = FALSE)
  rownames(raw1) <- make.unique(rownames(raw1)); colnames(raw1) <- make.unique(colnames(raw1))
  ## variance + PCA directly (the tiny bundled panel is too shallow to survive pagoda2's default QC/gene
  ## filter, which `run()` would apply; real-sized data passes it). Conos only needs the PCA reduction.
  p2 <- pagoda2::Pagoda2$from(raw1, verbose = FALSE)
  p2$runVariance(verbose = FALSE); p2$runReduction(method = "pca", verbose = FALSE)
  seu <- make_seurat_sample(small_panel.preprocessed[[2]])

  con <- Conos$new(list(bm = p2, cb = seu), n.cores = 1)
  expect_setequal(vapply(con$samples, function(s) class(s)[1], character(1)), c("Pagoda2", "Seurat"))
  con$runGraph(k = 15, k.self = 5, space = "PCA", ncomps = 20, n.odgenes = 1000, verbose = FALSE)
  con$runClustering(verbose = FALSE)
  expect_gt(nlevels(con$clusters$leiden$groups), 1)

  de <- con$runMarkers(z.threshold = 1, verbose = FALSE)
  expect_gt(length(Filter(function(d) is.data.frame(d) && nrow(d) > 0, de)), 0)
})
