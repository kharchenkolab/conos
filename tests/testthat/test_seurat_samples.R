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
