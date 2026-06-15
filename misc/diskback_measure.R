## conos §5 — evidence that pagoda2.1 disk-backing lowers a sample's resident RAM and streams through
## the exact accessors conos reads (getExpressionBlock / viewColMeanVar), with bit-exact values.
## Run with the dev pagoda2 + lstar (.Rlib) on a real 10x sample. Tiny-data variant: swap in
## conos::small_panel.preprocessed[[1]]$misc$rawCounts for `cxg`.
suppressMessages({
  .libPaths(c("/home/pkharchenko/p21/lstar/.Rlib", .libPaths()))
  library(Matrix)
  pkgload::load_all("/home/pkharchenko/p21/pagoda2", quiet = TRUE)
  library(lstar)
})

d   <- "/home/pkharchenko/p21/data"
pre <- "GSM5746259_MGI0369_1_SLAB-145-0"
m   <- Matrix::readMM(file.path(d, paste0(pre, ".matrix.mtx.gz")))         # genes x cells
rownames(m) <- make.unique(read.delim(gzfile(file.path(d, paste0(pre, ".features.tsv.gz"))),
                                      header = FALSE)$V2)
colnames(m) <- readLines(gzfile(file.path(d, paste0(pre, ".barcodes.tsv.gz"))))
cxg <- as(Matrix::t(m), "CsparseMatrix")                                    # cells x genes
cxg <- cxg[Matrix::rowSums(cxg) >= 500, ]
cxg <- as(cxg[, Matrix::colSums(cxg > 0) >= 10], "CsparseMatrix")           # light, realistic filter
cat(sprintf("real sample %s: %d cells x %d genes; counts = %.1f MB\n",
            pre, nrow(cxg), ncol(cxg), as.numeric(object.size(cxg)) / 1e6))

mk <- function(backend, dir = NULL) {
  p <- Pagoda2$new(Matrix::t(cxg), modelType = "plain", n.cores = 1, verbose = FALSE)
  p$addFacet("RNA2", cxg, modelType = "plain", featureType = "gene",
             backend = backend, backend.dir = dir)
  p
}
p.mem  <- mk("memory")
store  <- tempfile("rnadisk_", fileext = ".lstar.zarr")
p.disk <- mk("lstar", store)

sz.mem  <- as.numeric(object.size(p.mem$misc$facetStore[["RNA2"]]$rawCounts))
sz.disk <- as.numeric(object.size(p.disk$misc$facetStore[["RNA2"]]$rawCounts))   # NULL -> ~0
store.b <- sum(file.info(list.files(store, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE)
cat(sprintf("facet rawCounts RESIDENT:  memory = %.1f MB  |  lstar = %.4f MB (NULL)  |  on-disk = %.1f MB\n",
            sz.mem / 1e6, sz.disk / 1e6, store.b / 1e6))

## --- streaming equivalence on the conos read path (getExpressionBlock) ---
g  <- head(p.disk$misc$facetStore[["RNA2"]]$featureNames, 50)
xm <- as.matrix(p.mem$getExpressionBlock(genes = g, facet = "RNA2", orientation = "cell_by_gene"))
xd <- as.matrix(p.disk$getExpressionBlock(genes = g, facet = "RNA2", orientation = "cell_by_gene"))
xm <- xm[rownames(xd), colnames(xd), drop = FALSE]
cat(sprintf("getExpressionBlock (conos read path): max|mem-disk| = %.3e  (dim %s)\n",
            max(abs(xm - xd)), paste(dim(xd), collapse = "x")))

## --- per-feature mean/var (compared positionally; in-memory facet labels rows by index) ---
mvm <- p.mem$viewColMeanVar(facet = "RNA2"); mvd <- p.disk$viewColMeanVar(facet = "RNA2")
n <- min(nrow(mvm), nrow(mvd))
cat(sprintf("viewColMeanVar (streamed): max|dmean|=%.3e max|dvar|=%.3e\n",
            max(abs(mvm$m[1:n] - mvd$m[1:n])), max(abs(mvm$v[1:n] - mvd$v[1:n]))))
