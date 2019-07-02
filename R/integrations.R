extendMatrix <- function(mtx, col.names) {
  new.names <- setdiff(col.names, colnames(mtx))
  ext.mtx <- matrix(0, nrow=nrow(mtx), ncol=length(new.names))
  colnames(ext.mtx) <- new.names
  return(cbind(mtx, ext.mtx)[,col.names])
}

mergeCountMatrices <- function(cms, transposed=F) {
  if (!transposed) {
    cms %<>% lapply(Matrix::t)
  }

  gene.union <- lapply(cms, colnames) %>% Reduce(union, .)

  res <- lapply(cms, extendMatrix, gene.union) %>% Reduce(rbind, .)
  if (!transposed) {
    res %<>% Matrix::t()
  }
  return(res)
}

checkSeuratV3 <- function() {
  if (!requireNamespace('Seurat', quietly = TRUE) || packageVersion('Seurat') < package_version(x = '3.0.0')) {
    stop("Use of Seurat v3-backed Conos objects requires Seurat v3.X installed")
  }
}

seuratProcV2 <- function(count.matrix, vars.to.regress=NULL, verbose=TRUE, do.par=TRUE, n.pcs=100, cluster=TRUE, tsne=TRUE, umap=FALSE) {
  if (verbose) {
    message("Running Seurat v2 workflow")
  }
  rownames(count.matrix) <- make.unique(rownames(count.matrix))
  max.n.pcs <- min(nrow(count.matrix) - 1, ncol(count.matrix) - 1, n.pcs)
  so <- Seurat::CreateSeuratObject(count.matrix, display.progress=verbose) %>%
    Seurat::NormalizeData(display.progress=verbose) %>%
    Seurat::ScaleData(vars.to.regress=vars.to.regress, display.progress=verbose, do.par=do.par) %>%
    Seurat::FindVariableGenes(do.plot = FALSE, display.progress=verbose) %>%
    Seurat::RunPCA(pcs.compute=max.n.pcs, do.print=FALSE)
  if (cluster) {
    so <- Seurat::FindClusters(so, n.iter=500, n.start=10, dims.use=1:n.pcs, print.output = F)
  }
  if (tsne) {
    so <- Seurat::RunTSNE(so, dims.use=1:n.pcs)
  }
  if (umap) {
    if (packageVersion('Seurat') < package_version(x = '2.3.1')) {
      warning("UMAP support in Seurat came in v2.3.1, please update to a newer version of Seurt to enable UMAP functionality", immediate. = TRUE)
    } else {
      so <- Seurat::RunUMAP(object = so, dims.use = 1:n.pcs)
    }
  }
  return(so)
}

seuratProcV3 <- function(count.matrix, vars.to.regress=NULL, verbose=TRUE, n.pcs=100, cluster=TRUE, tsne=TRUE, umap=FALSE, ...) {
  if (verbose) {
    message("Running Seurat v3 workflow")
  }
  so <- Seurat::CreateSeuratObject(counts = count.matrix)
  so <- Seurat::NormalizeData(object = so, verbose = verbose)
  so <- Seurat::FindVariableFeatures(object = so, verbose = verbose)
  so <- Seurat::ScaleData(object = so, vars.to.regress = vars.to.regress, verbose = verbose)
  so <- Seurat::RunPCA(object = so, npcs = n.pcs, verbose = verbose)
  if (cluster) {
    so <- Seurat::FindNeighbors(object = so, dims = 1:n.pcs, verbose = verbose)
    so <- Seurat::FindClusters(object = so, n.iter = 500, n.start = 10, verbose = verbose)
  }
  if (tsne) {
    so <- Seurat::RunTSNE(object = so, dims = 1:n.pcs)
  }
  if (umap) {
    so <- Seurat::RunUMAP(object = so, dims = 1:n.pcs, verbose = verbose)
  }
  return(so)
}

#' Save Conos object on disk to read it from ScanPy
#'
#' @param con conos object
#' @param output.path path to a folder, where intermediate files will be saved
#' @param metadata.df data.frame with additional metadata with rownames corresponding to cell ids, which should be passed to ScanPy.
#' If NULL, only information about cell ids and origin dataset will be saved.
#' @param n.dims to emulate PCA, Conos embed the graph to a space with `n.dims` dimensions and saves it as PCA
#' @param verbose print more messages
#'
#' @export
saveConosForScanPy <- function(con, output.path, metadata.df=NULL, n.dims=100, verbose=F) {
  if (!dir.exists(output.path))
    stop("Path", output.path, "doesn't exist")

  if (verbose) cat("Merge count matrices... ")
  count.matrix.merged <- lapply(con$samples, getRawCountMatrix) %>% mergeCountMatrices()
  if (verbose) cat("Done.\n")

  cell.ids <- colnames(count.matrix.merged)

  if (!is.null(metadata.df)) {
    metadata.df %<>% .[cell.ids, , drop=F] %>% dplyr::mutate(CellId=cell.ids)
  } else {
    metadata.df <- tibble::tibble(CellId=cell.ids)
  }

  metadata.df$Dataset <- as.character(con$getDatasetPerCell()[cell.ids])

  if (verbose) cat("Create psudo-PCA space\n")
  gene.df <- data.frame(gene=rownames(count.matrix.merged))
  pseudopca.df <- con$embedGraph(target.dims=n.dims, method="largeVis", verbose=verbose)[cell.ids, ]

  if (verbose) cat("Done\n")

  graph.conn <- igraph::as_adjacency_matrix(con$graph, attr="weight")[cell.ids, cell.ids]
  graph.dist <- graph.conn
  graph.dist@x <- 1 - graph.dist@x

  if (verbose) cat("Write data to disk... ")
  Matrix::writeMM(count.matrix.merged, paste0(output.path, "/count_matrix.mtx"))
  data.table::fwrite(metadata.df, paste0(output.path, "/metadata.csv"))
  data.table::fwrite(gene.df, paste0(output.path, "/genes.csv"))
  data.table::fwrite(pseudopca.df, paste0(output.path, "/pca.csv"))
  Matrix::writeMM(graph.conn, paste0(output.path, "/graph_connectivities.mtx"))
  Matrix::writeMM(graph.dist, paste0(output.path, "/graph_distances.mtx"))
  if (verbose) cat("Done\n")
}

#' Create and preprocess a Seurat object
#'
#' @description Create Seurat object from gene count matrix
#'
#' @param count.matrix gene count matrix
#' @param vars.to.regress variables to regress with Seurat
#' @param verbose verbose mode
#' @param do.par use parallel processing for regressing out variables faster, for
#' Seurat v3, use \code{\link[future]{plan}} instead
#' @param n.pcs number of principal components
#' @param cluster do clustering
#' @param tsne do tSNE embedding
#' @param umap do UMAP embedding, works only for Seurat v2.3.1 or higher
#'
#' @return Seurat object
#'
#'
#' @export
#'
basicSeuratProc <- function(count.matrix, vars.to.regress=NULL, verbose=TRUE, do.par=TRUE, n.pcs=100, cluster=TRUE, tsne=TRUE, umap=FALSE) {
  if (!requireNamespace("Seurat")) {
    stop("You need to install 'Seurat' package to be able to use this function")
  }
  proc.fxn <- ifelse(
    test = packageVersion('Seurat') < package_version(x = '3.0.0'),
    yes = seuratProcV2,
    no = seuratProcV3
  )
  so <- proc.fxn(
    count.matrix = count.matrix,
    vars.to.regress = vars.to.regress,
    verbose = verbose,
    do.par = do.par,
    n.pcs = n.pcs,
    cluster = cluster,
    tsne = tsne,
    umap = umap
  )
}
