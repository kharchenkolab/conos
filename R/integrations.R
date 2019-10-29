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
#' @param norm logical, whether to include the matrix of normalised counts (only raw counts saved by default)
#' @param pseudo.pca logical, to produce an emulated PCA, Conos embeds the graph to a space with `n.dims` dimensions and saves it as a pseudoPCA
#' @param pca logical, whether to include PCA of all the samples (not batch corrected)
#' @param n.dims number of dimensions for calculating PCA and/or pseudoPCA
#' @param embed logical, whether to include the current conos embedding
#' @param connect logical, whether to include graph connectivities and distances
#' @param verbose print more messages
#'
#' @export
saveConosForScanPy <- function(con, output.path, metadata.df=NULL, norm=FALSE, pseudo.pca=FALSE, pca=FALSE, n.dims=100, embed=FALSE, connect=FALSE, verbose=FALSE) {
  if (!dir.exists(output.path))
    stop("Path", output.path, "doesn't exist")

  if (verbose) cat("Merge raw count matrices...\t")
  raw.count.matrix.merged <- con$getJointCountMatrix(raw=TRUE)
  if (verbose) cat("Done.\n")

  cell.ids <- colnames(raw.count.matrix.merged)
  gene.df <- data.frame(gene=rownames(raw.count.matrix.merged))

  if (!is.null(metadata.df)) {
    metadata.df %<>% .[cell.ids, , drop=F] %>% dplyr::mutate(CellId=cell.ids)
  } else {
    metadata.df <- tibble::tibble(CellId=cell.ids)
  }
  metadata.df$Dataset <- as.character(con$getDatasetPerCell()[cell.ids])

  if (norm) {
    if (verbose) cat("Merge count matrices...\t\t")
    count.matrix.merged <- con$getJointCountMatrix(raw=FALSE)
    if (verbose) cat("Done.\n")
  }

  # Create a batch-free embedding that can be used instead of PCA space
  if (pseudo.pca) {
    if (verbose) cat("Create psudo-PCA space\t")
    pseudopca.df <- con$embedGraph(target.dims=n.dims, method="largeVis", verbose=verbose)[cell.ids, ] %>% as.data.frame()
    if (verbose) cat("Done\n")
  }
  
  if (pca){
    if (verbose) cat("Save PCA space\t")
    pca.df <- pcaFromConos(con$samples, ncomps=n.dims, n.odgenes=2000, verbose=verbose) %>% as.data.frame()
    if (verbose) cat("Done\n")
  }

  if (embed){
    if (verbose) cat("Save the embedding\t")
    embed.df <- con$embedding[cell.ids,] %>% as.data.frame()
    if (verbose) cat("Done\n")
  }
  
  if (connect){
    if (verbose) cat("Save graph matrices\t")
    graph.conn <- igraph::as_adjacency_matrix(con$graph, attr="weight")[cell.ids, cell.ids]
    graph.dist <- graph.conn
    graph.dist@x <- 1 - graph.dist@x
    if (verbose) cat("Done\n")
  }

  if (verbose) cat("Write data to disk...\t\t")
  Matrix::writeMM(raw.count.matrix.merged, paste0(output.path, "/raw_count_matrix.mtx"))
  data.table::fwrite(metadata.df, paste0(output.path, "/metadata.csv"))
  data.table::fwrite(gene.df, paste0(output.path, "/genes.csv"))
  if (norm) Matrix::writeMM(count.matrix.merged, paste0(output.path, "/count_matrix.mtx"))
  if (pseudo.pca) data.table::fwrite(pseudopca.df, paste0(output.path, "/pseudopca.csv"))
  if (pca) data.table::fwrite(pseudopca.df, paste0(output.path, "/pca.csv"))
  if (embed) data.table::fwrite(embed.df, paste0(output.path, "/embed.csv"))
  if (connect) {
    Matrix::writeMM(graph.conn, paste0(output.path, "/graph_connectivities.mtx"))
    Matrix::writeMM(graph.dist, paste0(output.path, "/graph_distances.mtx"))
  }
  if (verbose) cat("Done.\n")
  if (verbose) cat("All Done!")
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

#' RNA velocity analysis on samples integrated with conos
#'
#' @description Create a list of objects to pass into gene.relative.velocity.estimates function from the velocyto.R package
#'
#' @param cms.list list of velocity files written out as cell.counts.matrices.rds files by running dropest with -V option
#' @param con conos object (after creating an embedding and running leiden clustering)
#' @param n.odgenes number of overdispersed genes to use for PCA
#' @param verbose verbose mode
#'
#'
#' @return List with cell distances, combined spliced expression matrix, combined unspliced expression matrix, combined matrix of spanning reads, cell colors for clusters and embedding (taken from conos)
#'
#'
#'
#' @export
#'
velocityInfoConos <- function(cms.list, con, n.odgenes=2e3, verbose=TRUE, min.max.cluster.average.emat=0.2, min.max.cluster.average.nmat=0.05, min.max.cluster.average.smat=0.01) {
  if (!requireNamespace("velocyto.R")) {
    stop("You need to install 'velocyto.R' package to be able to use this function")
  }

  if (verbose) cat("Merging raw count matrices...\n")
  # Merge samples to get names of relevant cells and genes 
  raw.count.matrix.merged <- con$getJointCountMatrix(raw=TRUE)
  
  if (verbose) cat("Merging velocity files...\n")
  # Intersect genes and cells between the conos object and all the velocity files
  cms.list <- lapply(cms.list, prepareVelocity, rownames(raw.count.matrix.merged), colnames(raw.count.matrix.merged))
  # Keep only genes present in velocity files from all the samples
  common.genes <-  Reduce(intersect, lapply(cms.list, function(x) {rownames(x[[1]])}))
  cms.list <- lapply(cms.list, function(x) {lapply(x, function(y) {y[row.names(y) %in% common.genes,]} )} )
  
  # Merge velocity files from different samples
  emat <- do.call(cbind, lapply(cms.list, function(x) {x[[1]]}))
  nmat <- do.call(cbind, lapply(cms.list, function(x) {x[[2]]}))
  smat <- do.call(cbind, lapply(cms.list, function(x) {x[[3]]}))
  
  cluster.label <- con$clusters$leiden$groups
  cell.colors <- fac2col(cluster.label)
  emb <- con$embedding
  
  # Keep the order of cells consistent between velocity matrices and the embedding (not really sure whether it's necessary...)
  emat <- emat[,order(match(colnames(emat), rownames(emb)))]
  nmat <- nmat[,order(match(colnames(nmat), rownames(emb)))]
  smat <- smat[,order(match(colnames(smat), rownames(emb)))]
  
  if (verbose) cat("Calculating cell distances...\n")
  # Get PCA results for all the samples from the conos object
  pcs <- pcaFromConos(con$samples, n.odgenes=n.odgenes)
  # Again, keep the order of cells consistent
  pcs <- pcs[order(match(rownames(pcs), rownames(emb))),]
  # Calculate the cell distances based on correlation
  cell.dist <- as.dist(1 - velocyto.R::armaCor(t(pcs)))
  
  if (verbose) cat("Filtering velocity...\n")
  emat %<>% velocyto.R::filter.genes.by.cluster.expression(cluster.label, min.max.cluster.average=min.max.cluster.average.emat)
  nmat %<>% velocyto.R::filter.genes.by.cluster.expression(cluster.label, min.max.cluster.average=min.max.cluster.average.nmat)
  smat %<>% velocyto.R::filter.genes.by.cluster.expression(cluster.label, min.max.cluster.average=min.max.cluster.average.smat)
  
  if (verbose) cat("All Done!")
  return(list(cell.dist=cell.dist, emat=emat, nmat=nmat, smat=smat, cell.colors=cell.colors, emb=emb))
}

# Intersect genes and cells between all the velocity files and the conos object
prepareVelocity <- function(cms.file, genes, cells) {
  exon.genes <- rownames(cms.file$exon)
  intron.genes <- rownames(cms.file$intron)
  spanning.genes <- rownames(cms.file$spanning)
  # Only common genes between the 3 files and conos object
  common.genes <- intersect(exon.genes, intron.genes) %>% intersect(spanning.genes) %>% intersect(genes)
  
  exon.cells <- colnames(cms.file$exon)
  intron.cells <- colnames(cms.file$intron)
  spanning.cells <- colnames(cms.file$spanning)
  # Only common cells between the 3 files and conos object
  common.cells <- intersect(exon.cells, intron.cells) %>% intersect(spanning.cells) %>% intersect(cells)
  
  cms.file$exon <- cms.file$exon[common.genes,common.cells]
  cms.file$intron <- cms.file$intron[common.genes,common.cells]
  cms.file$spanning <- cms.file$spanning[common.genes,common.cells]
  return(cms.file)
}

# Get PCA results for all the samples from the conos object
# This is a modification of the quickPlainPCA function
pcaFromConos <- function(p2.list, data.type='counts', k=30, ncomps=100, n.odgenes=NULL, verbose=TRUE) {
  od.genes <- commonOverdispersedGenes(p2.list, n.odgenes, verbose = FALSE)
  if(length(od.genes)<5) return(NULL)
  
  if(verbose) cat('Calculating PCs for',length(p2.list),' datasets ...')
  
  # Get scaled matrices from a list of pagoda2 objects
  sm <- scaledMatrices(p2.list, data.type=data.type, od.genes=od.genes, var.scale=TRUE, neighborhood.average=FALSE);
  # Transpose the scaled matrices since we want to run PCA on cells and not genes (like in quickPlainPCA)
  sm <- lapply(sm, t)
  # Get the names of all the cells
  nms <- Reduce(union, lapply(sm, colnames))
  
  pcs <- lapply(sm, function(x) {
    cm <- Matrix::colMeans(x);
    ncomps <- min(c(nrow(cm)-1,ncol(cm)-1,ncomps))
    res <- irlba::irlba(x, nv=ncomps, nu=0, center=cm, right_only=F, reorth=T)
    res
  })
  
  pcj <- do.call(rbind,lapply(pcs,function(x) x$v))
  rownames(pcj) <- nms
  res <- pcj
  
  return(res)
}