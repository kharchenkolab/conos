
#' @keywords internal
extendMatrix <- function(mtx, col.names) {
  new.names <- setdiff(col.names, colnames(mtx))
  ext.mtx <- matrix(0, nrow=nrow(mtx), ncol=length(new.names))
  colnames(ext.mtx) <- new.names
  return(cbind(mtx, ext.mtx)[,col.names])
}

#' @keywords internal
mergeCountMatrices <- function(cms, transposed=FALSE) {
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

#' @keywords internal
checkSeuratV3 <- function() {
  if (!requireNamespace('Seurat', quietly = TRUE) || packageVersion('Seurat') < package_version(x = '3.0.0')) {
    stop("Use of Seurat v3-backed Conos objects requires Seurat v3.X installed")
  }
}

#' @keywords internal
seuratProcV2 <- function(count.matrix, vars.to.regress=NULL, verbose=TRUE, do.par=TRUE, n.pcs=100, cluster=TRUE, tsne=TRUE, umap=FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it, as described here: <https://satijalab.org/seurat/articles/install.html>.", call. = FALSE)
  }
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

#' @keywords internal
seuratProcV3 <- function(count.matrix, vars.to.regress=NULL, verbose=TRUE, n.pcs=100, cluster=TRUE, tsne=TRUE, umap=FALSE, ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it, as described here: <https://satijalab.org/seurat/articles/install.html>.", call. = FALSE)
  }
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
#' @param hdf5_filename name of HDF5 written with ScanPy files. Note: the \pkg{\link{rhdf5}} package is required
#' @param metadata.df data.frame with additional metadata with rownames corresponding to cell ids, which should be passed to ScanPy (default=NULL)
#'     If NULL, only information about cell ids and origin dataset will be saved.
#' @param cm.norm boolean Whether to include the matrix of normalised counts (default=FALSE).
#' @param embedding boolean Whether to include the current conos embedding (default=TRUE).
#' @param pseudo.pca boolean Whether to produce an emulated PCA by embedding the graph to a space with `n.dims` dimensions and save it as a pseudoPCA (default=FALSE).
#' @param pca boolean Whether to include PCA of all the samples (not batch corrected) (default=FALSE).
#' @param n.dims numeric Number of dimensions for calculating PCA and/or pseudoPCA (default=100).
#' @param alignment.graph boolean Whether to include graph of connectivities and distances (default=TRUE).
#' @param verbose boolean Whether to use verbose mode (default=FALSE)
#' @seealso The rhdf5 package documentation here: <https://www.bioconductor.org/packages/release/bioc/html/rhdf5.html>
#' @return AnnData object for ScanPy, saved to disk
#' @export
saveConosForScanPy <- function(con, output.path, hdf5_filename, metadata.df=NULL, cm.norm=FALSE, pseudo.pca=FALSE, 
  pca=FALSE, n.dims=100, embedding=TRUE, alignment.graph=TRUE, verbose=FALSE) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) {
    stop("The package rhdf5 is required for saveConosForScanPy(). Please install.")
  }

  if (!dir.exists(output.path)){
    stop("Path", output.path, "doesn't exist")
  }

  if (tools::file_ext(hdf5_filename) != "h5"){
    stop("File", hdf5_filename, "must have the file extension *.h5")
  }

  if (verbose) message("Merge raw count matrices...\t")
  raw.count.matrix.merged <- con$getJointCountMatrix(raw=TRUE)
  if (verbose) message("Done.\n")

  cell.ids <- rownames(raw.count.matrix.merged)
  gene.df <- data.frame(gene=colnames(raw.count.matrix.merged))

  if (!is.null(metadata.df)) {
    metadata.df %<>% .[cell.ids, , drop=FALSE] %>% dplyr::mutate(CellId=cell.ids)
  } else {
    metadata.df <- tibble::tibble(CellId=cell.ids)
  }
  metadata.df$Dataset <- as.character(con$getDatasetPerCell()[cell.ids])

  if (cm.norm) {
    if (verbose) message("Merge count matrices...\t\t")
    count.matrix.merged <- con$getJointCountMatrix(raw=FALSE)
    if (verbose) message("Done.\n")
  }

  if (embedding){
    if (verbose) message("Save the embedding...\t\t")
    if (length(con$embedding)>1) {
      embedding.df <- con$embedding[cell.ids,] %>% as.data.frame()
      if (verbose) message("Done.\n")
    } else {
      warning("\n No embedding found in the conos object. Skipping... \n")
      embedding <- FALSE
    }

  }

  # Create a batch-free embedding that can be used instead of PCA space
  if (pseudo.pca) {
    if (verbose) message("Create psudo-PCA space...\t")
    pseudopca.df <- con$embedGraph(target.dims=n.dims, method="largeVis", verbose=FALSE)[cell.ids, ] %>% as.data.frame()
    if (verbose) message("Done.\n")
  }

  if (pca){
    if (verbose) message("Save PCA space...\t\t")
    pca.df <- pcaFromConos(con$samples, ncomps=n.dims, n.odgenes=2000, verbose=FALSE) %>% as.data.frame()
    if (verbose) message("Done.\n")
  }

  if (alignment.graph){
    if (verbose) message("Save graph matrices...\t\t")
    if (!is.null(con$graph)) {
      graph.conn <- igraph::as_adjacency_matrix(con$graph, attr="weight")[cell.ids, cell.ids]
      graph.dist <- graph.conn
      graph.dist@x <- 1 - graph.dist@x
      if (verbose) message("Done.\n")
    } else {
      warning("\n No graph found in the conos object. Skipping... \n")
      alignment.graph <- FALSE
    }
  }

  if (verbose) message("Write data to disk...\t\t")
  ## create HDF5 file
  total_hdf5file_path = paste0(output.path, "/", hdf5_filename)
  rhdf5::h5createFile(total_hdf5file_path)
  ## raw.count.matrix.merged
  rhdf5::h5createGroup(total_hdf5file_path, "raw_count_matrix")
  rhdf5::h5write(raw.count.matrix.merged@x, total_hdf5file_path, "raw_count_matrix/data")
  rhdf5::h5write(dim(raw.count.matrix.merged), total_hdf5file_path, "raw_count_matrix/shape")
  rhdf5::h5write(raw.count.matrix.merged@i, total_hdf5file_path, "raw_count_matrix/indices")
  rhdf5::h5write(raw.count.matrix.merged@p, total_hdf5file_path, "raw_count_matrix/indptr")
  ## metadata
  rhdf5::h5createGroup(total_hdf5file_path, "metadata")
  rhdf5::h5write(metadata.df, total_hdf5file_path, "metadata/metadata.df")
  ## genes
  rhdf5::h5createGroup(total_hdf5file_path, "genes")
  rhdf5::h5write(gene.df, total_hdf5file_path, "genes/genes.df")
  ## count_matrix
  if (cm.norm) {
    rhdf5::h5createGroup(total_hdf5file_path, "count_matrix")
    rhdf5::h5write(count.matrix.merged@x, total_hdf5file_path, "count_matrix/data")
    rhdf5::h5write(dim(count.matrix.merged), total_hdf5file_path, "count_matrix/shape")
    rhdf5::h5write(count.matrix.merged@i, total_hdf5file_path, "count_matrix/indices")
    rhdf5::h5write(count.matrix.merged@p, total_hdf5file_path, "count_matrix/indptr")
  }
  if (embedding) {
    rhdf5::h5createGroup(total_hdf5file_path, "embedding")
    rhdf5::h5write(embedding.df, total_hdf5file_path, "embedding/embedding.df")
  }
  if (pseudo.pca) {
    rhdf5::h5createGroup(total_hdf5file_path, "pseudopca")
    rhdf5::h5write(pseudopca.df, total_hdf5file_path, "pseudopca/pseudopca.df")
  }
  if (pca) {
    rhdf5::h5createGroup(total_hdf5file_path, "pca")
    rhdf5::h5write(pca.df, total_hdf5file_path, "pca/pca.df")
  }
  if (alignment.graph) {
    ## graph_connectivities
    rhdf5::h5createGroup(total_hdf5file_path, "graph_connectivities")
    rhdf5::h5write(graph.conn@x, total_hdf5file_path, "graph_connectivities/data")  
    rhdf5::h5write(dim(graph.conn), total_hdf5file_path, "graph_connectivities/shape") 
    rhdf5::h5write(graph.conn@i, total_hdf5file_path, "graph_connectivities/indices") 
    rhdf5::h5write(graph.conn@p, total_hdf5file_path, "graph_connectivities/indptr") 
    ## graph_distances
    rhdf5::h5createGroup(total_hdf5file_path, "graph_distances")
    rhdf5::h5write(graph.dist@x, total_hdf5file_path, "graph_distances/data")  
    rhdf5::h5write(dim(graph.dist), total_hdf5file_path, "graph_distances/shape") 
    rhdf5::h5write(graph.dist@i, total_hdf5file_path, "graph_distances/indices") 
    rhdf5::h5write(graph.dist@p, total_hdf5file_path, "graph_distances/indptr") 
  }
  if (verbose) message("All Done!")
}

#' Create and preprocess a Seurat object
#'
#' @param count.matrix gene count matrix
#' @param vars.to.regress variables to regress with Seurat (default=NULL)
#' @param verbose boolean Verbose mode (default=TRUE)
#' @param do.par boolean Use parallel processing for regressing out variables faster, for
#'     Seurat v3, use \code{\link[future]{plan}} instead (default=TRUE)
#' @param n.pcs numeric Number of principal components (default=100)
#' @param cluster boolean Whether to perform clustering (default=TRUE)
#' @param tsne boolean Whether to construct tSNE embedding (default=TRUE)
#' @param umap boolean Whether to construct UMAP embedding, works only for Seurat v2.3.1 or higher (default=FALSE)
#' @return Seurat object
#' @export
basicSeuratProc <- function(count.matrix, vars.to.regress=NULL, verbose=TRUE, do.par=TRUE, n.pcs=100, cluster=TRUE, tsne=TRUE, umap=FALSE) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" needed for this function to work. Please install it, as described here: <https://satijalab.org/seurat/articles/install.html>.", call. = FALSE)
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
#' Create a list of objects to pass into gene.relative.velocity.estimates function from the velocyto.R package
#'
#' @param cms.list list of velocity files written out as cell.counts.matrices.rds files by running dropest with -V option
#' @param con conos object (after creating an embedding and running leiden clustering)
#' @param clustering name of clustering in the conos object to use (default=NULL). Either 'clustering' or 'groups' must be provided. 
#' @param groups set of clusters to use (default=NULL). Ignored if 'clustering' is not NULL. 
#' @param n.odgenes numeric Number of overdispersed genes to use for PCA (default=2000).
#' @param verbose boolean Whether to use verbose mode (default=TRUE)
#' @return List with cell distances, combined spliced expression matrix, combined unspliced expression matrix, combined matrix of spanning reads, cell colors for clusters and embedding (taken from conos)
#' @export
velocityInfoConos <- function(cms.list, con, clustering=NULL, groups=NULL, n.odgenes=2e3, verbose=TRUE, min.max.cluster.average.emat=0.2, min.max.cluster.average.nmat=0.05, min.max.cluster.average.smat=0.01) {
  if (!requireNamespace("velocyto.R")) {
    stop("You need to install 'velocyto.R' package to be able to use this function")
  }

  groups <- parseCellGroups(con, clustering, groups)
  cell.colors <- fac2col(groups)

  if (!is.null(con$embedding)){
    emb <- con$embedding
  } else {
    stop("No embedding found in the conos object. Run 'con$embedGraph()' before running this function.")
  }

  if (verbose) message("Merging raw count matrices...\n")
  # Merge samples to get names of relevant cells and genes
  raw.count.matrix.merged <- con$getJointCountMatrix(raw=TRUE)

  if (verbose) message("Merging velocity files...\n")
  # Intersect genes and cells between the conos object and all the velocity files
  cms.list <- lapply(cms.list, prepareVelocity, genes=colnames(raw.count.matrix.merged), cells=rownames(raw.count.matrix.merged))
  # Keep only genes present in velocity files from all the samples
  common.genes <-  Reduce(intersect, lapply(cms.list, function(x) {rownames(x[[1]])}))
  cms.list <- lapply(cms.list, function(x) {lapply(x, function(y) {y[row.names(y) %in% common.genes,]} )} )

  # Merge velocity files from different samples
  emat <- do.call(cbind, lapply(cms.list, function(x) {x[[1]]}))
  nmat <- do.call(cbind, lapply(cms.list, function(x) {x[[2]]}))
  smat <- do.call(cbind, lapply(cms.list, function(x) {x[[3]]}))

  # Keep the order of cells consistent between velocity matrices and the embedding (not really sure whether it's necessary...)
  emat <- emat[,order(match(colnames(emat), rownames(emb)))]
  nmat <- nmat[,order(match(colnames(nmat), rownames(emb)))]
  smat <- smat[,order(match(colnames(smat), rownames(emb)))]

  if (verbose) message("Calculating cell distances...\n")
  # Get PCA results for all the samples from the conos object
  pcs <- pcaFromConos(con$samples, n.odgenes=n.odgenes)
  # Again, keep the order of cells consistent
  pcs <- pcs[order(match(rownames(pcs), rownames(emb))),]
  # Calculate the cell distances based on correlation
  cell.dist <- as.dist(1 - velocyto.R::armaCor(t(pcs)))

  if (verbose) message("Filtering velocity...\n")
  emat %<>% velocyto.R::filter.genes.by.cluster.expression(groups, min.max.cluster.average=min.max.cluster.average.emat)
  nmat %<>% velocyto.R::filter.genes.by.cluster.expression(groups, min.max.cluster.average=min.max.cluster.average.nmat)
  smat %<>% velocyto.R::filter.genes.by.cluster.expression(groups, min.max.cluster.average=min.max.cluster.average.smat)

  if (verbose) message("All Done!")
  return(list(cell.dist=cell.dist, emat=emat, nmat=nmat, smat=smat, cell.colors=cell.colors, emb=emb))
}

# Intersect genes and cells between all the velocity files and the conos object
#' @keywords internal
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
#' @keywords internal
pcaFromConos <- function(p2.list, data.type='counts', k=30, ncomps=100, n.odgenes=NULL, verbose=TRUE) {

  od.genes <- commonOverdispersedGenes(p2.list, n.odgenes, verbose = FALSE)
  if (length(od.genes)<5){
    return(NULL)
  }

  if(verbose) message('Calculating PCs for',length(p2.list),' datasets...\n')

  # Get scaled matrices from a list of pagoda2 objects
  sm <- scaledMatrices(p2.list, data.type=data.type, od.genes=od.genes, var.scale=TRUE)
  # Transpose the scaled matrices since we want to run PCA on cells and not genes (like in quickPlainPCA)
  sm <- lapply(sm, t)
  # Get the names of all the cells
  nms <- Reduce(union, lapply(sm, colnames))

  pcs <- lapply(sm, function(x) {
    cm <- Matrix::colMeans(x);
    ncomps <- min(c(nrow(cm)-1,ncol(cm)-1,ncomps))
    res <- irlba::irlba(x, nv=ncomps, nu=0, center=cm, right_only=FALSE, reorth=TRUE)
    res
  })

  pcj <- do.call(rbind,lapply(pcs,function(x) x$v))
  rownames(pcj) <- nms
  res <- pcj

  return(res)
}

#' Convert Conos object to Pagoda2 object
#'
#' @param con Conos object
#' @param n.pcs numeric Number of principal components (default=100)
#' @param n.odgenes numeric Number of overdispersed genes (default=2000)
#' @param verbose boolean Whether to give verbose output (default=TRUE)
#' @param ... parameters passed to Pagoda2$new()
#' @return pagoda2 object 
#' @export
convertToPagoda2 <- function(con, n.pcs=100, n.odgenes=2000, verbose=TRUE, ...) {
  if (!requireNamespace('pagoda2', quietly=TRUE)) {
    stop("'pagoda2' must be installed to convert a Conos object to a Pagoda2 object. Please refer to <https://github.com/kharchenkolab/pagoda2>.")
  }

  p2 <- con$getJointCountMatrix(raw=TRUE) %>% Matrix::t() %>%
    Pagoda2$new(n.cores=con$n.cores, verbose=verbose, ...)

  if (n.pcs > 0) {
    if (verbose) message("Estimating PCA... ")
    p2$reductions$PCA <- pcaFromConos(con$samples, ncomps=n.pcs, n.odgenes=n.odgenes, verbose=FALSE)
    if (verbose) message("Done\n")
  }

  p2$graphs$conos <- con$graph

  if (!("uninitializedField" %in% class(con$embedding))) {
    p2$embeddings$PCA$conos <- con$embedding
  }

  if (length(con$clusters) > 0) {
    con.clusters <- con$clusters %>% setNames(paste0("conos_", names(.))) %>%
      lapply(`[[`, "groups") %>% lapply(as.factor)
    p2$clusters %<>% c(con.clusters)
  }

  p2$clusters$dataset <- con$getDatasetPerCell()

  return(p2)
}



#' Utility function to generate a pagoda2 app from a conos object
#' 
#' @param conos Conos object
#' @param cdl list Optional list of raw matrices (so that gene merging doesn't have to be redone) (default=NULL)
#' @param metadata list Optional list of (named) metadata factors (default=NULL)
#' @param filename string Name of the *.bin file to seralize for the pagoda2 application if save=TRUE (default='conos_app.bin')
#' @param save boolean Save serialized *bin file specified in filename (default=TRUE)
#' @param n.cores integer Number of cores (default=2)
#' @param go.env GO environment for the organism of interest (default=NULL)
#' @param cell.subset string Cells to subset with the conos embedding conos$embedding. If NULL, uses all cells via rownames(conos$embedding) (default-NULL)
#' @param max.cells numeric Limit to the cells that are included in the conos. If Inf, there is no limit (default=Inf)
#' @param additional.embeddings list Additional embeddings to add to conos for the pagoda2 app (default=NULL)
#' @param test.pathway.overdispersion boolean Find all IDs using GO category against either org.Hs.eg.db ('hs') or org.Mm.eg.db ('mm') (default=FALSE
#' @param organism string Organism of interest, either 'hs' (Homo sapiens) or 'mm' (Mus musculus, i.e. mouse) (default=NULL). Only used if test.pathway.overdispersion is TRUE.
#' @param return.details boolean If TRUE, return list of p2 application, pagoda2 object, list of raw matrices, and cell names. If FALSE, simply return pagoda2 app object. (default=FALSE)
#' @return pagoda2 app object
#' @export 
p2app4conos <- function(conos, cdl=NULL, metadata=NULL, filename='conos_app.bin', 
  save=TRUE, n.cores=2, go.env=NULL, cell.subset=NULL, max.cells=Inf, additional.embeddings=NULL, 
  test.pathway.overdispersion=FALSE, organism=NULL, return.details=FALSE) {
  
  if (!requireNamespace("pagoda2", quietly = TRUE)) {
    stop("You have to install the pagoda2 package to use p2app4conos(). Please refer to <https://github.com/kharchenkolab/pagoda2>.")
  }

  if(is.null(cdl)) {
    #cdl <- lapply(conos$samples,function(p) t(p$misc$rawCounts))
    cdl <- lapply(rawMatricesWithCommonGenes(conos),t)
  }
  samf <- lapply(conos$samples,function(x) rownames(x$counts))
  samf <- as.factor(setNames(rep(names(samf),unlist(lapply(samf,length))),unlist(samf)))

  if(!is.null(cell.subset)) {
    cell.subset <- intersect(cell.subset,rownames(conos$embedding))
  } else {
    cell.subset <- rownames(conos$embedding)  
  }
  
  samf <- droplevels(samf[names(samf) %in% cell.subset])
  cdl <- lapply(cdl,function(x) x[,colnames(x) %in% cell.subset,drop=F])
  
  # limit to the cells that are included in the conos
  vc <- unlist(lapply(cdl,colnames));
  if (length(vc)>max.cells) { # subsample
    cat("subsampling",length(vc),"cells down to",max.cells,'...');
    vc <- sample(vc,max.cells)
    cat('done\n');
  } 
  cdl <- lapply(cdl,function(d) d[,colnames(d) %in% intersect(vc,names(samf))])
  cm <- do.call(cbind,cdl);
  
   
  cp2 <- pagoda2::basicP2proc(cm, min.cells.per.gene=1, nPcs=50, get.tsne=TRUE, get.largevis=FALSE, make.geneknn=TRUE, n.cores=n.cores)
  
  if (test.pathway.overdispersion) {
    if (organism =='mm') { 
      suppressMessages(library(org.Mm.eg.db))
      # translate gene names to ids
      ids <- unlist(lapply(mget(colnames(cp2$counts),org.Mm.egALIAS2EG,ifnotfound=NA),function(x) x[1]))
      # reverse map
      rids <- names(ids); names(rids) <- ids;
      # list all the ids per GO category
      go.env <- list2env(eapply(org.Mm.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))
    } else if (organism =='hs') {
      suppressMessages(library(org.Hs.eg.db))
      ids <- unlist(lapply(mget(colnames(cp2$counts),org.Hs.egALIAS2EG,ifnotfound=NA),function(x) x[1]))
      rids <- names(ids); names(rids) <- ids;
      # list all the ids per GO category
      go.env <- list2env(eapply(org.Hs.egGO2ALLEGS,function(x) as.character(na.omit(rids[x]))))
    } else { 
      stop("Unknown organism")
    }
  }
    
  if(!is.null(go.env)) {
    #cp2$getHierarchicalDiffExpressionAspects(type='PCA',clusterName='community',z.threshold=3)
    # here we test for pathway overdispersion, but recursive diff expression, as shown above will work just as well
    cp2$testPathwayOverdispersion(go.env,verbose=TRUE,correlation.distance.threshold=0.95,recalculate.pca=F,top.aspects=15)
    
    library(GO.db)
    termDescriptions <- Term(GOTERM[names(go.env)]); # saves a good minute or so compared to individual lookups
    geneSets <- lapply(sn(names(go.env)),function(x) {
      list(properties=list(locked=T,genesetname=x,shortdescription=as.character(termDescriptions[x])),genes=c(go.env[[x]]))
    })
  } else {
    hdea <- cp2$getHierarchicalDiffExpressionAspects(type='PCA',z.threshold=3)
    geneSets <- hierDiffToGenesets(hdea);
  }
  

  # add all kinds of embeddings
  # various joint embeddding versions
  #cp2$embeddings$Conos <- list("All"=t(conO$embedding),"cochlea"=t(conC$embedding),"DRG"=t(conD$embedding),"Merged"=t(conO$embedding),"Refined"=t(con2$embedding))
  cn <- rownames(cp2$counts); # cell names
  cp2$embeddings$Conos <- list("All"=conos$embedding[cn,])
  # hack: add placeholder spaces, so that old p2 version picks up the new embeddings
  cp2$reductions$Conos <- list(); 

  cp2$embeddings$sample <- lapply(conos$samples,function(x) { em <- x$embeddings$PCA[[1]]; em[rownames(em) %in% cn,] }) # embeddings of the individual samples
  cp2$embeddings$sample <- cp2$embeddings$sample[!unlist(lapply(cp2$embeddings$sample,is.null))]
  if(length(cp2$embeddings$sample)>1) {
    cp2$reductions$sample <- list();
  } else {
    cp2$embeddings$sample <- NULL;
    
  }
  
  if(!is.null(additional.embeddings)) {
    cp2$embeddings$Other <- lapply(additional.embeddings, function(em) { em[rownames(em) %in% cn,] })
    cp2$reductions$Other <- list();
  }
  
  
  # additional metadata with different factors .. you probably want to include something like sample or tissue (or patient type)
  metadata <- c(metadata,list(sample=samf));
  metadata <- lapply(metadata, function(d) d[cn])
  meta <- lapply(sn(names(metadata)),function(n) p2.metadata.from.factor(droplevels(as.factor(metadata[[n]])),displayname=n))
  
  p2app <- make.p2.app(cp2, dendrogramCellGroups = as.factor(conos$clusters[[1]]$groups[cn]), additionalMetadata = meta, geneSets = geneSets,innerOrder='odPCA');
  
  # Optional showing of app
  #show.app(p2app, name='newPagoda',browse=F)
  # Save serialised web object, RDS app and session image
  if(save){
    p2app$serializeToStaticFast(binary.filename = filename,verbose=TRUE)
  } 
  if(return.details) {
    return(list(app=p2app,p2=cp2,cdl=cdl,cn=cn))
  } else {
    invisible(p2app)
  }
}


