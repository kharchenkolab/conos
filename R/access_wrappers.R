setGeneric("getPca", function(sample) standardGeneric("getPca"))
setMethod("getPca", signature("Pagoda2"), function(sample) sample$reductions$PCA)
setMethod("getPca", signature("seurat"), function(sample) sample@dr$pca@cell.embeddings)
setMethod(
  f = 'getPca',
  signature = signature('Seurat'),
  definition = function(sample) {
    checkSeuratV3()
    return(Seurat::Embeddings(object = sample))
  }
)

setGeneric("getOverdispersedGenes", function(sample, n.odgenes=1000) standardGeneric("getOverdispersedGenes"))
setMethod("getOverdispersedGenes", signature("Pagoda2"), function(sample, n.odgenes=NULL) sample$getOdGenes(n.odgenes))
setMethod("getOverdispersedGenes", signature("seurat"), function(sample, n.odgenes=NULL)
  if (is.null(NULL)) sample@var.genes else head(rownames(sample@hvg.info), n.odgenes))

#' @importFrom rlang %||%
#'
setMethod(
  f = 'getOverdispersedGenes',
  signature = signature('Seurat'),
  definition = function(sample, n.odgenes = NULL) {
    checkSeuratV3()
    vf <- Seurat::VariableFeatures(object = sample) %||% rownames(x = sample)
    n.odgenes <- n.odgenes %||% length(x = vf)
    return(head(x = vf, n = n.odgenes))
  }
)

setGeneric("getCellNames", function(sample) standardGeneric("getCellNames"))
setMethod("getCellNames", signature("Pagoda2"), function(sample) rownames(sample$counts))
setMethod("getCellNames", signature("seurat"), function(sample) colnames(sample@data))
setMethod(f = 'getCellNames', signature = signature('Seurat'), definition = function(sample) return(colnames(x = sample)))

setGeneric("getGenes", function(sample) standardGeneric("getGenes"))
setMethod("getGenes", signature("Pagoda2"), function(sample) colnames(sample$counts))
setMethod("getGenes", signature("seurat"), function(sample) rownames(sample@data))
setMethod(f = 'getGenes', signature = signature('Seurat'), definition = function(sample) return(rownames(x = sample)))

setGeneric("edgeMat<-", function(sample, value) standardGeneric("edgeMat<-"))
setMethod("edgeMat<-", signature("Pagoda2"), function(sample, value) {sample$misc$edgeMat <- value; sample})
setMethod("edgeMat<-", signature("seurat"), function(sample, value) {sample@misc$edgeMat <- value; sample})
setMethod(
  f = 'edgeMat<-',
  signature = signature('Seurat'),
  definition = function(sample, value) {
    checkSeuratV3()
    Seurat::Misc(object = sample, slot = 'edgeMat') <- value
    return(sample)
  }
)

setGeneric("edgeMat", function(sample, value) standardGeneric("edgeMat"))
setMethod("edgeMat", signature("Pagoda2"), function(sample) sample$misc$edgeMat)
setMethod("edgeMat", signature("seurat"), function(sample) sample@misc$edgeMat)
setMethod(
  f = 'edgeMat',
  signature = signature('Seurat'),
  definition = function(sample) {
    checkSeuratV3()
    return(Seurat::Misc(object = sample, slot = 'edgeMat'))
  }
)

setGeneric("getCountMatrix", function(sample) standardGeneric("getCountMatrix"))
setMethod("getCountMatrix", signature("Pagoda2"), function(sample) t(sample$counts))
setMethod("getCountMatrix", signature("seurat"), function(sample) if (is.null(sample@scale.data)) sample@data else sample@scale.data)
setMethod(
  f = 'getCountMatrix',
  signature = signature('Seurat'),
  definition = function(sample) {
    checkSeuratV3()
    dat <- Seurat::GetAssayData(object = sample, slot = 'scale.data')
    dims <- dim(x = dat)
    dat.na <- all(dims == 1) && all(is.na(x = dat))
    if (all(dims == 0) || dat.na) {
      dat <- Seurat::GetAssayData(object = sample, slot = 'data')
    }
    return(dat)
  }
)

setGeneric("getRawCountMatrix", function(sample, transposed=F) standardGeneric("getRawCountMatrix"))
setMethod("getRawCountMatrix", signature("Pagoda2"), function(sample, transposed=F) if (transposed) sample$misc$rawCounts else t(sample$misc$rawCounts))
setMethod(
  f = "getRawCountMatrix",
  signature = signature("seurat"),
  definition = function(sample, transposed = F) {
    mi <- match(x = sample@cell.names, table = colnames(sample@raw.data))
    x <- sample@raw.data[, mi, drop = FALSE]
    if (transposed) {
      return(t(x = x))
    } else {
      return(x)
    }
  }
)
setMethod(
  f = 'getRawCountMatrix',
  signature = signature('Seurat'),
  definition = function(sample, transposed = FALSE) {
    checkSeuratV3()
    rd <- Seurat::GetAssayData(object = sample, slot = 'counts')
    # Raw data can be empty in Seurat v3
    # If it is, use data instead
    dims <- dim(x = rd)
    rd.na <- all(dims == 1) && all(is.na(x = rd))
    if (all(dims == 0) || rd.na) {
      rd <- Seurat::GetAssayData(object = sample, slot = 'data')
    }
    mi <- match(x = colnames(x = sample), table = colnames(x = rd))
    rd <- rd[, mi, drop = FALSE]
    if (transposed) {
      rd <- t(x = rd)
    }
    return(rd)
  }
)

setGeneric("getEmbedding", function(sample, type) standardGeneric("getEmbedding"))
setMethod("getEmbedding", signature("Pagoda2"), function(sample, type) sample$embeddings$PCA[[type]])
setMethod("getEmbedding", signature("seurat"), function(sample, type) if (is.null(sample@dr[[type]])) NULL else as.data.frame(sample@dr[[type]]@cell.embeddings))
setMethod(
  f = 'getEmbedding',
  signature = signature('Seurat'),
  definition = function(sample, type) {
    checkSeuratV3()
    emb <- tryCatch(
      expr = Seurat::Embeddings(object = sample, reduction = type),
      error = function(...) {
        return(NULL)
      }
    )
    return(emb)
  }
)

setGeneric("getClustering", function(sample, type) standardGeneric("getClustering"))
setMethod("getClustering", signature("Pagoda2"), function(sample, type) sample$clusters$PCA[[type]])
setMethod("getClustering", signature("seurat"), function(sample, type) {if (!is.null(type)) warning("Seurat support only single type of clustering"); sample@ident})
setMethod(
  f = 'getClustering',
  signature = signature('Seurat'),
  definition = function(sample, type) {
    checkSeuratV3()
    if (missing(x = type)) {
      type <- NULL
    } else if (!is.null(x = type) && !type %in% colnames(x = sample[[]])) {
      warning(
        "Cannot find ",
        type,
        " in sample metadata, using normal identities",
        call. = FALSE,
        immediate. = TRUE
      )
      type <- NULL
    }
    idents <- if (is.null(x = type)) {
      Seurat::Idents(object = sample)
    } else {
      ids <- sample[[type]]
      if (!is.factor(x = ids)) {
        ids <- factor(x = ids)
      }
      ids
    }
    return(idents)
  }
)
