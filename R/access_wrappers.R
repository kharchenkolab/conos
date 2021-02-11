#' Access PCA from sample
#' 
#' @param sample sample from which to access PCA
#' @rdname getPca
#' @export 
setGeneric("getPca", function(sample) standardGeneric("getPca"))

#' @rdname getPca
setMethod("getPca", signature("Pagoda2"), function(sample) sample$reductions$PCA)

#' @rdname getPca
setMethod("getPca", signature("seurat"), function(sample) sample@dr$pca@cell.embeddings)

#' @rdname getPca
setMethod(
  f = 'getPca',
  signature = signature('Seurat'),
  definition = function(sample) {
    checkSeuratV3()
    return(Seurat::Embeddings(object = sample))
  }
)


#' Access overdispersed genes from sample
#' 
#' @param sample sample from which to overdispereed genes
#' @param n.odgenes numeric Number of overdisperesed genes to get
#' @rdname getOverdispersedGenes
#' @export 
setGeneric("getOverdispersedGenes", function(sample, n.odgenes=1000) standardGeneric("getOverdispersedGenes"))

#' @rdname getOverdispersedGenes
setMethod("getOverdispersedGenes", signature("Pagoda2"), function(sample, n.odgenes=NULL) sample$getOdGenes(n.odgenes))

#' @rdname getOverdispersedGenes
setMethod("getOverdispersedGenes", signature("seurat"), function(sample, n.odgenes=NULL)
  if (is.null(NULL)) sample@var.genes else head(rownames(sample@hvg.info), n.odgenes))

#' @importFrom rlang %||%
#' @rdname getOverdispersedGenes
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

#' @rdname getOverdispersedGenes
setMethod("getOverdispersedGenes", signature("Conos"), function(sample, n.odgenes=NULL) commonOverdispersedGene(sample$samples,n.odgenes, verbose=FALSE))


#' Access cell names from sample
#' 
#' @param sample sample from which to cell names
#' @rdname getCellNames
#' @export 
setGeneric("getCellNames", function(sample) standardGeneric("getCellNames"))

#' @rdname getCellNames
setMethod("getCellNames", signature("Pagoda2"), function(sample) rownames(sample$counts))

#' @rdname getCellNames
setMethod("getCellNames", signature("seurat"), function(sample) colnames(sample@data))

#' @rdname getCellNames
setMethod(f = 'getCellNames', signature = signature('Seurat'), definition = function(sample) return(colnames(x = sample)))

#' @rdname getCellNames
setMethod("getCellNames", signature("Conos"), function(sample) unlist(lapply(sample$samples,getCellNames)))


#' Access genes from sample
#' 
#' @param sample sample from which to get genes
#' @rdname getGenes
#' @export 
setGeneric("getGenes", function(sample) standardGeneric("getGenes"))

#' @rdname getGenes
setMethod("getGenes", signature("Pagoda2"), function(sample) colnames(sample$counts))

#' @rdname getGenes
setMethod("getGenes", signature("seurat"), function(sample) rownames(sample@data))

#' @rdname getGenes
setMethod(f = 'getGenes', signature = signature('Seurat'), definition = function(sample) return(rownames(x = sample)))

#' @rdname getGenes
setMethod("getGenes", signature("Conos"), function(sample) unique(unlist(lapply(sample$samples, getGenes))))


#' Set edge matrix edgeMat with certain values on sample
#' 
#' @param sample sample from which to set edge matrix edgeMat with certain values
#' @param value values to set with edgeMat<-
#' @rdname edgeMat
#' @export 
setGeneric("edgeMat<-", function(sample, value) standardGeneric("edgeMat<-"))

#' @rdname edgeMat
setMethod("edgeMat<-", signature("Pagoda2"), function(sample, value) {sample$misc$edgeMat <- value; sample})

#' @rdname edgeMat
setMethod("edgeMat<-", signature("seurat"), function(sample, value) {sample@misc$edgeMat <- value; sample})

#' @rdname edgeMat
setMethod(
  f = 'edgeMat<-',
  signature = signature('Seurat'),
  definition = function(sample, value) {
    checkSeuratV3()
    Seurat::Misc(object = sample, slot = 'edgeMat') <- value
    return(sample)
  }
)



#' Access edgeMat from sample
#' 
#' @param sample sample from which to access edge matrix edgeMat
#' @rdname edgeMat
#' @export 
setGeneric("edgeMat", function(sample) standardGeneric("edgeMat"))

#' @rdname edgeMat
setMethod("edgeMat", signature("Pagoda2"), function(sample) sample$misc$edgeMat)

#' @rdname edgeMat
setMethod("edgeMat", signature("seurat"), function(sample) sample@misc$edgeMat)

#' @rdname edgeMat
setMethod(
  f = 'edgeMat',
  signature = signature('Seurat'),
  definition = function(sample) {
    checkSeuratV3()
    return(Seurat::Misc(object = sample, slot = 'edgeMat'))
  }
)


#' Access count matrix from sample
#' 
#' @param sample sample from which to get the count matrix
#' @param transposed boolean Whether the count matrix should be transposed (default=FALSE)
#' @rdname getCountMatrix
#' @export 
setGeneric("getCountMatrix", function(sample, transposed=FALSE) standardGeneric("getCountMatrix"))

#' @rdname getCountMatrix
setMethod("getCountMatrix", signature("Pagoda2"), function(sample, transposed=FALSE) if (transposed) sample$counts else Matrix::t(sample$counts))

#' @rdname getCountMatrix
setMethod("getCountMatrix", signature("seurat"), function(sample, transposed=FALSE) {
  cm <- if (is.null(sample@scale.data)) sample@data else sample@scale.data
  if (transposed)
    return(Matrix::t(cm))

  return(cm)
})

#' @rdname getCountMatrix
setMethod('getCountMatrix', signature('Seurat'), function(sample, transposed=FALSE) {
    checkSeuratV3()
    dat <- Seurat::GetAssayData(object = sample, slot = 'scale.data')
    dims <- dim(x = dat)
    dat.na <- all(dims == 1) && all(is.na(x = dat))
    if (all(dims == 0) || dat.na) {
      dat <- Seurat::GetAssayData(object = sample, slot = 'data')
    }

    if (transposed)
      return(Matrix::t(dat))

    return(dat)
  }
)



#' Access gene expression from sample
#' 
#' @param sample sample from which to access gene expression
#' @param gene character vector Genes to access
#' @rdname getGeneExpression
#' @export 
setGeneric("getGeneExpression", function(sample, gene) standardGeneric("getGeneExpression"))

#' @rdname getGeneExpression
setMethod("getGeneExpression", signature("Pagoda2"), function(sample, gene) {
  if (gene %in% colnames(sample$counts)) {
    return(sample$counts[, gene])
  }

  return(stats::setNames(rep(NA, nrow(sample$counts)), rownames(sample$counts)))
})

#' @rdname getGeneExpression
setMethod("getGeneExpression", signature("Conos"), function(sample, gene) {
  lapply(sample$samples, getGeneExpression, gene) %>% Reduce(c, .)
})

getGeneExpression.default <- function(sample, gene) {
  count.matrix <- getCountMatrix(sample)
  if(gene %in% rownames(count.matrix)) {
    return(count.matrix[gene,])
  }

  return(stats::setNames(rep(NA, ncol(count.matrix)), colnames(count.matrix)))
}


#' Access raw count matrix from sample
#' 
#' @param sample sample from which to get the raw count matrix
#' @param transposed boolean Whether the raw count matrix should be transposed (default=FALSE)
#' @rdname getRawCountMatrix
#' @export 
setGeneric("getRawCountMatrix", function(sample, transposed=FALSE) standardGeneric("getRawCountMatrix"))

#' @rdname getRawCountMatrix
setMethod("getRawCountMatrix", signature("Pagoda2"), function(sample, transposed=FALSE) if (transposed) sample$misc$rawCounts else t(sample$misc$rawCounts))

#' @rdname getRawCountMatrix
setMethod(
  f = "getRawCountMatrix",
  signature = signature("seurat"),
  definition = function(sample, transposed=FALSE) {
    mi <- match(x = sample@cell.names, table = colnames(sample@raw.data))
    x <- sample@raw.data[, mi, drop = FALSE]
    if (transposed) {
      return(t(x = x))
    } else {
      return(x)
    }
  }
)

#' @rdname getRawCountMatrix
setMethod(
  f = 'getRawCountMatrix',
  signature = signature('Seurat'),
  definition = function(sample, transposed=FALSE) {
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

#' @rdname getRawCountMatrix
setMethod("getRawCountMatrix", signature("Conos"), function(sample, transposed=FALSE) { m <- sample$getJointCountMatrix(raw=TRUE); if (transposed) t(m) else m })



#' Access embedding from sample
#' 
#' @param sample sample from which to get the embedding
#' @param type character Type of embedding to get
#' @rdname getEmbedding
#' @export 
setGeneric("getEmbedding", function(sample, type) standardGeneric("getEmbedding"))

#' @rdname getEmbedding
setMethod("getEmbedding", signature("Pagoda2"), function(sample, type) sample$embeddings$PCA[[type]])

#' @rdname getEmbedding
setMethod("getEmbedding", signature("seurat"), function(sample, type) if (is.null(sample@dr[[type]])) NULL else as.data.frame(sample@dr[[type]]@cell.embeddings))

#' @rdname getEmbedding
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

#' @rdname getEmbedding
setMethod("getEmbedding", signature("Conos"), function(sample, type) sample$embedding)



#' Access clustering from sample
#' 
#' @param sample sample from which to get the clustering
#' @param type character Type of clustering to get
#' @rdname getClustering
#' @export 
setGeneric("getClustering", function(sample, type) standardGeneric("getClustering"))

#' @rdname getClustering
setMethod("getClustering", signature("Pagoda2"), function(sample, type) sample$clusters$PCA[[type]])

#' @rdname getClustering
setMethod("getClustering", signature("seurat"), function(sample, type) {if (!is.null(type)) warning("Seurat support only single type of clustering"); sample@ident})

#' @rdname getClustering
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

#' @rdname getClustering
setMethod("getClustering", signature("Conos"), function(sample, type) { cl <- sample$clusters[[type]]; if(is.null(cl)) NULL else cl$groups })
