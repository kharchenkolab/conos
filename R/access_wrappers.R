setGeneric("getPca", function(sample) standardGeneric("getPca"))
setMethod("getPca", signature("Pagoda2"), function(sample) sample$reductions$PCA)
setMethod("getPca", signature("seurat"), function(sample) sample@dr$pca@cell.embeddings)

setGeneric("getOverdispersedGenes", function(sample, n.odgenes=1000) standardGeneric("getOverdispersedGenes"))
setMethod("getOverdispersedGenes", signature("Pagoda2"), function(sample, n.odgenes=NULL) sample$getOdGenes(n.odgenes))
setMethod("getOverdispersedGenes", signature("seurat"), function(sample, n.odgenes=NULL)
  if (is.null(NULL)) sample@var.genes else head(rownames(sample@hvg.info), n.odgenes))

setGeneric("getCellNames", function(sample) standardGeneric("getCellNames"))
setMethod("getCellNames", signature("Pagoda2"), function(sample) rownames(sample$counts))
setMethod("getCellNames", signature("seurat"), function(sample) colnames(sample@data))

setGeneric("getGenes", function(sample) standardGeneric("getGenes"))
setMethod("getGenes", signature("Pagoda2"), function(sample) colnames(sample$counts))
setMethod("getGenes", signature("seurat"), function(sample) rownames(sample@data))

setGeneric("edgeMat<-", function(sample, value) standardGeneric("edgeMat<-"))
setMethod("edgeMat<-", signature("Pagoda2"), function(sample, value) {sample$misc$edgeMat <- value; sample})
setMethod("edgeMat<-", signature("seurat"), function(sample, value) {sample@misc$edgeMat <- value; sample})

setGeneric("edgeMat", function(sample, value) standardGeneric("edgeMat"))
setMethod("edgeMat", signature("Pagoda2"), function(sample) sample$misc$edgeMat)
setMethod("edgeMat", signature("seurat"), function(sample) sample@misc$edgeMat)

setGeneric("getCountMatrix", function(sample) standardGeneric("getCountMatrix"))
setMethod("getCountMatrix", signature("Pagoda2"), function(sample) t(sample$counts))
setMethod("getCountMatrix", signature("seurat"), function(sample) sample@data)

setGeneric("getRawCountMatrix", function(sample, transposed=F) standardGeneric("getRawCountMatrix"))
setMethod("getRawCountMatrix", signature("Pagoda2"), function(sample, transposed=F) if (transposed) sample$misc$rawCounts else t(sample$misc$rawCounts))
setMethod("getRawCountMatrix", signature("seurat"), function(sample, transposed=F) { mi <- match(sample@cell.names,colnames(sample@raw.data)); x <- sample@raw.data[,mi,drop=FALSE]; if(transposed) return(t(x)) else return(x) })

setGeneric("getEmbedding", function(sample, type) standardGeneric("getEmbedding"))
setMethod("getEmbedding", signature("Pagoda2"), function(sample, type) sample$embeddings$PCA[[type]])
setMethod("getEmbedding", signature("seurat"), function(sample, type) if (is.null(sample@dr[[type]])) NULL else as.data.frame(sample@dr[[type]]@cell.embeddings))

setGeneric("getClustering", function(sample, type) standardGeneric("getClustering"))
setMethod("getClustering", signature("Pagoda2"), function(sample, type) sample$clusters$PCA[[type]])
setMethod("getClustering", signature("seurat"), function(sample, type) {if (!is.null(type)) warning("Seurat support only single type of clustering"); sample@ident})
