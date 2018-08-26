# setGeneric("getCounts", function(sample) {
#   standardGeneric("getCounts")
# })

# setMethod("getCounts", signature("Pagoda2"), function(sample) sample$counts)
# setMethod("getCounts", signature("seurat"), function(sample) t(sample@data))

setGeneric("getPca", function(sample) standardGeneric("getPca"))
setMethod("getPca", signature("Pagoda2"), function(sample) sample$reductions$PCA)
setMethod("getPca", signature("seurat"), function(sample) sample@dr$pca@cell.embeddings)

setGeneric("getScaledCounts", function(sample) standardGeneric("getScaledCounts"))
setMethod("getScaledCounts", signature("Pagoda2"), function(sample) stop("Not implemented"))
setMethod("getScaledCounts", signature("seurat"), function(sample) t(sample@scale.data))

setGeneric("getOverdispersedGenes", function(sample, n.odgenes=1000) standardGeneric("getOverdispersedGenes"))
setMethod("getOverdispersedGenes", signature("Pagoda2"), function(sample, n.odgenes=NULL) sample$getOdGenes(n.odgenes))
setMethod("getOverdispersedGenes", signature("seurat"), function(sample, n.odgenes=NULL) if (is.null(NULL)) sample@var.genes else sample@var.genes[1:min(n.odgenes, length(sample@var.genes))])

setGeneric("getGenes", function(sample) standardGeneric("getGenes"))
setMethod("getGenes", signature("Pagoda2"), function(sample) colnames(sample$counts))
setMethod("getGenes", signature("seurat"), function(sample) rownames(sample@data))

setGeneric("edgeMat<-", function(sample, value) standardGeneric("edgeMat<-"))
setMethod("edgeMat<-", signature("Pagoda2"), function(sample, value) {sample$misc$edgeMat <- value; sample})
setMethod("edgeMat<-", signature("seurat"), function(sample, value) {sample@misc$edgeMat <- value; sample})

setGeneric("edgeMat", function(sample, value) standardGeneric("edgeMat"))
setMethod("edgeMat", signature("Pagoda2"), function(sample) sample$misc$edgeMat)
setMethod("edgeMat", signature("seurat"), function(sample) sample@misc$edgeMat)