# setGeneric("getCounts", function(sample) {
#   standardGeneric("getCounts")
# })

setGeneric("getPca", function(sample) {
  standardGeneric("getPca")
})

setGeneric("getScaledCounts", function(sample) {
  standardGeneric("getScaledCounts")
})

setGeneric("getOverdispersedGenes", function(sample, n.odgenes=1000) {
  standardGeneric("getOverdispersedGenes")
})

setGeneric("getGenes", function(sample) {
  standardGeneric("getGenes")
})

# setMethod("getCounts", signature("Pagoda2"), function(sample) sample$counts)
# setMethod("getCounts", signature("seurat"), function(sample) t(sample@data))

setMethod("getPca", signature("Pagoda2"), function(sample) sample$reductions$PCA)
setMethod("getPca", signature("seurat"), function(sample) sample@dr$pca@cell.embeddings)

setMethod("getScaledCounts", signature("Pagoda2"), function(sample) stop("Not implemented"))
setMethod("getScaledCounts", signature("seurat"), function(sample) t(sample@scale.data))

setMethod("getOverdispersedGenes", signature("Pagoda2"), function(sample, n.odgenes=NULL) sample$getOdGenes(n.odgenes))
setMethod("getOverdispersedGenes", signature("seurat"), function(sample, n.odgenes=NULL) if (is.null(NULL)) sample@var.genes else sample@var.genes[1:min(n.odgenes, length(sample@var.genes))])

setMethod("getGenes", signature("Pagoda2"), function(sample) colnames(sample$counts))
setMethod("getGenes", signature("seurat"), function(sample) rownames(sample@data))