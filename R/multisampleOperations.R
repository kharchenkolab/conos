#' Subset a list of pagoda2 application down to specified genes
#' @param r.n list of pagoda2 applications
#' @param genes list of genes
#' @param remove.null.apps logical, remove apps that end up with no cells?
#' @return a list of apps
#' @export subsetAllappsToGenes
subsetAllappsToGenes <- function(r.n, genes, remove.null.apps = TRUE) {
  ret.apps <- lapply(r.n, function(o) {
    tryCatch({
      ret <- NULL
      app.genes.keep <- genes[genes %in% colnames(o$misc$rawCounts)]
      if (length(app.genes.keep) > 10) {
        p2 <- Pagoda2$new(t(o$misc$rawCounts[,app.genes.keep]), n.cores=32)
        p2$adjustVariance(plot=F,gam.k=20)
        p2$calculatePcaReduction(nPcs=100,n.odgenes=2000,maxit=3000)
        p2$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T);
        p2$makeKnnGraph(k=30, type='PCA', center=T, weight.type='none', n.cores=32, distance='cosine')
        p2$getKnnClusters(method = infomap.community, type = 'PCA' ,name = 'infomap')
        ret <- p2
      }
      ret
    }, error = function(e) {
      NULL
    })
  })
  if(remove.null.apps) ret.apps <- removeNullapps(ret.apps)
}

#' Subset all apps to clusters
#' @export subsetAllappsToClusters
subsetAllappsToClusters <- function(r.n, cl, cl.keep) {
  cells.keep <- names(cl)[cl %in% cl.keep]
  cat(length(cells.keep),'\n')
  lapply(r.n, function(o) {
    tryCatch({
      rn <- rownames(o$misc$rawCounts)
      app.cells.keep <- rn %in% cells.keep
      if (length(app.cells.keep) < 100) { NULL  }
      p2 <- Pagoda2$new(t(o$misc$rawCounts[app.cells.keep,]), n.cores=20)
      p2$adjustVariance(plot=F,gam.k=20)
      p2$calculatePcaReduction(nPcs=100,n.odgenes=1000,maxit=3000)
      p2$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T);
      p2$makeKnnGraph(k=30, type='PCA', center=T, weight.type='none', n.cores=20, distance='cosine')
      p2$getKnnClusters(method = infomap.community, type = 'PCA' ,name = 'infomap')
      p2
    }, warning = function(w) {
      NULL
    }, error = function(e) {
      NULL
    })
  })
}

#' Get common genes between multiple pagoda apps
#' @param r.n list of pagoda2 apps
#' @param cutoff minimum number of apps that must have genes
#' @return character vector of common genes
#' @export getCommonGenesCutoff
getCommonGenesCutoff <- function(r.n, cutoff = 3) {
  gl <- lapply(r.n, function(r) colnames(r$counts))
  all.genes <- unique(unlist(gl))
  gc <- do.call(rbind, lapply(gl, function(x) all.genes %in% x))
  common.genes <- all.genes[apply(gc,2,sum) > cutoff]
  common.genes
}

#' Get the genes that the apps have in common
#' @param r.n a list of pagoda2 apps
#' @return a character vector of common genes
#' @export getCommonGenes
getCommonGenes <- function(r.n) {
  Reduce(intersect, lapply(r.n, function(o) { colnames(o$counts) }))
}

#' From a list of pagoda2 application remove any that are NULL
#' @param os list of pagoda2 applications
#' @return list of pagoda2 application filtered for NULLs
#' @export removeNullapps
removeNullapps <- function(os) {
  os[!unlist(lapply(os,FUN=is.null))]
}
