### Functions for working with multiple pagoda2 objects

#' Perform gsea on a pagoda2 table of results
#' @param tbl pagoda2 table of results
#' @param mc.cores number of CPUs to use
#' @export gsea.on.p2.table
gsea.on.p2.table <- function(tbl, mc.cores =32) {
    require('liger')
    zsrs <- tbl$Z
    names(zsrs) <- rownames(tbl)
    go.l <- lapply(org.Hs.GO2Symbol.list, length)
    go.keep <- go.l > 10 & go.l < 1000
    go.sets <- org.Hs.GO2Symbol.list[go.keep]
    r.gsea <- bulk.gsea(zsrs, set.list=go.sets, mc.cores=mc.cores)
    ## Filter
    r.gsea <- r.gsea[r.gsea$q.val < 0.001,]
    ## Annotate with description
    descriptions <- lapply(rownames(r.gsea), function(x) {
        tryCatch({
            gt <- GOTERM[[x]]
            ret <- NULL
            if (!is.null(gt)) ret <- gt@Term
            ret
        }, warning = function(w) {
            NULL
        }, error = function(e) {
            NULL
        })
    })
    r.gsea$desc <- descriptions
    r.gsea[order(r.gsea$q.val,decreasing=F),]
}

#' Get number of differential genes at specified cutoff
#' @export getDEcountAtCutoff
getDEcountAtCutoff <- function(de.res, z.cutoff=6, M.cutoff=2) {
    de.c <- lapply(cl.merged.clean.renamed_wilcox, function(cl.list) {
        lapply(cl.list, function(x) {
            ret <- 0
            if(!is.null(x$s1 )) {
                s1 <- x$s1
                s1 <- abs(s1$M) > M.cutoff & abs(s1$Z) > z.cutoff
                ret <- dim(x$s1)[1]
            } else {
                0
            }
        })
    })
    de.counts <- do.call("cbind", de.c)
    de.counts
}

#' Get hierarchy of cell groups across multiple apps
#' @param r.n list of pagoda2 objects
#' @param clusters factor of clusters/cell grousp to summarise to
#' @export multiSampleClusterHierarchy
multiSampleClusterHierarchy <- function(r.n, cls) {
    cms <- lapply(r.n, function(o) { o$counts })
    common.genes <- Reduce(intersect, lapply(r.n, function(o) { colnames(o$counts) }))
    ## genes are rows now
    bcm <- do.call(rbind, lapply(cms, function(o) {(o[,common.genes])}))
    cls <- cls[rownames(bcm)]
    clsums <- pagoda2:::colSumByFac(bcm, cls)
    ## Set the gene names
    colnames(clsums) <- colnames(bcm)
    ## Remove NA sum
    clsums <- clsums[-1,]
    rownames(clsums) <- levels(cls)
    ## Get numbers of cells in each cluster
    cl.counts <- table(cls)[levels(cls)]
    ## Get normalised cluster centers
    clsums.norm <- sweep(clsums, 1, cl.counts, FUN='/')
    ## Get correlation distance dendrogram
    hc.cor.euclidean <- hclust(dist(clsums.norm))
    hc.cor.pearson <- hclust(as.dist(1-cor(t(clsums.norm),method='pearson')))
    hc.cor.spearman <- hclust(as.dist(1-cor(t(clsums.norm),method='spearman')))
    hc <- list(hc.cor.euclidean=hc.cor.euclidean,hc.cor.pearson=hc.cor.pearson,hc.cor.spearman=hc.cor.spearman)
    hc
}

#' Get differential expression markers from a pagoda2 differntial expression result
#' @description return the cluster-specific upregulated genes per cluster
#' @param de.res pagoda2 differential expression result
#' @return a list of markers
#' @export getMarkersFromDE
getMarkersFromDE <- function(de.res) {
    lapply(de.res, function(x) {
        (rownames(subset(x,highest==TRUE,Z>0)))
    })
}

#' Get cluster specific markers from a pagoda2 differential expression results
#' @param de.res pagoda2 differential expression result
#' @return list of character vectors of marker genes
#' @export getMarkersFromDE
getMarkersFromDE <- function(de.res) {
    lapply(de.res, function(x) {
        (rownames(subset(x,highest==TRUE,Z>0)))
    })
}

