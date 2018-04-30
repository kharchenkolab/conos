#' Postprocess joint clusters generated with cpcaJCp3
#' @param jfac result of cpcaJCp3 with return.details TRUE and run with walktrap (or other
#' hierarchical method)
#' @param p2ens a p2ensembl object containing the p2 objects in the p2objs slot as a list
#' @param no.cl number of clusters to request from walktrap
#' @param size.cutoff minumum number of cells in a joint cluster to use to reassigned
#' @param n.cores number of cores to use
#' @return a named factor of adjusted memberships
postProcessWalktrapClusters <- function(jfac, p2ens, no.cl=200, size.cutoff=10, n.cores=4) {
    cls <- jfac$cls
    p2objs <- p2ens$p2objs
    ##
    library(Matrix)
    factorBreakdown <- function(f) {tapply(names(f),f, identity) }
    ## Get joint global clusters at requested number
    global.cluster <- cut_at(cls, no=no.cl)
    names(global.cluster) <-  names(membership(jcl3$cls))
    ## identify clusters to merge
    fqs <- as.data.frame(table(global.cluster))
    cl.to.merge <- fqs[fqs$Freq < size.cutoff,]$global.cluster
    cl.to.keep <- fqs[fqs$Freq >= size.cutoff,]$global.cluster
    ## Memberships to keep
    global.cluster.filtered <- as.factor(global.cluster[global.cluster %in% cl.to.keep])
    ## Get new assignments for all the cells
    new.assign <- unlist(unname(parallel::mclapply(p2objs, function(p2o) {
        try({
            ## get global cluster filter centroids for cells in this app
            global.cluster.filtered.bd <- factorBreakdown(global.cluster.filtered)
            global.cl.centers <- do.call(rbind,lapply(global.cluster.filtered.bd ,function(cells) {
                cells <- cells[cells %in% rownames(p2o$counts)]
                cells
                if(length(cells) > 1) {
                    colSums(p2o$counts[cells,])
                } else {
                    NULL
                }
            }))
            ## cells to reassign in this app
            cells.reassign <- names(global.cluster[global.cluster %in% cl.to.merge])
            cells.reassign <- cells.reassign[cells.reassign %in% rownames(p2o$counts)]
            xcor <- cor(t(as.matrix(p2o$counts[cells.reassign,,drop=FALSE])), t(as.matrix(global.cl.centers)))
            ## Get new cluster assignments
            new.cluster.assign <- apply(xcor,1, function(x) {colnames(xcor)[which.max(x)]})
            new.cluster.assign
        })
    },mc.cores=n.cores)))
    ## Merge
    x <- as.character(global.cluster.filtered)
    names(x) <- names(global.cluster.filtered)
    new.clusters <- as.factor(c(x,new.assign))
    new.clusters
}
