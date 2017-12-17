
## From the fossil package
rand.index <- function (group1, group2)
{
    x <- abs(sapply(group1, function(x) x - group1))
    x[x > 1] <- 1
    y <- abs(sapply(group2, function(x) x - group2))
    y[y > 1] <- 1
    sg <- sum(abs(x - y))/2
    bc <- choose(dim(x)[1], 2)
    ri <- 1 - sg/bc
    return(ri)
}

#' Assess cross sample clustering stability by dropping
#' each sample sequentially and checking for consistent classification
#' of remaining samples
#' @param p2.objects the pagoda2 samples
#' @param n.samples.drop number of samples to drop
assessDropOneStability <- function(p2.objects, clusteringFunction,  n.samples.drop = NULL, drop.order = NULL, clusteringFunctionParams=list() ) {


    if (mode(clusteringFunction) != 'function') stop('clusteringFunction is not a function')
    
    if (is.null(drop.order)) {
        ## drop all one by one by default
        if (is.null(n.samples.drop)) n.samples.drop <- length(p2.objects);
        if (n.samples.drop > length(p2.objects)) {
            warning("can't drop more sampels that objects provided");
        }
        drop.order <- sample(length(p2.objects), n.samples.drop);
    }

    ## Calculate with all samples present
    cl.fn.param <- clusteringFunctionParams;
    cl.fn.param$p2.objects <- p2.objects;
    all.cls <- do.call(clusteringFunction, cl.fn.param)

    ## Calculate with samles absent
    subset.cls <- lapply(drop.order, function(drop.sample) {
        ## Get clusters in the absence of this sample
        samples.include <- seq_along(p2.objects)
        samples.include <- samples.include[samples.include != drop.sample]
        ## Prepare the parameters
        cl.fn.param <- clusteringFunctionParams;
        cl.fn.param$p2.objects <- p2.objects[samples.include];
        ## Run the function and record time
        iter.cl <- do.call(clusteringFunction, cl.fn.param);
        ## Return the list
        iter.cl
    })

    ## Get Rand index for clustering of cell in both
    rand.idxs <- lapply(subset.cls, function(cl) {
        ## compare cl to all.cls
        cell.use <- intersect(names(cl)[!is.na(cl)], names(all.cls)[!is.na(all.cls)])
        ri <- rand.index(as.numeric(cl[cell.use]),as.numeric(all.cls[cell.use]))
        ri
    })

    rand.idxs <- unlist(rand.idxs)

    list( m = mean(rand.idxs), sd = sd(rand.idxs), idx=rand.idxs)
}

##' Assess stability by removing one joint cluster from all samples
##' at a time and reclustering
assessStabilityRemoveClusters <- function(p2.objects, clusteringFunction, cluster.keep.pc.cutoff = 0.05,
                                          clusteringFunctionParams = list()) {
    if (mode(clusteringFunction) != 'function') stop('clusteringFunction is not a function')
    ## Calculate with all samples present
    cl.fn.param <- clusteringFunctionParams;
    cl.fn.param$p2.objects <- p2.objects;
    all.cls <- do.call(clusteringFunction, cl.fn.param)
    cl.mem.counts <- table(all.cls)
    n.cells <- length(all.cls)
    cluster.count.cutoff <- floor(n.cells * cluster.keep.pc.cutoff)
    ## Drop the clusters that over at least cluster.keep.pc.cutoff
    clusters.drop <- names(which(cl.mem.counts > cluster.count.cutoff))
    ## Outer lapply for dropping clusters one by one
    subset.cls <- lapply(clusters.drop, function(cl.drop) {
        cells.drop <- names(all.cls)[which(all.cls == cl.drop)]
        ## Regenerate the pagoda2 apps with out the cells in the designated global cluster
        subsetp2.objs <- lapply(p2.objects, function(p2o) {            
            ## Get matrix in correct orientation with cells removed
            nm <- p2o$misc$raw
            nm <- t(nm[!rownames(nm) %in% cells.drop,])
            ## Make a new pagoda object (Object New)
            on <- Pagoda2$new(nm);
            on$adjustVariance(plot=F,gam.k=10)
            on$calculatePcaReduction(nPcs=100, n.odgenes=3e3, maxit=3000)
            on
        })
        ## Perform new joint clustering on these subsetted objects
        cl.fn.param.1 <- clusteringFunctionParams;
        cl.fn.param.1$p2.objects <- subsetp2.objs;
        ## Generate new clustering
        new.cl <- do.call(clusteringFunction, cl.fn.param.1)
        new.cl
    })
    ## Get Rand index for clustering of cell in both
    rand.idxs <- lapply(subset.cls, function(cl) {
        ## compare cl to all.cls
        cell.use <- intersect(names(cl)[!is.na(cl)], names(all.cls)[!is.na(all.cls)])
        ri <- rand.index(as.numeric(cl[cell.use]),as.numeric(all.cls[cell.use]))
        ri
    })
    ## Make results
    rand.idxs <- unlist(rand.idxs)
    list( m = mean(rand.idxs), sd = sd(rand.idxs), idx=rand.idxs)    
}


##' Assess stability by shrinking one joint cluster from all samples
##' at a time and reclustering
assessStabilityShrinkClusters <- function(p2.objects, clusteringFunction, cluster.keep.pc.cutoff = 0.05,
                                          clusteringFunctionParams = list(), clusterShrinkFactor = 0.3) {
    if (mode(clusteringFunction) != 'function') stop('clusteringFunction is not a function')
    ## Calculate with all samples present
    cl.fn.param <- clusteringFunctionParams;
    cl.fn.param$p2.objects <- p2.objects;
    all.cls <- do.call(clusteringFunction, cl.fn.param)
    cl.mem.counts <- table(all.cls)
    n.cells <- length(all.cls)
    cluster.count.cutoff <- floor(n.cells * cluster.keep.pc.cutoff)
    ## Drop the clusters that over at least cluster.keep.pc.cutoff
    clusters.drop <- names(which(cl.mem.counts > cluster.count.cutoff))
    ## Outer lapply for dropping clusters one by one
    subset.cls <- lapply(clusters.drop, function(cl.drop) {
        cells.drop <- names(all.cls)[which(all.cls == cl.drop)]
        ## Sample the cells to drop according to clusterShrinkFactor
        cells.drop <- sample(cells.drop, floor(length(cells.drop) * (1 - clusterShrinkFactor)))
        ## Regenerate the pagoda2 apps with out the cells in the designated global cluster
        subsetp2.objs <- lapply(p2.objects, function(p2o) {            
            ## Get matrix in correct orientation with cells removed
            nm <- p2o$misc$raw
            nm <- t(nm[!rownames(nm) %in% cells.drop,])
            ## Make a new pagoda object (Object New)
            on <- Pagoda2$new(nm);
            on$adjustVariance(plot=F,gam.k=10)
            on$calculatePcaReduction(nPcs=100, n.odgenes=3e3, maxit=3000)
            on
        })
        ## Perform new joint clustering on these subsetted objects
        cl.fn.param.1 <- clusteringFunctionParams;
        cl.fn.param.1$p2.objects <- subsetp2.objs;
        ## Generate new clustering
        new.cl <- do.call(clusteringFunction, cl.fn.param.1)
        new.cl
    })
    ## Get Rand index for clustering of cell in both
    rand.idxs <- lapply(subset.cls, function(cl) {
        ## compare cl to all.cls
        cell.use <- intersect(names(cl)[!is.na(cl)], names(all.cls)[!is.na(all.cls)])
        ri <- rand.index(as.numeric(cl[cell.use]),as.numeric(all.cls[cell.use]))
        ri
    })
    ## Make results
    rand.idxs <- unlist(rand.idxs)
    list( m = mean(rand.idxs), sd = sd(rand.idxs), idx=rand.idxs)    
}

##' Assess stability by removing shifting the cluster composition
##' @param p2.objects list of pagoda2 objects
##' @param clusteringFunction the function that does the joint clusting
##' @param N the number of repeats to run
##' @param clusteringFunctionParams a list of arguments to pass to the clusting fucntion
assessStabilityShiftClusterComposition <- function(p2.objects, clusteringFunction, N = 3,
                                          clusteringFunctionParams = list()) {

    if (mode(clusteringFunction) != 'function') stop('clusteringFunction is not a function')

    ## Calculate jc with all samples present
    cl.fn.param <- clusteringFunctionParams;
    cl.fn.param$p2.objects <- p2.objects;
    all.cls <- do.call(clusteringFunction, cl.fn.param)

    subset.cls <- lapply(seq(1,N), function(p0) {
        
        ## Regenerate the pagoda2 apps with with shifted cluster compositions
        subsetp2.objs <- lapply(p2.objects, function(p2o) {
            ## DEVEL
            ## p2o <- p2.objects[[1]]
            
            ## Get matrix in correct orientation with cells removed
            nm <- p2o$misc$raw

            ## All the cells in this cluster in this app
            all.cls.app <- all.cls[rownames(nm)]

            ## For every cluster, keep a variable fraction
            keep.cells <- unlist(lapply(levels(all.cls.app), function(cls) {
                ## All the cells in this cluster in this app
                cl.cells.app <- names(all.cls.app)[all.cls.app == cls]
                ## The fraction of cells to keep in this app and this cluster
                fraction.keep <- runif(1)
                cl.cells.app <- sample(cl.cells.app, floor(length(cl.cells.app) * fraction.keep))
                cl.cells.app
            }))
            ## Not sure why we need this, but we end up with 7 NAs
            keep.cells <- keep.cells[!is.na(keep.cells)]
            
            ## Subset original clustering to the cells in this app
            nm <- t(nm[keep.cells,])
            
            ## Make a new pagoda object (Object New)
            on <- Pagoda2$new(nm);
            on$adjustVariance(plot=F,gam.k=10)
            on$calculatePcaReduction(nPcs=100, n.odgenes=3e3, maxit=3000)
            on
        })

        ## Perform new joint clustering on these subsetted objects
        cl.fn.param.1 <- clusteringFunctionParams;
        cl.fn.param.1$p2.objects <- subsetp2.objs;
        ## Generate new clustering
        new.cl <- do.call(clusteringFunction, cl.fn.param.1)

        new.cl
    })

    ## Get Rand index for clustering of cell in both
    rand.idxs <- lapply(subset.cls, function(cl) {
        ## compare cl to all.cls
        cell.use <- intersect(names(cl)[!is.na(cl)], names(all.cls)[!is.na(all.cls)])
        ri <- rand.index(as.numeric(cl[cell.use]),as.numeric(all.cls[cell.use]))
        ri
    })
    
    ## Make results
    rand.idxs <- unlist(rand.idxs)
    list( m = mean(rand.idxs), sd = sd(rand.idxs), idx=rand.idxs)    
}


####################################################
## RUN
###################################################

## library('nbHelpers')
## library('pagoda2')
## library('Rjnmf')

## ## Load a dataset
## bm.sampleset1 <- readRDS('bm.sampleset1.rds')

## ## Wrapper for clustering function, set maxiter to 20 for speed
## testWrapper <- function(p2.objects, ...) {
##     jnmfJCp(r.n = p2.objects, maxiter = 20, ...);
## }

## ## Run the tests here
## t0 <- Sys.time()
## test1 <- assessDropOneStability(p2.objects = bm.sampleset1[1:3], clusteringFunction = testWrapper)
## test2 <- assessStabilityRemoveClusters(p2.objects = bm.sampleset1[1:3], clusteringFunction = testWrapper)
## test3 <- assessStabilityShiftClusterComposition(bm.sampleset1[1:3], testWrapper,N=1)
## t1 <- Sys.time()

## test1
## test2
## test3

## t1 - t0
