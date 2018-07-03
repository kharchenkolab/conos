library('Rcpp')

Sys.setenv('NMSLIB_PATH'='/home/barkasn/lib/nmslib-1.6/')
Sys.setenv('PKG_CXXFLAGS'='-I"$(NMSLIB_PATH)/similarity_search/include" -I"../inst/include" -I"$(NMSLIB_PATH)/similarity_search/include" -I "$(NMSLIB_PATH)/similarity_search/lshkit/include" $(SHLIB_OPENMP_CXXFLAGS)')
Sys.setenv('PKG_LIBS'='-L/usr/lib/ -L$(NMSLIB_PATH)/similarity_search/release -lNonMetricSpaceLib -lgsl -lpthread -lboost_filesystem -lboost_system -lstdc++ `$(R_HOME)/bin/Rscript -e "Rcpp:::LdFlags()"` $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SHLIB_OPENMP_CXXFLAGS) -lm')
Sys.setenv('CXX_STD' = 'CXX11')

## compile CPP code
sourceCpp('panelClust.cpp')


#' Calculate sparse matrix covariance
sparse.cov <- function(x, cMeans = NULL) {
    n <- nrow(x);
    if(is.null(cMeans)) {  cMeans <- Matrix::colMeans(x) }
    covmat <- (as.matrix(Matrix::crossprod(x)) - n*Matrix::tcrossprod(cMeans))/(n-1)
}


quickNULL <- function(p2.objs = NULL, n.odgenes = NULL, var.scale = T,
                      verbose = TRUE, neighborhood.average=FALSE) {
    if(length(p2.objs) != 2) stop('quickNULL only supports pairwise alignment');
    ## Get common set of genes
    if(is.null(n.odgenes)) {
        odgenes <- table(unlist(lapply(p2.objs,function(x) x$misc$odgenes)))
    } else {
        odgenes <- table(unlist(lapply(p2.objs,function(x) rownames(x$misc$varinfo)[(order(x$misc$varinfo$lp,decreasing=F)[1:min(ncol(x$counts),n.odgenes)])])))
    }
    odgenes <- odgenes[names(odgenes) %in% Reduce(intersect,lapply(p2.objs,function(x) colnames(x$counts)))]
    odgenes <- names(odgenes)[1:min(length(odgenes),n.odgenes)]
    ## Common variance scaling
    if (var.scale) {
        cgsf <- do.call(cbind,lapply(p2.objs,function(x) x$misc$varinfo[odgenes,]$gsf))
        cgsf <- exp(rowMeans(log(cgsf)))
        names(cgsf) <- odgenes
    }
    ## Prepare the matrices
    cproj <- lapply(p2.objs,function(r) {
        x <- r$counts[,odgenes];
        if(var.scale) {
            x@x <- x@x*rep(cgsf,diff(x@p))
        }
        if(neighborhood.average) {
            ## use the averaged matrices
            xk <- r$misc$edgeMat$quickCPCA;
            x <- Matrix::t(xk) %*% x
        }
        x
    })
    list(genespace1=cproj[[1]], genespace2=cproj[[2]],cgsf=cgsf)
}

#' Perform pairwise JNMF
quickJNMF <- function(p2.objs = NULL, n.comps = 30, n.odgenes=NULL, var.scale=TRUE,
                      verbose =TRUE, max.iter=1000, neighborhood.average=FALSE) {
    ## Stop if more than 2 samples
    if (length(p2.objs) != 2) stop('quickJNMF only supports pairwise alignment');
    ## Get common set of genes
    if(is.null(n.odgenes)) {
        odgenes <- table(unlist(lapply(p2.objs,function(x) x$misc$odgenes)))
    } else {
        odgenes <- table(unlist(lapply(p2.objs,function(x) rownames(x$misc$varinfo)[(order(x$misc$varinfo$lp,decreasing=F)[1:min(ncol(x$counts),n.odgenes)])])))
    }
    odgenes <- odgenes[names(odgenes) %in% Reduce(intersect,lapply(p2.objs,function(x) colnames(x$counts)))]
    odgenes <- names(odgenes)[1:min(length(odgenes),n.odgenes)]
    ## Common variance scaling
    if (var.scale) {
        cgsf <- do.call(cbind,lapply(p2.objs,function(x) x$misc$varinfo[odgenes,]$gsf))
        cgsf <- exp(rowMeans(log(cgsf)))
        names(cgsf) <- odgenes
    }
    ## Prepare the matrices
    cproj <- lapply(p2.objs,function(r) {
        x <- r$counts[,odgenes];
        if(var.scale) {
            x@x <- x@x*rep(cgsf,diff(x@p))
        }
        if(neighborhood.average) {
            ## use the averaged matrices
            xk <- r$misc$edgeMat$quickCPCA;
            x <- Matrix::t(xk) %*% x
        }
        x
    })
    ## Convert to matrix
    cproj <- lapply(cproj, as.matrix)
    rjnmf.seed <- 12345
    ## Do JNMF
    z <- Rjnmf::Rjnmf(Xs=t(cproj[[1]]), Xu=t(cproj[[2]]), k=n.comps, alpha=0.5, lambda = 0.5, epsilon = 0.001,
                 maxiter= max.iter, verbose=F, seed=rjnmf.seed)
    rot1 <- cproj[[1]] %*% z$W
    rot2 <- cproj[[2]] %*% z$W
    ## return
    list(rot1=rot1, rot2=rot2,z=z,cgsf=cgsf)
}

#' Perform pairwise CPCA
quickCPCA <- function(p2.objs = NULL, n.comps = 30, n.odgenes = NULL, var.scale = T, verbose = T,
                      n.cores = 32, max.iter=1000, neighborhood.average = FALSE) {
    ## TODO: Check the input arguments
    ## Select a common set of genes
    if (is.null(n.odgenes)) {
        odgenes <- table(unlist(lapply(p2.objs, function(x) x$misc$odgenes)))
    } else {
        odgenes <- table(unlist(lapply(p2.objs,function(x) rownames(x$misc$varinfo)[(order(x$misc$varinfo$lp,decreasing=F)[1:min(ncol(x$counts),n.odgenes)])])))
    }
    odgenes <- odgenes[names(odgenes) %in% Reduce(intersect,lapply(p2.objs,function(x) colnames(x$counts)))]
    odgenes <- names(odgenes)[1:min(length(odgenes),n.odgenes)]
    ## Common variance scaling
    if (var.scale) {
        cgsf <- do.call(cbind,lapply(p2.objs,function(x) x$misc$varinfo[odgenes,]$gsf))
        cgsf <- exp(rowMeans(log(cgsf)))
        names(cgsf) <- odgenes
    }
    if (verbose) cat('calculating covariates for 2 datasets...');
    covl <- lapply(p2.objs, function(r) {
        x <- r$counts[,odgenes];
        if (var.scale) {
            x@x <- x@x * rep(cgsf, diff(x@p))
        }
        if (neighborhood.average) {
            xk <- r$misc$edgeMat$quickCPCA;
            x <- Matrix::t(xk) %*% x
        }
        sparse.cov(x)
    })
    if(verbose) cat('done\n');
    ## Perform CPCA
    if(verbose) cat('getting common PCs...');
    ncells <- unlist(lapply(p2.objs, function(x) nrow(x$counts)));
    xcp <- cpca::cpca(covl, ncells, ncomp=n.comps, maxit=max.iter)
    rownames(xcp$CPC) <- odgenes;
    if(verbose) cat('done\n');
    ## TODO: Check for convergence
    xcp
}


panelClust <- function(p2list,
                       var.scale = TRUE, n.odgenes = 1000,
                       reduction.method = 'CPCA', n.comps = 30, 
                       neighborhood.average = TRUE, neighborhood.average.k = 5,
                       k = 30, k.self = 10, k.self.weight = 0.01,
                       community.detection.method = "multilevel",
                       n.cores = 32, verbose=TRUE, xl=NULL) {
    ## TODO: Check that the p2list names are unique
    ## neighborhood averaging
    if (neighborhood.average) {
        if (verbose) cat("neighborhood averaging")
        p2list <- lapply(p2list, function(r) {
            ## TODO: check that the PCA reductions exist
            xk <- crossNN(mA = r$reductions$PCA, mB = r$reductions$PCA,
                          k = neighborhood.average.k, spaceType = 2,
                          lpSpaceP = 2.0, verbose = FALSE,
                          nThreads = n.cores)
            ## swap 0 to 1 scale, ensuring that nothing falls below 0
            xk@x <- pmax(1-xk@x,0)
            ## set diag
            diag(xk) <- 1;
            xk <- Matrix::t(Matrix::t(xk) / Matrix::colSums(xk))
            colnames(xk) <- colnames(xk) <- rownames(r$reductions$PCA)
            ## TODO: Change the name of this slot
            r$misc$edgeMat$quickCPCA <- xk;
            if (verbose) cat(".");
            r
        })
        if (verbose) cat(" done\n")
    }
    ## for every pair of samples
    cis <- combn(names(p2list),2)
    ## if pre-calculated matches are not provided run them
    if ( is.null(xl) ) {
        ## xl not provided
        cat('Calculating',ncol(cis),'pairwise', reduction.method, 'reductions...');
        ## TODO: n.cores is used twice here both in the outer loop and the call to the quick 
        xl <- BiocParallel::bplapply(1:ncol(cis), function(i) {
            if (reduction.method == 'CPCA') {
                cpca.max.iter <- 1000
                xcp <- quickCPCA(p2.objs = p2list[cis[,i]], n.comps = n.comps, n.odgenes = n.odgenes,
                                 var.scale = var.scale, verbose = FALSE,
                                 n.cores=n.cores, max.iter=cpca.max.iter,neighborhood.average=neighborhood.average);
            } else if (reduction.method == 'JNMF') {
                jnmf.max.iter <- 100
                xcp <- quickJNMF(p2.objs = p2list[cis[,i]], n.comps = n.comps, n.odgenes = n.odgenes,
                                 var.scale = var.scale, verbose = FALSE,
                                 max.iter = jnmf.max.iter, neighborhood.average=neighborhood.average);
            } else if (reduction.method == 'GeneSpace') {
                xcp <- quickNULL(p2.objs = p2list[cis[,i]], n.odgenes = n.odgenes,
                                 var.scale = var.scale, verbose = FALSE,
                                 neighborhood.average=neighborhood.average);
            } else {
                ## TODO: Add other methods
                stop(paste("unknown reduction method",reduction.method));
            }
            cat('.');
            return(xcp)
        }, BPPARAM = BiocParallel::MulticoreParam(workers = n.cores))
        names(xl) <- apply(cis,2,paste,collapse='.vs.')
        cat(' done\n')
    } else {
        ## precalculated xl has been provided
        ## TODO: check next few lines
        mi <- rep(NA,ncol(cis));
        nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        cat('matched',sum(!is.na(mi)),'out of',length(mi),'CPCA pairs ... ')
        ## if any rotations are missing run them here
        if(any(is.na(mi))) {
            cat('running',sum(is.na(mi)),'additional CPCAs ')
            xl2 <- BiocParallel::bplapply(which(is.na(mi)), function(i) {
                if(reduction.method=='CPCA') {
                    cpca.max.iter <- 1000
                    xcp <- quickCPCA(p2.objs = p2list[cis[,i]], n.comps = n.comps, n.odgenes = n.odgenes,
                                     var.scale = var.scale, verbose = FALSE,
                                     n.cores=n.cores, max.iter=cpca.max.iter, neighborhood.average);
                } else if (reduction.method=='JNMF') {
                    jnmf.max.iter <- 100
                    xcp <- quickJNMF(p2.objs = p2list[cis[,i]], n.comps = n.comps, n.odgenes = n.odgenes,
                                     var.scale = var.scale, verbose = FALSE,
                                     max.iter = jnmf.max.iter, neighborhood.average=neighborhood.average);
                } else if (reduction.method == 'GeneSpace') {
                    xcp <- quickNULL(p2.objs = p2list[cis[,i]], n.odgenes = n.odgenes,
                                     var.scale = var.scale, verbose = FALSE,
                                     neighborhood.average=neighborhood.average);
                } else {
                    ## TODO: Add other methods
                    stop(paste("unknown reduction method",reduction.method));
                }
                cat('.');
                return(xcp)
            }, BPPARAM = BiocParallel::MulticoreParam(workers = n.cores))
            names(xl2) <- apply(cis[,which(is.na(mi)),drop=F],2,paste,collapse='.vs.');
            ## append new matches to old 
            xl <- c(xl,xl2);
        }
        ## re-do the match and order
        mi <- rep(NA,ncol(cis));
        nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
        mi[which(!is.na(nm))] <- na.omit(nm);
        if(any(is.na(mi))) { stop("unable to get complete set of CPCA results") }
        xl <- xl[mi]
        cat(" done\n");
    } ## if is.null(xl)
    ## Get nearest neighbours between samples
    cat('Getting inter-sample nearest neighbors...')
    nnres <- lapply(1:ncol(cis), function(i) {
        r.ns <- p2list[cis[,i]]
        if (!is.null(xl[[i]]$rot1)) {
            ## JNMF
            n12 <- crossNN(xl[[i]]$rot1,xl[[i]]$rot2,k,2,2.0,FALSE,n.cores)
            n21 <- crossNN(xl[[i]]$rot2,xl[[i]]$rot1,k,2,2.0,FALSE,n.cores)
            ## mnn <- drop0(n21*t(n12))
            mnn <- as(n21*Matrix::t(n12),'dgTMatrix')
            ret.df <- data.frame(
                'mA.lab'=rownames(xl[[i]]$rot1)[mnn@i+1],
                'mB.lab'=rownames(xl[[i]]$rot2)[mnn@j+1],
                'w'=pmax(1-mnn@x,0),
                stringsAsFactors=F
            )
            ## TODO: recalculate distances
            ret.df
        } else if (!is.null(xl[[i]]$CPC)) {
            ## CPCA
            rot <- xl[[i]]$CPC;
            odgenes <- rownames(xl[[i]]$CPC);
            if (var.scale) {
                cgsf <- do.call(cbind,lapply(r.ns,function(x) x$misc$varinfo[odgenes,]$gsf));
                cgsf <- exp(rowMeans(log(cgsf)))
            }
            ## create matrices
            cproj <- lapply(r.ns,function(r) {
                x <- r$counts[,odgenes];
                if(var.scale) {
                    x@x <- x@x*rep(cgsf,diff(x@p))
                }
                if(neighborhood.average) {
                    xk <- r$misc$edgeMat$quickCPCA;
                    x <- Matrix::t(xk) %*% x
                }
                x
            })
            ## Perform centering
            ## TODO: make optional
            ncells <- unlist(lapply(cproj,nrow));
            centering <- colSums(do.call(rbind,lapply(cproj,Matrix::colMeans))*ncells)/sum(ncells)
            cpproj <- lapply(cproj,function(x) {
                x <- t(as.matrix(Matrix::t(x))-centering)
                x %*% rot;
            })
            n1 <- cis[1,i]; n2 <- cis[2,i]
            cat(".")
            ## get the nearest neighbours
            n12 <- crossNN(cpproj[[n1]],cpproj[[n2]],k,2,2.0,FALSE,n.cores)
            n21 <- crossNN(cpproj[[n2]],cpproj[[n1]],k,2,2.0,FALSE,n.cores)
            ## Calculate mutual neighbours
            mnn <- Matrix::drop0(n21*Matrix::t(n12))
            mnn <- as(n21*Matrix::t(n12),'dgTMatrix')
            ret.df <- data.frame(
                'mA.lab'=rownames(cpproj[[n1]])[mnn@i+1],
                'mB.lab'=rownames(cpproj[[n2]])[mnn@j+1],
                'w'=pmax(1-mnn@x,0),stringsAsFactors=F
            )
        } else if (!is.null(xl[[i]]$genespace1)) {
            ## Overdispersed Gene space
            is.null(xl[[i]]$genespace1)
            n12 <- crossNN(as.matrix(xl[[i]]$genespace1), as.matrix(xl[[i]]$genespace2),k,2,2.0,FALSE,n.cores)
            n21 <- crossNN(as.matrix(xl[[i]]$genespace2), as.matrix(xl[[i]]$genespace1),k,2,2.0,FALSE,n.cores)
            ##
            mnn <- Matrix::drop0(n21*Matrix::t(n12))
            mnn <- as(n21*Matrix::t(n12),'dgTMatrix')
            ## return
            ret.df <- data.frame(
                'mA.lab'=rownames(xl[[i]]$genespace1)[mnn@i+1],
                'mB.lab'=rownames(xl[[i]]$genespace2)[mnn@j+1],
                'w'=pmax(1-mnn@x,0),stringsAsFactors=F
            )
        } else {
            stop('unknown reduction provided');
        }
    }) ## lapply for nearest neighbours
    ## merge the knn networks
    el <- do.call(rbind, nnres)
    ## append some local edges
    if(k.self > 0) {
        cat('kNN pairs ')
        x <- data.frame(do.call(rbind,lapply(p2list, function(x) {
            xk <- crossNN(
                mA=x$reductions$PCA,
                mB=x$reductions$PCA,
                k = k.self,
                spaceType = 2,
                lpSpaceP = 2.0,
                verbose = FALSE,
                nThreads = n.cores)
            diag(xk) <- 0;
            xk <- as(xk,'dgTMatrix')
            cat(".")
            data.frame(
                'mA.lab'=rownames(x$reductions$PCA)[xk@i+1],
                'mB.lab'=rownames(x$reductions$PCA)[xk@j+1],
                'w'=pmax(1-xk@x,0),
                stringsAsFactors=F
            )
        })),stringsAsFactors = F)
        ## downweight the internal edges
        x$w <- x$w * k.self.weight
        cat(' done\n')
        ## merge internal and external edges
        el <- rbind(el,x)
    }
    ## Generate igraph object
    g <- igraph::graph_from_edgelist(as.matrix(el[,c(1,2)]),directed=FALSE)
    igraph::E(g)$weight <- el[,3]
    ## Do community detection
    ## TODO: add more methods
    cat('Detecting clusters...');
    if (community.detection.method == 'multilevel') {
        cls <- igraph::multilevel.community(g);
        cat('done\n')
    } else {
        stop(paste0('unknown community.detection.method',community.detection.method))
    }
    ## Extract cluster memberships
    cls.mem <- igraph::membership(cls)
    ## return
    list(xl=xl,el=el,cls=cls,cls.mem=cls.mem)
}
                       


###########################
## Tests

## Data for developmetn
datl.p2 <- readRDS('/home/barkasn/work/ninib/bmet.feb22/input_data/datl.p2.rds')


## Data subsets for tests
p2list <- datl.p2
p2list <- datl.p2[c( "BMM4-Whole","BMM6-Whole","BMM5-Whole")]
p2list <- datl.p2[grep('Whole',names(datl.p2))]
p2list <- datl.p2[grep('Niche',names(datl.p2))]

length(p2list)

pjc1 <- panelClust(p2list, reduction.method='CPCA')
pjc2 <- panelClust(p2list,reduction.method='JNMF')
pjc3 <- panelClust(p2list,reduction.method='GeneSpace')


nbHelpers::setDisplay('10.0')


X11()
par(mfrow=c(1,3))
lapply(names(p2list), function(n) {
    p2list[[n]]$plotEmbedding(type='PCA',embeddingType='tSNE',groups=pjc1$cls.mem,main=n)
})


X11()
par(mfrow=c(1,3))
lapply(names(p2list), function(n) {
    p2list[[n]]$plotEmbedding(type='PCA',embeddingType='tSNE',groups=pjc2$cls.mem,main=n)
})


X11()
par(mfrow=c(3,3))
lapply(names(p2list), function(n) {
    p2list[[n]]$plotEmbedding(type='PCA',embeddingType='tSNE',groups=pjc3$cls.mem,main=n,mark.clusters=TRUE)
})

##################################
## Benchmarks
##################################

#' Get % of clusters that are private to one sample
#' @param p2list list of pagoda2 objects on which the panelClust() was run
#' @param pjc result of panelClust()
#' @param priv.cutoff percent of total cells of a cluster that have to come from a single cluster
#' for it to be called private
getClusterPrivacy <- function(p2list, pjc, priv.cutoff= 0.99) {
    ## Get the clustering factor
    cl <- pjc$cls.mem
    ## Cell metadata
    meta <- do.call(rbind, lapply(names(p2list), function(n) {
        x <- p2list[[n]];
        data.frame(
            p2name = c(n),
            cellid = rownames(x$counts)
        )
    }))
    ## get sample / cluster counts
    meta$cl <- cl[meta$cellid]
    cl.sample.counts <- reshape2::acast(meta, p2name ~ cl, fun.aggregate=length,value.var='cl')
    ## Get clusters that are sample private
    private.clusters <- names(which(apply(sweep(cl.sample.counts, 2, apply(cl.sample.counts,2,sum), FUN='/') > priv.cutoff,2,sum) > 0))
    ## percent clusters that are private
    length(private.clusters) / length(unique(cl))
}

## Examples
getClusterPrivacy(p2list, pjc1)
getClusterPrivacy(p2list, pjc2)
getClusterPrivacy(p2list, pjc3)

namedNames <- function (g) {
    n <- names(g)
    names(n) <- n
    n
}



#' Evaluate consistency of cluster relationships
#' @description Using the clustering we are generating per-sample dendrograms
#' and we are examining their similarity between different samples
#' More information about similarity measures
#' https://www.rdocumentation.org/packages/dendextend/versions/1.8.0/topics/cor_cophenetic
#' https://www.rdocumentation.org/packages/dendextend/versions/1.8.0/topics/cor_bakers_gamma
#' @param p2list list of pagoda2 object
#' @param pjc a clustering factor
#' @return list of cophenetic and bakers_gama similarities of the dendrograms from each sample
getClusterRelationshipConsistency <- function(p2list, pjc) {
    hcs <- lapply(namedNames(p2list), function(n) {
        x <- p2list[[n]]
        app.cl <- pjc[names(pjc) %in% rownames(x$counts)]
        cpm <- sweep(rowsum(as.matrix(x$misc$rawCounts),
                            app.cl[rownames(x$misc$rawCounts)]),1, table(app.cl), FUN='/') * 1e6
        as.dendrogram(hclust(as.dist( 1 - cor(t(cpm)))))
    })
    ## Compare all dendrograms pairwise
    cis <- combn(names(hcs), 2)
    dend.comp <- lapply(1:ncol(cis), function(i) {
        s1 <- cis[1,i]
        s2 <- cis[2,i]
        dl1 <- dendextend::intersect_trees(hcs[[s1]],hcs[[s2]])
        list(
            cophenetic=dendextend::cor_cophenetic(dl1[[1]],dl1[[2]]),
            bakers_gamma=dendextend::cor_bakers_gamma(dl1[[1]],dl1[[2]])
        )
    })
    ## return mean and sd
    list(
        mean.cophenetic = mean(unlist(lapply(dend.comp, function(x) {x$cophenetic}))),
        sd.cophenetic = sd(unlist(lapply(dend.comp, function(x) {x$cophenetic}))),
        mean.bakers_gamma = mean(unlist(lapply(dend.comp, function(x) {x$bakers_gamma}))),
        sd.bakers_gamma = sd(unlist(lapply(dend.comp, function(x) {x$bakers_gamma})))
    )
}

## example
getClusterRelationshipConsistency(p2list, pjc3$cls.mem)

pjc <- pjc3

#' Evaluate how many clusters are global
#' @param p2list list of pagoda2 object on which clustering was generated
#' @param pjc the result of joint clustering
#' @param pc.samples.cutoff the percent of the number of the total samples that a cluster has to span to be considered global
#' @param min.cell.count.per.samples minimum number of cells of cluster in sample to be considered as represented in that sample
#' @return percent of clusters that are global given the above criteria
getPercentGlobalClusters <- function(p2list, pjc, pc.samples.cutoff = 0.9, min.cell.count.per.sample = 10) {
    ## get the cluster factor
    cl <- pjc$cls.mem
    ## get metadata table
    meta <- do.call(rbind, lapply(names(p2list), function(n) {
        x <- p2list[[n]];
        data.frame(
            p2name = c(n),
            cellid = rownames(x$counts)
        )
    }))
    ## get sample / cluster counts
    meta$cl <- cl[meta$cellid]
    cl.sample.counts <- reshape2::acast(meta, p2name ~ cl, fun.aggregate=length,value.var='cl')
    ## which clusters are global
    global.cluster <- apply(cl.sample.counts > min.cell.count.per.sample, 2, sum) >=  ceiling(nrow(cl.sample.counts) * pc.samples.cutoff)
    ## pc global clusters
    sum(global.cluster) / length(global.cluster)
}

