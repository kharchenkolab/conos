#' @useDynLib conos
#' @import Matrix
#' @import igraph
#' @import sccore
#' @importFrom dplyr %>%
#' @importFrom magrittr %<>%
#' @importFrom magrittr %$%
NULL

#' Wrapper to make is.label.fixed optional
smoothMatrixOnGraph <- function(edges, edge.weights, matrix, is.label.fixed=logical(), ...) {
  smooth_count_matrix(edges, edge.weights, matrix, is_label_fixed=is.label.fixed, ...)
}

scaledMatricesP2 <- function(p2.objs, data.type, od.genes, var.scale, neighborhood.average) {
  ## Common variance scaling
  if (var.scale) {
    # use geometric means for variance, as we're trying to focus on the common variance components
    cqv <- do.call(cbind, lapply(p2.objs,function(x) x$misc$varinfo[od.genes,]$qv)) %>%  log() %>% rowMeans() %>% exp()
  }
  ## Prepare the matrices
  cproj <- lapply(p2.objs,function(r) {
    if (data.type == 'counts') {
      x <- r$counts[,od.genes];
    } else if (data.type %in% names(r$reductions)){
      if (!all(od.genes %in% colnames(r$reductions[[data.type]]))) {
        stop("Reduction '", data.type, "' should have columns indexed by gene, with all overdispersed genes presented")
      }
      x <- r$reductions[[data.type]][,od.genes];
    } else {
      stop("No reduction named '", data.type, "' in pagoda")
    }

    if(var.scale) {
      cgsf <- sqrt(cqv/exp(r$misc$varinfo[od.genes,]$v))
      cgsf[is.na(cgsf) | !is.finite(cgsf)] <- 0;
      x@x <- x@x*rep(cgsf,diff(x@p))
    }
    if(neighborhood.average) {
      ## use the averaged matrices
      x <- Matrix::t(edgeMat(r)$mat) %*% x # TODO: looks like `$mat` doesn't exist anymore
    }
    x
  })

  return(cproj)
}

scaledMatricesSeurat <- function(so.objs, data.type, od.genes, var.scale, neighborhood.average) {
  if (var.scale) {
    warning("Seurat doesn't support variance scaling")
  }

  if (data.type == 'scaled') {
    x.data <- lapply(so.objs, function(so) t(so@scale.data)[,od.genes])
  } else if (data.type == 'counts') {
    x.data <- lapply(so.objs, function(so) t(so@data)[,od.genes])
  } else {

  res <- mapply(function(so, x) if(neighborhood.average) Matrix::t(edgeMat(so)$mat) %*% x else x,
                so.objs, x.data)

  return(res)
  }
}

scaledMatricesSeuratV3 <- function(so.objs, data.type, od.genes, var.scale, neighborhood.average) {
  checkSeuratV3()
  if (var.scale) {
    warning("Seurat doesn't support variance scaling")
  }
  slot <- switch(
    EXPR = data.type,
    'scaled' = 'scale.data',
    'counts' = 'data',
    stop("Unknown Seurat data type: ", data.type)
  )
  x.data <- lapply(
    X = so.objs,
    FUN = function(so) {
      return(t(x = Seurat::GetAssayData(object = so, slot = slot))[, od.genes])
    }
  )
  res <- mapply(
    FUN = function(so, x) {
      return(if (neighborhood.average) {
        Matrix::t(x = edgeMat(so)$mat) %*% x
      } else {
        x
      })
    },
    so.objs,
    x.data
  )
  return(res)
}

scaledMatrices <- function(samples, data.type, od.genes, var.scale, neighborhood.average) {
  if (class(samples[[1]]) == "Pagoda2") {
    return(scaledMatricesP2(samples, data.type = data.type, od.genes, var.scale, neighborhood.average))
  } else if (class(samples[[1]]) == "seurat") {
    return(scaledMatricesSeurat(samples, data.type = data.type, od.genes, var.scale, neighborhood.average))
  } else if (inherits(x = samples[[1]], what = 'Seurat')) {
    return(scaledMatricesSeuratV3(
      so.objs = samples,
      data.type = data.type,
      od.genes = od.genes,
      var.scale = var.scale,
      neighborhood.average = neighborhood.average
    ))
  }
  stop("Unknown class of sample: ", class(samples[[1]]))
}

commonOverdispersedGenes <- function(samples, n.odgenes, verbose) {
  od.genes <- sort(table(unlist(lapply(samples, getOverdispersedGenes, n.odgenes))),decreasing=T)
  common.genes <- Reduce(intersect, lapply(samples, getGenes));
  if(length(common.genes)==0) { warning(paste("samples",paste(names(samples),collapse=' and '),'do not share any common genes!')) }
  if(length(common.genes)<n.odgenes) { warning(paste("samples",paste(names(samples),collapse=' and '),'do not share enoguh common genes!')) }
  od.genes <- od.genes[names(od.genes) %in% common.genes]

  if(verbose) cat("using",length(od.genes),"od genes\n")

  return(names(od.genes)[1:min(length(od.genes),n.odgenes)])
}

quickNULL <- function(p2.objs, data.type='counts', n.odgenes = NULL, var.scale = T,
                      verbose = TRUE, neighborhood.average=FALSE) {
  if(length(p2.objs) != 2) stop('quickNULL only supports pairwise alignment');

  od.genes <- commonOverdispersedGenes(p2.objs, n.odgenes, verbose=verbose)
  if(length(od.genes)<5) return(NULL);

  cproj <- scaledMatrices(p2.objs, data.type=data.type, od.genes=od.genes, var.scale=var.scale,
                          neighborhood.average=neighborhood.average)

  return(list(genespace1=cproj[[1]], genespace2=cproj[[2]]))
}

#' Perform pairwise JNMF
quickJNMF <- function(p2.objs, data.type='counts', n.comps = 30, n.odgenes=NULL, var.scale=TRUE, verbose =TRUE, max.iter=1000, neighborhood.average=FALSE) {
  ## Stop if more than 2 samples
  if (length(p2.objs) != 2) stop('quickJNMF only supports pairwise alignment');

  od.genes <- commonOverdispersedGenes(p2.objs, n.odgenes, verbose=verbose)
  if(length(od.genes)<5) return(NULL);

  cproj <- scaledMatrices(p2.objs, data.type=data.type, od.genes=od.genes, var.scale=var.scale, neighborhood.average=neighborhood.average) %>%
    lapply(as.matrix)

  rjnmf.seed <- 12345
  ## Do JNMF
  z <- Rjnmf::Rjnmf(Xs=t(cproj[[1]]), Xu=t(cproj[[2]]), k=n.comps, alpha=0.5, lambda = 0.5, epsilon = 0.001, maxiter= max.iter, verbose=F, seed=rjnmf.seed)
  rot1 <- cproj[[1]] %*% z$W
  rot2 <- cproj[[2]] %*% z$W

  res <- list(rot1=rot1, rot2=rot2,z=z);


  return(res)
}

cpcaFast <- function(covl,ncells,ncomp=10,maxit=1000,tol=1e-6,use.irlba=TRUE,verbose=F) {
  ncomp <- min(c(nrow(covl)-1,ncol(covl)-1,ncomp));
  if(use.irlba) {
    # irlba initialization
    p <- nrow(covl[[1]]);
    S <- matrix(0, nrow = p, ncol = p)
    for(i in 1:length(covl)) {
      S <- S + (ncells[i] / sum(ncells)) * covl[[i]]
    }
    ncomp <- min(c(nrow(S)-1,ncol(S)-1,ncomp));
    ev <- irlba::irlba(S, ncomp, maxit=10000)
    cc <- abind::abind(covl,along=3)
    cpcaF(cc,ncells,ncomp,maxit,tol,eigenvR=ev$v,verbose)
  } else {
    cpcaF(cc,ncells,ncomp,maxit,tol,verbose=verbose)
  }
}

#' Perform cpca on two samples
#' @param r.n list of p2 objects
#' @param ncomps number of components to calculate (default=100)
#' @param n.odgenes number of overdispersed genes to take from each dataset
#' @param var.scale whether to scale variance (default=TRUE)
#' @param verbose whether to be verbose
#' @param neighborhood.average use neighborhood average values
#' @param n.cores number of cores to use
quickCPCA <- function(r.n,data.type='counts',ncomps=100,n.odgenes=NULL,var.scale=TRUE,verbose=TRUE,neighborhood.average=FALSE, score.component.variance=FALSE) {
  od.genes <- commonOverdispersedGenes(r.n, n.odgenes, verbose=verbose)
  if(length(od.genes)<5) return(NULL);

  ncomps <- min(ncomps, length(od.genes) - 1)

  if(verbose) cat('calculating covariances for',length(r.n),' datasets ...')

  ## # use internal C++ implementation
  ## sparse.cov <- function(x,cMeans=NULL){
  ##   if(is.null(cMeans)) {  cMeans <- Matrix::colMeans(x) }
  ##   covmat <- spcov(x,cMeans);
  ## }

  sm <- scaledMatrices(r.n, data.type=data.type, od.genes=od.genes, var.scale=var.scale, neighborhood.average=neighborhood.average)
  covl <- lapply(sm,function(x) spcov(as(x, "dgCMatrix"), Matrix::colMeans(x)))
  ## # centering
  ## if(common.centering) {
  ##   ncells <- unlist(lapply(covl,nrow));
  ##   centering <- colSums(do.call(rbind,lapply(covl,colMeans))*ncells)/sum(ncells)
  ## } else {
  ##   centering <- NULL;
  ## }

  ## covl <- lapply(covl,sparse.cov,cMeans=centering)

  if(verbose) cat(' done\n')

  #W: get counts
  ncells <- unlist(lapply(sm,nrow))
  if(verbose) cat('common PCs ...')
  #xcp <- cpca(covl,ncells,ncomp=ncomps)
  res <- cpcaFast(covl,ncells,ncomp=ncomps,verbose=verbose,maxit=500,tol=1e-5);
  #system.time(res <- cpca:::cpca_stepwise_base(covl,ncells,k=ncomps))
  #res <- cpc(abind(covl,along=3),k=ncomps)
  rownames(res$CPC) <- od.genes;
  if(score.component.variance) {
    v0 <- lapply(sm,function(x) sum(apply(x,2,var)))
    v1 <- lapply(1:length(sm),function(i) {
      x <- sm[[i]];
      cm <- Matrix::colMeans(x);
      rot <- as.matrix(t(t(x %*% res$CPC) - t(cm %*% res$CPC)))
      apply(rot,2,var)/v0[[i]]
    })
    # calculate projection
    res$nv <- v1;
  }
  if(verbose) cat(' done\n')
  return(res);
}

#' Use space of combined sample-specific PCAs as a space
#' @param r.n list of p2 objects
#' @param ncomps number of components to calculate (default=30)
#' @param n.odgenes number of overdispersed genes to take from each dataset
#' @param var.scale whether to scale variance (default=TRUE)
#' @param verbose whether to be verbose
#' @param cgsf an optional set of common genes to align on
#' @param neighborhood.average use neighborhood average values
#' @param n.cores number of cores to use
quickPlainPCA <- function(r.n,data.type='counts',ncomps=30,n.odgenes=NULL,var.scale=TRUE,verbose=TRUE,neighborhood.average=FALSE, score.component.variance=FALSE, n.cores=30) {
  od.genes <- commonOverdispersedGenes(r.n, n.odgenes, verbose=verbose)
  if(length(od.genes)<5) return(NULL);

  if(verbose) cat('calculating PCs for',length(r.n),' datasets ...')

  sm <- scaledMatrices(r.n, data.type=data.type, od.genes=od.genes, var.scale=var.scale, neighborhood.average=neighborhood.average);
  pcs <- lapply(sm, function(x) {
    cm <- Matrix::colMeans(x);
    ncomps <- min(c(nrow(cm)-1,ncol(cm)-1,round(ncomps/2)));
    res <- irlba::irlba(x, nv=ncomps, nu =0, center=cm, right_only = F, reorth = T);
    if(score.component.variance) {
      # calculate projection
      rot <- as.matrix(t(t(x %*% res$v) - as.vector(t(cm %*% res$v))))
      # note: this could be calculated a lot faster, but would need to account for the variable matrix format
      v0 <- apply(x,2,var)
      v1 <- apply(rot,2,var)/sum(v0)
      res$nv <- v1;
    }
    res
  })

  pcj <- do.call(cbind,lapply(pcs,function(x) x$v))
  # interleave the column order so that selecting top n components balances out the two datasets
  interleave <- function(n1,n2) { order(c((1:n1)-0.5,1:n2)) }
  ncols <- unlist(lapply(pcs,function(x) ncol(x$v)))
  pcj <- pcj[,interleave(ncols[1],ncols[2])]
  
  rownames(pcj) <- od.genes;
  res <- list(CPC=pcj);

  if(score.component.variance) {
    res$nv <- lapply(pcs,function(x) x$nv)
  }
  if(verbose) cat(' done\n')

  return(res);
}


#' Perform CCA (using PMA package or otherwise) on two samples
#' @param r.n list of p2 objects
#' @param ncomps number of components to calculate (default=100)
#' @param n.odgenes number of overdispersed genes to take from each dataset
#' @param var.scale whether to scale variance (default=TRUE)
#' @param verbose whether to be verbose
#' @param neighborhood.average use neighborhood average values
#' @param n.cores number of cores to use
quickCCA <- function(r.n,data.type='counts',ncomps=100,n.odgenes=NULL,var.scale=TRUE,verbose=TRUE,neighborhood.average=FALSE, PMA=FALSE, score.component.variance=FALSE) {

  od.genes <- commonOverdispersedGenes(r.n, n.odgenes, verbose=verbose)
  if(length(od.genes)<5) return(NULL);

  ncomps <- min(ncomps, length(od.genes) - 1)
  sm <- scaledMatrices(r.n, data.type=data.type, od.genes=od.genes, var.scale=var.scale, neighborhood.average=neighborhood.average)
  sm <- lapply(sm,function(m) m[rowSums(m)>0,])
  sm <- lapply(sm,scale,scale=F) # center
  if(PMA) {
    if (!requireNamespace("PMA", quietly=T))
      stop("You need to install package 'PMA' to use the PMA flag.")

    res <- PMA::CCA(t(sm[[1]]),t(sm[[2]]),K=ncomps,trace=FALSE,standardize=FALSE)
  } else {
    res <- irlba::irlba(sm[[1]] %*% t(sm[[2]]),ncomps)
  }
  rownames(res$u) <- rownames(sm[[1]])
  rownames(res$v) <- rownames(sm[[2]])

  res$ul <- t(sm[[1]]) %*% res$u # MASS::ginv(sm[[1]]) %*% res$u
  res$vl <- t(sm[[2]]) %*% res$v

  res$ul <- apply(res$ul,2,function(x) x/sqrt(sum(x*x)))
  res$vl <- apply(res$vl,2,function(x) x/sqrt(sum(x*x)))

  v0 <- lapply(sm,function(x) sum(apply(x,2,var)))
  res$nv <- list(apply(sm[[1]] %*% res$ul,2,var)/v0[[1]],
                 apply(sm[[2]] %*% res$vl,2,var)/v0[[2]]);
  names(res$nv) <- names(sm);

  #res$sm <- sm;
  # end DEBUG

  # adjust component weighting
  cw <- sqrt(res$d); cw <- cw/max(cw)
  res$u <- res$u %*% diag(cw)
  res$v <- res$v %*% diag(cw)

  return(res);
}


# dendrogram modification functions
#' Set dendrogram node width by breadth of the provided factor
#' @param d dendrogram
#' @param fac across cells
#' @param leafContent $leafContent output of greedy.modularity.cut() providing information about which cells map to which dendrogram leafs
#' @param min.width minimum line width
#' @param max.width maximum line width
#' @export
dendSetWidthByBreadth <- function(d,fac,leafContent,min.width=1,max.width=4) {
  cc2width <- function(cc) {
    ent <- entropy::entropy(cc[-1],method='MM',unit='log2')/log2(length(levels(fac)))
    min.width+ent*(max.width-min.width)
  }

  cbm <- function(d,fac) {
    if(is.leaf(d)) {
      lc <- fac[leafContent[[attr(d,'label')]]]
      cc <- c(sum(is.na(lc)),table(lc));
      lwd <- cc2width(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(lwd=lwd))
      attr(d,'cc') <- cc;
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm,fac=fac);
      attributes(d) <- oa;
      cc <- attr(d[[1]],'cc')+attr(d[[2]],'cc')
      lwd <- cc2width(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(lwd=lwd))
      attr(d,'cc') <- cc;
      return(d);
    }
  }
  cbm(d,fac);
}

#' Set dendrogram colors according to a 2- or 3-level factor mixture
#' @param d dendrogram
#' @param fac across cells
#' @param leafContent $leafContent output of greedy.modularity.cut() providing information about which cells map to which dendrogram leafs
#' @export
dendSetColorByMixture <- function(d,fac,leafContent) {
  fac <- as.factor(fac);
  if(length(levels(fac))>3) stop("factor with more than 3 levels are not supported")
  if(length(levels(fac))<2) stop("factor with less than 2 levels are not supported")

  cc2col <- function(cc,base=0.1) {
    if(sum(cc)==0) {
      cc <- rep(1,length(cc))
    } else {
      cc <- cc/sum(cc)
    }

    if(length(cc)==3) { # 2-color
      cv <- c(cc[2],0,cc[3])+base; cv <- cv/max(cv) * (1-base)
      #rgb(base+cc[2],base,base+cc[3],1)
      rgb(cv[1],cv[2],cv[3],1)
    } else if(length(cc)==4) { # 3-color
      cv <- c(cc[2],cc[3],cc[4])+base; cv <- cv/max(cv) * (1-base);
      #rgb(base+cc[2],base+cc[3],base+cc[4],1)
      rgb(cv[1],cv[2],cv[3],1)

    }
  }

  cbm <- function(d,fac) {
    if(is.leaf(d)) {
      lc <- fac[leafContent[[attr(d,'label')]]]
      cc <- c(sum(is.na(lc)),table(lc));
      col <- cc2col(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      attr(d,'cc') <- cc;
      return(d);
    } else {
      oa <- attributes(d);
      d <- lapply(d,cbm,fac=fac);
      attributes(d) <- oa;
      cc <- attr(d[[1]],'cc')+attr(d[[2]],'cc')
      col <- cc2col(cc)
      attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
      attr(d,'cc') <- cc;
      return(d);
    }
  }
  cbm(d,fac);
}

# other functions

# use mclapply if available, fall back on BiocParallel, but use regular
# lapply() when only one core is specified
papply <- function(...,n.cores=parallel::detectCores(), mc.preschedule=FALSE) {
  if(n.cores>1) {
    if(requireNamespace("parallel", quietly = TRUE)) {
      res <- parallel::mclapply(...,mc.cores=n.cores,mc.preschedule=mc.preschedule)
    }
    else if(requireNamespace("BiocParallel", quietly = TRUE)) {
      # It should never happen because parallel is specified in Imports
      res <- BiocParallel::bplapply(... , BPPARAM = BiocParallel::MulticoreParam(workers = n.cores))
    }
  } else {
    # fall back on lapply
    res <- lapply(...)
  }

  is.error <- (sapply(res, class) == "try-error")
  if (any(is.error)) {
    stop(paste("Errors in papply:", res[is.error]))
  }

  return(res)
}

##################################
## Benchmarks
##################################

#' Get percent of clusters that are private to one sample
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
            cellid = getCellNames(x)
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

sn <- function(x) { names(x) <- x; x }


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
    hcs <- lapply(sn(names(p2list)), function(n) {
        x <- p2list[[n]]
        app.cl <- pjc[names(pjc) %in% getCellNames(x)]
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
            cellid = getCellNames(x)
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


## helper function for breaking down a factor into a list
factorBreakdown <- function(f) {tapply(names(f),f, identity) }

#' Post process clusters generated with walktrap to control granularity
#' @param p2list list of pagoda2 objects
#' @param pjc joint clustering that was performed with walktrap
#' @param no.cl number of clusters to get from the walktrap dendrogram
#' @param size.cutoff cutoff below which to merge the clusters
#' @param n.cores number of cores to use
postProcessWalktrapClusters <- function(p2list, pjc, no.cl = 200, size.cutoff = 10, n.cores=4) {
    ##devel
    ## pjc <- pjc3
    ## no.cl <- 200
    ## size.cutoff <- 10
    ## n.cores <- 4
    ## rm(pjc, no.cl,size.cutoff, n.cores)
    ##
    global.cluster <- igraph::cut_at(cls, no=no.cl)
    names(global.cluster) <- names(igraph::membership(cls))
    ## identify clusters to merge
    fqs <- as.data.frame(table(global.cluster))
    cl.to.merge <- fqs[fqs$Freq < size.cutoff,]$global.cluster
    cl.to.keep <- fqs[fqs$Freq >= size.cutoff,]$global.cluster
    ## Memberships to keep
    global.cluster.filtered <- as.factor(global.cluster[global.cluster %in% cl.to.keep])
    ## Get new assignments for all the cells
    new.assign <- unlist(unname(papply(p2list, function(p2o) {
        try({
            ## get global cluster centroids for cells in this app
            global.cluster.filtered.bd <- factorBreakdown(global.cluster.filtered)

            # W: counts accessor
            global.cl.centers <- do.call(rbind, lapply(global.cluster.filtered.bd, function(cells) {
                cells <- cells[cells %in% getCellNames(p2o)]
                if (length(cells) > 1) {
                    Matrix::colSums(p2o$counts[cells,])
                } else {
                    NULL
                }
            }))
            ## cells to reassign in this app
            cells.reassign <- names(global.cluster[global.cluster %in% cl.to.merge])
            cells.reassign <- cells.reassign[cells.reassign %in% rownames(p2o$counts)]

            # W: counts accessor
            xcor <- cor(t(as.matrix(p2o$counts[cells.reassign,,drop=FALSE])), t(as.matrix(global.cl.centers)))
            ## Get new cluster assignments
            new.cluster.assign <- apply(xcor,1, function(x) {colnames(xcor)[which.max(x)]})
            new.cluster.assign
        })
    },n.cores=n.cores)))
    ## Merge
    x <- as.character(global.cluster.filtered)
    names(x) <- names(global.cluster.filtered)
    new.clusters <- as.factor(c(x,new.assign))
    new.clusters
}

#' Get top overdispersed genes across samples
#' @param samples list of pagoda2 objects
#' @param n.genes number of overdispersed genes to extract
getOdGenesUniformly <- function(samples, n.genes) {
  if (!("Pagoda2" %in% class(con$samples[[1]])))
    stop("This function is currently supported only for Pagoda2 objects")

  gene.info <- lapply(samples, function(s)
    tibble::rownames_to_column(s$misc[['varinfo']], "gene") %>%
      dplyr::mutate(rank=rank(lpa, ties.method="first"))
  )
  genes <- gene.info %>% dplyr::bind_rows() %>% dplyr::group_by(gene) %>%
    dplyr::summarise(rank=min(rank)) %>% dplyr::arrange(rank) %>% .$gene

  return(genes[1:min(n.genes, length(genes))])
}


projectSamplesOnGlobalAxes <- function(samples, cms.clust, data.type, neighborhood.average, verbose, n.cores) {
  if(verbose) cat('calculating global projections ');

  # calculate global eigenvectors

  gns <- Reduce(intersect,lapply(cms.clust,rownames))
  if(verbose) cat('.');
  if(length(gns) < length(cms.clust)) stop("insufficient number of common genes")
  tcc <- Reduce('+',lapply(cms.clust,function(x) x[gns,]))
  tcc <- t(tcc)/colSums(tcc)*1e6;
  gv <- apply(tcc,2,var);
  gns <- gns[is.finite(gv) & gv>0]
  tcc <- tcc[,gns,drop=F];

  if(verbose) cat('.');
  global.pca <- prcomp(log10(tcc+1),center=T,scale=T,retx=F)
  # project samples onto the global axes
  global.proj <- papply(samples,function(s) {
    smat <- as.matrix(scaledMatrices(list(s), data.type=data.type, od.genes=gns, var.scale=F, neighborhood.average=neighborhood.average)[[1]])
    if(verbose) cat('.')
    #smat <- as.matrix(conos:::getRawCountMatrix(s,transposed=TRUE)[,gns])
    #smat <- log10(smat/rowSums(smat)*1e3+1)
    smat <- scale(smat,scale=T,center=T); smat[is.nan(smat)] <- 0;
    sproj <- smat %*% global.pca$rotation
  },n.cores=n.cores)
  if(verbose) cat('. done\n');

  return(global.proj)
}

getDecoyProjections <- function(samples, samf, data.type, var.scale, cproj, neighborhood.average, base.groups, decoy.threshold, n.decoys) {
  cproj.decoys <- lapply(cproj, function(d) {
    tg <- tabulate(as.integer(base.groups[rownames(d)]),nbins=length(levels(base.groups)))
    nvi <- which(tg < decoy.threshold)
    if(length(nvi)>0) {
      # sample cells from other datasets
      decoy.cells <- names(base.groups)[unlist(lapply(nvi,function(i) {
        vc <- which(as.integer(base.groups)==i & (!samf[names(base.groups)] %in% names(cproj)))
        if(length(vc) > n.decoys) {
          vc <- sample(vc, n.decoys)
        }
      }))]
      if(length(decoy.cells)>0) {
        # get the matrices
        do.call(rbind,lapply(samples[unique(samf[decoy.cells])],function(s) {
          gn <- intersect(getGenes(s),colnames(d));
          m <- scaledMatrices(list(s),data.type=data.type, od.genes=gn, var.scale=var.scale, neighborhood.average=neighborhood.average)[[1]]
          m <- m[rownames(m) %in% decoy.cells,,drop=F]
          # append missing genes
          gd <- setdiff(colnames(d),gn)
          if(length(gd)>0) {
            m <- cbind(m,Matrix(0,nrow=nrow(m),ncol=length(gd),dimnames=list(rownames(m),gd),sparse=T))
            m <- m[,colnames(d),drop=F] # fix gene order
          }
        }))
      } else {
        # empty matrix
        Matrix(0,nrow=0,ncol=ncol(d),dimnames=list(c(),colnames(d)))
      }
    }
  })

  #if(verbose) cat(paste0("+",sum(unlist(lapply(cproj.decoys,nrow)))))

  return(cproj.decoys)
}

getLocalEdges <- function(samples, k.self, k.self.weight, metric, l2.sigma, verbose, n.cores) {
  if(verbose) cat('local pairs ')
  x <- data.frame(do.call(rbind, papply(samples, function(x) {
    pca <- getPca(x)
    if (is.null(pca)) {
      stop("PCA must be estimated for all samples")
    }

    xk <- n2Knn(pca, k.self + 1, 1, verbose=FALSE, indexType=metric) # +1 accounts for self-edges that will be removed in the next line
    diag(xk) <- 0; # no self-edges
    xk <- as(drop0(xk),'dgTMatrix')
    if(verbose) cat(".")

    data.frame(mA.lab=rownames(pca)[xk@i+1], mB.lab=rownames(pca)[xk@j+1],
               w=convertDistanceToSimilarity(xk@x, metric=metric, l2.sigma=l2.sigma), stringsAsFactors=F)
  }, n.cores=n.cores, mc.preschedule=TRUE)), stringsAsFactors=F)

  x$w <- x$w * k.self.weight
  x$type <- 0;
  if(verbose) cat(' done\n')

  return(x)
}


##' Find threshold of cluster detectability
##'
##' For a given clustering, walks the walktrap result tree to find
##' a subtree with max(min(sens,spec)) for each cluster, where sens is sensitivity, spec is specificity
##' @param res walktrap result object (igraph)
##' @param clusters cluster factor
##' @return a list of $thresholds - per cluster optimal detectability values, and $node - internal node id (merge row) where the optimum was found
##' @export
bestClusterThresholds <- function(res,clusters,clmerges=NULL) {
  clusters <- as.factor(clusters);
  # prepare cluster vectors
  cl <- as.integer(clusters[res$names]);
  clT <- tabulate(cl,nbins=length(levels(clusters)))
  # run
  res$merges <- igraph:::complete.dend(res,FALSE)
  #x <- conos:::findBestClusterThreshold(res$merges-1L,matrix(cl-1L,nrow=1),clT)
  if(is.null(clmerges)) {
    x <- conos:::treeJaccard(res$merges-1L,matrix(cl-1L,nrow=1),clT)
    names(x$threshold) <- levels(clusters);
  } else {
    x <- conos:::treeJaccard(res$merges-1L,matrix(cl-1L,nrow=1),clT,clmerges-1L)
  }
  x
}

##' Find threshold of cluster detectability in trees of clusters
##'
##' For a given clustering, walks the walktrap (of clusters) result tree to find
##' a subtree with max(min(sens,spec)) for each cluster, where sens is sensitivity, spec is specificity
##' @param res walktrap result object (igraph) where the nodes were clusters
##' @param leaf.factor a named factor describing cell assignments to the leaf nodes (in the same order as res$names)
##' @param clusters cluster factor
##' @return a list of $thresholds - per cluster optimal detectability values, and $node - internal node id (merge row) where the optimum was found
##' @export
bestClusterTreeThresholds <- function(res,leaf.factor,clusters,clmerges=NULL) {
  clusters <- as.factor(clusters);
  # prepare cluster vectors
  cl <- as.integer(clusters[names(leaf.factor)]);
  clT <- tabulate(cl,nbins=length(levels(clusters)))
  # prepare clusters matrix: cluster (rows) counts per leaf of the merge tree (column)
  mt <- table(cl,leaf.factor)
  # run
  merges <- igraph:::complete.dend(res,FALSE)
  #x <- conos:::findBestClusterThreshold(res$merges-1L,as.matrix(mt),clT)
  if(is.null(clmerges)) {
    x <- conos:::treeJaccard(res$merges-1L,as.matrix(mt),clT)
    names(x$threshold) <- levels(clusters);
  } else {
    x <- conos:::treeJaccard(res$merges-1L,as.matrix(mt),clT,clmerges-1L)
  }

  x
}


##' performs a greedy top-down selective cut to optmize modularity
##'
##' @param wt walktrap rsult
##' @param N number of top greedy splits to take
##' @param leaf.labels leaf sample label factor, for breadth calculations - must be a named factor containing all wt$names, or if wt$names is null, a factor listing cells in the same order as wt leafs
##' @param minsize minimum size of the branch (in number of leafs)
##' @param minbreadth minimum allowed breadth of a branch (measured as normalized entropy)
##' @param flat.cut whether to simply take a flat cut (i.e. follow provided tree; default=TRUE). Does no observe minsize/minbreadth restrictions
##' @return list(hclust - hclust structure of the derived tree, leafContent - binary matrix with rows corresponding to old leaves, columns to new ones, deltaM - modularity increments)
##' @export
greedyModularityCut <- function(wt,N,leaf.labels=NULL,minsize=0,minbreadth=0,flat.cut=TRUE) {
  # prepare labels
  nleafs <- nrow(wt$merges)+1;
  if(is.null(leaf.labels)) {
    ll <- integer(nleafs);
  } else {
    if(is.null(wt$names)) {
      # assume that leaf.labels are provided in the correct order
      if(length(leaf.labels)!=nleafs) stop("leaf.labels is of incorrct length and wt$names is NULL")
      ll <- as.integer(as.factor(leaf.labels))-1L;
    } else {
      if(!all(wt$names %in% names(leaf.labels))) { stop("leaf.labels do not cover all wt$names")}
      ll <- as.integer(as.factor(leaf.labels[wt$names]))-1L;
    }
  }
  x <- greedyModularityCutC(wt$merges-1L,-1*diff(wt$modularity),N,minsize,ll,minbreadth,flat.cut)
  if(length(x$splitsequence)<1) {
    stop("unable to make a single split using specified size/breadth restrictions")
  }
  # transfer cell names for the leaf content
  if(!is.null(wt$names)) { rownames(x$leafContent) <- wt$names; } else { rownames(x$leafContent) <- c(1:nrow(x$leafContent)) }
  m <- x$merges+1; nleafs <- nrow(m)+1; m[m<=nleafs] <- -1*m[m<=nleafs]; m[m>0] <- m[m>0]-nleafs;
  hc <- list(merge=m,height=1:nrow(m),labels=c(1:nleafs),order=c(1:nleafs)); class(hc) <- 'hclust'
  # fix the ordering so that edges don't intersects
  hc$order <- order.dendrogram(as.dendrogram(hc))
  leafContentCollapsed <- apply(x$leafContent,2,function(z)rownames(x$leafContent)[which(z>0)])
  clfac <- as.factor(apply(x$leafContent,1,which.max))
  return(list(hc=hc,groups=clfac,leafContentArray=x$leafContent,leafContent=leafContentCollapsed,deltaM=x$deltaM,breadth=as.vector(x$breadth),splits=x$splitsequence))
}

##' determine number of detectable clusters given a reference walktrap and a bunch of permuted walktraps
##'
##' @param refwt reference walktrap rsult
##' @param tests a list of permuted walktrap results
##' @param min.threshold min detectability threshold
##' @param min.size minimum cluster size (number of leafs)
##' @param average.thresholds report a single number of detectable clusters for averaged detected thresholds (a list of detected clusters for each element of the tests list is returned by default)
##' @return number of detectable stable clusters
##' @export
stableTreeClusters <- function(refwt,tests,min.threshold=0.8,min.size=10,n.cores=30,average.thresholds=FALSE) {
  # calculate detectability thresholds for each node against entire list of tests
  #i<- 0;
  refwt$merges <- igraph:::complete.dend(refwt,FALSE)
  for(i in 1:length(tests)) tests[[i]]$merges <- igraph:::complete.dend(tests[[i]],FALSE)
  thrs <- papply(tests,function(testwt) {
    #i<<- i+1; cat("i=",i,'\n');
    idmap <- match(refwt$names,testwt$names)-1L;
    idmap[is.na(idmap)] <- -1;
    x <- scoreTreeConsistency(testwt$merges-1L,refwt$merges-1L,idmap ,min.size)
    x$thresholds;
  },n.cores=n.cores)
  if(length(tests)==1) {
    x <- maxStableClusters(refwt$merges-1L,thrs[[1]],min.threshold,min.size)
    return(length(x$terminalnodes))
  } else {
    if(average.thresholds) {
      # calculate average detection threshold and
      x <- maxStableClusters(refwt$merges-1L,rowMeans(do.call(cbind,thrs)),min.threshold,min.size)
      return(length(x$terminalnodes))
    } else { # reporting the resulting numbers of clusters for each
      xl <- lapply(thrs,function(z) maxStableClusters(refwt$merges-1L,z,min.threshold,min.size))
      return(unlist(lapply(xl,function(x) length(x$terminalnodes))))
    }
  }
}

convertDistanceToSimilarity <- function(distances, metric, l2.sigma=1e5, cor.base=1) {
  if(metric=='angular') {
    return(pmax(0, cor.base - distances))
  }

  return(exp(-distances / l2.sigma))
}

getPcaBasedNeighborMatrix <- function(sample.pair, od.genes, rot, k, k1=k, data.type='counts', var.scale=T, neighborhood.average=F, common.centering=T,
                                      matching.method='mNN', metric='angular', l2.sigma=1e5, cor.base=1, subset.cells=NULL,
                                      base.groups=NULL, append.decoys=F, samples=NULL, samf=NULL, decoy.threshold=1, n.decoys=k*2, append.global.axes=T, global.proj=NULL) {
  # create matrices, adjust variance
  cproj <- scaledMatrices(sample.pair, data.type=data.type, od.genes=od.genes, var.scale=var.scale, neighborhood.average=neighborhood.average)

  # determine the centering
  if (common.centering) {
    ncells <- unlist(lapply(cproj, nrow));
    centering <- setNames(rep(colSums(do.call(rbind, lapply(cproj, colMeans)) * ncells) / sum(ncells), length(cproj)), names(cproj))
  } else {
    centering <- lapply(cproj,colMeans)
  }

  # append decoy cells if needed
  if(!is.null(base.groups) && append.decoys) {
    cproj.decoys <- getDecoyProjections(samples, samf, data.type, var.scale, cproj, neighborhood.average, base.groups, decoy.threshold, n.decoys)
    cproj <- lapply(sn(names(cproj)),function(n) rbind(cproj[[n]],cproj.decoys[[n]]))
  }

  cpproj <- lapply(sn(names(cproj)),function(n) {
    x <- cproj[[n]]
    x <- t(as.matrix(t(x))-centering[[n]])
    x %*% rot;
  })

  if(!is.null(base.groups) && append.global.axes) {
    #cpproj <- lapply(sn(names(cpproj)),function(n) cbind(cpproj[[n]],global.proj[[n]])) # case without decoys
    cpproj <- lapply(sn(names(cpproj)),function(n) {
      gm <- global.proj[[n]]
      if (append.decoys) {
        decoy.cells <- rownames(cproj.decoys[[n]])
        if (length(decoy.cells)>0) {
          gm <- rbind(gm, do.call(rbind, lapply(global.proj[unique(samf[decoy.cells])],
                                                function(m) m[rownames(m) %in% decoy.cells,,drop=F])))
        }
      }
      # append global axes
      cbind(cpproj[[n]],gm[rownames(cpproj[[n]]),])
    })
  }

  if (!is.null(subset.cells)) {
    cpproj <- lapply(cpproj, function(proj) proj[intersect(rownames(proj), subset.cells), ])
  }

  mnn <- getNeighborMatrix(cpproj[[names(sample.pair)[1]]], cpproj[[names(sample.pair)[2]]], k, k1=k1, matching=matching.method, metric=metric, l2.sigma=l2.sigma, cor.base=cor.base)

  if (!is.null(base.groups) && append.decoys) {
    # discard edges connecting to decoys
    decoy.cells <- unlist(lapply(cproj.decoys,rownames))
    mnn <- mnn[, !colnames(mnn) %in% decoy.cells, drop=F]
    mnn <- mnn[!rownames(mnn) %in% decoy.cells, , drop=F]
  }

  return(mnn)
}

##' Establish rough neighbor matching between samples given their projections in a common space
##'
##' @param p1 projection of sample 1
##' @param p2 projection of sample 2
##' @param k neighborhood radius
##' @param matching mNN (default) or NN
##' @param metric distance type (default: "angular", can also be 'L2')
##' @param l2.sigma L2 distances get transformed as exp(-d/sigma) using this value (default=1e5)
##' @param min.similarity minimal similarity between two cells, required to have an edge
##' @return matrix with the similarity (!) values corresponding to weight (1-d for angular, and exp(-d/l2.sigma) for L2)
getNeighborMatrix <- function(p1,p2,k,k1=k,matching='mNN',metric='angular',l2.sigma=1e5, cor.base=1, min.similarity=1e-5) {
  quiet.knn <- (k1 > k)
  if (is.null(p2)) {
    n12 <- n2CrossKnn(p1, p1,k1,1, FALSE, metric, quiet=quiet.knn)
    n21 <- n12
  } else {
    n12 <- n2CrossKnn(p1, p2, k1, 1, FALSE, metric, quiet=quiet.knn)
    n21 <- n2CrossKnn(p2, p1, k1, 1, FALSE, metric, quiet=quiet.knn)
  }


  n12@x <- convertDistanceToSimilarity(n12@x, metric=metric, l2.sigma=l2.sigma, cor.base=cor.base)
  n21@x <- convertDistanceToSimilarity(n21@x, metric=metric, l2.sigma=l2.sigma, cor.base=cor.base)

  if (matching=='NN') {
    adj.mtx <- drop0(n21+t(n12));
    adj.mtx@x <- adj.mtx@x/2;
  } else if (matching=='mNN') {
    adj.mtx <- drop0(n21*t(n12))
    adj.mtx@x <- sqrt(adj.mtx@x)
  } else {
    stop("Unrecognized type of NN matching:", matching)
  }

  rownames(adj.mtx) <- rownames(p1); colnames(adj.mtx) <- rownames(p2);
  adj.mtx@x[adj.mtx@x < min.similarity] <- 0
  adj.mtx <- drop0(adj.mtx);

  if(k1 > k) { # downsample edges
    adj.mtx <- reduceEdgesInGraphIteratively(adj.mtx,k)
  }

  return(as(drop0(adj.mtx),'dgTMatrix'))
}

# 1-step edge reduction
reduceEdgesInGraph <- function(adj.mtx,k,klow=k,preserve.order=TRUE) {
  if(preserve.order) { co <- colnames(adj.mtx); ro <- rownames(adj.mtx); }
  adj.mtx <- adj.mtx[,order(diff(adj.mtx@p),decreasing=T)]
  adj.mtx@x <- pareDownHubEdges(adj.mtx,tabulate(adj.mtx@i+1),k,klow)
  adj.mtx <- t(drop0(adj.mtx))
  adj.mtx <- adj.mtx[,order(diff(adj.mtx@p),decreasing=T)]
  adj.mtx@x <- pareDownHubEdges(adj.mtx,tabulate(adj.mtx@i+1),k,klow)
  adj.mtx <- t(drop0(adj.mtx));
  if(preserve.order) { adj.mtx <- adj.mtx[match(ro,rownames(adj.mtx)),match(co,colnames(adj.mtx))]; }

  return(adj.mtx)
}

# a simple multi-step strategy to smooth out remaining hubs
# max.kdiff gives approximate difference in the degree of the resulting nodes that is tolerable
reduceEdgesInGraphIteratively <- function(adj.mtx,k,preserve.order=TRUE,max.kdiff=5,n.steps=3) {
  cc <- diff(adj.mtx@p); rc <- tabulate(adj.mtx@i+1);
  maxd <- max(max(cc),max(rc));
  if(maxd<=k) return(adj.mtx); # nothing to be done - already below k
  klow <- max(min(k,3),k-max.kdiff); # allowable lower limit
  # set up a stepping strategy
  n.steps <- min(n.steps,round(maxd/k))
  if(n.steps>1) {
    ks <- round(exp(seq(log(maxd),log(k),length.out=n.steps+1))[-1])
  } else {
    ks <- c(k);
  }
  for(ki in ks) {
    adj.mtx <- reduceEdgesInGraph(adj.mtx,ki,preserve.order=preserve.order)
  }
  cc <- diff(adj.mtx@p); rc <- tabulate(adj.mtx@i+1);
  maxd <- max(max(cc),max(rc));
  if(maxd-k > max.kdiff) {
    # do a cleanup step
    adj.mtx <- reduceEdgesInGraph(adj.mtx,k,klow=klow,preserve.order=preserve.order)
  }

  return(adj.mtx)
}

adjustWeightsByCellBalancing <- function(adj.mtx, factor.per.cell, balance.weights, same.factor.downweight=1.0, n.iters=50, verbose=F) {
  adj.mtx %<>% .[colnames(.), colnames(.)] %>% as("dgTMatrix")
  factor.per.cell %<>% .[colnames(adj.mtx)] %>% as.factor() %>% droplevels()

  weights.adj <- adj.mtx@x

  if (balance.weights) {
    for (i in 0:(n.iters-1)) {
      factor.frac.per.cell <- getSumWeightMatrix(weights.adj, adj.mtx@i, adj.mtx@j, as.integer(factor.per.cell))
      w.dividers <- factor.frac.per.cell * rowSums(factor.frac.per.cell > 1e-10)
      weights.adj <- adjustWeightsByCellBalancingC(weights.adj, adj.mtx@i, adj.mtx@j, as.integer(factor.per.cell), w.dividers)

      if (verbose && i %% 10 == 0) {
        cat("Difference from balanced state:", sum(abs(w.dividers[w.dividers > 1e-10] - 1)), "\n")
      }
    }
  }

  if (abs(same.factor.downweight - 1) > 1e-5) {
    weights.adj[factor.per.cell[adj.mtx@i + 1] == factor.per.cell[adj.mtx@j + 1]] %<>% `*`(same.factor.downweight)
  }

  mtx.res <- adj.mtx
  mtx.res@x <- weights.adj

  return(mtx.res)
}

## Correct unloading of the library
.onUnload <- function (libpath) {
  library.dynam.unload("conos", libpath)
}


##' Scan joint graph modularity for a range of k (or k.self) values
##'
##' Builds graph with different values of k (or k.self if scan.k.self=TRUE), evaluating modularity of the resulting multilevel clustering
##' note: will run evaluations in parallel using con$n.cores (temporarily setting con$n.cores to 1 in the process)
##' @param con Conos object to test
##' @param min minimal value of k to test
##' @param max vlaue of k to test
##' @param by scan step (defaults to 1)
##' @param scan.k.self whether to test dependency on scan.k.self
##' @param ... other parameters will be passed to con$buildGraph()
##' @return a data frame with $k $m columns giving k and the corresponding modularity
##' @export
scanKModularity <- function(con, min=3, max=50, by=1, scan.k.self=FALSE, omit.internal.edges=TRUE, verbose=TRUE, plot=TRUE, ... ) {
  k.seq <- seq(min,max,by=by);
  n.cores <- con$n.cores;
  con$n.cores <- 1;
  if(verbose) cat(paste0(ifelse(scan.k.self,'k.self=(','k=('),min,', ',max,') ['))
  xl <- conos:::papply(k.seq,function(kv) {
    if(scan.k.self) {
      x <- con$buildGraph(k.self=kv, ..., verbose=FALSE)
    } else {
      x <- con$buildGraph(k=kv, ..., verbose=FALSE)
    }
    if(verbose) cat('.')
    if(omit.internal.edges) {
      x <- delete_edges(x,which(E(x)$type==0))
      #adj.mtx <- as_adj(x,attr='weight')
      #adj.mtx <- conos:::reduceEdgesInGraphIteratively(adj.mtx,kv)
      #adj.mtx <- drop0(adj.mtx*t(adj.mtx))
      #adj.mtx@x <- sqrt(adj.mtx@x)
      #x <- graph_from_adjacency_matrix(adj.mtx,mode = "undirected",weighted=TRUE)
    }
    xc <- multilevel.community(x)
    modularity(xc)
  },n.cores=30)
  if(verbose) cat(']\n')
  con$n.cores <- n.cores;

  k.sens <- data.frame(k=k.seq,m=as.numeric(unlist(xl)))
  if(plot) {
    ggplot2::ggplot(k.sens,aes(x=k,y=m))+theme_bw()+ggplot2::geom_point()+ggplot2::geom_smooth()+ggplot2::xlab('modularity')+ggplot2::ylab('k')
  }

  return(k.sens);
}

##' Merge into a common matrix, entering 0s for the missing ones
mergeCountMatrices <- function(cms, transposed=F) {
  extendMatrix <- function(mtx, col.names) {
    new.names <- setdiff(col.names, colnames(mtx))
    ext.mtx <- Matrix::Matrix(0, nrow=nrow(mtx), ncol=length(new.names), sparse=T) %>%
      as(class(mtx)) %>% `colnames<-`(new.names)
    return(cbind(mtx, ext.mtx)[,col.names])
  }

  if (!transposed) {
    cms %<>% lapply(Matrix::t)
  }

  gene.union <- lapply(cms, colnames) %>% Reduce(union, .)
  res <- lapply(cms, extendMatrix, gene.union) %>% Reduce(rbind, .)

  if (!transposed) {
    res %<>% Matrix::t()
  }

  return(res)
}

getSampleNamePerCell=function(samples) {
  cl <- lapply(samples, getCellNames)
  return(rep(names(cl), sapply(cl, length)) %>% stats::setNames(unlist(cl)) %>% as.factor())
}

#' Estimate labeling distribution for each vertex, based on provided labels using Random Walk
#' @param labels vector of factor or character labels, named by cell names
#' @param max.iters: maximal number of iterations. Default: 100.
#' @param tol: absolute tolerance as a stopping criteria. Default: 0.025
#' @param verbose: verbose mode. Default: TRUE.
#' @param fixed.initial.labels: prohibit changes of initial labels during diffusion. Default: TRUE.
propagateLabelsDiffusion <- function(graph, labels, max.iters=100, diffusion.fading=10.0, diffusion.fading.const=0.1, tol=0.025, verbose=TRUE, fixed.initial.labels=TRUE) {
  if (is.factor(labels)) {
    labels <- as.character(labels) %>% setNames(names(labels))
  }

  edges <- igraph::as_edgelist(graph)
  edge.weights <- igraph::edge.attributes(graph)$weight
  labels <- labels[intersect(names(labels), igraph::vertex.attributes(graph)$name)]
  label.distribution <- propagate_labels(edges, edge.weights, vert_labels=labels, max_n_iters=max.iters, verbose=verbose,
                                         diffusion_fading=diffusion.fading, diffusion_fading_const=diffusion.fading.const,
                                         tol=tol, fixed_initial_labels=fixed.initial.labels)
  return(label.distribution)
}

#' Propagate labels using Zhu, Ghahramani, Lafferty (2003) algorithm
#' http://mlg.eng.cam.ac.uk/zoubin/papers/zgl.pdf
#' TODO: change solver here for something designed for Laplacians. Need to look Daniel Spielman's research
propagateLabelsSolver <- function(graph, labels, solver="mumps") {
  if (!solver %in% c("mumps", "Matrix"))
    stop("Unknown solver: ", solver, ". Only 'mumps' and 'Matrix' are currently supported")

  if (!requireNamespace("rmumps", quietly=T)) {
    warning("Package 'rmumps' is required to use 'mumps' solver. Fall back to 'Matrix'")
    solver <- "Matrix"
  }

  adj.mat <- igraph::as_adjacency_matrix(graph, attr="weight")
  labeled.cbs <- intersect(colnames(adj.mat), names(labels))
  unlabeled.cbs <- setdiff(colnames(adj.mat), names(labels))

  labels <- as.factor(labels[labeled.cbs])

  weight.sum.mat <- Matrix::Diagonal(x=Matrix::colSums(adj.mat)) %>%
    `dimnames<-`(dimnames(adj.mat))

  laplasian.uu <- (weight.sum.mat[unlabeled.cbs, unlabeled.cbs] - adj.mat[unlabeled.cbs, unlabeled.cbs])

  type.scores <- Matrix::sparseMatrix(i=1:length(labels), j=as.integer(labels), x=1.0) %>%
    `colnames<-`(levels(labels)) %>% `rownames<-`(labeled.cbs)

  right.side <- Matrix::drop0(adj.mat[unlabeled.cbs, labeled.cbs] %*% type.scores)

  if (solver == "Matrix") {
    res <- Matrix::solve(laplasian.uu, right.side)
  } else {
    res <- rmumps::Rmumps$new(laplasian.uu, copy=F)$solve(right.side)
  }

  colnames(res) <- levels(labels)
  rownames(res) <- unlabeled.cbs
  return(rbind(res, type.scores))
}

#' Increase resolution for a specific set of clusters
#'
#' @param con conos object
#' @param target.clusters clusters for which the resolution should be increased
#' @param clustering name of clustering in the conos object to use. Either 'clustering' or 'groups' must be provided. Default: NULL
#' @param groups set of clusters to use. Ignored if 'clustering' is not NULL. Default: NULL
#' @param method function, used to find communities. Default: leiden.community
#' @param ... additional params passed to the community function
#' @export
findSubcommunities <- function(con, target.clusters, clustering=NULL, groups=NULL, method=leiden.community, ...) {
  groups <- parseCellGroups(con, clustering, groups)

  groups.raw <- as.character(groups) %>% setNames(names(groups))
  groups <- groups[intersect(names(groups), V(con$graph)$name)]

  if(length(groups) == 0) {
    stop("'groups' not defined for graph object.")
  }

  groups <- droplevels(as.factor(groups)[groups %in% target.clusters])
  if(length(groups) == 0) {
    stop("None of 'target.clusters' can be found in 'groups'.")
  }

  subgroups <- split(names(groups), groups)
  for (n in names(subgroups)) {
    if (length(subgroups[[n]]) < 2)
      next

    new.clusts <- method(induced_subgraph(con$graph, subgroups[[n]]), ...)
    groups.raw[new.clusts$names] <- paste0(n, "_", new.clusts$membership)
  }

  return(groups.raw)
}

parseCellGroups <- function(con, clustering, groups, parse.clusters=T) {
  if (!parse.clusters)
    return(groups)

  if (!is.null(groups)) {
    if (!any(names(groups) %in% names(con$getDatasetPerCell())))
      stop("'groups' aren't defined for any of the cells.")

    return(groups)
  }

  if (is.null(clustering)) {
    if (length(con$clusters) > 0)
      return(con$clusters[[1]]$groups)

    stop("Either 'groups' must be provided or the conos object must have some clustering estimated")
  }
  if(is.null(clusters[[clustering]]))
    stop(paste("clustering",clustering,"doesn't exist, run findCommunity() first"))

  return(con$clusters[[clustering]]$groups)
}

#' Estimate entropy of edge weights per cell according to the specified factor.
#' Can be used to visualize alignment quality according to this factor.
#'
#' @param con conos object
#' @param factor.per.cell some factor, which group cells, such as sample or a specific condition
#' @export
estimteWeightEntropyPerCell <- function(con, factor.per.cell) {
  adj.mat <- igraph::as_adjacency_matrix(con$graph, attr="weight") %>% as("dgTMatrix")
  factor.per.cell %<>% as.factor() %>% .[rownames(adj.mat)]
  weight.sum.per.fac.cell <- getSumWeightMatrix(adj.mat@x, adj.mat@i, adj.mat@j, as.integer(factor.per.cell)) %>%
    `colnames<-`(levels(adj.mat)) %>% `rownames<-`(rownames(adj.mat))

  xt <- table(factor.per.cell)

  entropy.per.cell <- apply(weight.sum.per.fac.cell, 1, entropy::KL.empirical, xt, unit=c('log2')) / log2(length(xt))

  return(entropy.per.cell)
}


