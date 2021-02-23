#' @useDynLib conos
#' @import Matrix
#' @import igraph
#' @import sccore
#' @importFrom dplyr %>%
#' @importFrom magrittr %<>%
#' @importFrom magrittr %$%
#' @importFrom rlang .data
NULL

## for magrittr and dplyr functions below
if (getRversion() >= "2.15.1"){
  utils::globalVariables(c(".", "x", "y"))
}


#' @keywords internal
scaledMatricesP2 <- function(p2.objs, data.type, od.genes, var.scale) {
  ## Common variance scaling
  if (var.scale) {
    # use geometric means for variance, as we're trying to focus on the common variance components
    cqv <- do.call(cbind, lapply(p2.objs,function(x) x$misc$varinfo[od.genes,]$qv)) %>%  log() %>% rowMeans() %>% exp()
  }
  ## Prepare the matrices
  cproj <- lapply(p2.objs, function(r) {
    if (data.type == 'counts') {
      x <- r$counts[,od.genes]
    } else if (data.type %in% names(r$reductions)){
      if (!all(od.genes %in% colnames(r$reductions[[data.type]]))) {
        stop("Reduction '", data.type, "' should have columns indexed by gene, with all overdispersed genes presented")
      }
      x <- r$reductions[[data.type]][,od.genes]
    } else {
      stop("No reduction named '", data.type, "' in pagoda2")
    }

    if(var.scale) {
      cgsf <- sqrt(cqv/exp(r$misc$varinfo[od.genes,]$v))
      cgsf[is.na(cgsf) | !is.finite(cgsf)] <- 0
      x@x <- x@x*rep(cgsf,diff(x@p))
    }
    return(x)
  })

  return(cproj)
}

#' @keywords internal
scaledMatricesSeurat <- function(so.objs, data.type, od.genes, var.scale) {
  if (var.scale) {
    warning("Seurat doesn't support variance scaling")
  }

  if (data.type == 'scaled') {
    x.data <- lapply(so.objs, function(so) t(so@scale.data)[,od.genes])
  } else if (data.type == 'counts') {
    x.data <- lapply(so.objs, function(so) t(so@data)[,od.genes])
  } else {

  res <- mapply(FUN = function(so, x) { return(x) }, so.objs, x.data )

  return(res)
  }
}

#' @keywords internal
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
  res <- mapply(FUN = function(so, x) { return(x) }, so.objs, x.data )

  return(res)
}

#' @keywords internal
scaledMatrices <- function(samples, data.type, od.genes, var.scale) {
  if ("Pagoda2" %in% class(samples[[1]])) {
    return(scaledMatricesP2(samples, data.type = data.type, od.genes, var.scale))
  } else if ("seurat" %in% class(samples[[1]])) {
    return(scaledMatricesSeurat(samples, data.type = data.type, od.genes, var.scale))
  } else if (inherits(x = samples[[1]], what = 'Seurat')) {
    return(scaledMatricesSeuratV3(
      so.objs = samples,
      data.type = data.type,
      od.genes = od.genes,
      var.scale = var.scale
    ))
  }
  stop("Unknown class of sample: ", class(samples[[1]]))
}

#' @keywords internal
commonOverdispersedGenes <- function(samples, n.odgenes, verbose) {
  od.genes <- sort(table(unlist(lapply(samples, getOverdispersedGenes, n.odgenes))),decreasing=TRUE)
  common.genes <- Reduce(intersect, lapply(samples, getGenes))
  if(length(common.genes)==0) { warning(paste("samples",paste(names(samples),collapse=' and '),'do not share any common genes!')) }
  if(length(common.genes)<n.odgenes) { warning(paste("samples",paste(names(samples),collapse=' and '),'do not share enough common genes!')) }
  od.genes <- od.genes[names(od.genes) %in% common.genes]

  if(verbose) message("using ",length(od.genes)," od genes\n")

  return(names(od.genes)[1:min(length(od.genes),n.odgenes)])
}

#' @keywords internal
quickNULL <- function(p2.objs, data.type='counts', n.odgenes = NULL, var.scale = T,
                      verbose = TRUE) {
  if (length(p2.objs) != 2){
    stop('quickNULL only supports pairwise alignment')
  }

  od.genes <- commonOverdispersedGenes(p2.objs, n.odgenes, verbose=verbose)
  if (length(od.genes)<5){
    return(NULL)
  }

  cproj <- scaledMatrices(p2.objs, data.type=data.type, od.genes=od.genes, var.scale=var.scale)

  return(list(genespace1=cproj[[1]], genespace2=cproj[[2]]))
}


## Performs Joint NMF for two matrices
#' @keywords internal
Rjnmf <- function(Xs, Xu, k, alpha, lambda, epsilon, maxiter, verbose, seed = 42) {
  set.seed(seed)

  ret <- RjnmfC(Xs, Xu, k, alpha, lambda, epsilon, maxiter, verbose)
  names(ret) <- c('Hs','Hu','W')
  ret
}

## Perform pairwise JNMF
#' @keywords internal
quickJNMF <- function(p2.objs, data.type='counts', n.comps=30, n.odgenes=NULL, var.scale=TRUE, verbose=TRUE, max.iter=1000, rjnmf.seed=12345) {

  ## Stop if more than 2 samples
  if (length(p2.objs) != 2){
    stop('quickJNMF only supports pairwise alignment')
  }

  od.genes <- commonOverdispersedGenes(p2.objs, n.odgenes, verbose=verbose)
  if (length(od.genes)<5){
    return(NULL)
  }

  cproj <- scaledMatrices(p2.objs, data.type=data.type, od.genes=od.genes, var.scale=var.scale) %>%
    lapply(as.matrix)

  ## Do JNMF
  z <- Rjnmf(Xs=t(cproj[[1]]), Xu=t(cproj[[2]]), k=n.comps, alpha=0.5, lambda = 0.5, epsilon = 0.001, maxiter= max.iter, verbose=FALSE, seed=rjnmf.seed)
  rot1 <- cproj[[1]] %*% z$W
  rot2 <- cproj[[2]] %*% z$W

  res <- list(rot1=rot1, rot2=rot2,z=z)

  return(res)
}

#' @keywords internal
cpcaFast <- function(covl,ncells,ncomp=10,maxit=1000,tol=1e-6,use.irlba=TRUE,verbose=FALSE) {
  ncomp <- min(c(nrow(covl)-1,ncol(covl)-1,ncomp))
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
#'
#' @param r.n list of pagoda2 objects
#' @param data.type character Type of data type in the input pagoda2 objects within r.n (default='counts')
#' @param ncomps numeric Number of components to calculate (default=100)
#' @param n.odgenes numeric Number of overdispersed genes to take from each dataset (default=NULL)
#' @param var.scale boolean Whether to scale variance (default=TRUE)
#' @param verbose boolean Whether to be verbose (default=TRUE)
#' @param score.component.variance boolean Whether to score component variance (default=FALSE)
#' @return cpca projection on two samples
#' @keywords internal
quickCPCA <- function(r.n, data.type='counts', ncomps=100, n.odgenes=NULL, var.scale=TRUE, verbose=TRUE, score.component.variance=FALSE) {
  od.genes <- commonOverdispersedGenes(r.n, n.odgenes, verbose=verbose)
  if(length(od.genes)<5){
    return(NULL)
  }

  ncomps <- min(ncomps, length(od.genes) - 1)

  if(verbose) message('calculating covariances for ',length(r.n),' datasets ...')

  ## # use internal C++ implementation
  ## sparse.cov <- function(x,cMeans=NULL){
  ##   if(is.null(cMeans)) {  cMeans <- Matrix::colMeans(x) }
  ##   covmat <- spcov(x,cMeans);
  ## }

  sm <- scaledMatrices(r.n, data.type=data.type, od.genes=od.genes, var.scale=var.scale)
  covl <- lapply(sm,function(x) spcov(as(x, "dgCMatrix"), Matrix::colMeans(x)))
  ## # centering
  ## if(common.centering) {
  ##   ncells <- unlist(lapply(covl,nrow));
  ##   centering <- colSums(do.call(rbind,lapply(covl,colMeans))*ncells)/sum(ncells)
  ## } else {
  ##   centering <- NULL;
  ## }

  ## covl <- lapply(covl,sparse.cov,cMeans=centering)

  if(verbose) message(' done\n')

  #W: get counts
  ncells <- unlist(lapply(sm,nrow))
  if(verbose) message('common PCs ...')
  #xcp <- cpca(covl,ncells,ncomp=ncomps)
  res <- cpcaFast(covl,ncells,ncomp=ncomps,verbose=verbose,maxit=500,tol=1e-5)
  #system.time(res <- cpca:::cpca_stepwise_base(covl,ncells,k=ncomps))
  #res <- cpc(abind(covl,along=3),k=ncomps)
  rownames(res$CPC) <- od.genes
  if (score.component.variance) {
    v0 <- lapply(sm,function(x) sum(apply(x,2,var)))
    v1 <- lapply(1:length(sm),function(i) {
      x <- sm[[i]];
      cm <- Matrix::colMeans(x);
      rot <- as.matrix(t(t(x %*% res$CPC) - t(cm %*% res$CPC)))
      apply(rot,2,var)/v0[[i]]
    })
    # calculate projection
    res$nv <- v1
  }
  if(verbose) message(' done\n')
  return(res)
}

#' Use space of combined sample-specific PCAs as a space
#'
#' @param r.n list of pagoda2 objects
#' @param data.type character Type of data type in the input pagoda2 objects within r.n (default='counts')
#' @param ncomps numeric Number of components to calculate (default=30)
#' @param n.odgenes numeric Number of overdispersed genes to take from each dataset (default=NULL)
#' @param var.scale boolean Whether to scale variance (default=TRUE)
#' @param verbose boolean Whether to be verbose (default=TRUE)
#' @param score.component.variance boolean Whether to score component variance (default=FALSE)
#' @param n.cores numeric Number of cores to use (default=1)
#' @return PCA projection, using space of combined sample-specific PCAs
#' @keywords internal
quickPlainPCA <- function(r.n, data.type='counts', ncomps=30, n.odgenes=NULL, var.scale=TRUE,
  verbose=TRUE, score.component.variance=FALSE, n.cores=1) {

  od.genes <- commonOverdispersedGenes(r.n, n.odgenes, verbose=verbose)

  if (length(od.genes)<5){
    return(NULL)
  }

  if(verbose) message('calculating PCs for ',length(r.n),' datasets ...')

  sm <- scaledMatrices(r.n, data.type=data.type, od.genes=od.genes, var.scale=var.scale)
  pcs <- lapply(sm, function(x) {
    cm <- Matrix::colMeans(x)
    ncomps <- min(c(nrow(cm)-1,ncol(cm)-1,round(ncomps/2)))
    res <- irlba::irlba(x, nv=ncomps, nu =0, center=cm, right_only = FALSE, reorth = TRUE)
    if(score.component.variance) {
      # calculate projection
      rot <- as.matrix(t(t(x %*% res$v) - as.vector(t(cm %*% res$v))))
      # note: this could be calculated a lot faster, but would need to account for the variable matrix format
      v0 <- apply(x,2,var)
      v1 <- apply(rot,2,var)/sum(v0)
      res$nv <- v1
    }
    res
  })

  pcj <- do.call(cbind,lapply(pcs,function(x) x$v))
  # interleave the column order so that selecting top n components balances out the two datasets
  interleave <- function(n1,n2) { order(c((1:n1)-0.5,1:n2)) }
  ncols <- unlist(lapply(pcs,function(x) ncol(x$v)))
  pcj <- pcj[,interleave(ncols[1],ncols[2])]

  rownames(pcj) <- od.genes
  res <- list(CPC=pcj)

  if(score.component.variance) {
    res$nv <- lapply(pcs,function(x) x$nv)
  }
  if(verbose) message(' done\n')

  return(res)
}


#' Perform CCA (using PMA package or otherwise) on two samples
#' @param r.n list of pagoda2 objects
#' @param data.type character Type of data type in the input pagoda2 objects within r.n (default='counts')
#' @param ncomps numeric Number of components to calculate (default=100)
#' @param n.odgenes numeric Number of overdispersed genes to take from each dataset
#' @param var.scale boolean Whether to scale variance (default=TRUE)
#' @param verbose boolean Whether to be verbose (default=TRUE)
#' @param score.component.variance boolean Whether to score component variance (default=FALSE)
#' @keywords internal
quickCCA <- function(r.n, data.type='counts', ncomps=100, n.odgenes=NULL, var.scale=TRUE, verbose=TRUE, PMA=FALSE, score.component.variance=FALSE) {

  od.genes <- commonOverdispersedGenes(r.n, n.odgenes, verbose=verbose)
  if(length(od.genes)<5){
    return(NULL)
  }

  ncomps <- min(ncomps, length(od.genes) - 1)
  sm <- scaledMatrices(r.n, data.type=data.type, od.genes=od.genes, var.scale=var.scale)
  sm <- lapply(sm,function(m) m[rowSums(m)>0,])
  sm <- lapply(sm,scale, scale=FALSE) # center
  if(PMA) {
    if (!requireNamespace("PMA", quietly=TRUE)){
      stop("You need to install package 'PMA' to use the PMA flag.")
    }
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
  res$nv <- list(apply(sm[[1]] %*% res$ul,2,var)/v0[[1]], apply(sm[[2]] %*% res$vl,2,var)/v0[[2]])
  names(res$nv) <- names(sm)

  #res$sm <- sm;
  # end DEBUG

  # adjust component weighting
  cw <- sqrt(res$d)
  cw <- cw/max(cw)
  res$u <- res$u %*% diag(cw)
  res$v <- res$v %*% diag(cw)

  return(res)
}


#### other functions

## complete.dend from igraph:::complete.dend
#' @keywords internal
complete.dend <- function(comm, use.modularity) {
  merges <- comm$merges
  if (nrow(merges) < comm$vcount-1) {
    if (use.modularity) {
      stop(paste("`use.modularity' requires a full dendrogram,",
                 "i.e. a connected graph"))
    }
    miss <- seq_len(comm$vcount + nrow(merges))[-as.vector(merges)]
    miss <- c(miss, seq_len(length(miss)-2) + comm$vcount+nrow(merges))
    miss <- matrix(miss, byrow=TRUE, ncol=2)
    merges <- rbind(merges, miss)
  }
  storage.mode(merges) <- "integer"

  merges
}

# use mclapply if available, fall back on BiocParallel, but use regular
# lapply() when only one core is specified
#' @keywords internal
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


#' @keywords internal
sn <- function(x) { names(x) <- x; x }


#' Evaluate consistency of cluster relationships
#' Using the clustering we are generating per-sample dendrograms
#' and we are examining their similarity between different samples
#' More information about similarity measures
#' <https://www.rdocumentation.org/packages/dendextend/versions/1.8.0/topics/cor_cophenetic>
#' <https://www.rdocumentation.org/packages/dendextend/versions/1.8.0/topics/cor_bakers_gamma>
#'
#' @param p2list list of pagoda2 object
#' @param pjc a clustering factor
#' @return list of cophenetic and bakers_gama similarities of the dendrograms from each sample
#' @keywords internal
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
#'
#' @param p2list list of pagoda2 object on which clustering was generated
#' @param pjc the result of joint clustering
#' @param pc.samples.cutoff numeric The percent of the number of the total samples that a cluster has to span to be considered global (default=0.9)
#' @param min.cell.count.per.samples numeric The minimum number of cells of cluster in sample to be considered as represented in that sample (default=10)
#' @return percent of clusters that are global given the above criteria
#' @keywords internal
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
#' @keywords internal
factorBreakdown <- function(f) {tapply(names(f),f, identity) }


#' Get top overdispersed genes across samples
#'
#' @param samples list of pagoda2 objects
#' @param n.genes number of overdispersed genes to extract
#' @keywords internal
getOdGenesUniformly <- function(samples, n.genes) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package \"tibble\" is needed for this function to work. Please install it.", call. = FALSE)
  }

  if (!("Pagoda2" %in% class(samples[[1]]))){
    stop("This function is currently supported only for Pagoda2 objects")
  }

  gene.info <- lapply(samples, function(s)
    tibble::rownames_to_column(s$misc[['varinfo']], "gene") %>%
      dplyr::mutate(rank=rank(.data$lpa, ties.method="first"))
  )
  genes <- gene.info %>% dplyr::bind_rows() %>% dplyr::group_by(.data$gene) %>%
    dplyr::summarise(rank=min(rank)) %>% dplyr::arrange(rank) %>% .$gene

  return(genes[1:min(n.genes, length(genes))])
}

#' @keywords internal
projectSamplesOnGlobalAxes <- function(samples, cms.clust, data.type, verbose, n.cores) {
  if(verbose) message('calculating global projections ')

  # calculate global eigenvectors

  gns <- Reduce(intersect,lapply(cms.clust,rownames))
  if(verbose) message('.')
  if(length(gns) < length(cms.clust)){
    stop("Insufficient number of common genes")
  }
  tcc <- Reduce('+',lapply(cms.clust,function(x) x[gns,]))
  tcc <- t(tcc)/colSums(tcc)*1e6
  gv <- apply(tcc,2,var)
  gns <- gns[is.finite(gv) & gv>0]
  tcc <- tcc[,gns,drop=FALSE]

  if(verbose) message('.')
  global.pca <- prcomp(log10(tcc+1),center=TRUE,scale=TRUE,retx=FALSE)
  # project samples onto the global axes
  global.proj <- papply(samples,function(s) {
    smat <- as.matrix(scaledMatrices(list(s), data.type=data.type, od.genes=gns, var.scale=FALSE)[[1]])
    if(verbose) message('.')
    #smat <- as.matrix(conos:::getRawCountMatrix(s,transposed=TRUE)[,gns])
    #smat <- log10(smat/rowSums(smat)*1e3+1)
    smat <- scale(smat,scale=TRUE,center=TRUE)
    smat[is.nan(smat)] <- 0
    sproj <- smat %*% global.pca$rotation
  },n.cores=n.cores)
  if(verbose) message('. done\n')

  return(global.proj)
}

#' @keywords internal
getDecoyProjections <- function(samples, samf, data.type, var.scale, cproj, base.groups, decoy.threshold, n.decoys) {
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
          gn <- intersect(getGenes(s),colnames(d))
          m <- scaledMatrices(list(s),data.type=data.type, od.genes=gn, var.scale=var.scale)[[1]]
          m <- m[rownames(m) %in% decoy.cells,,drop=FALSE]
          # append missing genes
          gd <- setdiff(colnames(d),gn)
          if(length(gd)>0) {
            m <- cbind(m,Matrix(0,nrow=nrow(m),ncol=length(gd),dimnames=list(rownames(m),gd),sparse=T))
            m <- m[,colnames(d),drop=FALSE] # fix gene order
          }
          m
        }))
      } else {
        # empty matrix
        Matrix(0,nrow=0,ncol=ncol(d),dimnames=list(c(),colnames(d)))
      }
    }
  })

  return(cproj.decoys)
}

#' @keywords internal
getLocalNeighbors <- function(samples, k.self, k.self.weight, metric, l2.sigma, verbose, n.cores) {
  if(verbose) message('local pairs ')
  x <- papply(samples, function(x) {
    pca <- getPca(x)
    if (is.null(pca)) {
      stop("PCA must be estimated for all samples")
    }

    xk <- N2R::Knn(pca, k.self + 1, 1, verbose=FALSE, indexType=metric) # +1 accounts for self-edges that will be removed in the next line
    diag(xk) <- 0 # no self-edges
    xk <- as(drop0(xk),'dgTMatrix')
    xk@x <- convertDistanceToSimilarity(xk@x, metric=metric, l2.sigma=l2.sigma) * k.self.weight
    rownames(xk) <- colnames(xk) <- rownames(pca)
    if(verbose) message(".")
    xk
  }, n.cores=n.cores, mc.preschedule=TRUE)
  if(verbose) message(' done\n')
  x
}

#' @keywords internal
getLocalEdges <- function(local.neighbors) {
  do.call(rbind,lapply(local.neighbors,function(xk) {
    data.frame(mA.lab=rownames(xk)[xk@i+1], mB.lab=colnames(xk)[xk@j+1],w=xk@x, type=0, stringsAsFactors=FALSE)
  }))
}



#' Find threshold of cluster detectability
#'
#' For a given clustering, walks the walktrap result tree to find
#' a subtree with max(min(sens,spec)) for each cluster, where sens is sensitivity, spec is specificity
#'
#' @param res walktrap result object (igraph)
#' @param clusters cluster factor
#' @param clmerges integer matrix of cluster merges (default=NULL). If NULL, the function treeJaccard() performs calculation without it.
#' @return a list of $thresholds - per cluster optimal detectability values, and $node - internal node id (merge row) where the optimum was found
#' @export
bestClusterThresholds <- function(res, clusters, clmerges=NULL) {
  clusters <- as.factor(clusters)
  # prepare cluster vectors
  cl <- as.integer(clusters[res$names])
  clT <- tabulate(cl,nbins=length(levels(clusters)))
  # run
  res$merges <- complete.dend(res,FALSE)
  #x <- conos:::findBestClusterThreshold(res$merges-1L,matrix(cl-1L,nrow=1),clT)
  if(is.null(clmerges)) {
    x <- treeJaccard(res$merges-1L,matrix(cl-1L,nrow=1),clT)
    names(x$threshold) <- levels(clusters)
  } else {
    x <- treeJaccard(res$merges-1L,matrix(cl-1L,nrow=1),clT,clmerges-1L)
  }
  x
}

#' Find threshold of cluster detectability in trees of clusters
#'
#' For a given clustering, walks the walktrap (of clusters) result tree to find
#' a subtree with max(min(sens,spec)) for each cluster, where sens is sensitivity, spec is specificity
#'
#' @param res walktrap result object (igraph) where the nodes were clusters
#' @param leaf.factor a named factor describing cell assignments to the leaf nodes (in the same order as res$names)
#' @param clusters cluster factor
#' @param clmerges integer matrix of cluster merges (default=NULL). If NULL, the function treeJaccard() performs calculation without it.
#' @return a list of $thresholds - per cluster optimal detectability values, and $node - internal node id (merge row) where the optimum was found
#' @export
bestClusterTreeThresholds <- function(res, leaf.factor, clusters, clmerges=NULL) {
  clusters <- as.factor(clusters)
  # prepare cluster vectors
  cl <- as.integer(clusters[names(leaf.factor)])
  clT <- tabulate(cl,nbins=length(levels(clusters)))
  # prepare clusters matrix: cluster (rows) counts per leaf of the merge tree (column)
  mt <- table(cl,leaf.factor)
  # run
  merges <- complete.dend(res,FALSE)
  #x <- conos:::findBestClusterThreshold(res$merges-1L,as.matrix(mt),clT)
  if(is.null(clmerges)) {
    x <- treeJaccard(res$merges-1L,as.matrix(mt),clT)
    names(x$threshold) <- levels(clusters)
  } else {
    x <- treeJaccard(res$merges-1L,as.matrix(mt),clT,clmerges-1L)
  }

  invisible(x)
}


#' Performs a greedy top-down selective cut to optmize modularity
#'
#' @param wt walktrap result
#' @param N numeric Number of top greedy splits to take
#' @param leaf.labels leaf sample label factor, for breadth calculations - must be a named factor containing all wt$names, or if wt$names is null, a factor listing cells in the same order as wt leafs (default=NULL)
#' @param minsize numeric Minimum size of the branch (in number of leafs) (default=0)
#' @param minbreadth numeric Minimum allowed breadth of a branch (measured as normalized entropy) (default=0)
#' @param flat.cut boolean Whether to simply take a flat cut (i.e. follow provided tree; default=TRUE). Does no observe minsize/minbreadth restrictions
#' @return list(hclust - hclust structure of the derived tree, leafContent - binary matrix with rows corresponding to old leaves, columns to new ones, deltaM - modularity increments)
#' @export
greedyModularityCut <- function(wt, N, leaf.labels=NULL, minsize=0, minbreadth=0, flat.cut=TRUE) {
  # prepare labels
  nleafs <- nrow(wt$merges)+1
  if(is.null(leaf.labels)) {
    ll <- integer(nleafs)
  } else {
    if(is.null(wt$names)) {
      # assume that leaf.labels are provided in the correct order
      if(length(leaf.labels)!=nleafs) stop("leaf.labels is of incorrct length and wt$names is NULL")
      ll <- as.integer(as.factor(leaf.labels))-1L
    } else {
      if(!all(wt$names %in% names(leaf.labels))) { stop("leaf.labels do not cover all wt$names")}
      ll <- as.integer(as.factor(leaf.labels[wt$names]))-1L
    }
  }
  x <- greedyModularityCutC(wt$merges-1L,-1*diff(wt$modularity),N,minsize,ll,minbreadth,flat.cut)
  if (length(x$splitsequence)<1) {
    stop("unable to make a single split using specified size/breadth restrictions")
  }
  # transfer cell names for the leaf content
  if (!is.null(wt$names)) {
    rownames(x$leafContent) <- wt$names
  } else {
    rownames(x$leafContent) <- c(1:nrow(x$leafContent))
  }
  m <- x$merges+1
  nleafs <- nrow(m)+1
  m[m<=nleafs] <- -1*m[m<=nleafs]
  m[m>0] <- m[m>0]-nleafs
  hc <- list(merge=m,height=1:nrow(m),labels=c(1:nleafs),order=c(1:nleafs))
  class(hc) <- 'hclust'
  # fix the ordering so that edges don't intersects
  hc$order <- order.dendrogram(as.dendrogram(hc))
  leafContentCollapsed <- apply(x$leafContent,2,function(z)rownames(x$leafContent)[which(z>0)])
  clfac <- as.factor(apply(x$leafContent,1,which.max))
  return(list(hc=hc,groups=clfac,leafContentArray=x$leafContent,leafContent=leafContentCollapsed,deltaM=x$deltaM,breadth=as.vector(x$breadth),splits=x$splitsequence))
}

#' Determine number of detectable clusters given a reference walktrap and a bunch of permuted walktraps
#'
#' @param refwt reference walktrap result
#' @param tests a list of permuted walktrap results
#' @param min.threshold numeric Min detectability threshold (default=0.8)
#' @param min.size numeric Minimum cluster size (number of leafs) (default=10)
#' @param n.cores numeric Number of cores (default=30)
#' @param average.thresholds boolean Report a single number of detectable clusters for averaged detected thresholds (default=FALSE) (a list of detected clusters for each element of the tests list is returned by default)
#' @return number of detectable stable clusters
#' @export
stableTreeClusters <- function(refwt, tests, min.threshold=0.8, min.size=10, n.cores=30, average.thresholds=FALSE) {
  # calculate detectability thresholds for each node against entire list of tests
  #i<- 0;
  refwt$merges <- complete.dend(refwt,FALSE)
  for (i in 1:length(tests)){
    tests[[i]]$merges <- complete.dend(tests[[i]],FALSE)
  }
  thrs <- papply(tests,function(testwt) {
    #i<<- i+1; cat("i=",i,'\n');
    idmap <- match(refwt$names,testwt$names)-1L
    idmap[is.na(idmap)] <- -1
    x <- scoreTreeConsistency(testwt$merges-1L,refwt$merges-1L,idmap ,min.size)
    x$thresholds
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

#' @keywords internal
convertDistanceToSimilarity <- function(distances, metric, l2.sigma=1e5, cor.base=1) {
  if (metric=='angular') {
    return(pmax(0, cor.base - distances))
  }

  return(exp(-distances / l2.sigma))
}

#' @keywords internal
getPcaBasedNeighborMatrix <- function(sample.pair, od.genes, rot, k, k1=k, data.type='counts', var.scale=TRUE, common.centering=TRUE,
                                      matching.method='mNN', metric='angular', l2.sigma=1e5, cor.base=1, subset.cells=NULL,
                                      base.groups=NULL, append.decoys=FALSE, samples=NULL, samf=NULL, decoy.threshold=1, n.decoys=k*2,
                                      append.global.axes=TRUE, global.proj=NULL) {
  # create matrices, adjust variance
  cproj <- scaledMatrices(sample.pair, data.type=data.type, od.genes=od.genes, var.scale=var.scale)

  # determine the centering
  if (common.centering) {
    ncells <- unlist(lapply(cproj, nrow));
    centering <- setNames(rep(colSums(do.call(rbind, lapply(cproj, colMeans)) * ncells) / sum(ncells), length(cproj)), names(cproj))
  } else {
    centering <- lapply(cproj,colMeans)
  }

  # append decoy cells if needed
  if(!is.null(base.groups) && append.decoys) {
    cproj.decoys <- getDecoyProjections(samples, samf, data.type, var.scale, cproj, base.groups, decoy.threshold, n.decoys)
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
                                                function(m) m[rownames(m) %in% decoy.cells,,drop=FALSE])))
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
    mnn <- mnn[, !colnames(mnn) %in% decoy.cells, drop=FALSE]
    mnn <- mnn[!rownames(mnn) %in% decoy.cells, , drop=FALSE]
  }

  return(mnn)
}

#' Establish rough neighbor matching between samples given their projections in a common space
#'
#' @param p1 projection of sample 1
#' @param p2 projection of sample 2
#' @param k neighborhood radius
#' @param k1 neighborhood radius
#' @param matching string mNN (default) or NN (default='mNN')
#' @param metric string Distance type (default: "angular", can also be 'L2')
#' @param l2.sigma numeric L2 distances get transformed as exp(-d/sigma) using this value (default=1e5)
#' @param cor.base numeric (default=1)
#' @param min.similarity minimal similarity between two cells, required to have an edge
#' @return matrix with the similarity (!) values corresponding to weight (1-d for angular, and exp(-d/l2.sigma) for L2)
#' @keywords internal
getNeighborMatrix <- function(p1, p2, k, k1=k, matching='mNN', metric='angular', l2.sigma=1e5, cor.base=1, min.similarity=1e-5) {
  quiet.knn <- (k1 > k)
  if (is.null(p2)) {
    n12 <- N2R::crossKnn(p1, p1,k1,1, FALSE, metric, quiet=quiet.knn)
    n21 <- n12
  } else {
    n12 <- N2R::crossKnn(p1, p2, k1, 1, FALSE, metric, quiet=quiet.knn)
    n21 <- N2R::crossKnn(p2, p1, k1, 1, FALSE, metric, quiet=quiet.knn)
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

  if (k1 > k) { # downsample edges
    adj.mtx <- reduceEdgesInGraphIteratively(adj.mtx,k)
  }

  return(as(drop0(adj.mtx),'dgTMatrix'))
}

## 1-step edge reduction
#' @keywords internal
reduceEdgesInGraph <- function(adj.mtx,k,klow=k,preserve.order=TRUE) {
  if(preserve.order) { co <- colnames(adj.mtx); ro <- rownames(adj.mtx); }
  adj.mtx <- adj.mtx[,order(diff(adj.mtx@p),decreasing=TRUE)]
  adj.mtx@x <- pareDownHubEdges(adj.mtx,tabulate(adj.mtx@i+1),k,klow)
  adj.mtx <- t(drop0(adj.mtx))
  adj.mtx <- adj.mtx[,order(diff(adj.mtx@p),decreasing=TRUE)]
  adj.mtx@x <- pareDownHubEdges(adj.mtx,tabulate(adj.mtx@i+1),k,klow)
  adj.mtx <- t(drop0(adj.mtx));
  if(preserve.order) { adj.mtx <- adj.mtx[match(ro,rownames(adj.mtx)),match(co,colnames(adj.mtx))]; }

  return(adj.mtx)
}

## a simple multi-step strategy to smooth out remaining hubs
## max.kdiff gives approximate difference in the degree of the resulting nodes that is tolerable
#' @keywords internal
reduceEdgesInGraphIteratively <- function(adj.mtx,k,preserve.order=TRUE,max.kdiff=5,n.steps=3) {
  cc <- diff(adj.mtx@p)
  rc <- tabulate(adj.mtx@i+1)
  maxd <- max(max(cc),max(rc))
  if (maxd<=k){
    return(adj.mtx) # nothing to be done - already below k
  }
  klow <- max(min(k,3),k-max.kdiff) # allowable lower limit
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
  cc <- diff(adj.mtx@p)
  rc <- tabulate(adj.mtx@i+1)
  maxd <- max(max(cc),max(rc))
  if(maxd-k > max.kdiff) {
    # do a cleanup step
    adj.mtx <- reduceEdgesInGraph(adj.mtx,k,klow=klow,preserve.order=preserve.order)
  }

  return(adj.mtx)
}

#' @keywords internal
adjustWeightsByCellBalancing <- function(adj.mtx, factor.per.cell, balance.weights, same.factor.downweight=1.0, n.iters=50, verbose=FALSE) {
  adj.mtx %<>% .[colnames(.), colnames(.)] %>% as("dgTMatrix")
  factor.per.cell %<>% .[colnames(adj.mtx)] %>% as.factor() %>% droplevels()

  weights.adj <- adj.mtx@x

  if (balance.weights) {
    for (i in 0:(n.iters-1)) {
      factor.frac.per.cell <- getSumWeightMatrix(weights.adj, adj.mtx@i, adj.mtx@j, as.integer(factor.per.cell))
      w.dividers <- factor.frac.per.cell * rowSums(factor.frac.per.cell > 1e-10)
      weights.adj <- adjustWeightsByCellBalancingC(weights.adj, adj.mtx@i, adj.mtx@j, as.integer(factor.per.cell), w.dividers)

      if (verbose && i %% 10 == 0) {
        message("Difference from balanced state: ", sum(abs(w.dividers[w.dividers > 1e-10] - 1)), "\n")
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
#' @keywords internal
.onUnload <- function (libpath) {
  library.dynam.unload("conos", libpath)
}


#' Scan joint graph modularity for a range of k (or k.self) values
#' Builds graph with different values of k (or k.self if scan.k.self=TRUE), evaluating modularity of the resulting multilevel clustering
#' NOTE: will run evaluations in parallel using con$n.cores (temporarily setting con$n.cores to 1 in the process)
#'
#' @param con Conos object to test
#' @param min numeric Minimal value of k to test (default=3)
#' @param max numeric Value of k to test (default=50)
#' @param by numeric Scan step (default=1)
#' @param scan.k.self boolean Whether to test dependency on scan.k.self (default=FALSE)
#' @param omit.internal.edges boolean Whether to omit internal edges of the graph (default=TRUE)
#' @param verbose boolean Whether to provide verbose output (default=TRUE)
#' @param plot boolean Whether to plot the output (default=TRUE)
#' @param ... other parameters will be passed to con$buildGraph()
#' @return a data frame with $k $m columns giving k and the corresponding modularity
#' @examples
#' \donttest{
#' library(pagoda2)
#' panel.preprocessed <- lapply(conosPanel::panel, basicP2proc, n.cores=1, min.cells.per.gene=0,
#'     n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)
#' con <- Conos$new(panel.preprocessed, n.cores=1)
#' scanKModularity(con)
#' }
#'
#' @export
scanKModularity <- function(con, min=3, max=50, by=1, scan.k.self=FALSE, omit.internal.edges=TRUE, verbose=TRUE, plot=TRUE, ... ) {
  k.seq <- seq(min,max, by=by)
  ncores <- con$n.cores
  con$n.cores <- 1
  if(verbose) message(paste0(ifelse(scan.k.self,'k.self=(','k=('),min,', ',max,') ['))
  xl <- papply(k.seq,function(kv) {
    if(scan.k.self) {
      x <- con$buildGraph(k.self=kv, ..., verbose=FALSE)
    } else {
      x <- con$buildGraph(k=kv, ..., verbose=FALSE)
    }
    if(verbose) message('.')
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
  },n.cores=ncores)
  if(verbose) message(']\n')
  con$n.cores <- ncores

  k.sens <- data.frame(k=k.seq,m=as.numeric(unlist(xl)))
  if(plot) {
    ggplot2::ggplot(k.sens,ggplot2::aes(x=.data$k,y=.data$m))+ggplot2::theme_bw()+ggplot2::geom_point()+ggplot2::geom_smooth()+ggplot2::xlab('modularity')+ggplot2::ylab('k')
  }

  return(k.sens)
}

## Merge into a common matrix, entering 0s for the missing ones
#' @keywords internal
mergeCountMatrices <- function(cms, transposed=FALSE) {
  extendMatrix <- function(mtx, col.names) {
    new.names <- setdiff(col.names, colnames(mtx))
    ext.mtx <- Matrix::Matrix(0, nrow=nrow(mtx), ncol=length(new.names), sparse=TRUE) %>%
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

#' Retrieve sample names per cell
#'
#' @param samples list of samples
#' @return list of sample names
#' getSampleNamePerCell(small_panel.preprocessed)
#'
#' @export
getSampleNamePerCell=function(samples) {
  cl <- lapply(samples, getCellNames)
  return(rep(names(cl), sapply(cl, length)) %>% stats::setNames(unlist(cl)) %>% as.factor())
}

#' Estimate labeling distribution for each vertex, based on provided labels using Random Walk
#'
#' @param graph input graph
#' @param labels vector of factor or character labels, named by cell names
#' @param max.iters maximal number of iterations (default=100)
#' @param tol numeric Absolute tolerance as a stopping criteria (default=0.025)
#' @param verbose boolean Verbose mode (default=TRUE)
#' @param fixed.initial.labels boolean Prohibit changes of initial labels during diffusion (default=TRUE)
#' @keywords internal
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

## Propagate labels using Zhu, Ghahramani, Lafferty (2003) algorithm
## http://mlg.eng.cam.ac.uk/zoubin/papers/zgl.pdf
## TODO: change solver here for something designed for Laplacians. Need to look Daniel Spielman's research

#' @keywords internal
propagateLabelsSolver <- function(graph, labels, solver="mumps") {
  if (!solver %in% c("mumps", "Matrix")){
    stop("Unknown solver: ", solver, ". Only 'mumps' and 'Matrix' are currently supported")
  }

  if (!requireNamespace("rmumps", quietly=TRUE)) {
    warning("Package 'rmumps' is required to use 'mumps' solver, which is the default option. Falliing back to solver='Matrix'")
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
    res <- rmumps::Rmumps$new(laplasian.uu, copy=FALSE)$solve(right.side)
  }

  colnames(res) <- levels(labels)
  rownames(res) <- unlabeled.cbs
  return(rbind(res, type.scores))
}

## exporting to inherit parameters below, leiden.community
#' @export
leidenAlg::leiden.community

#' Increase resolution for a specific set of clusters
#'
#' @param con conos object
#' @param target.clusters clusters for which the resolution should be increased
#' @param clustering name of clustering in the conos object to use. Either 'clustering' or 'groups' must be provided (default=NULL).
#' @param groups set of clusters to use. Ignored if 'clustering' is not NULL (default=NULL).
#' @param method function, used to find communities (default=leiden.community).
#' @param ... additional params passed to the community function
#' @return set of clusters with increased resolution
#' @export
findSubcommunities <- function(con, target.clusters, clustering=NULL, groups=NULL, method=leiden.community, ...) {
  groups <- parseCellGroups(con, clustering, groups)

  groups.raw <- as.character(groups) %>% setNames(names(groups))
  groups <- groups[intersect(names(groups), V(con$graph)$name)]

  if (length(groups) == 0) {
    stop("'groups' not defined for graph object.")
  }

  groups <- droplevels(as.factor(groups)[groups %in% target.clusters])

  if (length(groups) == 0) {
    stop("None of 'target.clusters' can be found in 'groups'.")
  }

  subgroups <- split(names(groups), groups)

  for (n in names(subgroups)) {
    if (length(subgroups[[n]]) < 2){
      next
    }

    new.clusts <- method(induced_subgraph(con$graph, subgroups[[n]]), ...)
    groups.raw[new.clusts$names] <- paste0(n, "_", new.clusts$membership)
  }

  return(groups.raw)
}

#' @keywords internal
parseCellGroups <- function(con, clustering, groups, parse.clusters=TRUE) {
  if (!parse.clusters){
    return(groups)
  }

  if (!is.null(groups)) {
    if (!any(names(groups) %in% names(con$getDatasetPerCell()))){
      stop("'groups' aren't defined for any of the cells.")
    }
    return(groups)
  }

  if (is.null(clustering)) {
    if (length(con$clusters) > 0){
      return(con$clusters[[1]]$groups)
    }
    stop("Either 'groups' must be provided or the conos object must have some clustering estimated")
  }
  if (is.null(clusters[[clustering]])){
    stop(paste("Clustering",clustering,"doesn't exist, run findCommunity() first"))
  }

  return(con$clusters[[clustering]]$groups)
}

#' Estimate entropy of edge weights per cell according to the specified factor.
#' Can be used to visualize alignment quality according to this factor.
#'
#' @param con conos object
#' @param factor.per.cell some factor, which group cells, such as sample or a specific condition
#' @return entropy of edge weights per cell
#' @export
estimateWeightEntropyPerCell <- function(con, factor.per.cell) {
  adj.mat <- igraph::as_adjacency_matrix(con$graph, attr="weight") %>% as("dgTMatrix")
  factor.per.cell %<>% as.factor() %>% .[rownames(adj.mat)]
  weight.sum.per.fac.cell <- getSumWeightMatrix(adj.mat@x, adj.mat@i, adj.mat@j, as.integer(factor.per.cell)) %>%
    `colnames<-`(levels(adj.mat)) %>% `rownames<-`(rownames(adj.mat))

  xt <- table(factor.per.cell)

  entropy.per.cell <- apply(weight.sum.per.fac.cell, 1, entropy::KL.empirical, xt, unit=c('log2')) / log2(length(xt))

  return(entropy.per.cell)
}


