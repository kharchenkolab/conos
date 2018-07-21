#' @import Matrix
#' @import igraph
#' @importFrom parallel mclapply
NULL

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

cpcaFast <- function(covl,ncells,ncomp=10,maxit=1000,tol=1e-6,use.irlba=TRUE,verbose=F) {
  if(use.irlba) {
    # irlba initialization
    p <- nrow(covl[[1]]);
    S <- matrix(0, nrow = p, ncol = p)  
    for(i in 1:length(covl)) {
      S <- S + (ncells[i] / sum(ncells)) * covl[[i]]
    }
    ev <- irlba::irlba(S,ncomp)
    cc <- abind::abind(covl,along=3)
    cpcaF(cc,ncells,ncomp,maxit,tol,eigenvR=ev$v,verbose)
  } else {
    cpcaF(cc,ncells,ncomp,maxit,tol,verbose=verbose)
  }
}


#' Perform cpca on two samples
#' @param r.n list of p2 objects
#' @param k neighborhood size to use
#' @param ncomps number of components to calculate (default=100)
#' @param n.odgenes number of overdispersed genes to take from each dataset
#' @param var.scale whether to scale variance (default=TRUE)
#' @param verbose whether to be verbose
#' @param cgsf an optional set of common genes to align on 
#' @param neighborhood.average use neighborhood average values
#' @param n.cores number of cores to use 
#' @export quickCPCA
quickCPCA <- function(r.n,k=30,ncomps=100,n.odgenes=NULL,var.scale=TRUE,verbose=T,cgsf=NULL,neighborhood.average=FALSE,n.cores=30) {
  #require(parallel)
  #require(cpca)
  #require(Matrix)
  
  # select a common set of genes
  if(is.null(cgsf)) {
    if(is.null(n.odgenes)) {
      odgenes <- table(unlist(lapply(r.n,function(x) x$misc$odgenes)))
    } else {
      odgenes <- table(unlist(lapply(r.n,function(x) rownames(x$misc$varinfo)[(order(x$misc$varinfo$lp,decreasing=F)[1:min(ncol(x$counts),n.odgenes)])])))
    }
    odgenes <- odgenes[names(odgenes) %in% Reduce(intersect,lapply(r.n,function(x) colnames(x$counts)))]
    odgenes <- names(odgenes)[1:min(length(odgenes),n.odgenes)]
  } else {
    odgenes <- names(cgsf)
  }
  cat("using",length(odgenes),"odgenes\n") 
  # common variance scaling
  if (var.scale) {
    if(is.null(cgsf)) {
      cgsf <- do.call(cbind,lapply(r.n,function(x) x$misc$varinfo[odgenes,]$gsf))
      cgsf <- exp(rowMeans(log(cgsf)))
    }
  }

  
  if(verbose) cat('calculating covariances for',length(r.n),' datasets ...')

  # use internal C++ implementation
  sparse.cov <- function(x,cMeans=NULL){
    if(is.null(cMeans)) {  cMeans <- Matrix::colMeans(x) }
    covmat <- spcov(x,cMeans);
  }

  
  covl <- lapply(r.n,function(r) {
    x <- r$counts[,odgenes];
    if(var.scale) {
      x@x <- x@x*rep(cgsf,diff(x@p))
    }
    if(neighborhood.average) {
      xk <- r$misc$edgeMat$quickCPCA;
      x <- t(xk) %*% x
    }
    sparse.cov(x)
  })
  
  ## # centering
  ## if(common.centering) {
  ##   ncells <- unlist(lapply(covl,nrow));
  ##   centering <- colSums(do.call(rbind,lapply(covl,colMeans))*ncells)/sum(ncells)
  ## } else {
  ##   centering <- NULL;
  ## }
  
  ## covl <- lapply(covl,sparse.cov,cMeans=centering)
  
  if(verbose) cat(' done\n')
  
  ncells <- unlist(lapply(r.n,function(x) nrow(x$counts)));
  if(verbose) cat('common PCs ...')
  #xcp <- cpca(covl,ncells,ncomp=ncomps)
  xcp <- cpcaFast(covl,ncells,ncomp=ncomps,verbose=TRUE,maxit=500,tol=1e-5);
  #system.time(xcp <- cpca:::cpca_stepwise_base(covl,ncells,k=ncomps))
  #xcp <- cpc(abind(covl,along=3),k=ncomps)
  rownames(xcp$CPC) <- odgenes;
  #xcp$rot <- xcp$CPC*cgsf;
  if(verbose) cat(' done\n')
  return(xcp);
}

##' Construct joint graph and detect communities
##'
##' @param r.n 
##' @param k 
##' @param k.self 
##' @param k.self.weight 
##' @param community.detection.method 
##' @param reduction.method 
##' @param matching.method 
##' @param var.scale 
##' @param min.group.size 
##' @param ncomps 
##' @param n.odgenes 
##' @param n.cores 
##' @param return.details 
##' @param xl 
##' @param neighborhood.average 
##' @param neighborhood.average.k 
##' @param groups 
##' @param common.centering 
##' @param ... 
##' @return 
##' @export
conosCluster <- function(r.n, k=30, k.self=10, k.self.weight=0.1,community.detection.method = multilevel.community, reduction.method='CPCA', matching.method='mNN', var.scale =TRUE, min.group.size = 0, ncomps=50, n.odgenes=1000, n.cores=30, return.details=T,xl=NULL,neighborhood.average=FALSE,neighborhood.average.k=10,groups=NULL, common.centering=TRUE, ...) {

  # keep k symmetric
  k1 <- k2 <- k;

  if(neighborhood.average) {
    cat("neighborhood averaging ")
    r.n <- lapply(r.n,function(r) {
      xk <- pagoda2:::n2Knn(r$reductions$PCA[rownames(r$counts),],neighborhood.average.k,n.cores,FALSE)
      xk@x <- pmax(1-xk@x,0);
      diag(xk) <- 1;
      xk <- t(t(xk)/colSums(xk))
      colnames(xk) <- rownames(xk) <- rownames(r$counts)
      r$misc$edgeMat$quickCPCA <- xk;
      cat(".")
      r
    })
    cat(" done\n")
  }

  cis <- combn(names(r.n),2)
  if(is.null(xl)) {
    cat('pairwise',reduction.method,'(',ncol(cis),'pairs) ')
    xl <- papply(1:ncol(cis), function(i) {
      if(reduction.method=='CPCA') {
        xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average)
      } else if(reduction.method=='JNMF') {
        xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average,maxiter=3e3)
      } else if (reduction.method == 'GeneSpace') {
                xcp <- quickNULL(p2.objs = p2list[cis[,i]], n.odgenes = n.odgenes, var.scale = var.scale, verbose = FALSE, neighborhood.average=neighborhood.average);
      } else {
        # TODO: add canonical correlation
        # TODO: add union of individual PCs
        stop(paste("unknown reduction method",reduction.method))
      }
      cat('.')
      xcp
    },n.cores=n.cores,mc.preschedule=T);
    names(xl) <- apply(cis,2,paste,collapse='.vs.');
    cat("done\n")
  } else {
    # match for all the pairs
    mi <- rep(NA,ncol(cis));
    nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    # try reverse match as well
    nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    cat('matched',sum(!is.na(mi)),'out of',length(mi),' pairs ... ')
    if(any(is.na(mi))) { # some pairs are missing
      cat('running',sum(is.na(mi)),'additional pairs ')
      xl2 <- papply(which(is.na(mi)), function(i) {
        if(reduction.method=='CPCA') {
          xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average)
        } else if(reduction.method=='JNMF') {
          xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average,maxiter=3e3)
        } else if (reduction.method == 'GeneSpace') {
          xcp <- quickNULL(p2.objs = p2list[cis[,i]], n.odgenes = n.odgenes, var.scale = var.scale, verbose = FALSE, neighborhood.average=neighborhood.average);
        }
           cat('.')
        xcp
      },n.cores=n.cores);
      names(xl2) <- apply(cis[,which(is.na(mi)),drop=F],2,paste,collapse='.vs.');
      xl <- c(xl,xl2);
    }
    # re-do the match and order
    mi <- rep(NA,ncol(cis));
    nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    if(any(is.na(mi))) { stop("unable to get complete set of pair comparison results") }
    xl <- xl[mi]
    cat(" done\n");
  }
  
  # run mNN separatly as it can't deal with multithreading
  cat('mNN ')
  mnnres <- papply(1:ncol(cis), function(i) {
    r.ns <- r.n[cis[,i]]
    if(!is.null(xl[[i]]$rot1)) {
      # JNMF
      n12 <- pagoda2:::n2CrossKnn(xl[[i]]$rot1,xl[[i]]$rot2,k,1,FALSE)
      n21 <- pagoda2:::n2CrossKnn(xl[[i]]$rot2,xl[[i]]$rot1,k,1,FALSE)
      mnn <- drop0(n21*t(n12))
      mnn <- as(n21*t(n12),'dgTMatrix')
      cat(".")
      
      return(data.frame('mA.lab'=rownames(xl[[i]]$rot1)[mnn@i+1],'mB.lab'=rownames(xl[[i]]$rot2)[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
      #return(data.frame('mA.lab'=rownames(xl[[i]]$rot1)[mnn@i+1],'mB.lab'=rownames(xl[[i]]$rot2)[mnn@j+1],'w'=1/pmax(1,log(mnn@x)),stringsAsFactors=F))
      
    } else if (!is.null(xl[[i]]$CPC)) { # CPCA or GSVD
      common.genes <- Reduce(intersect,lapply(r.ns,function(x) colnames(x$counts)))
      if(!is.null(xl[[i]]$CPC)) {
        # CPCA
        odgenes <- intersect(rownames(xl[[i]]$CPC),common.genes)
        rot <- xl[[i]]$CPC[odgenes,];
      } else if(!is.null(xl[[i]]$o$Q)) {
        # GSVD
        rot <- xl[[i]]$o$Q;
        odgenes <- rownames(rot) <- colnames(xl[[i]]$o$A);
      } else {
        stop("unknown reduction provided")
      }
      if (var.scale) {
        cgsf <- do.call(cbind,lapply(r.ns,function(x) x$misc$varinfo[odgenes,]$gsf))
        cgsf <- exp(rowMeans(log(cgsf)))
      }
      # create matrices, adjust variance
      cproj <- lapply(r.ns,function(r) {
        x <- r$counts[,odgenes];
        if(var.scale) {
          x@x <- x@x*rep(cgsf,diff(x@p))
        }
        if(neighborhood.average) {
          xk <- r$misc$edgeMat$quickCPCA;
          x <- t(xk) %*% x
        }
        x
      })
      if(common.centering) {
        ncells <- unlist(lapply(cproj,nrow));
        centering <- colSums(do.call(rbind,lapply(cproj,colMeans))*ncells)/sum(ncells)
      } else {
        centering <- NULL;
      }
      cpproj <- lapply(cproj,function(x) {
        if(is.null(centering)) { centering <- colMeans(x) }
        x <- t(as.matrix(t(x))-centering)
        x %*% rot;
      })
      n1 <- cis[1,i]; n2 <- cis[2,i]
      cat(".")
      
      n12 <- pagoda2:::n2CrossKnn(cpproj[[n1]],cpproj[[n2]],k,1,FALSE)
      n21 <- pagoda2:::n2CrossKnn(cpproj[[n2]],cpproj[[n1]],k,1,FALSE)
      #colnames(n12) <- rownames(n21) <- rownames(cpproj[[n1]])
      #colnames(n21) <- rownames(n12) <- rownames(cpproj[[n2]])
      mnn <- drop0(n21*t(n12))
      mnn <- as(n21*t(n12),'dgTMatrix')
      return(data.frame('mA.lab'=rownames(cpproj[[n1]])[mnn@i+1],'mB.lab'=rownames(cpproj[[n2]])[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
      

    } else if (!is.null(xl[[i]]$genespace1)) {
      ## Overdispersed Gene space
      n12 <- pagoda2:::n2CrossKnn(as.matrix(xl[[i]]$genespace1), as.matrix(xl[[i]]$genespace2),k,1,FALSE)
      n21 <- pagoda2:::n2CrossKnn(as.matrix(xl[[i]]$genespace2), as.matrix(xl[[i]]$genespace1),k,1,FALSE)

      ##
      mnn <- drop0(n21*t(n12))
      mnn <- as(n21*t(n12),'dgTMatrix')

      ## return
      ret.df <- data.frame(
        'mA.lab'=rownames(xl[[i]]$genespace1)[mnn@i+1],
        'mB.lab'=rownames(xl[[i]]$genespace2)[mnn@j+1],
        'w'=pmax(1-mnn@x,0),stringsAsFactors=F
      )
    } else {
      stop('unknown reduction provided');
    }
    mnnres
  },n.cores=n.cores)
  cat(" done\n")
  ## Merge the results into a edge table
  #el <- do.call(rbind, mnnres)[,c('mA.lab','mB.lab')]
  #el <- do.call(rbind, mnnres)[,c('mA.lab','mB.lab','dist')]
  #colnames(el) <- c("mA.lab","mB.lab","w")
  el <- do.call(rbind,mnnres)

  #el$w <- 1
  
  # append some local edges
  if(k.self>0) {
    cat('kNN pairs ')
    x <- data.frame(do.call(rbind,lapply(r.n,function(x) {
      xk <- pagoda2:::n2Knn(x$reductions$PCA,k.self,1,FALSE)
      diag(xk) <- 0;
      xk <- as(xk,'dgTMatrix')
      cat(".")
      return(data.frame('mA.lab'=rownames(x$reductions$PCA)[xk@i+1],'mB.lab'=rownames(x$reductions$PCA)[xk@j+1],'w'=pmax(1-xk@x,0),stringsAsFactors=F))
    })),stringsAsFactors = F)
    x$w <- k.self.weight
    cat(' done\n')
    el <- rbind(el,x)
  }
  
  g  <- graph_from_edgelist(as.matrix(el[,c(1,2)]), directed =FALSE)
  E(g)$weight <- el[,3]

  ## Do community detection on this graph
  cat('detecting clusters ...');
  cls <- community.detection.method(g, ...)
  cat(' done\n')

  ## Extract groups from this graph
  cls.mem <- membership(cls)
  cls.groups <- as.character(cls.mem)
  names(cls.groups) <- names(cls.mem)
  
  ## Filter groups
  lvls.keep <- names(which(table(cls.groups)  > min.group.size))
  cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
  cls.groups <- as.factor(cls.groups)
  
  if(return.details) {
    return(list(groups=cls.groups,xl=xl,cls=cls,mnnres=mnnres,g=g))
  } else {
    cls.groups
  }
  
}



# use mclapply if available, fall back on BiocParallel, but use regular
# lapply() when only one core is specified
papply <- function(...,n.cores=detectCores(), mc.preschedule=FALSE) {
  if(n.cores>1) {
    if(requireNamespace("parallel", quietly = TRUE)) {
      return(mclapply(...,mc.cores=n.cores,mc.preschedule=mc.preschedule))
    }

    if(requireNamespace("BiocParallel", quietly = TRUE)) {
      # It should never happen because parallel is specified in Imports
      return(BiocParallel::bplapply(... , BPPARAM = BiocParallel::MulticoreParam(workers = n.cores)))
    }
  }

  # fall back on lapply
  lapply(...)
}

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
    new.assign <- unlist(unname(parallel::mclapply(p2list, function(p2o) {
        try({
            ## get global cluster centroids for cells in this app
            global.cluster.filtered.bd <- factorBreakdown(global.cluster.filtered)
            global.cl.centers <- do.call(rbind, lapply(global.cluster.filtered.bd, function(cells) {
                cells <- cells[cells %in% rownames(p2o$counts)]
                if (length(cells) > 1) {
                    Matrix::colSums(p2o$counts[cells,])
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
