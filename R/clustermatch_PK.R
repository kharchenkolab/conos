
#' @export cpcaJC
cpcaJC <- function(r.n, k=30, k.self=0, k.self.weight=1,community.detection.method = multilevel.community, var.scale =TRUE, min.group.size = 10,ncomps=100, n.odgenes=1000,n.cores=30,return.details=F,verbose=T,neighborhood.average=FALSE,neighborhood.average.k=10,xcp=NULL, ...) {
  require(parallel)
  require(cpca)
  require(Matrix)
  require(abind)
  require(igraph)
  
  k1 <- k2 <- k
  
  if(is.null(xcp)) {
    xcp <- quickCPCA(r.n,k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=verbose,var.scale=var.scale)
  } 
  odgenes <- rownames(xcp$CPC)  
  # common variance scaling
  if (var.scale) {
    cgsf <- do.call(cbind,lapply(r.n,function(x) x$misc$varinfo[odgenes,]$gsf))
    cgsf <- exp(rowMeans(log(cgsf)))
  }
  
  # determine common centering
  cat('centering ...')
  cproj <- lapply(r.n,function(r) {
    x <- r$counts[,odgenes];
    if(var.scale) {
      x@x <- x@x*rep(cgsf,diff(x@p))
    }
    x
  })
  ncells <- unlist(lapply(cproj,nrow));
  centering <- colSums(do.call(rbind,lapply(cproj,colMeans))*ncells)/sum(ncells)
  cat(' done\n')
  
  cat('projecting ...')
  cpproj <- mclapply(cproj,function(x) {
    x <- t(as.matrix(t(x))-centering)
    x %*% xcp$CPC
  },mc.cores=n.cores)
  cat(' done\n')

    ## Get all non-redundant pair of apps
  comb <- combn(names(r.n),2)
  cat('mNN pairs ')
  mnnres <- lapply(1:ncol(comb),function(i) {
    n1 <- comb[1,i]; n2 <- comb[2,i]
    cat(".")
    mnnres <- pagoda2:::interNN(cpproj[[n1]], cpproj[[n2]], k1, k2, 2, verbose=F,neighbourhoodAverage=neighborhood.average,neighbourAvgKA=neighborhood.average.k,neighbourAvgKB=neighborhood.average.k,TRUE)
    mnnres$mA.lab <- rownames(cpproj[[n1]])[mnnres$mA.id]
    mnnres$mB.lab <- rownames(cpproj[[n2]])[mnnres$mB.id]
    mnnres
  })
  cat(' done\n')
  
  ## Merge the results into a edge table
  el <- do.call(rbind, mnnres)[,c('mA.lab','mB.lab')]
  el$w <- 1
  
  # append some local edges
  if(k.self>0) {
    cat('kNN pairs ')
    x <- data.frame(do.call(rbind,lapply(cpproj,function(x) {
      xk <- pagoda2:::hnswKnn2(x,k.self,n.cores,verbose=F)
      xk <- xk[xk$s!=xk$e,]
      cat(".")
      cbind("mA.lab"=rownames(x)[xk$s+1],"mB.lab"=rownames(x)[xk$e+1])
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
  cat('done\n')
  ## Extract groups from this graph
  cls.mem <- membership(cls)
  cls.groups <- as.character(cls.mem)
  names(cls.groups) <- names(cls.mem)
  
  ## Filter groups
  lvls.keep <- names(which(table(cls.groups)  > min.group.size))
  cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
  cls.groups <- as.factor(cls.groups)
  
  if(return.details) {
    return(list(groups=cls.groups,cpproj=cpproj,xcp=xcp,g=g))
  } else {
    cls.groups
  }
  
}

#' @export quickCCA
quickCCA <- function(r.n,k=30,ncomps=100,n.odgenes=NULL,var.scale=T,verbose=T,cgsf=NULL) {
  require(RGCCA)
  require(Matrix)
  
  if(length(r.n)!=2) stop("quickCCA supports only pair alignment")
  
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
  
  # common variance scaling
  if (var.scale) {
    if(is.null(cgsf)) {
      cgsf <- do.call(cbind,lapply(r.n,function(x) x$misc$varinfo[odgenes,]$gsf))
      cgsf <- exp(rowMeans(log(cgsf)))
    }
  }

  # determine common centering
  cproj <- lapply(r.n,function(r) {
    x <- r$counts[,odgenes];
    if(var.scale) {
      x@x <- x@x*rep(cgsf,diff(x@p))
    }
    x
  })
  ncells <- unlist(lapply(cproj,nrow));
  centering <- colSums(do.call(rbind,lapply(cproj,colMeans))*ncells)/sum(ncells)

  cproj <- lapply(cproj,function(x) {
    x <- t(as.matrix(t(x))-centering)
  })
  
  #x <- matrix(rnorm(150), 50, 3)
  #y <- matrix(rnorm(250), 50, 5)
  #z <- rgcca(A=list(x,y),C=matrix(c(0,1,1,0),2,2),tau=c(1,1),ncomp=c(2,2))
  #str(z)
  
  z <- rgcca(A=list(t(cproj[[1]]),t(cproj[[2]])),C=matrix(c(0,1,1,0),2,2),tau=c(0.1,0.1),ncomp=c(ncomps,ncomps),scale=FALSE,verbose=verbose)
  #z$Ys <- lapply(z$Y,function(x) t(t(x)/sqrt(colSums(x^2))))
  #z$P <- list(cproj[[1]] %*% z$Ys[[1]],cproj[[2]] %*% z$Ys[[2]])
  return(z);
}

#' @export ccaJCp
ccaJCp <- function(r.n, k=30, k.self=0, k.self.weight=1,community.detection.method = multilevel.community, var.scale =TRUE, min.group.size = 10,ncomps=100, n.odgenes=1000, n.cores=30, return.details=F,xl=NULL,neighborhood.average=FALSE,neighborhood.average.k=10,verbose=TRUE, ...) {
  require(parallel)
  require(Matrix)
  require(igraph)

  cis <- combn(names(r.n),2)
  if(is.null(xl)) {
    cat('pairwise CCA ')
    xl <- pagoda2:::papply(1:ncol(cis), function(i) {
      xcp <- quickCCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=ifelse(n.cores==1,verbose,FALSE),var.scale=var.scale)
      cat('.')
      xcp
    },n.cores=n.cores);
    names(xl) <- apply(cis,2,paste,collapse='.vs.');
    cat(" done\n")
  } else {
    # match for all the pairs
    mi <- rep(NA,ncol(cis));
    cat('matched ',sum(!is.na(mi)),' out of ',length(mi),' CCA results ... ')
    mi[which(!is.na(match(apply(cis,2,paste,collapse='.vs.'),names(xl))))] <- T;
    mi[which(!is.na(match(apply(cis[c(2,1),],2,paste,collapse='.vs.'),names(xl))))] <- T;
    if(any(is.na(mi))) {
      cat('running ',sum(is.na(mi)),' additional CCAs ')
    }
    xl2 <- pagoda2:::papply(which(is.na(mi)), function(i) {
      xcp <- quickCCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=ifelse(n.cores==1,verbose,FALSE),var.scale=var.scale)
      cat('.')
      xcp
    },n.cores=n.cores);
    names(xl2) <- apply(cis[,which(is.na(mi)),drop=F],2,paste,collapse='.vs.');
    xl <- c(xl,xl2);
    cat(" done\n");
  }
  
  # run mNN separatly as it can't deal with multithreading
  cat('mNN ')
  mnnres <- lapply(1:ncol(cis), function(i) {
    cat(".")
    mnnres <- pagoda2:::interNN(xl[[i]]$a[[1]],xl[[i]]$a[[2]], k, k, 2, verbose=F,neighbourhoodAverage=neighborhood.average,neighbourAvgKA=neighborhood.average.k,neighbourAvgKB=neighborhood.average.k,TRUE)
    #mnnres <- pagoda2:::interNN(cpproj[[n1]], cpproj[[n2]], k1, k2, 2, verbose=F,neighborhood.average,neighborhood.average.k,neighborhood.average.k,TRUE)
    mnnres$mA.lab <- rownames(xl[[i]]$a[[1]])[mnnres$mA.id]
    mnnres$mB.lab <- rownames(xl[[i]]$a[[2]])[mnnres$mB.id]
    mnnres
  })
  cat("done\n")
  ## Merge the results into a edge table
  el <- do.call(rbind, mnnres)[,c('mA.lab','mB.lab','dist')]
  colnames(el) <- c("mA.lab","mB.lab","w")
  
  # append some local edges
  if(k.self>0) {
    cat('kNN pairs ')
    x <- data.frame(do.call(rbind,lapply(r.n,function(x) {
      xk <- pagoda2:::hnswKnn2(x$reductions$PCA,k.self,n.cores,verbose=F)
      xk <- xk[xk$s!=xk$e,]
      cat(".")
      cbind("mA.lab"=rownames(x$reductions$PCA)[xk$s+1],"mB.lab"=rownames(x$reductions$PCA)[xk$e+1])
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
  cat('done\n')
  ## Extract groups from this graph
  cls.mem <- membership(cls)
  cls.groups <- as.character(cls.mem)
  names(cls.groups) <- names(cls.mem)
  
  ## Filter groups
  lvls.keep <- names(which(table(cls.groups)  > min.group.size))
  cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
  cls.groups <- as.factor(cls.groups)
  
  if(return.details) {
    return(list(groups=cls.groups,xl=xl,cls=cls,g=g))
  } else {
    cls.groups
  }
  
}

#' @export quickJNMF
quickJNMF <- function(r.n, k = 30, ncomps =100, n.odgenes=NULL, var.scale=T, verbose =T, cgsf=NULL, maxiter=1000) {
    require(Matrix)
    require(Rjnmf)
    if(length(r.n)!=2) stop('quickJNMF only supports pair alignment')
    
    ## select a common set of genes
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
    ## common variance scaling
    if (var.scale) {
      if(is.null(cgsf)) {
        cgsf <- do.call(cbind,lapply(r.n,function(x) x$misc$varinfo[odgenes,]$gsf))
        cgsf <- exp(rowMeans(log(cgsf)))
      }
    }

    cproj <- lapply(r.n,function(r) {
      x <- r$counts[,odgenes];
      if(var.scale) {
        x@x <- x@x*rep(cgsf,diff(x@p))
      }
      x
    })

    ## cproj <- lapply(r.n,function(r) {
    ##   x <- r$misc$rawCounts[,odgenes];
    ##   if(var.scale) {
    ##     x@x <- x@x*rep(cgsf,diff(x@p))
    ##   }
    ##   x
    ## })
    
     ## cproj <- lapply(r.n,function(r) {
     ##  x <- r$misc$rawCounts[,odgenes];
     ##  if(var.scale) {
     ##    x@x <- log10(x@x*rep(cgsf,diff(x@p))+1)
     ##  }
     ##  x
     ## })
    
    #ncells <- unlist(lapply(cproj,nrow));
    #centering <- colSums(do.call(rbind,lapply(cproj,colMeans))*ncells)/sum(ncells)

    ## cproj <- lapply(cproj,function(x) {
    ##     x <- t(as.matrix(t(x))-centering)
    ## })

    ## # Make sure all values are > 0
    ## cproj <- lapply(cproj, function(x) {
    ##     x <- x - min(x) + 1e-6
    ## })
    cproj <- lapply(cproj,function(x) as.matrix(x))
    
    z <- Rjnmf(t(cproj[[1]]),t(cproj[[2]]),k=ncomps, alpha=0.5,lambda=0.5, epsilon=0.001, maxiter=maxiter, verbose=F, seed=12345)

    rot1 <- cproj[[1]] %*% z$W
    rot2 <- cproj[[2]] %*% z$W

    list(rot1=rot1, rot2=rot2,z=z,cgsf=cgsf)

}

#' @export quickCPCA
quickCPCA <- function(r.n,k=30,ncomps=100,n.odgenes=NULL,var.scale=T,verbose=T,cgsf=NULL,neighborhood.average=FALSE,n.cores=30) {
  require(parallel)
  require(cpca)
  require(Matrix)
  
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
  
  # common variance scaling
  if (var.scale) {
    if(is.null(cgsf)) {
      cgsf <- do.call(cbind,lapply(r.n,function(x) x$misc$varinfo[odgenes,]$gsf))
      cgsf <- exp(rowMeans(log(cgsf)))
    }
  }

  
  if(verbose) cat('calculating covariances for',length(r.n),' datasets ...')
  
  sparse.cov <- function(x){
    n <- nrow(x)
    cMeans <- Matrix::colMeans(x)
    covmat <- (as.matrix(Matrix::crossprod(x)) - n*Matrix::tcrossprod(cMeans))/(n-1)
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
  if(verbose) cat(' done\n')
  
  ncells <- unlist(lapply(r.n,function(x) nrow(x$counts)));
  if(verbose) cat('common PCs ...')
  xcp <- cpca(covl,ncells,ncomp=ncomps)
  #system.time(xcp <- cpca:::cpca_stepwise_base(covl,ncells,k=ncomps))
  #xcp <- cpc(abind(covl,along=3),k=ncomps)
  rownames(xcp$CPC) <- odgenes;
  #xcp$rot <- xcp$CPC*cgsf;
  if(verbose) cat(' done\n')
  return(xcp);
}

#' @export cpcaJCp
cpcaJCp <- function(r.n, k=30, k.self=0, k.self.weight=1,community.detection.method = multilevel.community, reduction.method='CPCA', var.scale =TRUE, min.group.size = 10,ncomps=100, n.odgenes=1000, n.cores=30, return.details=F,xl=NULL,neighborhood.average=FALSE,neighborhood.average.k=10,groups=NULL, ...) {
  require(parallel)
  require(cpca)
  require(Matrix)
  require(abind)
  require(igraph)
  k1 <- k2 <- k;

  
  cis <- combn(names(r.n),2)
  if(is.null(xl)) {
    cat('pairwise',reduction.method)
    xl <- pagoda2:::papply(1:ncol(cis), function(i) {
      if(reduction.method=='CPCA') {
        xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
      } else if(reduction.method=='JNMF') {
        xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
      }
      cat('.')
      xcp
    },n.cores=n.cores);
    names(xl) <- apply(cis,2,paste,collapse='.vs.');
    cat("done\n")
  } else {
    # match for all the pairs
    mi <- rep(NA,ncol(cis));
    nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    cat('matched',sum(!is.na(mi)),'out of',length(mi),'CPCA pairs ... ')
    if(any(is.na(mi))) {
      cat('running',sum(is.na(mi)),'additional CPCAs ')
      xl2 <- pagoda2:::papply(which(is.na(mi)), function(i) {
        if(reduction.method=='CPCA') {
          xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
        } else if(reduction.method=='JNMF') {
          xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
        }
        #xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
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
    if(any(is.na(mi))) { stop("unable to get complete set of CPCA results") }
    xl <- xl[mi]
    cat(" done\n");
  }
  
  # run mNN separatly as it can't deal with multithreading
  cat('mNN ')
  mnnres <- lapply(1:ncol(cis), function(i) {
    r.ns <- r.n[cis[,i]]
    if(!is.null(xl[[i]]$rot1)) {
      # JNMF
      #system.time(yn <- pagoda2:::crossNN(x,x,k,2,2.0,verbose,n.cores))
      mnnres <- pagoda2:::interNN(xl[[i]]$rot1, xl[[i]]$rot2, k, k, 2, verbose=F,neighbourhoodAverage=neighborhood.average,neighbourAvgKA=neighborhood.average.k,neighbourAvgKB=neighborhood.average.k,TRUE)
      mnnres$mA.lab <- rownames(xl[[i]]$rot1)[mnnres$mA.id]
      mnnres$mB.lab <- rownames(xl[[i]]$rot2)[mnnres$mB.id]
    } else {
      if(!is.null(xl[[i]]$CPC)) {
        # CPCA
        rot <- xl[[i]]$CPC;
        odgenes <- rownames(xl[[i]]$CPC)
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
        x
      })
      ncells <- unlist(lapply(cproj,nrow));
      centering <- colSums(do.call(rbind,lapply(cproj,colMeans))*ncells)/sum(ncells)
      cpproj <- lapply(cproj,function(x) {
        x <- t(as.matrix(t(x))-centering)
        x %*% rot;
      })
      n1 <- cis[1,i]; n2 <- cis[2,i]
      cat(".")
      mnnres <- pagoda2:::interNN(cpproj[[n1]], cpproj[[n2]], k, k, 2, verbose=F,neighbourhoodAverage=neighborhood.average,neighbourAvgKA=neighborhood.average.k,neighbourAvgKB=neighborhood.average.k,TRUE)
      #mnnres <- pagoda2:::interNN(cpproj[[n1]], cpproj[[n2]], k1, k2, 2, verbose=F,neighborhood.average,neighborhood.average.k,neighborhood.average.k,TRUE)
      mnnres$mA.lab <- rownames(cpproj[[n1]])[mnnres$mA.id]
      mnnres$mB.lab <- rownames(cpproj[[n2]])[mnnres$mB.id]
    }
    mnnres
  })
  cat("done\n")
  ## Merge the results into a edge table
  #el <- do.call(rbind, mnnres)[,c('mA.lab','mB.lab')]
  el <- do.call(rbind, mnnres)[,c('mA.lab','mB.lab','dist')]
  colnames(el) <- c("mA.lab","mB.lab","w")
    

  el$w <- 1
  
  # append some local edges
  if(k.self>0) {
    cat('kNN pairs ')
    x <- data.frame(do.call(rbind,lapply(r.n,function(x) {
      xk <- pagoda2:::hnswKnn2(x$reductions$PCA,k.self,n.cores,verbose=F)
      xk <- xk[xk$s!=xk$e,]
      cat(".")
      cbind("mA.lab"=rownames(x$reductions$PCA)[xk$s+1],"mB.lab"=rownames(x$reductions$PCA)[xk$e+1])
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
  cat('done\n')
  ## Extract groups from this graph
  cls.mem <- membership(cls)
  cls.groups <- as.character(cls.mem)
  names(cls.groups) <- names(cls.mem)
  
  ## Filter groups
  lvls.keep <- names(which(table(cls.groups)  > min.group.size))
  cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
  cls.groups <- as.factor(cls.groups)
  
  if(return.details) {
    return(list(groups=cls.groups,xl=xl,cls=cls,g=g))
  } else {
    cls.groups
  }
  
}

#' @export cpcaJCp2
cpcaJCp2 <- function(r.n, k=30, k.self=0, k.self.weight=1,community.detection.method = multilevel.community, reduction.method='CPCA', var.scale =TRUE, min.group.size = 10,ncomps=100, n.odgenes=1000, n.cores=30, return.details=F,xl=NULL,neighborhood.average=FALSE,neighborhood.average.k=5,groups=NULL, ...) {
  require(parallel)
  require(cpca)
  require(Matrix)
  require(abind)
  require(igraph)
  k1 <- k2 <- k;

  if(neighborhood.average) {
    cat("neighborhood averaging ")
    r.n <- lapply(r.n,function(r) {
      xk <- pagoda2:::crossNN(r$reductions$PCA,r$reductions$PCA,neighborhood.average.k,2,2.0,FALSE,n.cores)
      xk@x <- pmax(1-xk@x,0);
      diag(xk) <- 1;
      xk <- t(t(xk)/colSums(xk))
      colnames(xk) <- rownames(xk) <- rownames(r$reductions$PCA)
      r$misc$edgeMat$quickCPCA <- xk;
      cat(".")
      r
    })
    cat(" done\n")
  }
    

  cis <- combn(names(r.n),2)
  if(is.null(xl)) {
    cat('pairwise',reduction.method,' ')
    xl <- pagoda2:::papply(1:ncol(cis), function(i) {
      if(reduction.method=='CPCA') {
        xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average)
      } else if(reduction.method=='JNMF') {
        xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
      }
      cat('.')
      xcp
    },n.cores=n.cores);
    names(xl) <- apply(cis,2,paste,collapse='.vs.');
    cat("done\n")
  } else {
    # match for all the pairs
    mi <- rep(NA,ncol(cis));
    nm <- match(apply(cis,2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(xl));
    mi[which(!is.na(nm))] <- na.omit(nm);
    cat('matched',sum(!is.na(mi)),'out of',length(mi),'CPCA pairs ... ')
    if(any(is.na(mi))) {
      cat('running',sum(is.na(mi)),'additional CPCAs ')
      xl2 <- pagoda2:::papply(which(is.na(mi)), function(i) {
        if(reduction.method=='CPCA') {
          xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average)
        } else if(reduction.method=='JNMF') {
          xcp <- quickJNMF(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
        }
        #xcp <- quickCPCA(r.n[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale)
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
    if(any(is.na(mi))) { stop("unable to get complete set of CPCA results") }
    xl <- xl[mi]
    cat(" done\n");
  }
  
  # run mNN separatly as it can't deal with multithreading
  cat('mNN ')
  mnnres <- lapply(1:ncol(cis), function(i) {
    r.ns <- r.n[cis[,i]]
    if(!is.null(xl[[i]]$rot1)) {
      # JNMF
      
      #mnnres <- pagoda2:::interNN(xl[[i]]$rot1, xl[[i]]$rot2, k, k, 2, verbose=F,neighbourhoodAverage=neighborhood.average,neighbourAvgKA=neighborhood.average.k,neighbourAvgKB=neighborhood.average.k,TRUE)
      #mnnres$mA.lab <- rownames(xl[[i]]$rot1)[mnnres$mA.id]
      #mnnres$mB.lab <- rownames(xl[[i]]$rot2)[mnnres$mB.id]

      n12 <- pagoda2:::crossNN(xl[[i]]$rot1,xl[[i]]$rot2,k,2,2.0,FALSE,n.cores)
      n21 <- pagoda2:::crossNN(xl[[i]]$rot2,xl[[i]]$rot1,k,2,2.0,FALSE,n.cores)
      mnn <- drop0(n21*t(n12))
      mnn <- as(n21*t(n12),'dgTMatrix')
      return(data.frame('mA.lab'=rownames(cpproj[[n1]])[mnn@i+1],'mB.lab'=rownames(cpproj[[n2]])[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
    } else {
      if(!is.null(xl[[i]]$CPC)) {
        # CPCA
        rot <- xl[[i]]$CPC;
        odgenes <- rownames(xl[[i]]$CPC)
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
      ncells <- unlist(lapply(cproj,nrow));
      centering <- colSums(do.call(rbind,lapply(cproj,colMeans))*ncells)/sum(ncells)
      cpproj <- lapply(cproj,function(x) {
        x <- t(as.matrix(t(x))-centering)
        x %*% rot;
      })
      n1 <- cis[1,i]; n2 <- cis[2,i]
      cat(".")

      #mnnres <- pagoda2:::interNN(cpproj[[n1]], cpproj[[n2]], k, k, 2, verbose=F,neighbourhoodAverage=neighborhood.average,neighbourAvgKA=neighborhood.average.k,neighbourAvgKB=neighborhood.average.k,TRUE)
      #mnnres <- pagoda2:::interNN(cpproj[[n1]], cpproj[[n2]], k1, k2, 2, verbose=F,neighborhood.average,neighborhood.average.k,neighborhood.average.k,TRUE)
      #mnnres$mA.lab <- rownames(cpproj[[n1]])[mnnres$mA.id]
      #mnnres$mB.lab <- rownames(cpproj[[n2]])[mnnres$mB.id]
      
      n12 <- pagoda2:::crossNN(cpproj[[n1]],cpproj[[n2]],k,2,2.0,FALSE,n.cores)
      n21 <- pagoda2:::crossNN(cpproj[[n2]],cpproj[[n1]],k,2,2.0,FALSE,n.cores)
      #colnames(n12) <- rownames(n21) <- rownames(cpproj[[n1]])
      #colnames(n21) <- rownames(n12) <- rownames(cpproj[[n2]])
      mnn <- drop0(n21*t(n12))
      mnn <- as(n21*t(n12),'dgTMatrix')
      return(data.frame('mA.lab'=rownames(cpproj[[n1]])[mnn@i+1],'mB.lab'=rownames(cpproj[[n2]])[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
      
    }
    mnnres
  })
  cat("done\n")
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
      xk <- pagoda2:::crossNN(x$reductions$PCA,x$reductions$PCA,k.self,2,2.0,FALSE,n.cores)
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
  cat('done\n')
  ## Extract groups from this graph
  cls.mem <- membership(cls)
  cls.groups <- as.character(cls.mem)
  names(cls.groups) <- names(cls.mem)
  
  ## Filter groups
  lvls.keep <- names(which(table(cls.groups)  > min.group.size))
  cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
  cls.groups <- as.factor(cls.groups)
  
  if(return.details) {
    return(list(groups=cls.groups,xl=xl,cls=cls,g=g))
  } else {
    cls.groups
  }
  
}

