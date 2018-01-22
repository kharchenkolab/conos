

#' Perform CCA alignment for all pairs of apps
#' @param r.n list of pagoda2 apps
#' @param xl list contains previouls performed alignments
#' @param n.cores number of cores to use, limited by the number of available pairs
#' @param k k parameter for quickCCA
#' @param ncomps ncomps for quickCCA
#' @param n.odgenes n.odgenes for quickCCA
#' @param verbose logical verbose
#' @param var.scale var.scae for quickCCA
#' @return a list of pairwise alignments
ccaPairs <- function(r.n, xl = NULL, n.cores = 30, k = 30,ncomps=100,n.odgenes=1000,verbose=TRUE,var.scale=TRUE) {
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
  xl
}


ccaJCp <- function(r.n, k=30, k.self=0, k.self.weight=1,community.detection.method = multilevel.community, var.scale =TRUE, min.group.size = 10,ncomps=100, n.odgenes=1000, n.cores=30, return.details=F,xl=NULL,neighborhood.average=FALSE,neighborhood.average.k=10,verbose=TRUE, ...) {
  #require(parallel)
  #require(Matrix)
  #require(igraph)

  xl <- ccaPairs(r.n = r.n, xl= xl, n.cores=n.cores, k = k, ncomps=ncomps,n.odgenes=n.odgenes,verbose=verbose, var.scale=var.scale)
  
  # run mNN separatly as it can't deal with multithreading
  cat('mNN ')
  mnnres <- lapply(1:ncol(cis), function(i) {
    cat(".")
    mnnres <- clusterMatch:::interNN(xl[[i]]$a[[1]],xl[[i]]$a[[2]], k, k, 2, verbose=F,neighbourhoodAverage=neighborhood.average,neighbourAvgKA=neighborhood.average.k,neighbourAvgKB=neighborhood.average.k,TRUE)
    #mnnres <- clusterMatch:::interNN(cpproj[[n1]], cpproj[[n2]], k1, k2, 2, verbose=F,neighborhood.average,neighborhood.average.k,neighborhood.average.k,TRUE)
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
