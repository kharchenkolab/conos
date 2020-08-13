# TODO: multitrap method
##' mutlilevel+walktrap communities
##'
##' Constructrs a two-step clustering, first running multilevel.communities, and then walktrap.communities within each
##' These are combined into an overall hierarchy
##' @param graph graph
##' @param n.cores number of cores to use
##' @param hclust.link link function to use when clustering multilevel communities (based on collapsed graph connectivity)
##' @param min.community.size minimal community size parameter for the walktrap communities .. communities smaller than that will be merged
##' @param verbose whether to output progress messages
##' @param level what level of multitrap clustering to use in the starting step. By default, uses the top level. An integer can be specified for a lower level (i.e. 1).
##' @param ... passed to walktrap
##' @return a fakeCommunities object that has methods membership() and as.dendrogram() to mimic regular igraph returns
multitrap.community <- function(graph, n.cores=parallel::detectCores(logical=FALSE), hclust.link='single', min.community.size=10, verbose=FALSE, level=NULL, ...) {
  .Deprecated()

  if(verbose) cat("running multilevel ... ");
  mt <- multilevel.community(graph);

  if(is.null(level)) {
    # get higest level (to avoid oversplitting at the initial step)
    mem <- membership(mt);
  } else {
    # get the specified level
    mem <- mt$memberships[level,]; names(mem) <- mt$names;
  }

  if(verbose) cat("found",length(unique(mem)),"communities\nrunning walktraps ... ")

  # calculate hierarchy on the multilevel clusters
  cgraph <- getClusterGraph(graph,mem)
  chwt <- walktrap.community(cgraph,steps=8)
  d <- as.dendrogram(chwt);



  wtl <- conos:::papply(sn(unique(mem)), function(cluster) {
    cn <- names(mem)[which(mem==cluster)]
    sg <- induced.subgraph(graph,cn)
    walktrap.community(induced.subgraph(graph,cn))
  },n.cores=n.cores)

  mbl <- lapply(wtl,membership);
  # correct small communities
  mbl <- lapply(mbl,function(x) {
    tx <- table(x)
    ivn <- names(tx)[tx<min.community.size]
    if(length(ivn)>1) {
      x[x %in% ivn] <- as.integer(ivn[1]); # collapse into one group
    }
    x
  })

  if(verbose) cat("found",sum(unlist(lapply(mbl,function(x) length(unique(x))))),"communities\nmerging dendrograms ... ")


  wtld <- lapply(wtl,as.dendrogram)
  max.height <- max(unlist(lapply(wtld,attr,'height')))

  # shift leaf ids to fill in 1..N range
  mn <- unlist(lapply(wtld,attr,'members'))
  shift.leaf.ids <- function(l,v) { if(is.leaf(l)) { la <- attributes(l); l <- as.integer(l)+v; attributes(l) <- la; }; l  }
  nshift <- cumsum(c(0,mn))[-(length(mn)+1)]; names(nshift) <- names(mn); # how much to shift ids in each tree

  get.heights <- function(l) {
    if(is.leaf(l)) {
      return(attr(l,'height'))
    } else {
      return(c(attr(l,'height'),unlist(lapply(l,get.heights))))
    }
  }
  min.d.height <- min(get.heights(d))
  height.scale <- length(wtld)*2
  height.shift <- 2

  shift.heights <- function(l,s) { attr(l,'height') <- attr(l,'height')+s; l }

  glue.dends <- function(l) {
    if(is.leaf(l)) {
      nam <- as.character(attr(l,'label'));
      id <- dendrapply(wtld[[nam]], shift.leaf.ids, v=nshift[nam])
      return(dendrapply(id,shift.heights,s=max.height-attr(id,'height')))

    }
    attr(l,'height') <- (attr(l,'height')-min.d.height)*height.scale + max.height + height.shift;
    l[[1]] <- glue.dends(l[[1]]); l[[2]] <- glue.dends(l[[2]])
    attr(l,'members') <- attr(l[[1]],'members') + attr(l[[2]],'members')
    return(l)
  }
  combd <- glue.dends(d)
  if(verbose) cat("done\n");

  # combined clustering factor
  fv <- unlist(lapply(sn(names(wtl)),function(cn) {
    paste(cn,as.character(mbl[[cn]]),sep='-')
  }))
  names(fv) <- unlist(lapply(mbl,names))

  # enclose in a masquerading class
  res <- list(membership=fv, dendrogram=combd, algorithm='multitrap', names=names(fv));
  class(res) <- rev("fakeCommunities");
  return(res);

}


##' mutlilevel+multilevel communities
##'
##' Constructrs a two-step clustering, first running multilevel.communities, and then walktrap.communities within each
##' These are combined into an overall hierarchy
##' @param graph graph
##' @param n.cores number of cores to use
##' @param hclust.link link function to use when clustering multilevel communities (based on collapsed graph connectivity)
##' @param min.community.size minimal community size parameter for the walktrap communities .. communities smaller than that will be merged
##' @param verbose whether to output progress messages
##' @param level what level of multitrap clustering to use in the starting step. By default, uses the top level. An integer can be specified for a lower level (i.e. 1).
##' @param ... passed to walktrap
##' @return a fakeCommunities object that has methods membership() and as.dendrogram() to mimic regular igraph returns
multimulti.community <- function(graph, n.cores=parallel::detectCores(logical=FALSE), hclust.link='single', min.community.size=10, verbose=FALSE, level=NULL, ...) {
  .Deprecated()

  if(verbose) cat("running multilevel 1 ... ");
  mt <- multilevel.community(graph);

  if(is.null(level)) {
    # get higest level (to avoid oversplitting at the initial step)
    mem <- membership(mt);
  } else {
    # get the specified level
    mem <- mt$memberships[level,]; names(mem) <- mt$names;
  }

  if(verbose) cat("found",length(unique(mem)),"communities\nrunning multilevel 2 ... ")

  # calculate hierarchy on the multilevel clusters
  cgraph <- getClusterGraph(graph,mem)
  chwt <- walktrap.community(cgraph,steps=8)
  d <- as.dendrogram(chwt);


  wtl <- conos:::papply(sn(unique(mem)), function(cluster) {
    cn <- names(mem)[which(mem==cluster)]
    sg <- induced.subgraph(graph,cn)
    multilevel.community(induced.subgraph(graph,cn))
  },n.cores=n.cores)

  mbl <- lapply(wtl,membership);
  # correct small communities
  mbl <- lapply(mbl,function(x) {
    tx <- table(x)
    ivn <- names(tx)[tx<min.community.size]
    if(length(ivn)>1) {
      x[x %in% ivn] <- as.integer(ivn[1]); # collapse into one group
    }
    x
  })

  if(verbose) cat("found",sum(unlist(lapply(mbl,function(x) length(unique(x))))),"communities\nmerging ... ")

  # combined clustering factor
  fv <- unlist(lapply(sn(names(wtl)),function(cn) {
    paste(cn,as.character(mbl[[cn]]),sep='-')
  }))
  names(fv) <- unlist(lapply(mbl,names))

  # enclose in a masquerading class
  res <- list(membership=fv, dendrogram=NULL, algorithm='multimulti', names=names(fv));
  class(res) <- rev("fakeCommunities");
  return(res);

}

##' Leiden algorithm community detection
##'
##' Detect communities using Leiden algorithm (implementation copied from https://github.com/vtraag/leidenalg)
##' @param graph graph on which communities should be detected
##' @param resolution resolution parameter (default=1.0) - higher numbers lead to more communities
##' @param n.iterations number of iterations that the algorithm should be run for(default =2)
##' @return community object
##' @export
leiden.community <- function(graph, resolution=1.0, n.iterations=2) {

  x <- leiden_community(graph,E(graph)$weight,resolution,n.iterations);

  # enclose in a masquerading class
  fv <- as.factor(setNames(x,V(graph)$name))
  res <- list(membership=fv, dendrogram=NULL, algorithm='leiden', resolution=resolution,
              n.iter=n.iterations, names=names(fv));
  class(res) <- rev("fakeCommunities");
  return(res);
}

##' recursive leiden communities
##'
##' Constructrs a n-step recursive clustering, using leiden.communities
##' @param graph graph
##' @param n.cores number of cores to use
##' @param max.depth recursive depth
##' @param min.community.size minimal community size parameter for the walktrap communities .. communities smaller than that will be merged
##' @param verbose whether to output progress messages
##' @param resolution resolution parameter passed to leiden.communities (either a single value, or a value equivalent to max.depth)
##' @param ... passed to leiden.communities
##' @return a fakeCommunities object that has methods membership() ... does not return a dendrogram ... see cltrap.community() to constructo that
##' @export
rleiden.community <- function(graph, max.depth=2, n.cores=parallel::detectCores(logical=FALSE), min.community.size=10, verbose=FALSE, resolution=1, cur.depth=1, hierarchical=TRUE, ...) {

  if(verbose & cur.depth==1) cat(paste0("running ",max.depth,"-recursive Leiden clustering: "));
  if(length(resolution)>1) {
    if(length(resolution)!=max.depth) { stop("resolution value must be either a single number or a vector of length max.depth")}
    res <- resolution[cur.depth]
  } else { res <- resolution }
  mt <- leiden.community(graph, resolution=res, ...);

  mem <- membership(mt);
  tx <- table(mem)
  ivn <- names(tx)[tx<min.community.size]
  if(length(ivn)>1) {
    mem[mem %in% ivn] <- as.integer(ivn[1]); # collapse into one group
  }
  if(verbose) cat(length(unique(mem)),' ');

  if(cur.depth<max.depth) {
    # start recursive run
    wtl <- papply(conos:::sn(unique(mem)), function(cluster) {
      cn <- names(mem)[which(mem==cluster)]
      sg <- induced.subgraph(graph,cn)
      rleiden.community(induced.subgraph(graph,cn), max.depth=max.depth, resolution=resolution, cur.depth=cur.depth+1, min.community.size=min.community.size, hierarchical=hierarchical, verbose=verbose, n.cores=1, ...)
    },n.cores=n.cores)

    # merge clusters, cleanup
    mbl <- lapply(wtl,membership);
    # combined clustering factor
    fv <- unlist(lapply(conos:::sn(names(wtl)),function(cn) {
      paste(cn,as.character(mbl[[cn]]),sep='-')
    }))
    names(fv) <- unlist(lapply(mbl,names))
  } else {
    fv <- mem;
    if(hierarchical) {
      # use walktrap on the last level
      wtl <- conos:::papply(sn(unique(mem)), function(cluster) {
        cn <- names(mem)[which(mem==cluster)]
        sg <- induced.subgraph(graph,cn)
        res <- walktrap.community(induced.subgraph(graph,cn))
        res$merges <- igraph:::complete.dend(res,FALSE)
        res
      },n.cores=n.cores)
    }
  }

  if(hierarchical) {
    # calculate hierarchy on the multilevel clusters
    if(length(wtl)>1) {
      cgraph <- getClusterGraph(graph,mem)
      chwt <- walktrap.community(cgraph,steps=8)
      d <- as.dendrogram(chwt);

      # merge hierarchical portions
      wtld <- lapply(wtl,as.dendrogram)
      max.height <- max(unlist(lapply(wtld,attr,'height')))

      # shift leaf ids to fill in 1..N range
      mn <- unlist(lapply(wtld,attr,'members'))
      shift.leaf.ids <- function(l,v) { if(is.leaf(l)) { la <- attributes(l); l <- as.integer(l)+v; attributes(l) <- la; }; l  }
      nshift <- cumsum(c(0,mn))[-(length(mn)+1)]; names(nshift) <- names(mn); # how much to shift ids in each tree

      get.heights <- function(l) {
        if(is.leaf(l)) {
          return(attr(l,'height'))
        } else {
          return(c(attr(l,'height'),unlist(lapply(l,get.heights))))
        }
      }
      min.d.height <- min(get.heights(d))
      height.scale <- length(wtld)*2
      height.shift <- 2

      shift.heights <- function(l,s) { attr(l,'height') <- attr(l,'height')+s; l }

      glue.dends <- function(l) {
        if(is.leaf(l)) {
          nam <- as.character(attr(l,'label'));
          id <- dendrapply(wtld[[nam]], shift.leaf.ids, v=nshift[nam])
          return(dendrapply(id,shift.heights,s=max.height-attr(id,'height')))

        }
        attr(l,'height') <- (attr(l,'height')-min.d.height)*height.scale + max.height + height.shift;
        l[[1]] <- glue.dends(l[[1]]); l[[2]] <- glue.dends(l[[2]])
        attr(l,'members') <- attr(l[[1]],'members') + attr(l[[2]],'members')
        return(l)
      }
      combd <- glue.dends(d)
    } else {
      combd <- as.dendrogram(wtl[[1]]);
    }
  } else {
    combd <- NULL;
  }

  if(cur.depth==1) {
    if(verbose) {
      cat(paste0(' detected a total of ',length(unique(fv)),' clusters '));
      cat("done\n");
    }
  }

  # enclose in a masquerading class
  res <- list(membership=fv, dendrogram=combd, algorithm='rleiden', names=names(fv));
  if(hierarchical & cur.depth==max.depth) {
    # reconstruct merges matrix
    hcm <- as.hclust(as.dendrogram(combd))$merge
    # translate hclust $merge to walktrap-like $merges
    res$merges <- hcm + nrow(hcm) + 1
    res$merges[hcm < 0] <- -hcm[hcm < 0] - 1
  }
  class(res) <- rev("fakeCommunities");
  return(res);
}


##' returns pre-calculated dendrogram
##'
##' @param obj fakeCommunities object
##' @param ... dropped
##' @return dendrogram
##' @export
as.dendrogram.fakeCommunities <- function(obj, ...) {
  return(obj$dendrogram)
}
##' returns pre-calculated membership factor
##'
##' @param obj fakeCommunities object
##' @return membership factor
##' @export
membership.fakeCommunities <- function(obj) {
  return(obj$membership)
}
