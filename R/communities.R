#' @importFrom graphics par
#' @importFrom grDevices adjustcolor colorRampPalette dev.size rainbow
#' @importFrom igraph E induced.subgraph membership V walktrap.community
#' @import Matrix.utils
#' @import Matrix
#' @import parallel
#' @import sccore
#' @importFrom stats as.dendrogram as.dist as.hclust cor dendrapply dnorm hclust is.leaf median 
#' @importFrom stats na.omit order.dendrogram p.adjust prcomp qnorm quantile relevel sd setNames var
#' @importFrom utils combn packageVersion write.table
NULL

#' Constructs a two-step clustering, first running multilevel.communities, and then walktrap.communities within each
#' These are combined into an overall hierarchy
#'
#' @param graph graph
#' @param n.cores numeric Number of cores to use (default=parallel::detectCores(logical=FALSE))
#' @param hclust.link character Link function to use when clustering multilevel communities (based on collapsed graph connectivity) (default='single')
#' @param min.community.size numeric Minimal community size parameter for the walktrap communities .. communities smaller than that will be merged (default=10)
#' @param verbose boolean Whether to output progress messages (default=FALSE)
#' @param level numeric What level of multitrap clustering to use in the starting step. By default, uses the top level. An integer can be specified for a lower level (i.e. 1) (default=NULL)
#' @param ... passed to walktrap
#' @return a fakeCommunities object that has methods membership() and as.dendrogram() to mimic regular igraph returns
#' @keywords internal
multitrap.community <- function(graph, n.cores=parallel::detectCores(logical=FALSE), hclust.link='single', min.community.size=10, verbose=FALSE, level=NULL, ...) {
  .Deprecated()

  if(verbose) message("running multilevel ... ");
  mt <- multilevel.community(graph);

  if(is.null(level)) {
    # get higest level (to avoid oversplitting at the initial step)
    mem <- membership(mt);
  } else {
    # get the specified level
    mem <- mt$memberships[level,]; names(mem) <- mt$names;
  }

  if(verbose) message("found ",length(unique(mem))," communities\nrunning walktraps ... ")

  # calculate hierarchy on the multilevel clusters
  cgraph <- getClusterGraph(graph,mem)
  chwt <- walktrap.community(cgraph,steps=8)
  d <- as.dendrogram(chwt);

  wtl <- papply(sn(unique(mem)), function(cluster) {
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

  if(verbose) message("found ",sum(unlist(lapply(mbl,function(x) length(unique(x)))))," communities\nmerging dendrograms ... ")

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
  if(verbose) message("done\n");

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


### mutlilevel+multilevel communities

#' Constructrs a two-step clustering, first running multilevel.communities, and then walktrap.communities within each
#' These are combined into an overall hierarchy
#'
#' @param graph graph
#' @param n.cores numeric Number of cores to use (default=parallel::detectCores(logical=FALSE))
#' @param hclust.link character Link function to use when clustering multilevel communities (based on collapsed graph connectivity) (default='single')
#' @param min.community.size numeric Minimal community size parameter for the walktrap communities .. communities smaller than that will be merged (default=10)
#' @param verbose boolean Whether to output progress messages (default=FALSE)
#' @param level numeric What level of multitrap clustering to use in the starting step. By default, uses the top level. An integer can be specified for a lower level (i.e. 1) (default=NULL)
#' @param ... arguments passed to walktrap
#' @return a fakeCommunities object that has methods membership() and as.dendrogram() to mimic regular igraph returns
#' @keywords internal
multimulti.community <- function(graph, n.cores=parallel::detectCores(logical=FALSE), hclust.link='single', min.community.size=10, verbose=FALSE, level=NULL, ...) {
  .Deprecated()

  if(verbose) message("running multilevel 1 ... ");
  mt <- multilevel.community(graph);

  if(is.null(level)) {
    # get higest level (to avoid oversplitting at the initial step)
    mem <- membership(mt);
  } else {
    # get the specified level
    mem <- mt$memberships[level,]; names(mem) <- mt$names;
  }

  if(verbose) message("found ",length(unique(mem))," communities\nrunning multilevel 2 ... ")

  # calculate hierarchy on the multilevel clusters
  cgraph <- getClusterGraph(graph,mem)
  chwt <- walktrap.community(cgraph,steps=8)
  d <- as.dendrogram(chwt);


  wtl <- papply(sn(unique(mem)), function(cluster) {
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

  if(verbose) message("found ",sum(unlist(lapply(mbl,function(x) length(unique(x)))))," communities\nmerging ... ")

  # combined clustering factor
  fv <- unlist(lapply(sn(names(wtl)),function(cn) {
    paste(cn,as.character(mbl[[cn]]),sep='-')
  }))
  names(fv) <- unlist(lapply(mbl,names))

  # enclose in a masquerading class
  res <- list(membership=fv, dendrogram=NULL, algorithm='multimulti', names=names(fv));
  class(res) <- rev("fakeCommunities")
  return(res)

}

