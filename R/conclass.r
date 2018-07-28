

##' A function for quickly plotting collections and joint clustering
##'
##' @name Conos_plotPanel
##'
##' @inheritParams getClusteringGroups
##' @inherit plotPagodas params return
NULL

# initialize or append samples to the panel
##' initialize or add a set of samples to the conos panel
##'
##' note: this will simply add samples, but will not update graph, clustering, etc.
##' @name Conos_addSamples
##' @param x a named list of pagoda2 objects (one per sample)
##' @param replace whether the existing samples should be purged before adding new ones
##' @return invisible view of the full sample list
NULL

##' find joint communities
##'
##' @name Conos_findCommunities
##' @param method community detection method (igraph syntax)
##' @param min.group.size minimal allowed community size
##' @param name optional name of the clustering result (will default to the algorithm name)
##' @param ... extra parameters are passed to the specified community detection method
##' @return invisible list containing identified communities (groups) and the full community detection result (result)
NULL

#' Conos reference class
#'
#' The class encompasses sample collections, providing methods for calculating and visualizing joint graph and communities.
#' @title Conos reference class
#' @import methods
#' @export Conos
#' @exportClass Conos
Conos <- setRefClass(
  "Conos",

  # overall data handling approach:
  #   - samples list will be used to store sample specific-results (either in p2$misc$ locations or in a future wrapper of Seurat objects
  #   - pairs will contain pairwise alignment results (potentially for multiple reduction methods)
  #   - graph - we will support only one graph for now
  #   - clusters will contain potentially multiple clustering results/data
  #   - embedding contains embedding of the joint graph

  fields=c('samples','pairs','graph','clusters','embedding','n.cores','misc', 'override.conos.plot.theme'),

  methods = list(
    initialize=function(x, ..., n.cores=parallel::detectCores(logical=F), verbose=TRUE, override.conos.plot.theme=FALSE) {
      # # init all the output lists
      samples <<- list();
      pairs <<- list();
      graph <<- NULL;
      clusters <<- list();
      misc <<-list();
      n.cores <<- n.cores;
      override.conos.plot.theme <<- override.conos.plot.theme;

      if(!missing(x) && class(x)=='Conos') { # copy constructor
        callSuper(x, ..., n.cores=n.cores);
      } else {
        callSuper(..., n.cores=n.cores);
        if(!missing(x)) { # interpret x as a list of samples
          if(!is.list(x)) {
            stop("x is not a list of pagoda2 or Seurat objects")
          }

          if (class(x[[1]]) == 'Pagoda2') { # Pagoda2
            addSamples(x);
          } else {
            stop("only pagoda2 result lists are currently supported");
          }
        }
      }
    },

    adjustTheme=function(theme) {
      if (is.null(theme)) {
        theme <- ggplot2::theme()
      }
      main.theme <- ggplot2::theme_bw() + ggplot2::theme(
        legend.background=ggplot2::element_rect(fill=ggplot2::alpha("white", 0.6)),
        plot.margin=ggplot2::margin()
      )

      if (override.conos.plot.theme) {
        return(main.theme + ggplot2::theme_get() + theme)
      }

      return(main.theme + theme)
    },

    addSamples=function(x, replace=FALSE, verbose=FALSE) {
      # check names
      if(is.null(names(x))) {
        stop("the sample list must be named")
      }
      if(replace || length(samples)==0) {
        if(length(x)<2) {
          stop("provided list contains less than 2 samples; 2 required, >3 recommended")
        }
      }
      if(any(names(x) %in% names(samples))) {
        stop("some of the names in the provided sample list are identical to already existing samples")
      }

      # TODO: package-independent wrapper
      samples <<- c(samples,x);
    },

    updatePairs=function(space='CPCA',ncomps=50,n.odgenes=1e3,var.scale=TRUE,neighborhood.average=FALSE,neighborhood.average.k=10, verbose=FALSE) {
      if(neighborhood.average) {
        # pre-calculate averaging matrices for each sample
        if(verbose)  cat("calculating local averaging neighborhoods ")
        lapply(samples,function(r) {
          # W: get PCA reduction
          if(is.null(r$misc$edgeMat$quickCPCA) || r$misc$edgeMat$quickCPCAk != neighborhood.average.k) {
            xk <- n2Knn(r$reductions$PCA[rownames(r$counts),],neighborhood.average.k,n.cores,FALSE)
            xk@x <- pmax(1-xk@x,0);
            diag(xk) <- 1;
            xk <- t(t(xk)/colSums(xk))
            colnames(xk) <- rownames(xk) <- rownames(r$counts)
            # W: store averaging neighborhoods
            r$misc$edgeMat$quickCPCA <- xk;
            r$misc$edgeMat$quickCPCAk <- neighborhood.average.k;
            if(verbose) cat(".")
          }
        })
        if(verbose) cat(" done\n")
      }

      # make a list of all pairs
      cis <- combn(names(samples),2);
      # TODO: add random subsampling for very large panels

      # determine the pairs that need to be calculated
      if(is.null(pairs[[space]])) { pairs[[space]] <<- list() }
      mi <- rep(NA,ncol(cis));
      nm <- match(apply(cis,2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      # try reverse match as well
      nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      if(verbose) cat('found',sum(!is.na(mi)),'out of',length(mi),'cached',space,' space pairs ... ')
      if(any(is.na(mi))) { # some pairs are missing
        if(verbose) cat('running',sum(is.na(mi)),'additional',space,' space pairs ')
        xl2 <- papply(which(is.na(mi)), function(i) {
          if(space=='CPCA') {
            xcp <- quickCPCA(samples[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average)
          } else if(space=='JNMF') {
            xcp <- quickJNMF(samples[cis[,i]],n.comps=ncomps,n.odgenes=n.odgenes,var.scale=var.scale,verbose=FALSE,max.iter=3e3,neighborhood.average=neighborhood.average)
          } else if (space == 'genes') {
            xcp <- quickNULL(p2.objs = samlpes[cis[,i]], n.odgenes = n.odgenes, var.scale = var.scale, verbose = FALSE, neighborhood.average=neighborhood.average);
          }
          if(verbose) cat('.')
          xcp
        },n.cores=n.cores);

        names(xl2) <- apply(cis[,which(is.na(mi)),drop=F],2,paste,collapse='.vs.');
        pairs[[space]] <<- c(pairs[[space]],xl2);
      }

      # re-do the match and order
      mi <- rep(NA,ncol(cis));
      nm <- match(apply(cis,2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      nm <- match(apply(cis[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      if(any(is.na(mi))) { stop("unable to get complete set of pair comparison results") }
      pairs[[space]] <<- pairs[[space]][mi]
      if(verbose) cat(" done\n");
      return(invisible(cis))
    },

    buildGraph=function(k=30, k.self=10, k.self.weight=0.1, space='CPCA', matching.method='mNN', var.scale =TRUE, ncomps=50, n.odgenes=1000, return.details=T,neighborhood.average=FALSE,neighborhood.average.k=10, common.centering=TRUE , verbose=TRUE) {

      supported.spaces <- c("CPCA","JNMF","genes")
      if(!space %in% supported.spaces) {
        stop(paste0("only the following spaces are currently supported: [",paste(supported.spaces),"]"))
      }

      supported.matching.methods <- c("mNN","NN");
      if(!matching.method %in% supported.matching.methods) {
        stop(paste0("only the following matching methods are currently supported: [",paste(supported.matching.methods),"]"))
      }
      

      # calculate or update pairwise alignments
      cis <- updatePairs(space=space,ncomps=ncomps,n.odgenes=n.odgenes,verbose=verbose,var.scale=var.scale,neighborhood.average=neighborhood.average,neighborhood.average.k=10)


      # determine inter-sample mapping
      if(verbose) cat('inter-sample links using ',matching.method,' ');
      xl <- pairs[[space]]
      mnnres <- papply(1:ncol(cis), function(i) {
        r.ns <- samples[cis[,i]]
        if(space=='JNMF') {
          n12 <- n2CrossKnn(xl[[i]]$rot1,xl[[i]]$rot2,k,1,FALSE)
          n21 <- n2CrossKnn(xl[[i]]$rot2,xl[[i]]$rot1,k,1,FALSE)
          mnn <- drop0(n21*t(n12))
          mnn <- as(n21*t(n12),'dgTMatrix')
          if(verbose) cat(".")

          return(data.frame('mA.lab'=rownames(xl[[i]]$rot1)[mnn@i+1],'mB.lab'=rownames(xl[[i]]$rot2)[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
          #return(data.frame('mA.lab'=rownames(xl[[i]]$rot1)[mnn@i+1],'mB.lab'=rownames(xl[[i]]$rot2)[mnn@j+1],'w'=1/pmax(1,log(mnn@x)),stringsAsFactors=F))

        } else if (space %in% c("CPCA","GSVD")) {
          # W: get counts
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
            # W: get variance scaling
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
              # W: get neighborhood averaging matrix
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
          if(verbose) cat(".")

          mnn <- get.neighbor.matrix(cpproj[[n1]],cpproj[[n2]],k,matching=matching.method)
          return(data.frame('mA.lab'=rownames(mnn)[mnn@i+1],'mB.lab'=colnames(mnn)[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
          
        } else if (space=='genes') {
          ## Overdispersed Gene space
          mnn <- get.neighbor.matrix(as.matrix(xl[[i]]$genespace1), as.matrix(xl[[i]]$genespace2),k,matching=matching.method)
          return(data.frame('mA.lab'=rownames(mnn)[mnn@i+1],'mB.lab'=colnames(mnn)[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
        }
        mnnres
      },n.cores=n.cores)
      if(verbose) cat(" done\n")
      ## Merge the results into a edge table
      el <- do.call(rbind,mnnres)
      #el <- el$type <- 1; # encode connection type 1- intersample, 0- intrasample

      # append some local edges
      if(k.self>0) {
        if(verbose) cat('local pairs ')
        x <- data.frame(do.call(rbind,papply(samples,function(x) {
          # W: get PCA reduction
          xk <- n2Knn(x$reductions$PCA,k.self,1,FALSE)
          diag(xk) <- 0;
          xk <- as(xk,'dgTMatrix')
          cat(".")
          return(data.frame('mA.lab'=rownames(x$reductions$PCA)[xk@i+1],'mB.lab'=rownames(x$reductions$PCA)[xk@j+1],'w'=pmax(1-xk@x,0),stringsAsFactors=F))
        },n.cores=n.cores)),stringsAsFactors = F)
        x$w <- k.self.weight
        #x$type <- 0;
        cat(' done\n')
        el <- rbind(el,x)
      }
      # TODO: add edge labels to distinguish inter and intra-sample edges
      g  <- graph_from_edgelist(as.matrix(el[,c(1,2)]), directed =FALSE)
      E(g)$weight <- el[,3]

      graph <<- g;
      return(invisible(g))
    },



    findCommunities=function(method=multilevel.community, min.group.size=0, name=NULL, ...) {

      cls <- method(graph, ...)
      if(is.null(name)) {
        name <- cls$algorithm;
        if(is.null(name)) {
          name <- "community";
        }
      }

      ## Extract groups from this graph
      cls.mem <- membership(cls)
      cls.groups <- as.character(cls.mem)
      names(cls.groups) <- names(cls.mem)

      ## Filter groups
      if(min.group.size>0) {
        lvls.keep <- names(which(table(cls.groups)  > min.group.size))
        cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
        cls.groups <- as.factor(cls.groups)
      }
      res <- list(groups=cls.groups,result=cls)
      clusters[[name]] <<- res;
      return(invisible(res))

    },

    plotPanel=function(clustering=NULL, groups=NULL, colors=NULL, gene=NULL, embedding.type='tSNE', ncol=NULL, nrow=NULL, raster=FALSE, panel.size=NULL,
                       adjust.func=NULL, use.local.clusters=FALSE, plot.theme=NULL, ...) {
      if (use.local.clusters) {
        if (is.null(clustering)) {
          stop("You have to provide 'clustering' parameter to be able to use local clusters")
        }

        groups <- Reduce(c, lapply(samples, function(x) x$clusters$PCA[[clustering]]))
        if (is.null(groups)) {
          stop(paste0("No clustering '", clustering, "' presented in the samples"))
        }
      }
      else if (is.null(groups) && is.null(colors) && is.null(gene)) {
        groups <- getClusteringGroups(clusters, clustering)
      }

      gg <- plotPagodas(samples, groups=groups, colors=colors, gene=gene, embedding.type=embedding.type, ncol=ncol, nrow=nrow, raster=raster,
                        panel.size=panel.size, adjust.func=adjust.func, plot.theme=adjustTheme(plot.theme), ...)
      return(gg)
    },

    embedGraph=function(method='largeVis',M=1,gamma=1,alpha=0.1,perplexity=50,sgd_batches=1e8,seed=1,verbose=TRUE) {
      if(method!='largeVis') { stop("currently, only largeVis embeddings are supported") }
      wij <- largeVis:::buildWijMatrix(as_adj(graph,attr='weight'),perplexity=perplexity,threads=n.cores)
      coords <- largeVis:::projectKNNs(wij = wij, dim=2, verbose = verbose,sgd_batches = sgd_batches,gamma=gamma, M=M, seed=seed, alpha=alpha, rho=1, threads=n.cores)
      colnames(coords) <- V(graph)$name
      embedding <<- coords;
      return(invisible(embedding))
    },

    plotGraph=function(conos, color.by='cluster', clustering=NULL, groups=NULL, colors=NULL, plot.theme=NULL, ...) {
      if(class(embedding)[1] == "uninitializedField") {
        embedGraph();
      }

      if(is.null(groups) && is.null(colors)) {
        if(color.by == 'cluster') {
          groups <- getClusteringGroups(clusters, clustering)
        } else if(color.by == 'sample') {
          cl <- lapply(samples, function(x) rownames(x$counts))
          groups <- rep(names(cl), sapply(cl, length)) %>% setNames(unlist(cl)) %>% as.factor()
        } else {
          stop('supported values of color.by are ("cluster" and "sample")')
        }
      }

      return(embeddingPlot(t(con$embedding), groups=groups, colors=colors, plot.theme=adjustTheme(plot.theme), ...))
    }
  )

);


## helper functions
##' Establish rough neighbor matching between samples given their projections in a common space
##'
##' @param p1 projection of sample 1
##' @param p2 projection of sample 2
##' @param k neighborhood radius
##' @param matching mNN (default) or NN
##' @param indexType distance type (default: "angular")
##' @return matrix 
get.neighbor.matrix <- function(p1,p2,k,matching='mNN',indexType='angular') {
  n12 <- n2CrossKnn(p1,p2,k,1,FALSE)
  n21 <- n2CrossKnn(p2,p1,k,1,FALSE)
  # Viktor's solution
  n12@x[n12@x<0] <- 0
  n21@x[n21@x<0] <- 0
  
  if(matching=='NN') {
    mnn <- n21+t(n12);
    mnn@x <- mnn@x/2;
  } else { # mNN
    mnn <- drop0(n21*t(n12))
    mnn@x <- sqrt(mnn@x)
  }
  mnn <- as(mnn,'dgTMatrix')
  rownames(mnn) <- rownames(p1); colnames(mnn) <- rownames(p2);
  mnn
}


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
##' @param ... passed to walktrap
##' @return a fakeCommunities object that has methods membership() and as.dendrogram() to mimic regular igraph returns
##' @export
multitrap.community <- function(graph, n.cores=parallel::detectCores(logical=F), hclust.link='single', min.community.size=10, verbose=FALSE, ...) {
  graph <- con$graph
  if(verbose) cat("running multilevel ... ");
  mt <- multilevel.community(graph);

  # get the highest level
  mem <- mt$memberships[1,]; names(mem) <- mt$names;
  if(verbose) cat("found",length(unique(mem)),"communities\nrunning walktraps ... ")
  
  # calculate hierarchy on the multilevel clusters
  cgraph <- get.cluster.graph(graph,mem)
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

  if(verbose) cat("found",sum(unlist(lapply(wtl,function(x) length(unique(x))))),"communities\nmerging dendrograms ... ")
  
  
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
  res <- list(membership=fv,dendrogram=combd,algorithm='multitrap');
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



##' Collapse vertices belonging to each cluster in a graph
##'
##' @param graph graph to be collapsed
##' @param groups factor on vertives describing cluster assignment (can specify integer vertex ids, or character vertex names which will be matched)
##' @param plot whether to show collapsed graph plot
##' @return collapsed graph
##' @export
get.cluster.graph <- function(graph,groups,plot=FALSE,node.scale=50,edge.scale=50,edge.alpha=0.3) {
  if(plot) V(graph)$num <- 1;
  gcon <- contract.vertices(graph,groups,vertex.attr.comb=list('num'='sum',"ignore"))
  gcon <- simplify(gcon, edge.attr.comb=list(weight="sum","ignore"))
  gcon <- induced.subgraph(gcon, unique(groups))
  
  if(plot) {
    set.seed(1)
    par(mar = rep(0.1, 4)) 
    plot.igraph(gcon, layout=layout_with_fr(gcon), vertex.size=V(gcon)$num/(sum(V(gcon)$num)/node.scale), edge.width=E(gcon)$weight/sum(E(gcon)$weight/edge.scale), edge.color=adjustcolor('black',alpha=edge.alpha))
  }
  return(invisible(gcon))
}


# TODO: p2 abstraction

