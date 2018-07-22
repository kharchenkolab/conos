

##' A function for quickly plotting collections and joint clustering
##'
##' @name Conos_plotPanel
##' @param filename optional name of a PNG file where to save the plot
##' @param panel.size width/height of an individual sample panel
##' @param shuffle.colors whether cluster colors should be randomly shuffled
##' @param embeddingType name of the embedding to use for each datasets (default='tSNE')
##' @param groups optional groups (factor) to show with colors (instead of res$groups)
##' @param gene name of the gene to show with colors
##' @param group.level.colors explicitly specify colors for group factor levels
##' @param n number of columns to arrange the samples in
##' @param mark.cluster.cex cex parameter for cluster labels
##' @param colors per-cell colors to show 
##' @param res results of the conosCluster
##' @return nothing
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
  
  fields=c('samples','pairs','graph','clusters','embedding','n.cores','misc'),

  methods = list(
    initialize=function(x, ..., n.cores=parallel::detectCores(logical=F), verbose=TRUE) {
      # # init all the output lists
      samples <<- list();
      pairs <<- list();
      graph <<- NULL;
      clusters <<- list();
      misc <<-list();
      n.cores <<- n.cores;

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
            xcp <- quickJNMF(samples[cis[,i]],k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average,maxiter=3e3)
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
          
          n12 <- n2CrossKnn(cpproj[[n1]],cpproj[[n2]],k,1,FALSE)
          n21 <- n2CrossKnn(cpproj[[n2]],cpproj[[n1]],k,1,FALSE)
          #colnames(n12) <- rownames(n21) <- rownames(cpproj[[n1]])
          #colnames(n21) <- rownames(n12) <- rownames(cpproj[[n2]])
          mnn <- drop0(n21*t(n12))
          mnn <- as(n21*t(n12),'dgTMatrix')
          return(data.frame('mA.lab'=rownames(cpproj[[n1]])[mnn@i+1],'mB.lab'=rownames(cpproj[[n2]])[mnn@j+1],'w'=pmax(1-mnn@x,0),stringsAsFactors=F))
          

        } else if (space=='genes') { 
          ## Overdispersed Gene space
          n12 <- n2CrossKnn(as.matrix(xl[[i]]$genespace1), as.matrix(xl[[i]]$genespace2),k,1,FALSE)
          n21 <- n2CrossKnn(as.matrix(xl[[i]]$genespace2), as.matrix(xl[[i]]$genespace1),k,1,FALSE)
          
          ##
          mnn <- drop0(n21*t(n12))
          mnn <- as(n21*t(n12),'dgTMatrix')

          ## return
          ret.df <- data.frame(
            'mA.lab'=rownames(xl[[i]]$genespace1)[mnn@i+1],
            'mB.lab'=rownames(xl[[i]]$genespace2)[mnn@j+1],
            'w'=pmax(1-mnn@x,0),stringsAsFactors=F
          )
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
      misc$el <<- el; # cache edge list, though we coudl get it back from the igraph object as well
      
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


    plotPanel=function(clustering=NULL, filename=NULL,panel.size=250,shuffle.colors=F,embeddingType='tSNE',groups=NULL,gene=NULL,group.level.colors=NULL,n=NULL,mark.cluster.cex=0.8,colors=NULL, alpha=0.2) {
      if(is.null(groups) && is.null(colors)) {
        if(length(clusters)<1) { stop("generate a joint clustering first") }
        if(is.null(clustering)) { # take the first one
          rg <- clusters[[1]]$groups
        } else {
          if(is.null(clusters[[clustering]])) { stop(paste("clustering",clustering,"hasn't been calculated")) }
          rg <- clusters[[clustering]]$groups
        }
      }
      conosPlot(samples,list(groups=rg),filename=filename,panel.size=panel.size,shuffle.colors=shuffle.colors,embeddingType=embeddingType,groups=groups,gene=gene,group.level.colors=group.level.colors,n=n,mark.cluster.cex=mark.cluster.cex,colors=colors)
    },

    embedGraph=function(method='largeVis',M=1,gamma=1,alpha=0.1,perplexity=50,sgd_batches=1e8,seed=1,verbose=TRUE) {
      if(method!='largeVis') { stop("currently, only largeVis embeddings are supported") }
      wij <- largeVis:::buildWijMatrix(as_adj(graph,attr='weight'),perplexity=perplexity,threads=n.cores)
      coords <- largeVis:::projectKNNs(wij = wij, dim=2, verbose = verbose,sgd_batches = sgd_batches,gamma=gamma, M=M, seed=seed, alpha=alpha, rho=1, threads=n.cores)
      colnames(coords) <- V(graph)$name
      embedding <<- coords;
      return(invisible(embedding))
    },

    plotGraph=function(color.by='cluster',clustering=NULL,groups=NULL,colors=NULL,alpha=0.2,do.par=TRUE) {
      if(is.null(embedding)) {
        embedGraph();
      }
      if(do.par) par(mfrow=c(1,1), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
      if(is.null(groups) && is.null(colors)) {
        if(color.by=='cluster') {
          if(length(clusters)<1) { stop("generate a joint clustering first") }
          if(is.null(clustering)) { # take the first one
            groups <- clusters[[1]]$groups
          } else {
            if(is.null(clusters[[clustering]])) { stop(paste("clustering",clustering,"hasn't been calculated")) }
            groups <- clusters[[clustering]]$groups
          }
        } else if(color.by=='sample') {
          # list all the cells in the samples
          cl <- lapply(samples,function(x) rownames(x$counts))
          groups <- rep(names(cl),unlist(lapply(cl,length)))
          names(groups) <- unlist(cl);
          groups <- as.factor(groups);
        } else {
          stop('supported values of color.by are ("cluster" and "sample")')
        }
      }
      coords <- embedding;
      spal <- pagoda2:::fac2col(groups,return.details=T)
      plot(t(coords),pch=19,col=adjustcolor(spal$colors[colnames(coords)],alpha=0.1),cex=0.5,axes=F,panel.first=grid()); box()
      cent.pos <- do.call(rbind,tapply(1:ncol(coords),groups,function(ii) apply(coords[,ii,drop=F],1,median)))
      #rownames(cent.pos) <- levels(groups);
      cent.pos <- na.omit(cent.pos);
      text(cent.pos[,1],cent.pos[,2],labels=rownames(cent.pos),cex=0.8)
    }
  )
    
);
