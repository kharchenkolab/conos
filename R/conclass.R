##' A function for quickly plotting collections and joint clustering
##'
##' @name Conos_plotPanel
##'
##' @inheritParams getClusteringGroups
##' @inherit plotSamples params return
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

  fields=c('samples','pairs','graph','clusters', 'expression.adj','embedding','n.cores','misc', 'override.conos.plot.theme'),

  methods = list(
    initialize=function(x, ..., n.cores=parallel::detectCores(logical=F), verbose=TRUE, override.conos.plot.theme=FALSE) {
      # # init all the output lists
      samples <<- list();
      pairs <<- list();
      graph <<- NULL;
      clusters <<- list();
      expression.adj <<- list();
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
          if (inherits(x = x[[1]], what = c('Pagoda2', 'seurat', 'Seurat'))) {
            addSamples(x);
          } else {
            stop("only Pagoda2 or Seurat result lists are currently supported");
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
      if(any(duplicated(names(samples)))) { stop("duplicate names found in the supplied samples") }
      if(any(names(x) %in% names(samples))) {
        stop("some of the names in the provided sample list are identical to already existing samples")
      }

      # TODO: package-independent wrapper
      samples <<- c(samples,x);
    },

    updatePairs=function(space='PCA', data.type='counts', ncomps=50, n.odgenes=1e3, var.scale=TRUE, neighborhood.average=FALSE, neighborhood.average.k=10, matching.mask=NULL, exclude.samples=NULL, score.component.variance=FALSE, verbose=FALSE) {
      if(neighborhood.average) {
        # pre-calculate averaging matrices for each sample
        if(verbose)  cat("calculating local averaging neighborhoods ")
        for (n in names(samples)) {
          r <- samples[[n]]
          if(!is.null(edgeMat(r)$mat) && edgeMat(r)$k != neighborhood.average.k)
            next

          xk <- n2Knn(getPca(r)[getCellNames(r),],neighborhood.average.k,n.cores,FALSE)
          xk@x <- pmax(1-xk@x,0);
          diag(xk) <- 1;
          xk <- t(t(xk)/colSums(xk))
          colnames(xk) <- rownames(xk) <- getCellNames(r)
          edgeMat(samples[[n]]) <<- list(mat=xk, k=neighborhood.average.k)
          if(verbose) cat(".")
        }
        if(verbose) cat(" done\n")
      }

      # make a list of all pairs
      sample.names <- names(samples);
      if(!is.null(exclude.samples)) {
        mi <- sample.names %in% exclude.samples;
        if(verbose) { cat("excluded", sum(mi), "out of", length(sample.names), "samples, based on supplied exclude.samples\n") }
        sample.names <- sample.names[!mi];
      }

      # TODO: add random subsampling for very large panels
      if(!is.null(matching.mask)) { # remove pairs that shouldn't be compared directly
        tryCatch(matching.mask <- matching.mask[sample.names, sample.names],
                 error=function(e) stop("matching.mask should have the same row- and colnames as provided samples. Error:", e))

        matching.mask <- matching.mask | t(matching.mask)
        selected.ids <- which(lower.tri(matching.mask) & matching.mask) - 1
        sn.pairs <- sample.names[selected.ids %/% length(sample.names) + 1] %>%
          cbind(sample.names[selected.ids %% length(sample.names) + 1]) %>%
          t()

        if(verbose) cat("Use", ncol(sn.pairs), "pairs, based on the passed exclude.pairs\n")
      } else {
        sn.pairs <- combn(sample.names, 2);
      }

      # determine the pairs that need to be calculated
      if(is.null(pairs[[space]])) { pairs[[space]] <<- list() }
      mi <- rep(NA,ncol(sn.pairs));
      nm <- match(apply(sn.pairs,2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      # try reverse match as well
      nm <- match(apply(sn.pairs[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      if(verbose) cat('found',sum(!is.na(mi)),'out of',length(mi),'cached',space,' space pairs ... ')
      if(any(is.na(mi))) { # some pairs are missing
        if(verbose) cat('running',sum(is.na(mi)),'additional',space,' space pairs ')
        xl2 <- papply(which(is.na(mi)), function(i) {
          if(space=='CPCA') {
            xcp <- quickCPCA(samples[sn.pairs[,i]],data.type=data.type,k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average, score.component.variance=score.component.variance)
          } else if(space=='JNMF') {
            xcp <- quickJNMF(samples[sn.pairs[,i]],data.type=data.type,n.comps=ncomps,n.odgenes=n.odgenes,var.scale=var.scale,verbose=FALSE,max.iter=3e3,neighborhood.average=neighborhood.average)
          } else if (space == 'genes') {
            xcp <- quickNULL(p2.objs = samples[sn.pairs[,i]], data.type=data.type, n.odgenes=n.odgenes, var.scale = var.scale, verbose = FALSE, neighborhood.average=neighborhood.average);
          } else if (space == 'PCA') {
            xcp <- quickPlainPCA(samples[sn.pairs[,i]], data.type=data.type, k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average, score.component.variance=score.component.variance)
          } else if (space == 'CCA' || space=='PMA') {
            xcp <- quickCCA(samples[sn.pairs[,i]],data.type=data.type,k=k,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale,neighborhood.average=neighborhood.average, score.component.variance=score.component.variance,PMA=(space=='PMA'))
          }
          if(verbose) cat('.')
          xcp
        },n.cores=n.cores,mc.preschedule=(space=='PCA'));

        names(xl2) <- apply(sn.pairs[,which(is.na(mi)),drop=F],2,paste,collapse='.vs.');
        xl2 <- xl2[!unlist(lapply(xl2,is.null))]
        pairs[[space]] <<- c(pairs[[space]],xl2);
      }

      # re-do the match and order
      mi <- rep(NA,ncol(sn.pairs));
      nm <- match(apply(sn.pairs,2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      nm <- match(apply(sn.pairs[c(2,1),,drop=F],2,paste,collapse='.vs.'),names(pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm);
      if(any(is.na(mi))) {
        warning("unable to get complete set of pair comparison results")
        sn.pairs <- sn.pairs[,!is.na(mi),drop=FALSE]
      }
      if(verbose) cat(" done\n");
      return(invisible(sn.pairs))
    },


    buildGraph=function(k=15, k.self=10, k.self.weight=0.1, alignment.strength=NULL, space='PCA', matching.method='mNN', metric='angular', k1=k, data.type='counts', l2.sigma=1e5, var.scale =TRUE, ncomps=40, n.odgenes=2000, neighborhood.average=FALSE, neighborhood.average.k=10, matching.mask=NULL, exclude.samples=NULL, common.centering=TRUE, verbose=TRUE, base.groups=NULL, append.global.axes=TRUE, append.decoys=TRUE, decoy.threshold=1, n.decoys=k*2, score.component.variance=FALSE, balance.edge.weights=FALSE, balancing.factor.per.cell=NULL, same.factor.downweight=1.0) {
      supported.spaces <- c("CPCA","JNMF","genes","PCA","PMA","CCA")
      if(!space %in% supported.spaces) {
        stop(paste0("only the following spaces are currently supported: [",paste(supported.spaces,collapse=' '),"]"))
      }
      supported.matching.methods <- c("mNN", "NN");
      if(!matching.method %in% supported.matching.methods) {
        stop(paste0("only the following matching methods are currently supported: ['",paste(supported.matching.methods,collapse="' '"),"']"))
      }
      supported.metrics <- c("L2","angular");
      if(!metric %in% supported.metrics) {
        stop(paste0("only the following distance metrics are currently supported: ['",paste(supported.metrics,collapse="' '"),"']"))
      }
      if (!is.null(alignment.strength)) {
        alignment.strength %<>% max(0) %>% min(1)
        k1 <- sapply(samples, function(sample) ncol(getCountMatrix(sample))) %>% max() %>%
          `*`(alignment.strength ^ 2) %>% round() %>% max(k)
      }
      if(k1<k) { stop("k1 must be >= k") }
      # calculate or update pairwise alignments
      sn.pairs <- updatePairs(
        space=space, ncomps=ncomps, n.odgenes=n.odgenes, verbose=verbose, var.scale=var.scale, neighborhood.average=neighborhood.average,
        neighborhood.average.k=10, matching.mask=matching.mask, exclude.samples=exclude.samples, score.component.variance=score.component.variance
      )
      if(ncol(sn.pairs)<1) { stop("insufficient number of comparable pairs") }
      if(!is.null(base.groups)) {
        samf <- lapply(samples,getCellNames)
        base.groups <- as.factor(base.groups[names(base.groups) %in% unlist(samf)]) # clean up the group factor
        if(length(base.groups) < 2) stop("provided base.gropus doesn't cover enough cells")
        # make a sample factor
        samf <- setNames(rep(names(samf),unlist(lapply(samf,length))),unlist(samf))
        if(append.global.axes) {
          cms.clust <- getClusterCountMatrices(groups=base.groups,common.genes=FALSE)
          global.proj <- projectSamplesOnGlobalAxes(samples, cms.clust, data.type, neighborhood.average, verbose, n.cores)
        }
      }
      # determine inter-sample mapping
      if(verbose) cat('inter-sample links using ',matching.method,' ');
      cached.pairs <- pairs[[space]]
      mnnres <- papply(1:ncol(sn.pairs), function(j) {
        # we'll look up the pair by name (possibly reversed), not to assume for the ordering of $pairs[[space]] to be the same
        i <- match(paste(sn.pairs[,j],collapse='.vs.'),names(cached.pairs));
        if(is.na(i)) { i <- match(paste(rev(sn.pairs[,j]),collapse='.vs.'),names(cached.pairs)) }
        if(is.na(i)) { stop(paste("unable to find alignment for pair",paste(sn.pairs[,j],collapse='.vs.'))) }

        if(space=='JNMF') {
          mnn <- getNeighborMatrix(cached.pairs[[i]]$rot1,cached.pairs[[i]]$rot2,k,matching=matching.method,metric=metric,l2.sigma=l2.sigma)
          if(verbose) cat(".")
          return(data.frame('mA.lab'=rownames(cached.pairs[[i]]$rot1)[mnn@i+1],'mB.lab'=rownames(cached.pairs[[i]]$rot2)[mnn@j+1],'w'=mnn@x,stringsAsFactors=F))
          #return(data.frame('mA.lab'=rownames(cached.pairs[[i]]$rot1)[mnn@i+1],'mB.lab'=rownames(cached.pairs[[i]]$rot2)[mnn@j+1],'w'=1/pmax(1,log(mnn@x)),stringsAsFactors=F))

        } else if (space %in% c("CPCA","GSVD","PCA")) {
          #common.genes <- Reduce(intersect,lapply(r.ns, getGenes))
          if(!is.null(cached.pairs[[i]]$CPC)) {
            # CPCA or PCA
            #od.genes <- intersect(rownames(cached.pairs[[i]]$CPC),common.genes)
            od.genes <- rownames(cached.pairs[[i]]$CPC);
            rot <- cached.pairs[[i]]$CPC[od.genes,];
          } else if(!is.null(cached.pairs[[i]]$o$Q)) {
            # GSVD
            rot <- cached.pairs[[i]]$o$Q;
            od.genes <- rownames(rot) <- colnames(cached.pairs[[i]]$o$A);
          } else {
            stop("unknown reduction provided")
          }

          # TODO: a more careful analysis of parameters used to calculate the cached version
          if(ncomps > ncol(rot)) {
            warning(paste0("specified ncomps (",ncomps,") is greater than the cached version (",ncol(rot),")"))
          } else {
            rot <- rot[,1:ncomps,drop=F]
          }

          mnn <- getPcaBasedNeighborMatrix(samples[sn.pairs[,j]], od.genes=od.genes, rot=rot, k=k, k1=k1, data.type=data.type,
                                           var.scale=var.scale, neighborhood.average=neighborhood.average, common.centering=common.centering,
                                           matching.method=matching.method, metric=metric, l2.sigma=l2.sigma, cor.base=1 + min(1, alignment.strength * 10),
                                           base.groups=base.groups, append.decoys=append.decoys, samples=samples, samf=samf, decoy.threshold=decoy.threshold,
                                           n.decoys=n.decoys, append.global.axes=append.global.axes, global.proj=global.proj)
          if(verbose) cat(".")

          return(data.frame('mA.lab'=rownames(mnn)[mnn@i+1],'mB.lab'=colnames(mnn)[mnn@j+1],'w'=mnn@x,stringsAsFactors=F))

        } else if (space=='genes') {
          ## Overdispersed Gene space
          mnn <- getNeighborMatrix(as.matrix(cached.pairs[[i]]$genespace1), as.matrix(cached.pairs[[i]]$genespace2),k,matching=matching.method,metric=metric,l2.sigma=l2.sigma)
          return(data.frame('mA.lab'=rownames(mnn)[mnn@i+1],'mB.lab'=colnames(mnn)[mnn@j+1],'w'=mnn@x,stringsAsFactors=F))
        } else if(space=='PMA' || space=='CCA') {
          mnn <- getNeighborMatrix(cached.pairs[[i]]$u,cached.pairs[[i]]$v,k,k1=k1,matching=matching.method,metric=metric,l2.sigma=l2.sigma,cor.base=1 + min(1, alignment.strength * 10));
          return(data.frame('mA.lab'=rownames(cached.pairs[[i]]$u)[mnn@i+1],'mB.lab'=rownames(cached.pairs[[i]]$v)[mnn@j+1],'w'=mnn@x,stringsAsFactors=F))
        }
        mnnres
      },n.cores=n.cores,mc.preschedule=TRUE)
      if(verbose) cat(" done\n")
      ## Merge the results into a edge table
      el <- do.call(rbind,mnnres)
      el$type <- 1; # encode connection type 1- intersample, 0- intrasample
      # append some local edges
      if(k.self>0) {
        if(verbose) cat('local pairs ')
        x <- getLocalEdges(samples, k.self, k.self.weight, metric, l2.sigma=l2.sigma, verbose, n.cores)
        el <- rbind(el,x)
      }
      if(verbose) cat('building graph .')
      g  <- graph_from_edgelist(as.matrix(el[,c(1,2)]), directed =FALSE)
      E(g)$weight <- el[,3]
      E(g)$type <- el[,4]
      if(verbose) cat('.')
      # collapse duplicate edges
      g <- simplify(g, edge.attr.comb=list(weight="sum", type = "first"))
      if(verbose) cat('done\n')
      if (balance.edge.weights || !is.null(balancing.factor.per.cell)) {
        if(verbose) cat('balancing edge weights ');

        if (is.null(balancing.factor.per.cell)) {
          balancing.factor.per.cell <- getDatasetPerCell()
        }

        g <- igraph::as_adjacency_matrix(g, attr="weight") %>%
          adjustWeightsByCellBalancing(factor.per.cell=balancing.factor.per.cell, balance.weights=balance.edge.weights,
                                       same.factor.downweight=same.factor.downweight) %>%
          igraph::graph_from_adjacency_matrix(mode="undirected", weighted=T)

        if(verbose) cat('done\n');
      }
      graph <<- g;
      return(invisible(g))
    },

    getDifferentialGenes=function(clustering=NULL, groups=NULL, z.threshold=3.0, upregulated.only=F, verbose=T, plot=FALSE, n.genes.to.show=10, inner.clustering=FALSE, n.cores=NULL) {
      if (!is.null(clustering)) {
        groups <- clusters[[clustering]]$groups
      }

      if (is.null(groups))
        stop("Either 'clustering' or 'groups' must be provided")

      if (class(samples[[1]]) != 'Pagoda2') # TODO: add Seurat
        stop("Only Pagoda onjects are supported for marker genes")

      n.jobs <- if (is.null(n.cores)) .self$n.cores else n.cores

      de.genes <- getDifferentialGenesP2(samples, groups=groups, z.threshold=z.threshold, upregulated.only=upregulated.only, verbose=verbose, n.cores=n.jobs)
      de.genes <- de.genes[levels(groups)]

      if(plot) {
        # genes to show
        vi <- unlist(lapply(de.genes,class))=='data.frame';
        x <- lapply(de.genes[vi],function(d) {  if(!is.null(d) && nrow(d)>0) { d[1:min(nrow(d),n.genes.to.show),] } else { NULL } })
        x <- lapply(x,rownames);
        genes <- unique(unlist(x))
        # make expression matrix
        cl <- lapply(samples,function(y) { m <- getCountMatrix(y); m[rownames(m) %in% genes,,drop=F] })
        # merge into a common matrix, entering 0s for the missing ones, convert to regular matrices
        ExtendMatrix <- function(mtx, col.names) {
          new.names <- setdiff(col.names, colnames(mtx))
          ext.mtx <- matrix(0, nrow=nrow(mtx), ncol=length(new.names))
          colnames(ext.mtx) <- new.names
          return(cbind(as.matrix(mtx), ext.mtx)[,col.names])
        }

        MergeCountMatrices <- function(cms) {
          cms <- lapply(cms, t)
          gene.union <- lapply(cms, colnames) %>% Reduce(union, .)

          res <- lapply(cms, ExtendMatrix, gene.union) %>% Reduce(rbind, .)
          return(Matrix::t(res))
        }

        em <- MergeCountMatrices(cl)
        # renormalize rows
        gradient.range.quantile <- 0.95; # make a parameter?
        if(all(sign(em)>=0)) {
          gradientPalette <- colorRampPalette(c('gray90','red'), space = "Lab")(1024)
          em <- t(apply(em,1,function(x) {
            zlim <- as.numeric(quantile(x,p=c(1-gradient.range.quantile,gradient.range.quantile)))
            if(diff(zlim)==0) {
              zlim <- as.numeric(range(x))
            }
            x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
            x <- (x-zlim[1])/(zlim[2]-zlim[1])
          }))
        } else {
          gradientPalette <- colorRampPalette(c("blue", "grey90", "red"), space = "Lab")(1024)
          em <- t(apply(em,1,function(x) {
            zlim <- c(-1,1)*as.numeric(quantile(abs(x),p=gradient.range.quantile))
            if(diff(zlim)==0) {
              zlim <- c(-1,1)*as.numeric(max(abs(x)))
            }
            x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
            x <- (x-zlim[1])/(zlim[2]-zlim[1])
          }))
        }

        # cluster cell types by averages
        gmap <- data.frame(gene=unlist(x),cl=rep(names(x),unlist(lapply(x,length))));
        rowfac <- factor(gmap$cl[match(rownames(em),gmap$gene)],levels=names(x))
        if(inner.clustering) {
          clclo <- hclust(as.dist(1-cor(do.call(cbind,tapply(1:nrow(em),rowfac,function(ii) Matrix::colMeans(em[ii,,drop=FALSE]))))),method='complete')$order
        } else {
          clclo <- 1:length(levels(rowfac))
        }

        if(inner.clustering) {
          # cluster genes within each cluster
          clgo <- tapply(1:nrow(em),rowfac,function(ii) {
            ii[hclust(as.dist(1-cor(t(em[ii,]))),method='complete')$order]
          })
        } else {
          clgo <- tapply(1:nrow(em),rowfac,I)
        }
        if(inner.clustering) {
          # cluster cells within each cluster
          clco <- tapply(1:ncol(em),groups[colnames(em)],function(ii) {
            if(length(ii)>3) {
              ii[hclust(as.dist(1-cor(em[,ii,drop=F])),method='complete')$order]
            } else {
              ii
            }
          })
        } else {
          clco <- tapply(1:ncol(em),groups[colnames(em)],I)
        }
        #clco <- clco[names(clgo)]
        # filter down to the clusters that are included
        #vic <- cols %in% clclo
        colors <- fac2col(groups[colnames(em)],v=0.95,s=0.95,return.details=TRUE)
        samf <- fac2col(getDatasetPerCell(),v=0.75,s=0.9,return.details=TRUE);
        cellcols <- colors$colors[unlist(clco[clclo])]
        samfcols <- samf$colors[unlist(clco[clclo])]
        genecols <- rev(rep(colors$palette,unlist(lapply(clgo,length)[clclo])))
        drawGroupNames <- FALSE;
        bottomMargin <- ifelse(drawGroupNames,4,0.5);

        browser()

        #pagoda2:::my.heatmap2(em[rev(unlist(clgo[clclo])),unlist(clco[clclo])],col=gradientPalette,Colv=NA,Rowv=NA,labRow=NA,labCol=NA,RowSideColors=genecols,ColSideColors=rbind(samfcols,cellcols),margins=c(bottomMargin,0.5),ColSideColors.unit.vsize=0.05,RowSideColors.hsize=0.05,useRaster=TRUE, box=TRUE)

        pagoda2:::my.heatmap2(em[rev(unlist(clgo[clclo])),unlist(clco[clclo])],col=gradientPalette,Colv=NA,Rowv=NA,labRow=NA,labCol=NA,RowSideColors=genecols,ColSideColors=rbind(samfcols,cellcols),margins=c(bottomMargin,0.5),ColSideColors.unit.vsize=0.05,RowSideColors.hsize=0.05,useRaster=TRUE, box=TRUE)
        abline(v=cumsum(unlist(lapply(clco[clclo],length))),col=1,lty=3)
        abline(h=cumsum(rev(unlist(lapply(clgo[clclo],length))))+0.5,col=1,lty=3)

    }


      return(de.genes)
    },

    findCommunities=function(method=leiden.community, min.group.size=0, name=NULL, test.stability=FALSE, stability.subsampling.fraction=0.95, stability.subsamples=100, verbose=TRUE, cls=NULL, sr=NULL, ...) {

      if(is.null(cls)) {
        cls <- method(graph, ...)
      }
      if(is.null(name)) {
        name <- cls$algorithm;
        if(is.null(name)) {
          name <- "community";
        }
      }

      ## Extract groups from this graph
      cls.mem <- membership(cls)
      cls.groups <- factor(setNames(as.character(cls.mem),names(cls.mem)),levels=1:max(cls.mem))
      cls.levs <- levels(cls.groups)
      res <- list(groups=cls.groups,result=cls)

      # test stability
      if(test.stability) {
        subset.clustering <- function(g,f=stability.subsampling.fraction,seed=NULL, ...) {
          if(!is.null(seed)) { set.seed(seed) }
          vi <- sample(1:length(V(g)),ceiling(length(V(g))*(f)))
          sg <- induced_subgraph(g,vi)
          method(sg,...)
        }
        if(verbose) { cat("running",stability.subsamples,"subsampling iterations ... ")}
        if(is.null(sr)) {
          sr <- papply(1:stability.subsamples,function(i) subset.clustering(graph,f=stability.subsampling.fraction,seed=i),n.cores=n.cores)
        }

        if(verbose) { cat("done\n")}

        if(verbose) cat("calculating flat stability stats ... ")
        # Jaccard coefficient for each cluster against all, plus random expecctation
        jc.stats <- do.call(rbind,conos:::papply(sr,function(o) {
          p1 <- membership(o);
          p2 <- cls.groups[names(p1)]; p1 <- as.character(p1)
          #x <- tapply(1:length(p2),factor(p2,levels=cls.levs),function(i1) {
          x <- tapply(1:length(p2),p2,function(i1) {
            i2 <- which(p1==p1[i1[[1]]])
            length(intersect(i1,i2))/length(unique(c(i1,i2)))
          })
        },n.cores=n.cores,mc.preschedule=T))

        # Adjusted rand index
        if(verbose) cat("adjusted Rand ... ")
        ari <- unlist(conos:::papply(sr,function(o) { ol <- membership(o); adjustedRand(as.integer(ol),as.integer(cls.groups[names(ol)]),randMethod='HA') },n.cores=n.cores))
        if(verbose) cat("done\n");

        res$stability <- list(flat=list(jc=jc.stats,ari=ari))

        # hierarchical measures
        if(verbose) cat("calculating hierarchical stability stats ... ")
        if(is.hierarchical(cls)) {
          # hierarchical to hierarchical stability analysis - cut reference
          # determine hierarchy of clusters (above the cut)
          t.get.walktrap.upper.merges <- function(res,n=length(unique(membership(res)))) {
            clm <- igraph:::complete.dend(res,FALSE)
            x <- tail(clm,n-1)
            x <- x - 2*nrow(res$merges) + nrow(x)-1
            # now all >=0 ids are cut leafs and need to be reassigned ids according to their rank
            xp <- x+nrow(x)+1
            xp[x<=0] <- rank(-x[x<=0])
            xp
          }

          clm <- t.get.walktrap.upper.merges(cls)
          res$stability$upper.tree <- clm

          if(verbose) cat("tree Jaccard ... ")
          jc.hstats <- do.call(rbind,mclapply(sr,function(z) conos:::bestClusterThresholds(z,cls.groups,clm)$threshold ,mc.cores=n.cores))

        } else {
          # compute cluster hierarchy based on cell mixing (and then something)
          # assess stability for that hierarchy (to visualize internal node stability)
          # for the original clustering and every subsample clustering,
          if(verbose) cat("upper clustering ... ")
          cgraph <- getClusterGraph(graph,cls.groups,plot=F,normalize=F)
          chwt <- walktrap.community(cgraph,steps=9)
          clm <- igraph:::complete.dend(chwt,FALSE)

          if(verbose) cat("clusterTree Jaccard ... ")
          jc.hstats <- do.call(rbind,mclapply(sr,function(st1) {
            mf <- membership(st1); mf <- as.factor(setNames(as.character(mf),names(mf)))
            st1g <- getClusterGraph(graph,mf,plot=F,normalize=T)
            st1w <- walktrap.community(st1g,steps=8)

            #merges <- st1w$merge; leaf.factor <- mf; clusters <- cls.groups
            x <- conos:::bestClusterTreeThresholds(st1w,mf,cls.groups,clm)
            x$threshold
          },mc.cores=n.cores))

        }
        res$stability$upper.tree <- clm
        res$stability$sr <- sr
        res$stability$hierarchical <- list(jc=jc.hstats);
        if(verbose) cat("done\n");

      }

      ## Filter groups
      if(min.group.size>0) {
        lvls.keep <- names(which(table(cls.groups)  > min.group.size))
        cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
        cls.groups <- as.factor(cls.groups)
      }
      res$groups <- cls.groups;
      clusters[[name]] <<- res;
      return(invisible(res))

    },

    plotPanel=function(clustering=NULL, groups=NULL, colors=NULL, gene=NULL, use.local.clusters=FALSE, plot.theme=NULL, ...) {
      if (use.local.clusters) {
        if (is.null(clustering) && !(inherits(x = samples[[1]], what = c('seurat', 'Seurat')))) {
          stop("You have to provide 'clustering' parameter to be able to use local clusters")
        }

        groups <- Reduce(c, lapply(samples, getClustering, clustering))
        if (is.null(groups)) {
          stop(paste0("No clustering '", clustering, "' presented in the samples"))
        }
      }
      else if (is.null(groups) && is.null(colors) && is.null(gene)) {
        groups <- getClusteringGroups(clusters, clustering)
      }

      gg <- plotSamples(samples, groups=groups, colors=colors, gene=gene, plot.theme=adjustTheme(plot.theme), ...)

      return(gg)
    },

    embedGraph=function(method='largeVis', M=1, gamma=1, alpha=0.1, perplexity=NA, sgd_batches=1e8, seed=1, verbose=TRUE, target.dims=2, n.cores=NULL, ...) {
      "Generate an embedding of a joint graph.\n
       Params:\n
       - method: embedding method. Currently largeVis and UMAP are supported\n
       - M, gamma, alpha, sgd__batched - largeVis parameters (defaults are 1, 1, 0.01, 1e8 respectively).\n
       - perplexity: perplexity passed to largeVis (defaults to NA).\n
       - seed: random seed for the largeVis algorithm. Default: 1.\n
       - target.dims: numer of dimensions for the reduction. Default: 2. Higher dimensions can be used to generate embeddings for subsequent reductions by other methods, such as tSNE\n
       - n.cores: number of cores, overrides class field
       - verbose: verbose mode. Default: TRUE.
       - ...: additional arguments, passed to UMAP embedding (run ?conos:::embedGraphUmap for more info)
      "

      n.jobs <- if (is.null(n.cores)) .self$n.cores else n.cores

      supported.methods <- c('largeVis', 'UMAP')
      if(!method %in% supported.methods) { stop(paste0("currently, only the following embeddings are supported: ",paste(supported.methods,collapse=' '))) }

      if (method == 'largeVis') {
        wij <- as_adj(graph,attr='weight');
        if(!is.na(perplexity)) {
          wij <- conos:::buildWijMatrix(wij,perplexity=perplexity,threads=n.jobs)
        }
        coords <- conos:::projectKNNs(wij = wij, dim=target.dims, verbose = verbose,sgd_batches = sgd_batches,gamma=gamma, M=M, seed=seed, alpha=alpha, rho=1, threads=n.jobs)
        colnames(coords) <- V(graph)$name
        embedding <<- t(coords);
      } else {
        if (!requireNamespace("uwot", quietly=T))
          stop("You need to install package 'uwot' to be able to use UMAP embedding.")

        embedding <<- embedGraphUmap(graph, verbose=verbose, return.all=F, n.cores=n.jobs, ...)
      }

      return(invisible(embedding))
    },

    plotClusterStability=function(clustering=NULL,what='all') {
      "Plot cluster stability statistics.\n
       Params:\n
       - clustering : name of the clustering result to show
       - what : show a specific plot (ari - adjusted rand index, fjc - flat Jaccard, hjc - hierarchical Jaccard, dend - cluster dendrogram)
      "

      if(is.null(clustering)) clustering <- names(clusters)[[1]]

      if(is.null(clusters[[clustering]]))
        stop(paste("clustering",clustering,"doesn't exist, run findCommunity() first"))

      if(is.null(clusters[[clustering]]$stability))
        stop(paste("clustering",clustering,"doesn't have stability info. Run findCommunity( ... , test.stability=TRUE) first"))

      st <- clusters[[clustering]]$stability
      nclusters <- ncol(st$flat$jc)
      jitter.alpha <- 0.1;

      if(what=='all' || what=='ari') {
        p.fai <- ggplot2::ggplot(data.frame(aRI=st$flat$ari), ggplot2::aes(x=1,y=aRI)) +
          ggplot2::geom_boxplot(notch=T,outlier.shape=NA) +
          ggplot2::geom_point(shape=16, position = ggplot2::position_jitter(), alpha=jitter.alpha) +
          ggplot2::guides(color=FALSE) +
          ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
          ggplot2::ylim(c(0,1)) + ggplot2::labs(x=" ", y="adjusted Rand Index") +
          ggplot2::theme(legend.position="none", axis.ticks.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank())

        if(what=='ari')
          return(p.fai)
      }

      if(what=='all' || what=='fjc') {
        df <- reshape2::melt(st$flat$jc);
        colnames(df) <- c('rep','cluster','jc')
        df$cluster <- factor(colnames(st$flat$jc)[df$cluster],levels=levels(clusters[[clustering]]$groups))

        p.fjc <- ggplot2::ggplot(df,aes(x=cluster,y=jc,color=cluster)) +
          ggplot2::geom_boxplot(aes(color=cluster),notch=T,outlier.shape=NA) +
          ggplot2::geom_jitter(shape=16, position=position_jitter(0.2),alpha=jitter.alpha) +
          ggplot2::guides(color=FALSE) +
          ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
          ggplot2::ylab("Jaccard coefficient (flat)") + ggplot2::ylim(c(0,1))

        if(what=='fjc') return(p.fjc)
      }

      if(what=='all' || what=='hjc') {
        # hierarchical
        df <- reshape2::melt(st$hierarchical$jc[,1:nclusters])
        colnames(df) <- c('rep','cluster','jc');
        df$cluster <- factor(colnames(st$flat$jc)[df$cluster],levels=levels(clusters[[clustering]]$groups))
        p.hjc <- ggplot2::ggplot(df,aes(x=cluster,y=jc,color=cluster)) +
          ggplot2::geom_boxplot(aes(color=cluster),notch=T,outlier.shape=NA) +
          ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha=jitter.alpha) +
          ggplot2::guides(color=FALSE) +
          ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
          ggplot2::ylab("Jaccard coefficient (hierarchical)") + ggplot2::ylim(c(0,1))

        if(what=='hjc') return(p.hjc)
      }

      if(what=='dend') {
        require(dendextend)
        m <- st$upper.tree; nleafs <- nrow(m)+1; m[m<=nleafs] <- -1*m[m<=nleafs]; m[m>0] <- m[m>0]-nleafs;
        hc <- list(merge=m,height=1:nrow(m),labels=levels(clusters[[clustering]]$groups),order=c(1:nleafs)); class(hc) <- 'hclust'
        # fix the ordering so that edges don't intersects
        hc$order <- order.dendrogram(as.dendrogram(hc))

        d <- as.dendrogram(hc) %>% hang.dendrogram()

        # depth-first traversal of a merge matrix
        t.dfirst <- function(m,i=nrow(m)) {
          rl <- m[i,1]; if(rl<0) { rl <- abs(rl) } else { rl <- t.dfirst(m,rl) }
          rr <- m[i,2]; if(rr<0) { rr <- abs(rr) } else { rr <- t.dfirst(m,rr) }
          c(i+nrow(m)+1,rl,rr)
        }
        xy <- get_nodes_xy(d)
        to <- t.dfirst(hc$merge)
        plot(d,las=2,axes=F)
        # flat on the left
        #x <- apply(st$flat$jc,2,median)
        #text(xy,labels=round(x[to],2),col='blue',adj=c(-0.1,-1.24),cex=0.8)
        x <- apply(st$hierarchical$jc,2,median)
        text(xy,labels=round(x[to],2),col='red',adj=c(-0.1,-0.12),cex=0.8)
        return(NULL)
      }

      cowplot::plot_grid(plotlist=list(p.fai,p.fjc,p.hjc),nrow=1,rel_widths=c(4,nclusters,nclusters))
    },


    plotGraph=function(color.by='cluster', clustering=NULL, groups=NULL, colors=NULL, gene=NULL, plot.theme=NULL, ...) {
      if(class(embedding)[1] == "uninitializedField") {
        embedGraph();
      }

      if (!is.null(gene)) {
        colors <- lapply(samples, getCountMatrix) %>% lapply(getGeneExpression, gene) %>% Reduce(c, .)
        if(all(is.na(colors))) warning(paste("gene",gene,"is not found in any of the samples"))
      }

      if(is.null(groups) && is.null(colors)) {
        if(color.by == 'cluster') {
          groups <- getClusteringGroups(clusters, clustering)
        } else if(color.by == 'sample') {
          groups <- getDatasetPerCell()
        } else {
          stop('supported values of color.by are ("cluster" and "sample")')
        }
      }

      return(embeddingPlot(embedding, groups=groups, colors=colors, plot.theme=adjustTheme(plot.theme), ...))
    },

    correctGenes=function(genes=NULL, n.od.genes=500, fading=10.0, fading.const=0.5, max.iters=15, tol=5e-3, name='diffusion', verbose=TRUE, count.matrix=NULL, normalize=TRUE) {
      "Smooth expression of genes, so they better represent structure of the graph.\n
       Use diffusion of expression on graph with the equation dv = exp(-a * (v + b))\n
       Params:\n
       - genes: list of genes for smoothing\n
       - n.od.genes: if 'genes' is NULL, top n.od.genes of overdispersed genes are taken across all samples. Default: 500.\n
       - fading: level of fading of expression change from distance on the graph (parameter 'a' of the equation). Default: 10.\n
       - fading.const: minimal penalty for each new edge during diffusion (parameter 'b' of the equation). Default: 0.5.\n
       - max.iters: maximal number of diffusion iterations. Default: 15.\n
       - tol: tolerance after which the diffusion stops. Default: 5e-3.\n
       - name: name to save the correction. Default: diffusion.\n
       - verbose: verbose mode. Default: TRUE.
       - count.matrix: alternative gene count matrix to correct. Default: joint count matrix for all datasets.
      "
      edges <- igraph::as_edgelist(graph)
      edge.weights <- igraph::edge.attributes(graph)$weight

      if (is.null(count.matrix)) {
        if (is.null(genes)) {
          genes <- getOdGenesUniformly(samples, n.genes=n.od.genes)
        }

        cms <- lapply(samples, `[[`, "counts")
        genes <- Reduce(intersect, lapply(cms, colnames)) %>% intersect(genes)
        count.matrix <- Reduce(rbind, lapply(cms, function(x) x[, genes])) %>% as.matrix()
      } else {
        count.matrix <- t(count.matrix)
      }

      cm <- smooth_count_matrix(edges, edge.weights, count.matrix, max_n_iters=max.iters, diffusion_fading=fading, diffusion_fading_const=fading.const, verbose=verbose, normalize=normalize)
      return(invisible(expression.adj[[name]] <<- cm))
    },

    propagateLabels=function(labels, max.iters=50, diffusion.fading=10.0, diffusion.fading.const=0.5, tol=5e-3, return.distribution=TRUE, verbose=TRUE, fixed.initial.labels=FALSE) {
    "Estimate labeling distribution for each vertex, based on provided labels.\n
     Params:\n
     - labels: vector of factor or character labels, named by cell names\n
     - max.iters: maximal number of iterations. Default: 50.\n
     - return.distribution: return distribution of labeling, but not single label for each vertex. Default: TRUE.\n
     - verbose: verbose mode. Default: TRUE.\n
     - fixed.initial.labels: prohibit changes of initial labels during diffusion.\n
     \n
     Return: matrix with distribution of label probabilities for each vertex by rows.
    "
      if (is.factor(labels)) {
        labels <- as.character(labels) %>% setNames(names(labels))
      }
      edges <- igraph::as_edgelist(graph)
      edge.weights <- igraph::edge.attributes(graph)$weight
      labels <- labels[intersect(names(labels), igraph::vertex.attributes(graph)$name)]

      label.distribution <- propagate_labels(edges, edge.weights, vert_labels=labels, max_n_iters=max.iters, verbose=verbose,
                                             diffusion_fading=diffusion.fading, diffusion_fading_const=diffusion.fading.const,
                                             tol=tol, fixed_initial_labels=fixed.initial.labels)
      if (return.distribution) {
        return(label.distribution)
      }

      label.distribution <- colnames(label.distribution)[apply(label.distribution, 1, which.max)] %>%
        setNames(rownames(label.distribution))
      return(label.distribution)
    },

    getClusterCountMatrices=function(clustering=NULL, groups=NULL,common.genes=TRUE,omit.na.cells=TRUE) {
      "Estimate per-cluster molecule count matrix by summing up the molecules of each gene for all of the cells in each cluster.\n\n
       Params:\n
       - clustering: the name of the clustering that should be used
       - groups: explicitly provided cell grouping\n
       - common.genes: bring individual sample matrices to a common gene list
       - omit.na.cells: if set to FALSE, the resulting matrices will include a first column named 'NA' that will report total molecule counts for all of the cells that were not covered by the provided factor.
       \n
       Return: a list of per-sample uniform dense matrices with rows being genes, and columns being clusters
      "
      if(is.null(groups)) {
        groups <- getClusteringGroups(clusters, clustering)
      }

      groups <- as.factor(groups)

      matl <- lapply(samples,function(s) {
        m <- conos:::getRawCountMatrix(s,trans=TRUE); # rows are cells
        cl <- factor(groups[match(rownames(m),names(groups))],levels=levels(groups));
        tc <- colSumByFactor(m,cl);
        if(omit.na.cells) { tc <- tc[-1,,drop=F] }
        t(tc);
      })
      # bring to a common gene space
      if(common.genes) {
        gs <- unique(unlist(lapply(matl,rownames)))
        matl <- lapply(matl,function(m) {
          nm <- matrix(0,nrow=length(gs),ncol=ncol(m))
          colnames(nm) <- colnames(m); rownames(nm) <- gs;
          mi <- match(rownames(m),gs)
          nm[mi,] <- m;
          nm
        })
      }
      matl
    },

    getDatasetPerCell=function() {
      cl <- lapply(samples, getCellNames)
      return(rep(names(cl), sapply(cl, length)) %>% stats::setNames(unlist(cl)) %>% as.factor())
    }
  )
)
