#' @import ggrepel
#' @import gridExtra
#' @import leidenAlg
#' @import R6
#' @import Rtsne
NULL

## exporting to inherit parameters below, leiden.community
#' @export
leidenAlg::leiden.community


#' @title Conos R6 class
#' @description The class encompasses sample collections, providing methods for calculating and visualizing joint graph and communities.
#' @import methods
#' @param x a named list of pagoda2 or Seurat objects (one per sample)
#' @param n.cores numeric Number of cores to use (default=parallel::detectCores(logical=FALSE))
#' @param verbose boolean Whether to provide verbose output (default=TRUE)
#' @param clustering string Name of the clustering to use
#' @param groups a factor on cells to use for coloring
#' @param colors a color factor (named with cell names) use for cell coloring
#' @param gene show expression of a gene
#' @export Conos
Conos <- R6::R6Class("Conos", lock_objects=FALSE,
  public = list(
    #' @field samples list of samples (Pagoda2 or Seurat objects)
    samples = list(),

    #' @field pairs pairwise alignment results
    pairs = list(),

    #' @field graph alignment graph
    graph = NULL,

    #' @field clusters list of clustering results named by clustering type
    clusters = list(),

    #' @field expression.adj adjusted expression values
    expression.adj = list(),

    #' @field embeddings list of joint embeddings
    embeddings = list(),

    #' @field embedding joint embedding
    embedding = NULL,

    #' @field n.cores number of cores
    n.cores = 1,

    #' @field misc list with unstructured additional info
    misc = list(),

    #' @field override.conos.plot.theme
    override.conos.plot.theme = FALSE,

    #' @description initialize Conos class
    #'
    #' @param override.conos.plot.theme boolean Whether to reset plot settings to the ggplot2 default (default=FALSE)
    #' @param ... additional parameters upon initializing Conos
    #' @return a new 'Conos' object
    #' @examples
    #' con <- Conos$new(small_panel.preprocessed, n.cores=1)
    #'
    initialize=function(x, ..., n.cores=parallel::detectCores(logical=FALSE), verbose=TRUE, override.conos.plot.theme=FALSE) {
      self$n.cores <- n.cores
      self$override.conos.plot.theme <- override.conos.plot.theme

      if (missing(x)){
        return()
      }

      if ('Conos' %in% class(x)) { # copy constructor
        for(n in ls(x)) {
          if (!is.function(get(n, x))) assign(n, get(n, x), self)
        }
      } else {
        if (!is.list(x)) {
          stop("x is not a list of pagoda2 or Seurat objects")
        }
        if (inherits(x = x[[1]], what = c('Pagoda2', 'seurat', 'Seurat'))) {
          self$addSamples(x)
        } else {
          stop("Only Pagoda2 or Seurat result lists are currently supported")
        }
      }
    },

    #' @description Initialize or add a set of samples to the conos panel. Note: this will simply add samples, but will not update graph, clustering, etc.
    #'
    #' @param replace boolean Whether the existing samples should be purged before adding new ones (default=FALSE)
    #' @param verbose boolean Whether to provide verbose output (default=FALSE)
    #' @return invisible view of the full sample list
    addSamples=function(x, replace=FALSE, verbose=FALSE) {
      # check names
      if(is.null(names(x))) {
        stop("The sample list must be named")
      }
      if(replace || length(self$samples)==0) {
        if(length(x)<2) {
          stop("The provided list contains less than 2 samples; 2 required, >3 recommended")
        }
      }
      if(any(duplicated(names(self$samples)))) {
        stop("duplicate names found in the supplied samples")
      }
      if(any(names(x) %in% names(self$samples))) {
        stop("Some of the names in the provided sample list are identical to already existing samples")
      }

      # TODO: package-independent wrapper
      self$samples <- c(self$samples, x)
    },

    #' @description Build the joint graph that encompasses all the samples, establishing weighted inter-sample cell-to-cell links
    #'
    #' @param k integer integer Size of the inter-sample neighborhood (default=15)
    #' @param k.self integer Size of the with-sample neighborhoods (default=10).
    #' @param k.self.weight numeric Weight multiplier on the intra-sample edges relative to inter-sample edges (default=0.1)
    #' @param alignment.strength numeric Alignment strength (default=NULL will result in alignment.strength=0)
    #' @param space character Reduced expression space used to establish putative alignments between pairs of samples (default='PCA'). Currently supported spaces are:
    #'     --- "CPCA" Common principal component analysis
    #'     --- "JNMF" Joint NMF
    #'     --- "genes" Gene expression space (log2 transformed)
    #'     --- "PCA" Principal component analysis
    #'     --- "CCA" Canonical correlation analysis
    #'     --- "PMA" (Penalized Multivariate Analysis <https://cran.r-project.org/web/packages/PMA/index.html>)
    #' @param matching.method character Matching method (default='mNN'). Currently supported methods are "NN" (nearest neighbors) or "mNN" (mututal nearest neighbors).
    #' @param metric character Distance metric to measure similarity (default='angular'). Currenlty supported metrics are "angular" and "L2".
    #' @param k1 numeric Neighborhood radius for identifying mutually-matching neighbors (default=k). Note that k1 must be greater than or equal to k, i.e. k1>=k. Increasing k1 beyond k will lead to more aggressive alignment of distinct subpopulations (i.e. increased alignment strengths).
    #' @param data.type character Type of data type in the input pagoda2 objects within r.n (default='counts').
    #' @param l2.sigma numeric L2 distances get transformed as exp(-d/sigma) using this value (default=1e5)
    #' @param var.scale boolean Whether to use common variance scaling (default=TRUE). If TRUE, use geometric means for variance, as we're trying to focus on the common variance components. See scaledMatricesP2() code.
    #' @param ncomps integer Number of components (default=40)
    #' @param n.odgenes integer Number of overdispersed genes to be used in each pairwise alignment (default=2000)
    #' @param matching.mask an optional matrix explicitly specifying which pairs of samples should be compared (a symmetrical matrix of logical values with row and column names corresponding to sample names). (default=NULL). By default, comparisons between all paris are allowed. The argument can be used to exclude comparisons across certain pairs of samples (e.g. techincal replicates, which are expected to show very high similarity).
    #' @param exclude.samples optional list of sample names that should be excluded from the alignment and the resulting graph (default=NULL)
    #' @param common.centering boolean When calculating reduced expression space for a given sample pair, whether the expression of genes should be centered using the mean from both samples (TRUE) or using the mean within each sample (FALSE) (default=TRUE)
    #' @param base.groups an optional factor on cells specifying previously-obtained cell grouping to be used for adjusting the sample alignment (default: NULL). Specifically, cell clusters specfiieid by the base.groups can be used to i) calculate global expression axes which are appended to the overall set of eigenvectors, ii) adding decoy cells.
    #' @param append.global.axes boolean Whether to project samples on global expression axes, as defined by pre-defined (typically crude) set of cell subpopulations as specified by the base.gruops parameter (default=TRUE, but works only if base.groups is specified)
    #' @param append.decoys boolean Whether to use pre-defined cell groups (specified by base.groups) to append decoy cells to the samples which are otherwise lacking any of the pre-specified cell groups (default=TRUE, but works only if base.groups is specified). The decoy cells can reduce the number of erroneous matches in highly heterogeneous sample collections, where some of the samples lack entire cell subpopulations which are found in other samples. The approach only works if the base.groups (typically a crude clustering of top-level cell types) can be established with a reasonable confidence.
    #' @param decoy.threshold integer Minimal number of cells of a given cell type that should exist in a given sample (according to base.groups) to avoid addition of decoy cells to that sample for the purposes of alignment (default=1)
    #' @param n.decoys integer Number of decoy cells that should be added to a sample that had less than decoy.threshold cells of a given cell type (default=k*2)
    #' @param score.component.variance boolean Whether to score the amount of total variance explained by different components (default=FALSE as it takes extra time to calculate)
    #' @param snn boolean Whether to transform the joint graph by computing a shared nearest neighborhood graph (analogous to Seurat 3), further weighting the edges between two matched cells based on the similarity (measured by Jaccard coefficient) of all of their predicted neighbors (across all of the samples) (default: FALSE)
    #' @param snn.quantile numeric Specifies how the shared neighborhood graph transformation will determine final edge weights. If snn.quantile=NULL, the edge weight will be simply equal to the Jaccard coefficient of the neighborhoods. If snn.quantile is a vector of two numeric values (p1, p2), they will be treated as quantile probabilities, and quantile values (q1,q2) on the set of all Jaccard coefficients (for all edges) will be determiend. The edge weights will then be reset, so that edges with Jaccard coefficients below or equal to q1 will be set to 0, and those with coefficients >=q2 will be set to 1. The rest of the weights will be mapped uniformly from [q1,q2]->[0,1] range. If a single numeric value is supplied, it will be treated as a symmetric quantile probability (i.e. snn.quantile=0.8 is equivalent to specifying snn.quantile=c(1-0.8,0.8)). (default: 0.9)
    #' @param min.snn.jaccard numeric Minimum Jaccard coefficient required for a shared neighborhood graph edge (default: 0). The edges with Jaccard coefficients below this threshold will be removed (i.e. weight set to 0)
    #' @param min.snn.weight numeric Shared nearest neighbor procedure will adjust the weights of the edges, and even eliminate some of the edges (by setting their weight to zero). The min.snn.weight parameter allows to set a minimal adjusted edge weight, so that the edge weight is never reduced beyond this level (and hence never deleted) (default: 0 - no adjustments)
    #' @param snn.k.self integer Size of the within-sample neighorhood to be used in shared nearest neighbor calculations (default=k.self)
    #' @param balance.edge.weights boolean Whether to balance edge weights to control for a cell- or sample- specific factor (default=FALSE)
    #' @param balancing.factor.per.sample A covariate factor per sample that should be controlled for by adjusting edge weights in the joint graph (default=NULL)
    #' @param balancing.factor.per.cell A per-cell factor (discrete factor, named with cell names) specifying a design difference should be controlled for by adjusting edge weights in the joint graph (default=NULL)
    #' @param same.factor.downweight numeric Optional weighting factor for edges connecting cells with the same cell factor level per cell balancing (default=1.0)
    #' @param k.same.factor integer An neighborhood size that should be used when aligning samples of the same balancing.factor.per.sample level. Setting a value smaller than k will lead to reduction of alingment strenth within the sample batches (default=k)
    #' @return joint graph to be used for downstream analysis
    #' @examples
    #' con <- Conos$new(small_panel.preprocessed, n.cores=1)
    #' con$buildGraph(k=10, k.self=5, space='PCA', ncomps=10, n.odgenes=20, matching.method='mNN', 
    #'     metric='angular', score.component.variance=TRUE, verbose=TRUE)
    #'
    #'
    buildGraph=function(k=15, k.self=10, k.self.weight=0.1, alignment.strength=NULL, space='PCA', matching.method='mNN', metric='angular', k1=k, data.type='counts', l2.sigma=1e5, var.scale=TRUE, ncomps=40,
                        n.odgenes=2000, matching.mask=NULL, exclude.samples=NULL, common.centering=TRUE, verbose=TRUE,
                        base.groups=NULL, append.global.axes=TRUE, append.decoys=TRUE, decoy.threshold=1, n.decoys=k*2, score.component.variance=FALSE,
                        snn=FALSE, snn.quantile=0.9, min.snn.jaccard=0, min.snn.weight=0, snn.k.self=k.self,
                        balance.edge.weights=FALSE, balancing.factor.per.cell=NULL, same.factor.downweight=1.0, k.same.factor=k, balancing.factor.per.sample=NULL) {

      supported.spaces <- c("CPCA","JNMF","genes","PCA","PMA","CCA")
      if (!space %in% supported.spaces) {
        stop(paste0("only the following spaces are currently supported: [",paste(supported.spaces,collapse=' '),"]"))
      }
      supported.matching.methods <- c("mNN", "NN")
      if (!matching.method %in% supported.matching.methods) {
        stop(paste0("only the following matching methods are currently supported: ['",paste(supported.matching.methods,collapse="' '"),"']"))
      }
      supported.metrics <- c("L2","angular")
      if (!metric %in% supported.metrics) {
        stop(paste0("only the following distance metrics are currently supported: ['",paste(supported.metrics,collapse="' '"),"']"))
      }

      if (!is.null(snn.quantile) && !is.na(snn.quantile)) {
        if(length(snn.quantile)==1)  {
          snn.quantile <- c(1-snn.quantile,snn.quantile)
        }
        snn.quantile <- sort(snn.quantile,decreasing=FALSE)
        if(snn.quantile[1]<0 | snn.quantile[2]>1) {
          stop("snn.quantile must be one or two numbers in the [0,1] range")
        }
      }

      if (!is.null(alignment.strength)) {
        alignment.strength %<>% max(0) %>% min(1)
        k1 <- sapply(self$samples, function(sample) ncol(getCountMatrix(sample))) %>% max() %>%
          `*`(alignment.strength ^ 2) %>% round() %>% max(k)
      } else {
        alignment.strength <- 0 # otherwise, estimation of cor.base uses NULL value
      }

      if (k1<k) { stop("k1 must be >= k") }
      # calculate or update pairwise alignments
      sn.pairs <- private$updatePairs(space=space, ncomps=ncomps, n.odgenes=n.odgenes, verbose=verbose, var.scale=var.scale, matching.mask=matching.mask, exclude.samples=exclude.samples, score.component.variance=score.component.variance)
      if (ncol(sn.pairs)<1) { stop("insufficient number of comparable pairs") }
      if (!is.null(base.groups)) {
        samf <- lapply(self$samples,getCellNames)
        base.groups <- as.factor(base.groups[names(base.groups) %in% unlist(samf)]) # clean up the group factor
        if (length(base.groups) < 2) stop("provided base.groups doesn't cover enough cells")
        # make a sample factor
        samf <- setNames(rep(names(samf),unlist(lapply(samf,length))),unlist(samf))
        if (append.global.axes) {
          cms.clust <- self$getClusterCountMatrices(groups=base.groups, common.genes=FALSE)
          global.proj <- projectSamplesOnGlobalAxes(self$samples, cms.clust, data.type, verbose, self$n.cores)
        }
      }

      if (snn){
        local.neighbors <- getLocalNeighbors(self$samples[! names(self$samples) %in% exclude.samples], snn.k.self, k.self.weight, metric, l2.sigma=l2.sigma, verbose, self$n.cores)
      } else {
        local.neighbors <- NULL
      }

      # determine inter-sample mapping
      if (verbose) message('inter-sample links using ',matching.method,' ')
      cached.pairs <- self$pairs[[space]]
      cor.base <- 1 + min(1, alignment.strength * 10)  ## see convertDistanceToSimilarity()
      mnnres <- papply(1:ncol(sn.pairs), function(j) {
        # we'll look up the pair by name (possibly reversed), not to assume for the ordering of $pairs[[space]] to be the same
        i <- match(paste(sn.pairs[,j],collapse='.vs.'),names(cached.pairs))
        if(is.na(i)) { i <- match(paste(rev(sn.pairs[,j]),collapse='.vs.'),names(cached.pairs)) }
        if(is.na(i)) { stop(paste("unable to find alignment for pair",paste(sn.pairs[,j],collapse='.vs.'))) }

        k.cur <- k
        if (!is.null(balancing.factor.per.sample) && (balancing.factor.per.sample[sn.pairs[1,j]] == balancing.factor.per.sample[sn.pairs[2,j]])) {
          k.cur <- min(k.same.factor, k1) # It always should be less then k1, though never supposed to be set higher
        }

        if(space=='JNMF') {
          mnn <- getNeighborMatrix(cached.pairs[[i]]$rot1, cached.pairs[[i]]$rot2,
                                   k=k.cur, k1=k1, matching=matching.method, metric=metric, l2.sigma=l2.sigma, cor.base=cor.base)
        } else if (space %in% c("CPCA","GSVD","PCA")) {
          #common.genes <- Reduce(intersect,lapply(r.ns, getGenes))
          if(!is.null(cached.pairs[[i]]$CPC)) {
            # CPCA or PCA
            #od.genes <- intersect(rownames(cached.pairs[[i]]$CPC),common.genes)
            od.genes <- rownames(cached.pairs[[i]]$CPC)
            rot <- cached.pairs[[i]]$CPC[od.genes,]
          } else if(!is.null(cached.pairs[[i]]$o$Q)) {
            # GSVD
            rot <- cached.pairs[[i]]$o$Q
            od.genes <- rownames(rot) <- colnames(cached.pairs[[i]]$o$A)
          } else {
            stop("unknown reduction provided")
          }

          # TODO: a more careful analysis of parameters used to calculate the cached version
          if(ncomps > ncol(rot)) {
            warning(paste0("specified ncomps (",ncomps,") is greater than the cached version (",ncol(rot),")"))
          } else {
            rot <- rot[,1:ncomps,drop=FALSE]
          }

          mnn <- getPcaBasedNeighborMatrix(self$samples[sn.pairs[,j]], od.genes=od.genes, rot=rot, data.type=data.type,
                                           k=k.cur, k1=k1, matching.method=matching.method, metric=metric, l2.sigma=l2.sigma, cor.base=cor.base,
                                           var.scale=var.scale, common.centering=common.centering,
                                           base.groups=base.groups, append.decoys=append.decoys, samples=self$samples, samf=samf, decoy.threshold=decoy.threshold,
                                           n.decoys=n.decoys, append.global.axes=append.global.axes, global.proj=global.proj)
        } else if (space=='genes') { ## Overdispersed Gene space
          mnn <- getNeighborMatrix(as.matrix(cached.pairs[[i]]$genespace1), as.matrix(cached.pairs[[i]]$genespace2),
                                   k=k.cur, k1=k1, matching=matching.method, metric=metric, l2.sigma=l2.sigma, cor.base=cor.base)
        } else if(space=='PMA' || space=='CCA') {
          mnn <- getNeighborMatrix(cached.pairs[[i]]$u, cached.pairs[[i]]$v,
                                   k=k.cur, k1=k1, matching=matching.method, metric=metric, l2.sigma=l2.sigma, cor.base=cor.base)
        } else {
          stop("Unknown space: ", space)
        }

        if(snn) { # optionally, perform shared neighbor weighting, a la SeuratV3, scran

          m1 <- cbind(mnn,t(local.neighbors[[sn.pairs[1,j] ]]))
          m2 <- rbind(local.neighbors[[sn.pairs[2,j] ]], mnn)

          mnn1 <- mnn
          mnn1@x <- rep(1,length(mnn1@x))
          m1@x <- rep(1,length(m1@x))
          m2@x <- rep(1,length(m2@x))

          x <- ((m1 %*% m2) * mnn1) / pmax(outer(rowSums(m1),colSums(m2),FUN=pmin),1)

          # scale by Jaccard coefficient

          if(min.snn.jaccard>0) {
            x@x[x@x<min.snn.jaccard] <- 0
          }

          if(!is.null(snn.quantile) && !is.na(snn.quantile)) {
            xq <- quantile(x@x,p=c(snn.quantile[1],snn.quantile[2]))
            x@x <- pmax(0,pmin(1,(x@x-xq[1])/pmax(1,diff(xq))))
          }

          x <- drop0(x)
          if (min.snn.weight>0) {
            mnn <- as(drop0(mnn*min.snn.weight + mnn*x),'dgTMatrix')
          } else {
            mnn <- as(drop0(mnn*x),'dgTMatrix')
          }

        }


        if(verbose) cat(".")
        return(data.frame('mA.lab'=rownames(mnn)[mnn@i+1],'mB.lab'=colnames(mnn)[mnn@j+1],'w'=mnn@x, stringsAsFactors=FALSE))
      }, n.cores=self$n.cores,mc.preschedule=TRUE)

      if (verbose) message(" done")
      ## Merge the results into a edge table
      el <- do.call(rbind,mnnres)
      if (nrow(el)==0) {
        el = data.frame('mA.lab'=0,'mB.lab'=0,'w'=0, 'type'=1, stringsAsFactors=FALSE)
      } else {
        el$type <- 1 # encode connection type 1- intersample, 0- intrasample
      }
      # append local edges
      if(k.self>0) {
        if(is.null(local.neighbors) || snn.k!=k.self) { # recalculate local neighbors
          local.neighbors <- getLocalNeighbors(self$samples[! names(self$samples) %in% exclude.samples], k.self, k.self.weight, metric, l2.sigma=l2.sigma, verbose, self$n.cores)
        }
        el <- rbind(el,getLocalEdges(local.neighbors))
      }
      if(verbose) message('building graph .')
      el <- el[el[,3]>0,]
      g  <- graph_from_edgelist(as.matrix(el[,c(1,2)]), directed =FALSE)
      E(g)$weight <- el[,3]
      E(g)$type <- el[,4]

      if(verbose) cat(".")

      # collapse duplicate edges
      g <- simplify(g, edge.attr.comb=list(weight="sum", type = "first"))
      if(verbose) message('done')

      if (!is.null(balancing.factor.per.sample)) {
        if (is.null(balancing.factor.per.cell)) {
          sf <- self$getDatasetPerCell()
          balancing.factor.per.cell <- setNames(balancing.factor.per.sample[as.character(sf)], names(sf))
        } else {
          warning("Both balancing.factor.per.cell and balancing.factor.per.sample are provided. Used the former for balancing edge weights")
        }
      }

      if (balance.edge.weights || !is.null(balancing.factor.per.cell)) {
        if(verbose) message('balancing edge weights ')

        if (is.null(balancing.factor.per.cell)) {
          balancing.factor.per.cell <- self$getDatasetPerCell()
        }

        g <- igraph::as_adjacency_matrix(g, attr="weight") %>%
          adjustWeightsByCellBalancing(factor.per.cell=balancing.factor.per.cell, balance.weights=balance.edge.weights,
                                       same.factor.downweight=same.factor.downweight) %>%
          igraph::graph_from_adjacency_matrix(mode="undirected", weighted=TRUE)

        if(verbose) message('done')
      }
      self$graph <- g
      return(invisible(g))
    },


    #' @description Calculate genes differentially expressed between cell clusters. Estimates base mean, z-score, p-values, specificity, precision, expressionFraction, AUC (if append.auc=TRUE)
    #'
    #' @param groups a cell factor (a factor named with cell names) specifying clusters of cells to be compared (one against all). To compare two cell clusters against each other, simply pass a factor containing only two levels (default: NULL, see clustering)
    #' @param clustering character Name of the clustering to use (see names(con$clusters)) for the value of the groups factor (default: NULL - if groups are not specified, the first clustering will be used)
    #' @param z.threshold numeric Minimum absolute value of a Z score for which the genes should be reported (default=3.0).
    #' @param upregulated.only boolean If TRUE, will report only genes significantly upregulated in each cluster; otherwise both up- and down-regulated genes will be reported (default=FALSE)
    #' @param append.specificity.metrics boolean Whether to append specificity metrics (default=TRUE)
    #' @param append.auc boolean Whether to append AUC scores (default=TRUE)
    #' @return list of DE results; each is a data frame with rows corresponding to the differentially expressed genes, and columns listing log2 fold change (M), signed Z scores (both raw and adjusted for mulitple hypothesis using BH correction), optional specificty/sensitivity and AUC metrics.
    #'
    getDifferentialGenes=function(clustering=NULL, groups=NULL, z.threshold=3.0, upregulated.only=FALSE, verbose=TRUE, append.specificity.metrics=TRUE, append.auc=TRUE) {

      groups <- parseCellGroups(self, clustering, groups)

      groups %<>% as.factor() %>% droplevels()
      # TODO: add Seurat
      '%ni%' <- Negate('%in%')
      if ('Pagoda2' %ni% class(self$samples[[1]])){
        stop("Only Pagoda2 objects are supported for marker genes")
      }

      de.genes <- getDifferentialGenesP2(self$samples, groups=groups, z.threshold=z.threshold, upregulated.only=upregulated.only, verbose=verbose, n.cores=self$n.cores)
      de.genes <- de.genes[levels(groups)]

      if (append.specificity.metrics) {
        if (verbose) message("Estimating specificity metrics")

        cm.merged <- self$getJointCountMatrix(raw=TRUE)
        groups.clean <- groups %>% .[!is.na(.)] %>% .[names(.) %in% rownames(cm.merged)]

        de.genes %<>% lapply(function(x) if ((length(x) > 0) && (nrow(x) > 0)) subset(x, complete.cases(x)) else x)
        de.genes %<>% names() %>% setNames(., .) %>%
          sccore::plapply(function(n) appendSpecificityMetricsToDE(de.genes[[n]], groups.clean, n, p2.counts=cm.merged, append.auc=append.auc), progress=verbose, n.cores=self$n.cores)
      }

      if (verbose) message("All done!")

      return(de.genes)
    },

    #' @description Find cell clusters (as communities on the joint graph)
    #'
    #' @param method community detection method (igraph syntax) (default=leiden.community)
    #' @param min.group.size numeric Minimal allowed community size (default=0)
    #' @param name character Optional name of the clustering result (will default to the algorithm name) (default=NULL will try to obtain the name from the community detection method, or will use 'community' as a default)
    #' @param test.stability boolean Whether to test stability of community detection (default=FALSE)
    #' @param stability.subsampling.fraction numeric Fraction of clusters to subset (default=0.95). Must be within range [0, 1].
    #' @param stability.subsamples integer Number of subsampling iterations (default=100)
    #' @param cls optional pre-calculated community result (may be useful for stability testing) (default: NULL)
    #' @param sr optional pre-calculated subsampled community results (useful for stability testing) (default: NULL)
    #' @param ... extra parameters are passed to the specified community detection method
    #' @return invisible list containing identified communities (groups) and the full community detection result (result); The results are stored in $clusters$name slot in the conos object. Each such slot contains an object with elements: $results which stores the raw output of the community detection method, and $groups which is a factor on cells describing the resulting clustering. The later can be used, for instance, in plotting: con$plotGraph(groups=con$clusters$leiden$groups). If test.stability==TRUE, then the result object will also contain a $stability slot.
    #' @examples
    #' con <- Conos$new(small_panel.preprocessed, n.cores=1)
    #' con$buildGraph(k=10, k.self=5, space='PCA', ncomps=10, n.odgenes=20, matching.method='mNN', 
    #'     metric='angular', score.component.variance=TRUE, verbose=TRUE)
    #' con$findCommunities(method = igraph::walktrap.community, steps=5)
    #'
    findCommunities=function(method=leiden.community, min.group.size=0, name=NULL, test.stability=FALSE, stability.subsampling.fraction=0.95, stability.subsamples=100, verbose=TRUE, cls=NULL, sr=NULL, ...) {

      if (is.null(cls)) {
        cls <- method(self$graph, ...)
      }
      if (is.null(name)) {
        name <- cls$algorithm
        if(is.null(name)) {
          name <- "community"
        }
      }

      ## Extract groups from this graph
      cls.mem <- membership(cls)
      if(suppressWarnings(any(is.na(as.numeric(cls.mem))))) {
        cls.groups <- as.factor(cls.mem)
      } else {
        cls.groups <- factor(setNames(as.character(cls.mem),names(cls.mem)),levels=sort(as.numeric(unique(as.character(cls.mem))),decreasing=FALSE))
      }
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
        if(verbose) { message("running ",stability.subsamples," subsampling iterations ... ")}
        if (is.null(sr)) {
          sr <- papply(1:stability.subsamples,function(i) subset.clustering(self$graph,f=stability.subsampling.fraction,seed=i),n.cores=self$n.cores)
        }

        if(verbose) { message("done")}

        if(verbose) message("calculating flat stability stats ... ")
        # Jaccard coefficient for each cluster against all, plus random expectation
        jc.stats <- do.call(rbind, papply(sr,function(o) {
          p1 <- membership(o)
          p2 <- cls.groups[names(p1)]
          p1 <- as.character(p1)
          #x <- tapply(1:length(p2),factor(p2,levels=cls.levs),function(i1) {
          x <- tapply(1:length(p2),p2,function(i1) {
            i2 <- which(p1==p1[i1[[1]]])
            length(intersect(i1,i2))/length(unique(c(i1,i2)))
          })

        }, n.cores=self$n.cores, mc.preschedule=TRUE))


        # v0.2.4 on Feb. 3, 2009 by Weiliang Qiu, <https://github.com/cran/clues/blob/master/R/adjustedRand.R>
        #  (1) moved some code out of the for loop
        #
        # cl1 --- partition 1 of the data set
        # cl2 --- partition 2 of the data set
        #
        # flag = 1 --- Rand index
        # flag = 2 --- Hubert and Arabie's adjusted Rand index
        # flag = 3 --- Morey and Agresti's adjusted Rand index
        # flag = 4 --- Fowlkes and Mallows's index
        # flag = 5 --- Jaccard index
        adjustedRand <- function(cl1, cl2, randMethod = c("Rand","HA", "MA", "FM", "Jaccard")){
            if(!is.vector(cl1)){
              stop("cl1 is not a vector!\n")
            }
            if(!is.vector(cl2)){
              stop("cl2 is not a vector!\n")
            }
            if(length(cl1) != length(cl2)){
              stop("Two vectors have different lengths!\n")
            }

            len <- length(randMethod)
            if(len == 0){
              stop("The argument 'randMethod' is empty!\n")
            }

            # unique values of elements in 'cl1'
            cl1u <- unique(cl1)
            # number of clusters in partition 1
            m1 <- length(cl1u)

            # unique values of elements in 'cl2'
            cl2u <- unique(cl2)
            # number of clusters in partition 2
            m2 <- length(cl2u)

            n <- length(cl1)
            randVec <- rep(0, len)
            names(randVec) <- randMethod
            for(i in 1:len){
                randMethod[i] <- match.arg(arg = randMethod[i],
                    choices = c("Rand","HA", "MA", "FM", "Jaccard"))

                flag <- match(randMethod[i],
                    c("Rand","HA", "MA", "FM", "Jaccard"))

                c.res <- .C("adjustedRand",
                    as.integer(cl1),
                    as.integer(cl1u),
                    as.integer(cl2),
                    as.integer(cl2u),
                    as.integer(m1),
                    as.integer(m2),
                    as.integer(n),
                    as.integer(flag),
                    r = as.double(0))
                randVec[i] <- c.res$r
            }
            return(randVec)
        }

        # Adjusted rand index
        if(verbose) message("adjusted Rand ... ")
        ari <- unlist(papply(sr,function(o) { ol <- membership(o); adjustedRand(as.integer(ol),as.integer(cls.groups[names(ol)]),randMethod='HA') }, n.cores=self$n.cores))
        if(verbose) message("done")

        res$stability <- list(flat=list(jc=jc.stats,ari=ari))

        # hierarchical measures
        if(verbose) message("calculating hierarchical stability stats ... ")
        if(is.hierarchical(cls)) {
          # hierarchical to hierarchical stability analysis - cut reference
          # determine hierarchy of clusters (above the cut)
          t.get.walktrap.upper.merges <- function(res,n=length(unique(membership(res)))) {
            clm <- complete.dend(res,FALSE)
            x <- tail(clm,n-1)
            x <- x - 2*nrow(res$merges) + nrow(x)-1
            # now all >=0 ids are cut leafs and need to be reassigned ids according to their rank
            xp <- x+nrow(x)+1
            xp[x<=0] <- rank(-x[x<=0])
            xp
          }

          clm <- t.get.walktrap.upper.merges(cls)
          res$stability$upper.tree <- clm

          if(verbose) message("tree Jaccard ... ")
          jc.hstats <- do.call(rbind, papply(sr,function(z) bestClusterThresholds(z,cls.groups,clm)$threshold, n.cores=self$n.cores))
        } else {
          # compute cluster hierarchy based on cell mixing (and then something)
          # assess stability for that hierarchy (to visualize internal node stability)
          # for the original clustering and every subsample clustering,
          if(verbose) message("upper clustering ... ")
          cgraph <- getClusterGraph(self$graph,cls.groups,plot=FALSE,normalize=FALSE)
          chwt <- walktrap.community(cgraph,steps=9)
          clm <- complete.dend(chwt,FALSE)

          if(verbose) message("clusterTree Jaccard ... ")
          jc.hstats <- do.call(rbind, papply(sr,function(st1) {
            mf <- membership(st1)
            mf <- as.factor(setNames(as.character(mf),names(mf)))

            st1g <- getClusterGraph(self$graph,mf, plot=FALSE, normalize=TRUE)

            st1w <- walktrap.community(st1g, steps=8)

            #merges <- st1w$merge; leaf.factor <- mf; clusters <- cls.groups
            x <- bestClusterTreeThresholds(st1w,mf,cls.groups,clm)
            x$threshold
          }, n.cores=self$n.cores))

        }
        res$stability$upper.tree <- clm
        res$stability$sr <- sr
        res$stability$hierarchical <- list(jc=jc.hstats)
        if(verbose) message("done")

      }

      ## Filter groups
      if(min.group.size>0) {
        lvls.keep <- names(which(table(cls.groups)  > min.group.size))
        cls.groups[! as.character(cls.groups) %in% as.character(lvls.keep)] <- NA
        cls.groups <- as.factor(cls.groups)
      }
      res$groups <- cls.groups
      self$clusters[[name]] <- res
      return(invisible(res))

    },

    #' @description Plot panel of individual embeddings per sample with joint coloring
    #'
    #' @param groups a cell factor (a factor named with cell names) specifying clusters of cells to be compared (one against all). To compare two cell clusters against each other, simply pass a factor containing only two levels (default=NULL, see clustering)
    #' @param clustering character Name of the clustering to use (see names(con$clusters)) for the value of the groups factor (default=NULL - if groups are not specified, the first clustering will be used)
    #' @param use.local.clusters boolean Whether clusters should be taken from the individual samples; otherwise joint clusters in the conos object will be used (see clustering) (default=FALSE).
    #' @param plot.theme string Theme for the plot, passed to plotSamples() (default=NULL)
    #' @param use.common.embedding boolean Whether a joint embedding in the conos object should be used (or embeddings determined for the individual samples) (default=FALSE)
    #' @param embedding (default=NULL) If a character value is passed, it is interpreted as an embedding name (a name of a joint embedding in conos when use.commmon.embedding=TRUE, or a name of an embedding within the individual objects when use.common.embedding=FALSE).
    #'     If a matrix is passed, it is interpreted as an actual embedding (then first two columns are interpreted as x/y coordinates, row names must be cell names). If NULL, the default embedding will be used.
    #' @param adj.list an optional list of additional ggplot2 directions to apply (default=NULL)
    #' @param ... Additional parameters passed to plotSamples(), plotEmbeddings(), sccore::embeddingPlot().
    #' @return cowplot grid object with the panel of plots
    #'
    plotPanel=function(clustering=NULL, groups=NULL, colors=NULL, gene=NULL, use.local.clusters=FALSE, plot.theme=NULL, use.common.embedding=FALSE, embedding=NULL, adj.list=NULL, ...) {

      if (use.local.clusters) {
        if (is.null(clustering) && !(inherits(x = self$samples[[1]], what = c('seurat', 'Seurat')))) {
          stop("You have to provide 'clustering' parameter to be able to use local clusters")
        }

        groups <- lapply(self$samples, getClustering, clustering) %>%
          lapply(function(cls) setNames(as.character(cls), names(cls))) %>% Reduce(c, .)

        if (is.null(groups)) {
          stop(paste0("No clustering '", clustering, "' presented in the samples"))
        }
      } else if (is.null(groups) && is.null(colors) && is.null(gene)) {
        groups <- getClusteringGroups(self$clusters, clustering)
      }

      if (use.common.embedding) { # look up the embedding within the conos object
        ## if use.common.embedding, pass the Conos embedding to plotSamples
        embedding <- private$getEmbedding(embedding)
        adj.list <- c(ggplot2::lims(x=range(embedding[,1]), y=range(embedding[,2])), adj.list)
      }

      plotSamples(self$samples, groups=groups, colors=colors, gene=gene, plot.theme=private$adjustTheme(plot.theme), embedding.type=embedding, adj.list=adj.list, ...)
    },

    #' @description Generate an embedding of a joint graph
    #'
    #' @param method Embedding method (default='largeVis'). Currently 'largeVis' and 'UMAP' are supported.
    #' @param embedding.name character Optional name of the name of the embedding set by user to store multiple embeddings (default: method name)
    #' @param M numeric (largeVis) The number of negative edges to sample for each positive edge to be used (default=1)
    #' @param gamma numeric (largeVis) The strength of the force pushing non-neighbor nodes apart (default=1)
    #' @param alpha numeric (largeVis) Hyperparameter used in the default distance function, \eqn{1 / (1 + \alpha \dot ||y_i - y_j||^2)} (default=0.1).  The function relates the distance
    #'     between points in the low-dimensional projection to the likelihood that the two points are nearest neighbors. Increasing \eqn{\alpha} tends
    #'     to push nodes and their neighbors closer together; decreasing \eqn{\alpha} produces a broader distribution. Setting \eqn{\alpha} to zero
    #'     enables the alternative distance function. \eqn{\alpha} below zero is meaningless.
    #' @param perplexity (largeVis) The perplexity passed to largeVis (default=NA)
    #' @param sgd_batches (largeVis) The number of edges to process during SGD (default=1e8). Defaults to a value set based on the size of the dataset. If the parameter given is
    #'     between \code{0} and \code{1}, the default value will be multiplied by the parameter.
    #' @param seed numeric Random seed for the largeVis algorithm (default=1)
    #' @param target.dims numeric Number of dimensions for the reduction (default=2). Higher dimensions can be used to generate embeddings for subsequent reductions by other methods, such as tSNE
    #' @param ... additional arguments, passed to UMAP embedding (run ?conos:::embedGraphUmap for more info)
    #'
    embedGraph=function(method='largeVis', embedding.name=method, M=1, gamma=1, alpha=0.1, perplexity=NA, sgd_batches=1e8, seed=1, verbose=TRUE, target.dims=2, ...) {
      supported.methods <- c('largeVis', 'UMAP')
      if(!method %in% supported.methods) {
        stop(paste0("Currently, only the following embeddings are supported: ",paste(supported.methods,collapse=' ')))
      }

      ## check if embedding.name already in list
      ## if so, throw warning
      if (length(self$embeddings)>0){
        ## check if embedding.name already created
        if (embedding.name %in% names(self$embeddings)){
          warning(paste0("Already created an embedding: ", embedding.name, ". Overwriting."))
        }
      }

      if (method == 'largeVis') {
        wij <- as_adj(self$graph,attr='weight')
        if(!is.na(perplexity)) {
          wij <- buildWijMatrix(wij,perplexity=perplexity, threads=self$n.cores)
        }
        coords <- projectKNNs(wij = wij, dim=target.dims, verbose = verbose,sgd_batches = sgd_batches,gamma=gamma, M=M, seed=seed, alpha=alpha, rho=1, threads=self$n.cores)
        colnames(coords) <- V(self$graph)$name
        self$embedding <- t(coords)
        embedding.result <- self$embedding
      } else {
        ## method == 'UMAP'
        if (!requireNamespace("uwot", quietly=TRUE)){
          stop("You need to install package 'uwot' to be able to use UMAP embedding. Please install it.")
        }

        self$embedding <- embedGraphUmap(self$graph, verbose=verbose, return.all=FALSE, n.cores=self$n.cores, target.dims=target.dims, ...)
        embedding.result <- self$embedding
      }

      self$embeddings[[embedding.name]] <- embedding.result
      self$embedding <- embedding.result # hang on to the latest embedding for backwards compatibility
      return(invisible(embedding.result))
    },

    #' @description Plot cluster stability statistics.
    #'
    #' @param clustering string Name of the clustering result to show (default=NULL)
    #' @param what string Show a specific plot (ari - adjusted rand index, fjc - flat Jaccard, hjc - hierarchical Jaccard, dend - cluster dendrogram, all - everything except 'dend') (default='all')
    #' @return cluster stability statistics
    plotClusterStability=function(clustering=NULL, what='all') {
      if(is.null(clustering)){
        clustering <- names(self$clusters)[[1]]
      }

      if(is.null(self$clusters[[clustering]])){
        stop(paste("clustering",clustering,"doesn't exist, run findCommunity() first"))
      }

      if(is.null(self$clusters[[clustering]]$stability)){
        stop(paste("clustering",clustering,"doesn't have stability info. Run findCommunity( ... , test.stability=TRUE) first"))
      }

      st <- self$clusters[[clustering]]$stability
      nclusters <- ncol(st$flat$jc)
      jitter.alpha <- 0.1

      if (what=='all' || what=='ari') {
        p.fai <- ggplot2::ggplot(data.frame(aRI=st$flat$ari), ggplot2::aes(x=1,y=aRI)) +
          ggplot2::geom_boxplot(notch=TRUE, outlier.shape=NA) +
          ggplot2::geom_point(shape=16, position = ggplot2::position_jitter(), alpha=jitter.alpha) +
          ggplot2::guides(color=FALSE) +
          ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
          ggplot2::ylim(c(0,1)) + ggplot2::labs(x=" ", y="adjusted Rand Index") +
          ggplot2::theme(legend.position="none", axis.ticks.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank())

        if(what=='ari')
          return(p.fai)
      }

      if (what=='all' || what=='fjc') {
        df <- reshape2::melt(st$flat$jc)
        colnames(df) <- c('rep','cluster','jc')
        df$cluster <- factor(colnames(st$flat$jc)[df$cluster],levels=levels(self$clusters[[clustering]]$groups))

        p.fjc <- ggplot2::ggplot(df,aes(x=cluster,y=jc,color=cluster)) +
          ggplot2::geom_boxplot(aes(color=cluster), notch=TRUE, outlier.shape=NA) +
          ggplot2::geom_jitter(shape=16, position=position_jitter(0.2), alpha=jitter.alpha) +
          ggplot2::guides(color=FALSE) +
          ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
          ggplot2::ylab("Jaccard coefficient (flat)") + ggplot2::ylim(c(0,1))

        if(what=='fjc') return(p.fjc)
      }

      if (what=='all' || what=='hjc') {
        # hierarchical
        df <- reshape2::melt(st$hierarchical$jc[,1:nclusters])
        colnames(df) <- c('rep','cluster','jc')
        df$cluster <- factor(colnames(st$flat$jc)[df$cluster],levels=levels(self$clusters[[clustering]]$groups))
        p.hjc <- ggplot2::ggplot(df,aes(x=cluster,y=jc,color=cluster)) +
          ggplot2::geom_boxplot(aes(color=cluster), notch=TRUE, outlier.shape=NA) +
          ggplot2::geom_jitter(shape=16, position=ggplot2::position_jitter(0.2), alpha=jitter.alpha) +
          ggplot2::guides(color=FALSE) +
          ggplot2::geom_hline(yintercept=1, linetype="dashed", alpha=0.2) +
          ggplot2::ylab("Jaccard coefficient (hierarchical)") + ggplot2::ylim(c(0,1))

        if(what=='hjc') return(p.hjc)
      }

      if (what=='dend') {
        m <- st$upper.tree
        nleafs <- nrow(m)+1
        m[m<=nleafs] <- -1*m[m<=nleafs]
        m[m>0] <- m[m>0]-nleafs
        hc <- list(merge=m, height=1:nrow(m), labels=levels(self$clusters[[clustering]]$groups), order=c(1:nleafs))
        class(hc) <- 'hclust'
        # fix the ordering so that edges don't intersects
        hc$order <- order.dendrogram(as.dendrogram(hc))

        d <- as.dendrogram(hc) %>% dendextend::hang.dendrogram()

        # depth-first traversal of a merge matrix
        t.dfirst <- function(m,i=nrow(m)) {
          rl <- m[i,1]; if(rl<0) { rl <- abs(rl) } else { rl <- t.dfirst(m,rl) }
          rr <- m[i,2]; if(rr<0) { rr <- abs(rr) } else { rr <- t.dfirst(m,rr) }
          c(i+nrow(m)+1,rl,rr)
        }
        xy <- dendextend::get_nodes_xy(d)
        to <- t.dfirst(hc$merge)
        plot(d,las=2,axes=FALSE)
        # flat on the left
        #x <- apply(st$flat$jc,2,median)
        #text(xy,labels=round(x[to],2),col='blue',adj=c(-0.1,-1.24),cex=0.8)
        x <- apply(st$hierarchical$jc,2,median)
        text(xy,labels=round(x[to],2),col='red',adj=c(-0.1,-0.12),cex=0.8)
        return(NULL)
      }

      cowplot::plot_grid(plotlist=list(p.fai,p.fjc,p.hjc),nrow=1,rel_widths=c(4,nclusters,nclusters))
    },

    #' @description Plot joint graph
    #'
    #' @param groups a cell factor (a factor named with cell names) specifying clusters of cells to be compared (one against all). To compare two cell clusters against each other, simply pass a factor containing only two levels (default: NULL, see clustering)
    #' @param clustering a character name of the clustering to use (see names(con$clusters)) for the value of the groups factor (default: NULL - if groups are not specified, the first clustering will be used)
    #' @param color.by character A shortcut to color the plot by 'cluster' or by 'sample' (default: 'cluster'). If any other string is input, an error is thrown.
    #' @param embedding A character name of an embedding, or a matrix of the actual embedding (rownames should correspond to cells, first to columns to x/y coordinates). If NULL (default: NULL), the latest generated embedding will be used
    #' @param colors a color factor (named with cell names) use for cell coloring (default=NULL)
    #' @param gene Show expression of a gene (default=NULL)
    #' @param plot.theme Theme for the plot, passed to sccore::embeddingPlot() (default=NULL)
    #' @param subset A subset of cells to show (default: NULL - shows all the cells)
    #' @param ... Additional parameters passed to sccore::embeddingPlot()
    #' @return ggplot2 plot of joint graph
    #'
    plotGraph=function(color.by='cluster', clustering=NULL, embedding=NULL, groups=NULL, colors=NULL, gene=NULL, plot.theme=NULL, subset=NULL, ...) {
      embedding <- private$getEmbedding(embedding)

      if (!is.null(subset)) {
        embedding <- embedding[rownames(embedding) %in% subset,,drop=FALSE]
      }

      if (!is.null(gene)) {
        colors <- lapply(self$samples, getGeneExpression, gene) %>% Reduce(c, .)
        if(all(is.na(colors))) stop(paste("Gene", gene,"is not found in any of the samples"))
      }

      if(is.null(groups) && is.null(colors)) {
        if(color.by == 'cluster') {
          groups <- getClusteringGroups(self$clusters, clustering)
        } else if(color.by == 'sample') {
          groups <- self$getDatasetPerCell()
        } else {
          stop('Supported values of color.by are ("cluster" and "sample"); Use groups/colors parameters to explicitly pass factor/numeric data for coloring')
        }
      }

      return(embeddingPlot(embedding, groups=groups, colors=colors, plot.theme=private$adjustTheme(plot.theme), ...))
    },

    #' @description Smooth expression of genes to minimize the batch effect between samples
    #'   Use diffusion of expression on graph with the equation dv = exp(-a * (v + b))
    #'
    #' @param genes List of genes to be smooothed smoothing (default=NULL will smooth top n.od.genes overdispersed genes)
    #' @param n.od.genes numeric If 'genes' is NULL, top n.od.genes of overdispersed genes are taken across all samples (default=500)
    #' @param fading numeric Level of fading of expression change from distance on the graph (parameter 'a' of the equation) (default=10)
    #' @param fading.const numeric Minimal penalty for each new edge during diffusion (parameter 'b' of the equation) (default=0.5)
    #' @param max.iters numeric Maximal number of diffusion iterations (default=15)
    #' @param tol numeric Tolerance after which the diffusion stops (default=5e-3)
    #' @param name string Name to save the correction (default='diffusion')
    #' @param verbose boolean Verbose mode (default=TRUE)
    #' @param count.matrix Alternative gene count matrix to correct (rows: genes, columns: cells; has to be dense matrix). Default: joint count matrix for all datasets.
    #' @param normalize boolean Whether to normalize values (default=TRUE)
    #' @return smoothed expression of the input genes
    #' 
    correctGenes=function(genes=NULL, n.od.genes=500, fading=10.0, fading.const=0.5, max.iters=15, tol=5e-3, name='diffusion', verbose=TRUE, count.matrix=NULL, normalize=TRUE) {
      edges <- igraph::as_edgelist(self$graph)
      edge.weights <- igraph::edge.attributes(self$graph)$weight

      if (is.null(count.matrix)) {
        if (is.null(genes)) {
          genes <- getOdGenesUniformly(self$samples, n.genes=n.od.genes)
        }

        cms <- lapply(self$samples, getCountMatrix, transposed=TRUE)
        genes <- Reduce(intersect, lapply(cms, colnames)) %>% intersect(genes)
        count.matrix <- Reduce(rbind, lapply(cms, function(x) x[, genes])) %>% as.matrix()
      } else {
        count.matrix <- t(count.matrix)
      }
      vn <- V(self$graph)$name
      if(!all(rownames(count.matrix)==vn)) { # subset to a common set of genes
        if(!all(vn %in% rownames(count.matrix))) {
          stop("count.matrix does not provide values for all the vertices in the alignment graph!")
        }
        count.matrix <- count.matrix[vn,]
      }

      ## Wrapper to make is.label.fixed optional
      smoothMatrixOnGraph <- function(edges, edge.weights, matrix, is.label.fixed=logical(), ...) {
        smooth_count_matrix(edges, edge.weights, matrix, is_label_fixed=is.label.fixed, ...)
      }

      cm <- smoothMatrixOnGraph(edges, edge.weights, count.matrix, max_n_iters=max.iters, diffusion_fading=fading,
                              diffusion_fading_const=fading.const, verbose=verbose, normalize=normalize)
      return(invisible(self$expression.adj[[name]] <<- cm))
    },

    #' @description Estimate labeling distribution for each vertex, based on a partial labeling of the cells.
    #' There are two methods used for the propagation to calculate the distribution of labels: "solver" and "diffusion".
    #' * "diffusion" (default) will estimate the labeling distribution for each vertex, based on provided labels using a random walk.
    #' * "solver" will propagate labels using the algorithm described by Zhu, Ghahramani, Lafferty (2003) <http://mlg.eng.cam.ac.uk/zoubin/papers/zgl.pdf>
    #' Confidence values are then calculated by taking the maximum value from this distribution of labels, for each cell.
    #'
    #' @param labels Input labels
    #' @param method type of propagation. Either 'diffusion' or 'solver'. 'solver' gives better result
    #'  but has bad asymptotics, so is inappropriate for datasets > 20k cells. (default='diffusion')
    #' @param ... additional arguments for conos:::propagateLabels* functions
    #' @return list with three fields:
    #' * labels = matrix with distribution of label probabilities for each vertex by rows.
    #' * uncertainty = 1 - confidence values
    #' * label.distribution = the distribution of labels calculated using either the methods "diffusion" or "solver"
    #'
    propagateLabels=function(labels, method="diffusion", ...) {
      if (method == "solver") {
        label.dist <- propagateLabelsSolver(self$graph, labels, ...)
      } else if (method == "diffusion") {
        label.dist <- propagateLabelsDiffusion(self$graph, labels, ...)
      } else {
        stop("Unknown method: ", method, ". Only 'solver' and 'diffusion' are supported.")
      }

      labels <- colnames(label.dist)[apply(label.dist, 1, which.max)] %>%
        setNames(rownames(label.dist))

      confidence <- apply(label.dist, 1, max) %>% setNames(rownames(label.dist))

      return(list(labels=labels, uncertainty=(1 - confidence), label.distribution=label.dist))
    },

    #' @description Calculate pseudo-bulk expression matrices for clusters (by adding up, for each gene, all of the molecules detected for all cells in a given cluster in a given sample)
    #'
    #' @param common.genes boolean Whether to bring individual sample matrices to a common gene list (default=TRUE)
    #' @param omit.na.cells boolean If set to FALSE, the resulting matrices will include a first column named 'NA' that will report total molecule counts for all of the cells that were not covered by the provided factor. (default=TRUE)
    #' @return a list of per-sample uniform dense matrices with rows being genes, and columns being clusters
    #'
    getClusterCountMatrices=function(clustering=NULL, groups=NULL, common.genes=TRUE, omit.na.cells=TRUE) {
      if(is.null(groups)) {
        groups <- getClusteringGroups(self$clusters, clustering)
      }

      groups <- as.factor(groups)

      matl <- lapply(self$samples,function(s) {
        m <- getRawCountMatrix(s,trans=TRUE) # rows are cells
        cl <- factor(groups[match(rownames(m),names(groups))],levels=levels(groups));
        tc <- colSumByFactor(m,cl)
        if(omit.na.cells) { tc <- tc[-1,,drop=FALSE] }
        t(tc)
      })
      # bring to a common gene space
      if(common.genes) {
        gs <- unique(unlist(lapply(matl,rownames)))
        matl <- lapply(matl,function(m) {
          nm <- matrix(0,nrow=length(gs),ncol=ncol(m))
          colnames(nm) <- colnames(m)
          rownames(nm) <- gs
          mi <- match(rownames(m),gs)
          nm[mi,] <- m
          nm
        })
      }
      matl
    },

    #' @description applies 'getCellNames()' on all samples
    #' @return list of cellnames for all samples
    #' @examples
    #' con <- Conos$new(small_panel.preprocessed, n.cores=1)
    #' con$getDatasetPerCell()
    #'
    getDatasetPerCell=function() {
      getSampleNamePerCell(self$samples)
    },

    #' @description Retrieve joint count matrices
    #'
    #' @param raw boolean If TRUE, return merged "raw" count matrices, using function getRawCountMatrix(). Otherwise, return the merged count matrices, using getCountMatrix(). (default=FALSE)
    #' @return list of merged count matrices
    #' @examples
    #' con <- Conos$new(small_panel.preprocessed, n.cores=1)
    #' con$getJointCountMatrix()
    #'
    getJointCountMatrix=function(raw=FALSE) {
      lapply(self$samples, (if (raw) getRawCountMatrix else getCountMatrix), transposed=TRUE) %>%
        mergeCountMatrices(transposed=TRUE)
    }
  ),



  ## Private functions
  private = list(
    adjustTheme=function(theme) {
      if (is.null(theme)) {
        theme <- ggplot2::theme()
      }
      main.theme <- ggplot2::theme_bw() + ggplot2::theme(
        legend.background=ggplot2::element_rect(fill=ggplot2::alpha("white", 0.6)),
        plot.margin=ggplot2::margin()
      )

      if (self$override.conos.plot.theme) {
        return(main.theme + ggplot2::theme_get() + theme)
      }

      return(main.theme + theme)
    },

    # a utility function to look up an embedding by name or accept an actual embedding data
    getEmbedding=function(embedding) {
      if (!is.null(embedding)) {
        if(class(embedding) %in% c('matrix')) { # actuall embedding was passed
          # check validity?
        } else if(class(embedding)=='character') {  # look up embedding by name
          ## check if embedding.name exists in list
          if (embedding %in% names(self$embeddings)) {
            ## embedding to plot
            embedding <- self$embeddings[[embedding]]
          } else {
            ## embedding.name not in list of self$embeddings, so the user is confused
            ## throw error
            stop(paste0("No embedding named '", embedding, "' found. Please generate this with embedGraph()."))
          }
        } else {
          stop('embedding must be either a character name of the embedding, or an actual matrix of embedding coordinates')
        }
      } else {
        if(!is.null(self$embedding)) { # use the latest
          embedding <- self$embedding
        }  else { # pick one from the list
          if(is.null(self$embeddings) || length(self$embeddings)<1) stop("no joint embeddings have been generated; use embedGraph() first")
          embedding <- self$embedding[length(self$embedding)] # by default, pick last-named embedding
        }
      }
      return(embedding)
    },

    updatePairs=function(space='PCA', data.type='counts', ncomps=50, n.odgenes=1e3, var.scale=TRUE, matching.mask=NULL, exclude.samples=NULL, score.component.variance=FALSE, verbose=FALSE) {

      # make a list of all pairs
      sample.names <- names(self$samples)
      if(!is.null(exclude.samples)) {
        mi <- sample.names %in% exclude.samples
        if(verbose) { message("excluded ", sum(mi), " out of ", length(sample.names), " samples, based on supplied exclude.samples") }
        sample.names <- sample.names[!mi]
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

        if(verbose) message("Use ", ncol(sn.pairs), " pairs, based on the passed exclude.pairs")
      } else {
        sn.pairs <- combn(sample.names, 2)
      }

      # determine the pairs that need to be calculated
      if (is.null(self$pairs[[space]])) {
        self$pairs[[space]] <- list()
      }
      mi <- rep(NA,ncol(sn.pairs))
      nm <- match(apply(sn.pairs,2,paste,collapse='.vs.'),names(self$pairs[[space]]))
      mi[which(!is.na(nm))] <- na.omit(nm)
      # try reverse match as well
      nm <- match(apply(sn.pairs[c(2,1),,drop=FALSE],2,paste,collapse='.vs.'),names(self$pairs[[space]]));
      mi[which(!is.na(nm))] <- na.omit(nm)
      if(verbose) message('found ',sum(!is.na(mi)),' out of ',length(mi),' cached ',space,' space pairs ... ')
      if(any(is.na(mi))) { # some pairs are missing
        if(verbose) message('running ',sum(is.na(mi)),' additional ',space,' space pairs ')
        xl2 <- sccore::plapply(which(is.na(mi)), function(i) {
          if(space=='CPCA') {
            xcp <- quickCPCA(self$samples[sn.pairs[,i]],data.type=data.type,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale, score.component.variance=score.component.variance)
          } else if(space=='JNMF') {
            xcp <- quickJNMF(self$samples[sn.pairs[,i]],data.type=data.type,n.comps=ncomps,n.odgenes=n.odgenes,var.scale=var.scale,verbose=FALSE,max.iter=3e3)
          } else if (space == 'genes') {
            xcp <- quickNULL(p2.objs = self$samples[sn.pairs[,i]], data.type=data.type, n.odgenes=n.odgenes, var.scale = var.scale, verbose = FALSE)
          } else if (space == 'PCA') {
            xcp <- quickPlainPCA(self$samples[sn.pairs[,i]], data.type=data.type,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale, score.component.variance=score.component.variance)
          } else if (space == 'CCA' || space=='PMA') {
            xcp <- quickCCA(self$samples[sn.pairs[,i]],data.type=data.type,ncomps=ncomps,n.odgenes=n.odgenes,verbose=FALSE,var.scale=var.scale, score.component.variance=score.component.variance,PMA=(space=='PMA'))
          }
          if(verbose) cat(".")
          xcp
        }, n.cores=self$n.cores, mc.preschedule=(space=='PCA'), progress=FALSE)

        names(xl2) <- apply(sn.pairs[,which(is.na(mi)),drop=FALSE],2,paste,collapse='.vs.')
        xl2 <- xl2[!unlist(lapply(xl2,is.null))]
        self$pairs[[space]] <- c(self$pairs[[space]],xl2)
      }

      # re-do the match and order
      mi <- rep(NA,ncol(sn.pairs))
      nm <- match(apply(sn.pairs,2,paste,collapse='.vs.'),names(self$pairs[[space]]))
      mi[which(!is.na(nm))] <- na.omit(nm)
      nm <- match(apply(sn.pairs[c(2,1),,drop=FALSE],2,paste,collapse='.vs.'),names(self$pairs[[space]]))
      mi[which(!is.na(nm))] <- na.omit(nm)
      if(any(is.na(mi))) {
        warning("unable to get complete set of pair comparison results")
        sn.pairs <- sn.pairs[,!is.na(mi),drop=FALSE]
      }
      if(verbose) message(" done")
      return(invisible(sn.pairs))
    }
  )
)
