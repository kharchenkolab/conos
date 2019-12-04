#' @importFrom dplyr %>%
NULL

#' @export
sccore::embeddingPlot

##' Extract specified clustering from list of conos clusterings
##'
##' @param clusters list of conos clusterings
##' @param clustering name of extracted clustering
##' @return vector of clusters, named with cell names
getClusteringGroups <- function(clusters, clustering) {
  if (length(clusters) < 1) {
    stop("generate a joint clustering first")
  }

  if (is.null(clustering)) { # take the first one
    return(clusters[[1]]$groups)
  }

  if (is.null(clusters[[clustering]])) {
    stop(paste("clustering", clustering, "hasn't been calculated"))
  }

  return(clusters[[clustering]]$groups)
}

##' Plot panel of specified embeddings
##'
##' @inheritParams sccore::embeddingPlot
##' @param embeddings list of two-column matrices with (x, y) coordinates of the embeddings. Each mutrix must have cell names in rownames.
##' @param ncol number of columns in the panel
##' @param nrow number of rows in the panel
##' @param panel.size vector with two numbers, which specified (width, height) of the panel in inches. Ignored if raster == FALSE.
##' @param adjust.func function to adjust plots before combining them to single panel. Can be used, for example, to provide color pallette of guides of the plots.
##' @return ggplot2 object with the panel of plots
plotEmbeddings <- function(embeddings, groups=NULL, colors=NULL, ncol=NULL, nrow=NULL, raster=FALSE, panel.size=NULL, adjust.func=NULL, title.size=6,
                           raster.width=NULL, raster.height=NULL, adj.list=NULL, ...) {
  if (is.null(panel.size)) {
    panel.size <- dev.size(units="in")
  } else if (length(panel.size) == 1) {
    panel.size <- c(panel.size, panel.size)
  }

  n.plots <- length(embeddings)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n.plots))
  }
  if (is.null(ncol)) {
    ncol <- ceiling(n.plots / nrow)
  }
  if (is.null(nrow)) {
    nrow <- ceiling(n.plots / ncol)
  }

  if (is.null(names(embeddings))) {
    names(embeddings) <- paste(1:length(embeddings))
  }

  if (is.null(raster.width)) {
    raster.width <- panel.size[1] / nrow
  }

  if (is.null(raster.height)) {
    raster.height <- panel.size[2] / ncol
  }

  plot.list <- lapply(names(embeddings), function(n)
    embeddingPlot(embeddings[[n]], groups=groups, colors=colors, raster=raster,
                  raster.width=raster.width, raster.height=raster.height, ...) +
      ggplot2::geom_label(data=data.frame(x=-Inf, y=Inf, label=n), mapping=ggplot2::aes(x=x, y=y, label=label),
                          fill=ggplot2::alpha("white", 0.6), hjust=0, vjust=1, size=title.size,
                          label.padding=ggplot2::unit(title.size / 4, "pt"), label.size = NA)
    )

  for (adj in adj.list) {
    plot.list %<>% lapply(`+`, adj)
  }

  if (!is.null(adjust.func)) {
    plot.list <- lapply(plot.list, adjust.func)
  }

  return(cowplot::plot_grid(plotlist=plot.list, ncol=ncol, nrow=nrow))
}

##' Plot panel of specified embeddings, extracting them from pagoda2 objects
##'
##' @inheritParams plotEmbeddings
##' @param samples list of pagoda2 or Seurat objects
##' @param gene gene name. If this parameter is provided, points are colored by expression of this gene.
##' @param embedding.type type of embedding. Default: tSNE. If a numeric matrix is passed, interpreted as an actual embedding.
##' @return ggplot2 object with the panel of plots
plotSamples <- function(samples, groups=NULL, colors=NULL, gene=NULL, embedding.type=NULL, ...) {
  if (!is.null(groups)) {
    groups <- as.factor(groups)
  }

  if (is.null(x = embedding.type)) {
    embedding.type <- if (inherits(x = samples[[1]], what = c('seurat', 'Seurat'))) {
      'tsne'
    } else {
      'tSNE'
    }
  }
  if(class(embedding.type)=='matrix') { # actual embedding was passed
    embeddings <- lapply(samples,function(r) embedding.type[rownames(embedding.type) %in% getCellNames(r),,drop=F])
    embeddings <- embeddings[unlist(lapply(embeddings,function(x) nrow(x)>0))]
  } else { # extract embeddings from samples
    embeddings <- lapply(samples, getEmbedding, embedding.type)
  }
  no.embedding <- sapply(embeddings, is.null)
  if (all(no.embedding)) {
    stop(paste0("No '", embedding.type, "' embedding presented in the samples"))
  }
  if (any(no.embedding)) {
    warning(paste0(sum(no.embedding), " of your samples doesn't have '", embedding.type, "' embedding"))
    embeddings <- embeddings[!no.embedding]
  }
  if (!is.null(gene)) {
    colors <- lapply(samples, getGeneExpression, gene) %>% Reduce(c, .)
  }
  return(plotEmbeddings(embeddings, groups = groups, colors = colors, ...))
}

#' Plots barplots per sample of composition of each pagoda2 application based on
#' selected clustering
#' @param conos.obj A conos object
#' @param clustering name of clustering in the current object
#' @param groups arbitrary grouping of cells (to use instead of the clustering)
#' @param sample.factor a factor describing cell membership in the samples (or some other category); will default to samples if not provided
#' @param show.entropy whether to include entropy barplot
#' @param show.size whether to include size barplot
#' @param show.composition whether to include composition barplot
#' @param legend.height relative hight of the legend panel
#' @return a ggplot object
#' @export
plotClusterBarplots <- function(conos.obj=NULL, clustering=NULL, groups=NULL,sample.factor=NULL,show.entropy=TRUE,show.size=TRUE,show.composition=TRUE,legend.height=0.2) {
  ## param checking
  if(!is.null(clustering)) {
    if(is.null(conos.obj)) stop('conos.obj must be passed if clustering name is specified');
    if(!clustering %in% names(conos.obj$clusters)) stop('specified clustering doesn\'t exist')
    groups <- conos.obj$clusters[[clustering]]$groups
  } else if (is.null(groups)) {
    if(is.null(conos.obj)) stop('either groups factor on the cells or a conos object needs to be specified')
    if(is.null(conos.obj$clusters[[1]])) stop('conos object lacks any clustering. run $findCommunities() first')
    groups <- conos.obj$clusters[[1]]$groups
  }

  groups <- as.factor(groups)
  if(is.null(sample.factor)) {
    sample.factor <- conos.obj$getDatasetPerCell(); # assignment to samples
  }

  xt <- table(sample.factor[match(names(groups),names(sample.factor))],groups)
  xt <- xt[rowSums(xt)>0,]; xt <- xt[,colSums(xt)>0]

  df <- reshape2::melt(xt); colnames(df) <- c("sample","cluster","f");  df$f <- df$f/colSums(xt)[as.character(df$cluster)]
  clp <- ggplot2::ggplot(df, ggplot2::aes(x=factor(cluster, levels=levels(groups)),y=f,fill=sample)) +
    ggplot2::geom_bar(stat='identity') + ggplot2::xlab('cluster') + ggplot2::ylab('fraction of cells') + ggplot2::theme_bw() +
    ggplot2::scale_y_continuous(expand=c(0, 0))

  if(!show.size && !show.entropy)
    return(clp);

  # extract legend
  leg <- cowplot::get_legend(clp + ggplot2::theme(legend.position="bottom"))
  pl <- list(clp + ggplot2::theme(legend.position="none"));

  if(show.entropy) {
    if (!requireNamespace("entropy", quietly=T))
      stop("You need to install 'entropy' package to use 'show.entropy=T'")

    n.samples <- nrow(xt);
    ne <- 1-apply(xt, 2, entropy::KL.empirical, y2=rowSums(xt), unit=c('log2')) / log2(n.samples) # relative entropy
    enp <- ggplot2::ggplot(data.frame(cluster=factor(colnames(xt),levels=levels(groups)),entropy=ne), ggplot2::aes(cluster, entropy)) +
      ggplot2::geom_bar(stat='identity',fill='grey65') + ggplot2::ylim(0,1) +
      ggplot2::geom_hline(yintercept=1, linetype="dashed", color = "grey30") + ggplot2::theme_bw()
    pl <- c(pl,list(enp))
  }

  if(show.size) {
    szp <- ggplot2::ggplot(data.frame(cluster=factor(colnames(xt),levels=levels(groups)), cells=colSums(xt)), ggplot2::aes(cluster,cells)) +
      ggplot2::geom_bar(stat='identity') + ggplot2::scale_y_continuous(trans='log10') + ggplot2::theme_bw() + ggplot2::ylab('number of cells')
    pl <- c(pl,list(szp))
  }

  pp <- cowplot::plot_grid(plotlist=pl,ncol=1,rel_heights=c(1,rep(0.3,length(pl)-1)))
  pp2 <- cowplot::plot_grid(leg,pp,ncol=1,rel_heights=c(legend.height,1))

  return(pp2)
}


#' Generate boxplot per cluster of the proportion of cells in each celltype
#' @param conos.obj conos object
#' @param clustering name of the clustering to use
#' @param apptypes a factor specifying how to group the samples
#' @param return.details if TRUE return a list with the plot and the summary data.frame
plotClusterBoxPlotsByAppType <- function(conos.obj, clustering=NULL, apptypes=NULL, return.details=FALSE) {
    type <- 'proportions'
    ## param checking
    if(is.null(clustering)) {
      clustering <- 'multi level'
    }
    if(is.null(apptypes)) stop('apptypes must be spectified')
    if(!is.factor(apptypes)) stop('apptypes must be a factor')
    if(!type %in% c('counts','proportions')) stop('argument type must be either counts or proportions')

    ## main function
    groups <- as.factor(conos.obj$clusters[[clustering]]$groups)
    plot.df <- do.call(rbind,lapply(names(conos.obj$samples), function(n) {
        o <- conos.obj$samples[[n]]
        grps1 <- groups[intersect(names(groups), rownames(o$counts))]
        tbl1 <- data.frame(
            clname=levels(grps1),
            val=tabulate(grps1),
            sample=c(n),
            stringsAsFactors = FALSE
        )
        if(type=='proportions') {
            tbl1$val <- tbl1$val / sum(tbl1$val)
        }
        tbl1
    }))
    ## append app type
    plot.df$apptype <- apptypes[plot.df$sample]
    ## Make the plot
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=apptype,y=val,fill=clname)) + ggplot2::facet_wrap(~clname) + ggplot2::geom_boxplot()
    if (type == 'counts') {
        gg <- gg + ggplot2::scale_y_continuous(name='counts')
    } else {
        gg <- gg + ggplot2::scale_y_continuous(name='% of sample')
    }
    gg <- gg + ggplot2::scale_x_discrete(name='cluster')
    ## return
    if(return.details)
        return(list(plot=gg,data=plot.df))

    return(gg)
}


#' Get markers for global clusters
#' @param conos.obj conos object
#' @param clustering name of the clustering to use
#' @param min.samples.expressing minimum number of samples that must have the genes upregulated in the respective cluster
#' @param min.percent.samples.expression minumum percent of samples that must have the gene upregulated
getGlobalClusterMarkers <- function(conos.obj, clustering='multi level',
                                    min.samples.expressing=0,min.percent.samples.expressing=0){
  .Deprecated("getDifferentialGenes")
    ## get the groups from the clusters
    groups <- as.factor(conos.obj$clusters[[clustering]]$groups)
    ## de lists
    delists <- lapply(conos.obj$samples, function(p2) {
        cells <- rownames(p2$counts)
        groups.p2 <- groups[cells]
        de <- p2$getDifferentialGenes(groups=groups.p2)
        de
    })
    ## get de genes per app
    z <- lapply(namedLevels(groups), function(l) {
        lapply(delists, function(x) {
            res <- x[[l]]
            rownames(res[res$Z > 0,])
        })
    })
    ## get consistent genes for each cluster
    zp <- lapply(z, function(k) {
        k <- lapply(k, unique)
        gns <- factor(unlist(unname(k)))
        t.gns <- tabulate(gns)
        names(t.gns) <- levels(gns)
        t.gns.pc <- t.gns / length(k)
        ## consistent genes
        names(t.gns[t.gns >= min.samples.expressing & t.gns.pc >= min.percent.samples.expressing])
    })
    ## return consistent genes
    zp
}

##' Plot fraction of variance explained by the successive reduced space components (PCA, CPCA)
##'
##' Requires buildGraph() or updatePairs() to be ran first with the argument score.component.variance=TRUE.
##'
##' @title Plot variance explained by the successive components
##' @param conos.obj conos object
##' @param space reduction space to be analyzed (currently, component variance scoring is only supported by PCA and CPCA)
##' @return ggplot
##' @export
plotComponentVariance <- function(conos.obj, space='PCA',plot.theme=theme_bw()) {
  pairs <- conos.obj$pairs[[space]]

  if(!is.null(pairs[[space]])) stop(paste("no pairs for space",space,"found. Please run buildGraph() or updatePairs() first, with score.component.variance=TRUE"))

  nvs <- lapply(pairs,'[[','nv'); nvs <- setNames(unlist(nvs,recursive=F,use.names=F),unlist(lapply(nvs,names)))
  if(length(nvs)<1) stop("no variance information found. Please run buildGraph() or updatePairs() with score.component.variance=TRUE")
  if(space=='PCA') { # omit duplicates
    nvs <- nvs[unique(names(nvs))]
  }

  df <- reshape2::melt(do.call(cbind,nvs))
  colnames(df) <- c('component','dataset','var')
  df$component <- factor(df$component,levels=sort(unique(df$component)))

  ggplot2::ggplot(df, ggplot2::aes(x=component,y=var)) +
    ggplot2::geom_point(shape=16, ggplot2::aes(color=dataset), position = ggplot2::position_jitter(), alpha=0.3) +
    ggplot2::geom_line(ggplot2::aes(group=dataset,color=dataset), alpha=0.2)+
    ggplot2::ylab('fraction of variance explained') + ggplot2::xlab('component number') +
    ggplot2::geom_boxplot(notch=F,outlier.shape=NA,fill=NA) + plot.theme + ggplot2::theme(legend.position='none')

}

plotDEGenes <- function(de.genes, samples, groups, n.genes.to.show=10, inner.clustering=FALSE, gradient.range.quantile=0.95) {
  # genes to show
  vi <- unlist(lapply(de.genes, class))=='data.frame';
  x <- lapply(de.genes[vi],function(d) {  if(!is.null(d) && nrow(d)>0) { d[1:min(nrow(d),n.genes.to.show),] } else { NULL } })
  x <- lapply(x,rownames);
  genes <- unique(unlist(x))
  # make expression matrix
  cl <- lapply(samples, function(y) getCountMatrix(y) %>% .[rownames(.) %in% genes,,drop=F] )

  em <- mergeCountMatrices(cl, transposed=T)

  # renormalize rows
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
    clgo <- tapply(1:nrow(em),rowfac,function(ii) ii[hclust(as.dist(1-cor(t(em[ii,]))),method='complete')$order])
    clco <- tapply(1:ncol(em),groups[colnames(em)],function(ii) {
      if(length(ii)>3) ii[hclust(as.dist(1-cor(em[,ii,drop=F])),method='complete')$order] else ii
    })
  } else {
    clclo <- 1:length(levels(rowfac))
    clgo <- tapply(1:nrow(em),rowfac,I)
    clco <- tapply(1:ncol(em),groups[colnames(em)],I)
  }

  #clco <- clco[names(clgo)]
  # filter down to the clusters that are included
  #vic <- cols %in% clclo
  colors <- fac2col(groups[colnames(em)],v=0.95,s=0.95,return.details=TRUE)
  samf <- fac2col(getSampleNamePerCell(samples),v=0.75,s=0.9,return.details=TRUE);
  cellcols <- colors$colors[unlist(clco[clclo])]
  samfcols <- samf$colors[unlist(clco[clclo])]
  genecols <- rev(rep(colors$palette,unlist(lapply(clgo,length)[clclo])))
  drawGroupNames <- FALSE;
  bottomMargin <- ifelse(drawGroupNames,4,0.5);

  # browser()

  #pagoda2:::my.heatmap2(em[rev(unlist(clgo[clclo])),unlist(clco[clclo])],col=gradientPalette,Colv=NA,Rowv=NA,labRow=NA,labCol=NA,RowSideColors=genecols,ColSideColors=rbind(samfcols,cellcols),margins=c(bottomMargin,0.5),ColSideColors.unit.vsize=0.05,RowSideColors.hsize=0.05,useRaster=TRUE, box=TRUE)

  pagoda2:::my.heatmap2(em[rev(unlist(clgo[clclo])),unlist(clco[clclo])],col=gradientPalette,Colv=NA,Rowv=NA,labRow=NA,labCol=NA,RowSideColors=genecols,ColSideColors=rbind(samfcols,cellcols),margins=c(bottomMargin,0.5),ColSideColors.unit.vsize=0.05,RowSideColors.hsize=0.05,useRaster=TRUE, box=TRUE)
  abline(v=cumsum(unlist(lapply(clco[clclo],length))),col=1,lty=3)
  abline(h=cumsum(rev(unlist(lapply(clgo[clclo],length))))+0.5,col=1,lty=3)
}
