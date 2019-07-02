#' @importFrom dplyr %>%
NULL

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

getGeneExpression <- function(count.matrix, gene) {
  if(gene %in% rownames(count.matrix)) {
    return(count.matrix[gene,])
  }

  return(stats::setNames(rep(NA, ncol(count.matrix)), colnames(count.matrix)))
}

##' Plot panel of specified embeddings
##'
##' @inheritParams embeddingPlot
##' @param embeddings list of two-column matrices with (x, y) coordinates of the embeddings. Each mutrix must have cell names in rownames.
##' @param ncol number of columns in the panel
##' @param nrow number of rows in the panel
##' @param panel.size vector with two numbers, which specified (width, height) of the panel in inches. Ignored if raster == FALSE.
##' @param adjust.func function to adjust plots before combining them to single panel. Can be used, for example, to provide color pallette of guides of the plots.
##' @return ggplot2 object with the panel of plots
plotEmbeddings <- function(embeddings, groups=NULL, colors=NULL, ncol=NULL, nrow=NULL, raster=FALSE, panel.size=NULL, adjust.func=NULL, title.size=6,
                           raster.width=NULL, raster.height=NULL, ...) {
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
##' @param embedding.type type of embedding. Default: tSNE.
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
  embeddings <- lapply(samples, getEmbedding, embedding.type)
  no.embedding <- sapply(embeddings, is.null)
  if (all(no.embedding)) {
    stop(paste0("No '", embedding.type, "' embedding presented in the samples"))
  }
  if (any(no.embedding)) {
    warning(paste0(sum(no.embedding), " of your samples doesn't have '", embedding.type, "' embedding"))
    embeddings <- embeddings[!no.embedding]
  }
  if (!is.null(gene)) {
    colors <- lapply(samples, getCountMatrix) %>% lapply(getGeneExpression, gene) %>% Reduce(c, .)
  }
  return(plotEmbeddings(embeddings, groups = groups, colors = colors, ...))
}

##' Plot embedding with provided labels / colors using ggplot2
##'
##' @inheritDotParams ggrepel::geom_label_repel
##' @param embedding two-column matrix with x and y coordinates of the embedding, rownames contain cell names and are used to match coordinates with groups or colors
##' @param groups vector of cluster labels, names contain cell names
##' @param colors vector of numbers, which must be shouwn with point colors, names contain cell names. This argument is ignored if groups are provided.
##' @param plot.na plot points, for which groups / colors are missed (TRUE / FALSE)
##' @param min.cluster.size labels for all groups with number of cells fewer than this parameter are considered as missed. This argument is ignored if groups aren't provided
##' @param mark.groups plot cluster labels above points
##' @param show.legend show legend
##' @param alpha opacity level [0; 1]
##' @param size point size
##' @param title plot title
##' @param plot.theme theme for the plot
##' @param palette function, which accepts number of colors and return list of colors (i.e. see colorRampPalette)
##' @param color.range controls range, in which colors are estimated. Pass "all" to estimate range based on all values of "colors", pass "data" to estimate it only based on colors, presented in the embedding. Alternatively you can pass vector of length 2 with (min, max) values.
##' @param font.size font size for cluster labels. It can either be single number for constant font size or pair (min, max) for font size depending on cluster size
##' @param show.ticks show ticks and tick labels
##' @param legend.position vector with (x, y) positions of the legend
##' @param legend.title legend title
##' @param gradient.range.quantile Winsorization quantile for the numeric colors and gene gradient
##' @param raster should layer with the points be rasterized (TRUE/ FALSE)? Setting of this argument to TRUE is useful when you need to export a plot with large number of points
##' @param raster.width width of the plot in inches. Ignored if raster == FALSE.
##' @param raster.height height of the plot in inches. Ignored if raster == FALSE.
##' @param raster.dpi dpi of the rasterized plot. Ignored if raster == FALSE.
##' @param shuffle.colors shuffle colors
##' @return ggplot2 object
##' @export
embeddingPlot <- function(embedding, groups=NULL, colors=NULL, plot.na=TRUE, min.cluster.size=0, mark.groups=TRUE,
                          show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, plot.theme=NULL, palette=NULL, color.range="all",
                          font.size=c(3, 7), show.ticks=FALSE, show.labels=FALSE, legend.position=NULL, legend.title=NULL,
                          gradient.range.quantile=1, raster=FALSE, raster.width=NULL, raster.height=NULL, raster.dpi=300,
                          shuffle.colors=FALSE,
                          ...) {
  if(is.numeric(colors) && gradient.range.quantile < 1) {
    x <- colors;
    zlim <- as.numeric(quantile(x, p=c(1 - gradient.range.quantile, gradient.range.quantile), na.rm=TRUE))
    if(diff(zlim)==0) {
      zlim <- as.numeric(range(x))
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    colors <- x;
  }

  labels <- ggplot2::labs(x='Component 1', y='Component 2')
  plot.df <- tibble::rownames_to_column(as.data.frame(embedding), "CellName")
  colnames(plot.df)[2:3] <- c("x", "y")

  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    if (packageVersion("ggrastr") <= "0.1.6") {
      geom_point_w <- function(...)
        ggrastr::geom_point_rast(..., width=raster.width, height=raster.height, dpi=raster.dpi)
    } else {
      geom_point_w <- function(...)
        ggrastr::geom_point_rast(..., raster.width=raster.width, raster.height=raster.height, raster.dpi=raster.dpi)
    }
  } else {
    if (raster) {
      warning("You have to install ggrastr package to be able to use 'raster' parameter")
    }
    geom_point_w <- ggplot2::geom_point
  }

  if (!is.null(groups)) {
    groups <- as.factor(groups)

    plot.df$Group <- factor(NA, levels=levels(groups))
    arr.ids <- match(names(groups), plot.df$CellName)
    plot.df$Group[arr.ids[!is.na(arr.ids)]] <- groups[!is.na(arr.ids)]

    big.clusts <- (plot.df %>% dplyr::group_by(Group) %>% dplyr::summarise(Size=n()) %>%
                     dplyr::filter(Size >= min.cluster.size))$Group %>% as.vector()

    plot.df$Group[!(plot.df$Group %in% big.clusts)] <- NA
    na.plot.df <- plot.df %>% dplyr::filter(is.na(Group))
    plot.df <- plot.df %>% dplyr::filter(!is.na(Group))

    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) +
      geom_point_w(ggplot2::aes(col=Group), alpha=alpha, size=size) +
      labels

    if (mark.groups) {
      labels.data <- plot.df %>% dplyr::group_by(Group) %>%
        dplyr::summarise(x=median(x), y=median(y), Size=n())

      if (length(font.size) == 1) {
        font.size <- c(font.size, font.size)
      }

      gg <- gg + ggrepel::geom_label_repel(
        data=labels.data, ggplot2::aes(label=Group, size=Size), color='black',
        fill=ggplot2::alpha('white', 0.7), label.size = NA,
        label.padding=ggplot2::unit(1, "pt"), seed=42, ...) +
        ggplot2::scale_size_continuous(range=font.size, trans='identity', guide='none')
    }

    if (is.null(legend.title)) {
      legend.title <- "Group"
    }

    if(is.null(palette)) {
      palette <- rainbow
    }

    color.vals <- palette(length(levels(groups)))
    if (shuffle.colors) {
      color.vals <- sample(color.vals)
    }
    gg <- gg + ggplot2::scale_color_manual(name=legend.title, values=color.vals, labels=levels(groups), drop=F) +
      ggplot2::guides(color=ggplot2::guide_legend(override.aes=list(alpha=1.0)))
  } else if (!is.null(colors)) {
    plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
    na.plot.df <- plot.df %>% dplyr::filter(is.na(Color))
    plot.df <- plot.df %>% dplyr::filter(!is.na(Color))

    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y))+labels;
    if(is.character(colors)) {
      gg <- gg + geom_point_w(color=plot.df$Color, alpha=alpha, size=size)
    } else {
      gg <- gg + geom_point_w(ggplot2::aes(col=Color), alpha=alpha, size=size)

      if (length(color.range) == 1) {
        if (color.range == "all") {
          color.range <- range(colors)
        } else if (color.range == "data") {
          color.range <- NULL
        }
      }

      if (!is.null(palette)) {
        gg <- gg + ggplot2::scale_colour_gradientn(colors=palette(100), limits=color.range)
      } else {
        if (prod(range(colors, na.rm=T)) < 0) {
          gg <- gg + ggplot2::scale_color_gradient2(low="#0000ff",mid="#d8d0d0", high="#ff0000", limits=color.range)
        } else {
          gg <- gg + ggplot2::scale_color_gradient(low="#d8d0d0", high="#ff0000", limits=color.range)
        }
      }
    }

    if (!is.null(legend.title)) {
      gg <- gg + ggplot2::guides(color=ggplot2::guide_colorbar(title=legend.title))
    }
  } else {
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) +
      geom_point_w(alpha=alpha, size=size) +
      labels
  }

  if (plot.na && exists("na.plot.df")) {
    gg <- gg + geom_point_w(data=na.plot.df, alpha=alpha, size=size, color='black', shape=4)
  }

  ## Styling
  if (!is.null(plot.theme)) {
    gg <- gg + plot.theme
  }

  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(title)
  }

  if (!is.null(legend.position)) {
    gg <- gg + ggplot2::theme(legend.position=legend.position,
                              legend.justification=legend.position)
  }

  if (!show.legend) {
    gg <- gg + ggplot2::theme(legend.position="none")
  }

  if (!show.ticks) {
    gg <- gg + ggplot2::theme(axis.ticks=ggplot2::element_blank(),
                              axis.text=ggplot2::element_blank())
  }

  if (!show.labels) {
    gg <- gg + ggplot2::theme(axis.title=ggplot2::element_blank())
  }

  return(gg)
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
    groups <- as.factor(conos.obj$clusters[[clustering]]$groups)
  } else if (is.null(groups)) {
    if(is.null(conos.obj)) stop('either groups factor on the cells or a conos object needs to be specified')
    if(is.null(conos.obj$clusters[[1]])) stop('conos object lacks any clustering. run $findCommunities() first')
    groups <- as.factor(conos.obj$clusters[[1]]$groups)
  }
  if(is.null(sample.factor)) {
    sample.factor <- conos.obj$getDatasetPerCell(); # assignment to samples
  }

  xt <- table(sample.factor[match(names(groups),names(sample.factor))],groups)
  xt <- xt[rowSums(xt)>0,]; xt <- xt[,colSums(xt)>0]

  df <- reshape2::melt(xt); colnames(df) <- c("sample","cluster","f");  df$f <- df$f/colSums(xt)[as.character(df$cluster)]
  clp <- ggplot2::ggplot(df, ggplot2::aes(x=factor(cluster, levels=levels(groups)),y=f,fill=sample)) +
    ggplot2::geom_bar(stat='identity') + ggplot2::xlab('cluster') + ggplot2::ylab('fraction of cells') + ggplot2::theme_bw()

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


# a utility function to translate factor into colors
fac2col <- function(x,s=1,v=1,shuffle=FALSE,min.group.size=1,return.details=F,unclassified.cell.color='gray50',level.colors=NULL) {
  nx <- names(x);
  x <- as.factor(x);
  if(min.group.size>1) {
    x <- factor(x,exclude=levels(x)[unlist(tapply(rep(1,length(x)),x,length))<min.group.size])
    x <- droplevels(x)
  }
  if(is.null(level.colors)) {
    col <- rainbow(length(levels(x)),s=s,v=v);
  } else {
    col <- level.colors[1:length(levels(x))];
  }
  names(col) <- levels(x);

  if(shuffle) col <- sample(col);

  y <- col[as.integer(x)]; names(y) <- names(x);
  y[is.na(y)] <- unclassified.cell.color;
  names(y) <- nx;
  if(return.details) {
    return(list(colors=y,palette=col))
  } else {
    return(y);
  }
}
