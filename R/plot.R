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

##' Plot panel of specified embeddings
##'
##' @inheritParams embeddingPlot
##' @param embeddings list of two-column matrices with (x, y) coordinates of the embeddings. Each mutrix must have cell names in rownames.
##' @param ncol number of columns in the panel
##' @param nrow number of rows in the panel
##' @param panel.size vector with two numbers, which specified (width, height) of the panel in inches. Ignored if raster == FALSE.
##' @param adjust.func function to adjust plots before combining them to single panel. Can be used, for example, to provide color pallette of guides of the plots.
##' @return ggplot2 object with the panel of plots
plotEmbeddings <- function(embeddings, groups=NULL, colors=NULL, ncol=NULL, nrow=NULL, raster=FALSE, panel.size=NULL, adjust.func=NULL, title.size=6, ...) {
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

  plot.list <- lapply(names(embeddings), function(n)
    embeddingPlot(embeddings[[n]], groups=groups, colors=colors, raster=raster,
                  raster.width=panel.size[1] / nrow, raster.height=panel.size[2] / ncol, ...) +
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
##' @param pagoda.samples list of pagoda2 objects
##' @param gene gene name. If this parameter is provided, points are colored by expression of this gene.
##' @param embedding.type type of pagoda2 embedding
##' @param gradient.range.quantile Winsorization quantile for the numeric colors and gene gradient
##' @return ggplot2 object with the panel of plots
plotPagodas <- function(pagoda.samples, groups=NULL, colors=NULL, gene=NULL, gradient.range.quantile=1, embedding.type='tSNE', ...) {
  if (!is.null(groups)) {
    groups <- as.factor(groups)
  }

  embeddings <- lapply(pagoda.samples, function(s) s$embeddings$PCA[[embedding.type]])
  no.embedding <- sapply(embeddings, is.null)
  if (all(no.embedding)) {
    stop(paste0("No '", embedding.type, "' embedding presented in the samples"))
  }

  if (any(no.embedding)) {
    warning(paste0(sum(no.embedding), " of your samples doesn't have '", embedding.type, "' embedding"))
    embeddings <- embeddings[!no.embedding]
  }

  if (!is.null(gene)) {
      colors <- stats::setNames(as.numeric(unlist(unname(lapply(pagoda.samples, function(d) {
          if(gene %in% colnames(d$counts)) {
            d$counts[,gene]
          }  else {
            rep(NA,nrow(d$counts))
          }
      })))),unlist(lapply(pagoda.samples,function(d) rownames(d$counts))))
  }
  
  if(is.numeric(colors) && gradient.range.quantile<1) {
    x <- colors;
    zlim <- as.numeric(quantile(x,p=c(1-gradient.range.quantile,gradient.range.quantile),na.rm=TRUE))
    if(diff(zlim)==0) {
      zlim <- as.numeric(range(x))
    }
    x[x<zlim[1]] <- zlim[1]; x[x>zlim[2]] <- zlim[2];
    colors <- x;
  }

  return(plotEmbeddings(embeddings, groups=groups, colors=colors, ...))
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
##' @param raster should layer with the points be rasterized (TRUE/ FALSE)? Setting of this argument to TRUE is useful when you need to export a plot with large number of points
##' @param raster.width width of the plot in inches. Ignored if raster == FALSE.
##' @param raster.height height of the plot in inches. Ignored if raster == FALSE.
##' @param raster.dpi dpi of the rasterized plot. Ignored if raster == FALSE.
##' @return ggplot2 object
##' @export
embeddingPlot <- function(embedding, groups=NULL, colors=NULL, plot.na=TRUE, min.cluster.size=0, mark.groups=TRUE,
                          show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, plot.theme=NULL, palette=NULL, color.range="all",
                          font.size=c(3, 7), show.ticks=FALSE, show.labels=FALSE, legend.position=NULL, legend.title=NULL,
                          raster=FALSE, raster.width=NULL, raster.height=NULL, raster.dpi=300, shuffle.colors=FALSE,
                         ...) {
  labels <- ggplot2::labs(x='Component 1', y='Component 2')
  plot.df <- tibble::rownames_to_column(as.data.frame(embedding), "CellName")
  colnames(plot.df)[2:3] <- c("x", "y")

  if (raster && requireNamespace("ggrastr", quietly = TRUE)) {
    geom_point_w <- function(...)
      ggrastr::geom_point_rast(..., width=raster.width, height=raster.height, dpi=raster.dpi)
  } else {
    if (raster) {
      warning("You have to install ggrastr package to be able to use 'raster' parameter")
    }
    geom_point_w <- ggplot2::geom_point
  }

  if (!is.null(groups)) {
    groups <- as.factor(groups)
    plot.df <- plot.df %>% dplyr::mutate(Group=groups[CellName])

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
        dplyr::summarise(x=mean(x, tirm=0.4), y=mean(y, trim=0.4), Size=n())

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
#' @param conosObjs A conos objects
#' @param type one of 'counts' or 'proportions' to select type of plot
#' @param clustering name of clustering in the current object
#' @return a ggplot object 
plotClusterBarplots <- function(conosObjs, type='counts',clustering=NULL, groups=NULL) {
    ## param checking
    #if(is.null(clustering)) clustering <- 'multi level'
    if(!type %in% c('counts','proportions')) stop('argument type must be either counts or proportions')
    ## main function
    if(!is.null(clustering)) {
        if(clustering %in% names(conosObjs$clusters)) stop('Specified clustering doesn\'t exist')
        groups <- as.factor(conosObjs$clusters[[clustering]]$groups)
    } else {
        if (is.null(groups)) stop('One of clustering or groups needs to be specified')
        groups <- as.factor(groups)
    }
    plot.df <- do.call(rbind,lapply(names(conosObjs$samples), function(n) {
        o <- conosObjs$samples[[n]]
        grps1 <- groups[intersect(names(groups), rownames(o$counts))]
        tbl1 <- data.frame(
            clname=levels(grps1),
            val=as.numeric(table(grps1)),
            sample=c(n),
            stringsAsFactors = FALSE
        )
        if(type=='proportions') {
            tbl1$val <- tbl1$val / sum(tbl1$val)
        }
        tbl1
    }))
    gg <- ggplot(plot.df, aes(x=clname,y=val,fill=clname)) + geom_bar(stat='identity') + facet_wrap(~sample)
    if (type == 'counts') {
        gg <- gg + scale_y_continuous(name='counts')
    } else {
        gg <- gg + scale_y_continuous(name='% of sample')
    }
    gg <- gg + scale_x_discrete(name='cluster')
    gg
}     


#' Generate boxplot per cluster of the proportion of cells in each celltype
#' @param conosObjs conos object
#' @param clustering name of the clustering to use
#' @param apptypes a factor specifying how to group the samples
#' @param return.details if TRUE return a list with the plot and the summary data.frame
plotClusterBoxPlotsByAppType <- function(conosObjs, clustering=NULL, apptypes=NULL, return.details=FALSE) {
    type <- 'proportions'
    ## param checking
    if(is.null(clustering)) clustering <- 'multi level'
    if(is.null(apptypes)) stop('apptypes must be spectified')
    if(!is.factor(apptypes)) stop('apptypes must be a factor')
    if(!type %in% c('counts','proportions')) stop('argument type must be either counts or proportions')
    ## main function
    groups <- as.factor(conosObjs$clusters[[clustering]]$groups)
    plot.df <- do.call(rbind,lapply(names(conosObjs$samples), function(n) {
        o <- conosObjs$samples[[n]]
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
    gg <- ggplot(plot.df, aes(x=apptype,y=val,fill=clname)) + facet <- wrap(~clname) + geom <- boxplot()
    if (type == 'counts') {
        gg <- gg + scale <- y <- continuous(name='counts')
    } else {
        gg <- gg + scale <- y <- continuous(name='% of sample')
    }
    gg <- gg + scale <- x <- discrete(name='cluster')
    ## return
    if(return.details) {
        list(plot=gg,data=plot.df)
    } else {
        gg
    }
}


#' Get markers for global clusters
#' @param conosObjs conos object
#' @param clustering name of the clustering to use
#' @param min.samples.expressing minimum number of samples that must have the genes upregulated in the respective cluster
#' @param min.percent.samples.expression minumum percent of samples that must have the gene upregulated
getGlobalClusterMarkers <- function(conosObjs, clustering='multi level',
                                    min.samples.expressing=0,min.percent.samples.expressing=0){
    ## get the groups from the clusters
    groups <- as.factor(conosObjs$clusters[[clustering]]$groups)
    ## de lists
    delists <- lapply(conosObjs$samples, function(p2) {
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
