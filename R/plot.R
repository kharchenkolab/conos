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
##' @inheritParams embedingPlot
##' @param embeddings list of two-column matrices with (x, y) coordinates of the embeddings. Each mutrix must have cell names in rownames.
##' @param ncol number of columns in the panel
##' @param nrow number of rows in the panel
##' @param panel.size vector with two numbers, which specified (width, height) of the panel in inches. Ignored if raster == FALSE.
##' @param adjust.func function to adjust plots before combining them to single panel. Can be used, for example, to provide color pallette of guides of the plots.
##' @return ggplot2 object with the panel of plots
plotEmbeddings <- function(embeddings, groups=NULL, colors=NULL, ncol=NULL, nrow=NULL, raster=FALSE, panel.size=NULL, adjust.func=NULL, ...) {
  if (is.null(panel.size)) {
    panel.size <- dev.size(units="in")
  } else if (length(panel.size) == 1) {
    panel.size <- c(panel.size, panel.size)
  }

  n.plots <- length(embeddings)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n.plots))
  }
  if (is.null(ncol))
    ncol <- ceiling(n.plots / nrow)
  if (is.null(nrow)) {
    nrow <- ceiling(n.plots / ncol)
  }

  plot.list <- lapply(names(embeddings), function(n)
    embedingPlot(embeddings[[n]], groups=groups, colors=colors, raster=raster, title=n,
                 raster.width=panel.size[1] / nrow, raster.height=panel.size[2] / ncol, ...))

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
##' @return ggplot2 object with the panel of plots
plotPagodas <- function(pagoda.samples, groups=NULL, colors=NULL, gene=NULL, embedding.type='tSNE', ...) {
  embeddings <- lapply(pagoda.samples, function(s) s$embeddings$PCA[[embedding.type]])
  if (!is.null(gene)) {
    colors <- Reduce(c, sapply(pagoda.samples, function(d)
      if(gene %in% colnames(d$counts)) d$counts[,gene] else setNames(rep(NA,nrow(d$counts)), rownames(d$counts))))
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
embedingPlot <- function(embedding, groups=NULL, colors=NULL, plot.na=TRUE, min.cluster.size=0, mark.groups=TRUE,
                         show.legend=FALSE, alpha=0.4, size=0.8, title=NULL, plot.theme=NULL,
                         font.size=c(3, 7), show.ticks=FALSE, show.labels=FALSE, legend.position=NULL, legend.title=NULL,
                         raster=FALSE, raster.width=NULL, raster.height=NULL, raster.dpi=300,
                         ...) {
  labels <- ggplot2::labs(x='Component 1', y='Component 2')
  plot.df <- tibble::as_tibble(embedding, rownames="CellName")
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
    plot.df <- plot.df %>% dplyr::mutate(Cluster=groups[CellName])

    plot.df$Cluster <- as.character(plot.df$Cluster)

    big.clusts <- (plot.df %>% dplyr::group_by(Cluster) %>% dplyr::summarise(Size=n()) %>%
                     dplyr::filter(Size >= min.cluster.size))$Cluster %>% as.vector()

    plot.df$Cluster[!(plot.df$Cluster %in% big.clusts)] <- NA
    na.plot.df <- plot.df %>% dplyr::filter(is.na(Cluster))
    plot.df <- plot.df %>% dplyr::filter(!is.na(Cluster))

    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) +
      geom_point_w(ggplot2::aes(col=Cluster), alpha=alpha, size=size) +
      labels

    if (plot.na) {
      gg <- gg + geom_point_w(data=na.plot.df, alpha=alpha, size=size, color='black', shape=4)
    }

    if (mark.groups) {
      labels.data <- plot.df %>% dplyr::group_by(Cluster) %>%
        dplyr::summarise(x=mean(x, tirm=0.4), y=mean(y, trim=0.4), Size=n())

      if (length(font.size) == 1) {
        font.size <- c(font.size, font.size)
      }

      gg <- gg + ggrepel::geom_label_repel(
        data=labels.data, ggplot2::aes(label=Cluster, size=Size), color='black',
        fill=ggplot2::alpha('white', 0.7), label.size = NA,
        label.padding=ggplot2::unit(1, "pt"), seed=42, ...) +
        ggplot2::scale_size_continuous(range=font.size, trans='identity', guide='none')
    }

    if (!is.null(legend.title)) {
      gg <- gg + ggplot2::guides(color=ggplot2::guide_legend(title=legend.title))
    }

  } else if (!is.null(colors)) {
    plot.df <- plot.df %>% dplyr::mutate(Color=colors[CellName])
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) +
      geom_point_w(ggplot2::aes(col=Color), alpha=alpha, size=size) +
      labels

    if (!is.null(legend.title)) {
      gg <- gg + ggplot2::guides(color=ggplot2::guide_colorbar(title=legend.title))
    }

  } else {
    gg <- ggplot2::ggplot(plot.df, ggplot2::aes(x=x, y=y)) +
      geom_point_w(alpha=alpha, size=size) +
      labels
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
