#' Get proportion plots
#' @export getClusterProportionPlots
getClusterProportionPlots <- function(p2o, cl, main = '', order.levels.numeric=FALSE, colors = NULL,
                                       only.clusters = NULL,ymax=NULL) {
  if (!is.null(only.clusters)) {
    cl <- factor2Char(cl)
    cl <- cl[cl %in% only.clusters]
    cl <- as.factor(cl)
  }
  require(ggplot2)
  this.app.cells <- rownames(p2o$counts)
  this.sample.annot <- rep('NA', length(this.app.cells))
  names(this.sample.annot) <- this.app.cells
  ## Cells from this app that are annotated
  this.app.cells.annot <- intersect(this.app.cells, names(cl))
  this.sample.annot[this.app.cells.annot] <- as.character(cl[this.app.cells.annot])
  this.sample.annot <- factor(this.sample.annot, levels=levels(cl))
  dftmp <- data.frame(table(this.sample.annot))
  ## sort levels
  if (order.levels.numeric) {
    lvls <- levels(dftmp$this.sample.annot)
    lvls <- lvls[order(as.numeric(as.character(lvls)))]
    dftmp$this.sample.annot <- factor(as.character(dftmp$this.sample.annot), levels=lvls)
  }
  if (is.null(colors)) {
    colors <- rainbow(nlevels(cl), s=1,v=0.8)
    names(colors) <- levels(cl)
  }
  dftmp$pc <- dftmp$Freq/sum(dftmp$Freq)
  if(is.null(ymax)) { limits <- NULL } else { limits <- c(0,ymax) }
  freq.plot <- ggplot(dftmp, aes(x=this.sample.annot, y= pc, fill=this.sample.annot)) +
    geom_bar(stat='identity')  + theme_bw() +
    theme(axis.text.x = element_text(angle = 33, hjust = 1)) + scale_x_discrete(name='') +
    scale_y_continuous(name='', limits=limits) + ggtitle(main) + 
    scale_fill_manual(values=colors) + guides(fill=FALSE) +
    theme(plot.title = element_text(size = 10, face = "bold")) 
  freq.plot
}

#' Get proportion plots for a list of pagoda2 objects
#' @param r.n list of pagoda2 objects
#' @param cl clusters
#' @param order.levels.numeric order leves as if they are numbers
#' @param only.clusters only include specified clusters
#' @export getProportionPlots
getProportionPlots <- function(r.n, cl, order.levels.numeric=FALSE, only.clusters = NULL,ymax=NULL) {
  require(cowplot)
  nms <- names(r.n)
  names(nms) <- nms
  clus.prop.plots <- lapply(nms, function(n) {
    o <- r.n[[n]]
    getClusterProportionPlots3(o, cl, n, order.levels.numeric, only.clusters = only.clusters,ymax=ymax)
  })
}

#' Plot proportion plots
#' @export plotProportionPlots
plotProportionPlots <- function(r.n, cl, order.levels.numeric=FALSE) {
  require(cowplot)
  clus.prop.plots <- lapply(names(r.n), function(n) {
    o <- r.n[[n]]
    getClusterProportionsPlots(o, cl, n, order.levels.numeric)
  })
  pls <- lapply(clus.prop.plots, function(x) {x$freq.plot})
  do.call(plot_grid, pls)      
}

#' Plot count plots
#' @export plotCountPlots
plotCountPlots <- function(r.n, cl, order.levels.numeric=FALSE) {
  require(cowplot)
  clus.prop.plots <- lapply(names(r.n), function(n) {
    o <- r.n[[n]]
    getClusterProportionsPlots(o, cl, n, order.levels.numeric)
  })
  pls <- lapply(clus.prop.plots, function(x) {x$count.plot})
  do.call(plot_grid, pls)      
}

