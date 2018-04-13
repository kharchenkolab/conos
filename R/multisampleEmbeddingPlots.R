
#' Get common set of genes from apps
#' @param p2objs list of pagoda2 objects
#' @export getCommonGenes
getCommonGenes <- function (p2objs) 
{
    app.genes <- lapply(p2objs, function(o) { colnames(o$counts) })
    Reduce(intersect, app.genes)
}

#' Return TRUE if the parameter is an error object
#' @param x the object to test
#' @export is.error
is.error <- function(x) {
    inherits(x, c("try-error", "error"))
}

#' Plot a set of apps with the distribution plot for a signature
#' @param p2.objs pagoda2 objects
#' @param signature set of genes
#' @param filename optional filename to save the signature in
#' @export plotAllWithSignatureDistribution
plotAllWithSignatureDistribution <- function (p2.objs, signature, filename = NULL, panel.size = 600, mark.cluster.cex = 0.8, verbose = F, breaks = 10, score.method='sum') 
{
    require(Cairo)
    mfrow = getParMfrow(length(p2.objs))
    ## Optionally open file
    if (!is.null(filename)) {
        CairoPNG(file = filename, height = mfrow[[1]] * panel.size, width = mfrow[[2]] * panel.size)
    }
    par(mfrow = mfrow) ##, mar = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.65, 0), cex = 0.85)
    common.genes <- getCommonGenes(p2.objs)
    if (verbose) cat(length(common.genes), " common genes identified\n");
    ## Gene genes in the signature
    genes <- signature[signature %in% common.genes]
    cat(length(genes), "of the ", length(signature), " genes were found\n")
    ## Plot them all
    lapply(names(p2.objs), function(dn) {
        plotOneWithSignatureDistribution(p2.objs[[dn]], genes  = genes, title=dn, breaks=breaks, score.method=score.method);
    })
    ## Close file if opened
    if (!is.null(filename)) {
        dev.off()
    }
    invisible(NULL);
}

#' @export plotOneWithSignatureDistribution
plotOneWithSignatureDistribution <- function(p2obj, genes = NULL, title = "", breaks=10, score.method = 'sum') {
    genes.in.app <- genes %in% colnames(p2obj$counts);
    if(!all(genes.in.app)) {
        warning('Some of the specified genes were not found, keeping found genes only.');
        genes <- genes[genes.in.app];
    }
    scores <- getCellSignatureScores(p2obj, genes, score.method=score.method)
    hist(scores, main= title, breaks=breaks)
    NULL
}

#' given a p2 app and a signature
#' score each cell
#' @param p2obj pagoda2 object
#' @param genes set of genes
#' @param score.method currently only sum supported
#' @export getCellSignatureScores
getCellSignatureScores <- function(p2obj, genes, score.method='sum') {
    if(score.method == 'sum') {
        Matrix::rowSums(p2obj$counts[,genes])
    } else {
        stop('Unsupported score method: ',score.method);
    }
}

#' get number of rows and cols to use when printing a n.items number of items
#' @param n.items number of items
#' @param square force number of columns and rows to be equal
#' @export getParMfrow
getParMfrow <- function(n.items, square = FALSE) {
    n <- ceiling(sqrt(n.items))
    if (square)  {
        c(n,n);
    } else {
        m <- ceiling(n.items/n)
        c(n,m)
    }
}

#' Plot multiple pagoda2 application with a specific genes and common zlim
#' @param p2.objs list of pagoda2 applications
#' @param gene name of genes to plot
#' @param filename if not NULL save to file
#' @param panel.size panel size for saving to file
#' @param mark.cluster.cex cex for marking clusters
#' @return NULL
#' @export plotAllWithGene
plotAllWithGene <- function(p2.objs, gene, filename=NULL,panel.size = 600,mark.cluster.cex=0.8,name.prefix='') {
  require(Cairo)
  n <- ceiling(sqrt(length(p2.objs)))
  if (!is.null(filename)) {
    CairoPNG(file=filename,height=n*panel.size,width=n*panel.size)
  }
  expr.vals <- unlist(lapply(p2.objs, function(o) {
    ret <- NA
    if (gene %in% colnames(o$counts)) {
      ret <- o$counts[,gene]
    }
    ret
  }))
  expr.vals <-expr.vals[!is.na(expr.vals)]
  zlim <- calcZlim(expr.vals)
  par(mfrow=c(n,n), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
  invisible(lapply(names(p2.objs),function(dn) {
    d <- p2.objs[[dn]];
    gene.names <- colnames(d$counts)
    if(!gene %in% gene.names) {
      colors <- rep('grey70', length(gene.names))
    } else {
      colors <- d$counts[,gene]
    }
    ## If no cells present fix
    d$plotEmbedding(type='PCA',embeddingType='tSNE',colors=colors,alpha=0.2,do.par=F,zlim=zlim);
    legend(x='topleft',bty='n',legend=paste0(name.prefix,dn))
  }))
  if(!is.null(filename)) {
    dev.off()
  }
  NULL
}


#' Plot an embedding with multiple gene sets
#' @param app pagoda2 app to plot
#' @param markers character vector of genes to lot
#' @param show.gene.count show how many genes were found
#' @param do.par run par()
#' @return NULL
#' @export p2PlotEmbeddingMultiGeneSets
p2PlotEmbeddingMultiGeneSets <- function(app, markers, show.gene.count=TRUE,do.par=TRUE) {
  if(do.par) {
    n <- ceiling(sqrt(length(markers)))
    par(mfrow=c(n,n))
  }
  lapply(names(markers), function(n) {
    m <- markers[[n]]
    p2PlotEmbeddingMultiGenes(app,m,main=n)
  })
  NULL
}

#' Plot the embedding of a pagoda2 app colors by the trend of expression
#' of multiple genes
#' @param app the pagoda2 app to plot
#' @param genes a list of genes to merge and plot
#' @param type specifies the embedding
#' @param embeddingType specified the embedding
#' @param main title of plot
#' @param show.gene.count show how many genes were found
#' @export p2PlotEmbeddingMultiGenes
p2PlotEmbeddingMultiGenes <- function(app, genes, type='PCA', embeddingType='tSNE', main='',
                                      show.gene.count=TRUE) {
  gns <- intersect(genes, colnames(app$counts))
  if (show.gene.count) {
    main <- paste0(main,' ',length(gns),'/',length(genes))
  }
  if(length(gns) > 1) {
    cns <- app$counts[,gns]
    cns <- apply(cns, 2, function(x) {
      pcntiles <- quantile(x,probs=seq(0,1,0.02))[c(1,99)]
      x[x < pcntiles[1]] <-  pcntiles[1]
      x[x > pcntiles[2]] <-  pcntiles[2]
      x
    })
    colors <- rowSums(scale(cns))
  } else if (length(gns) == 1) {
    colors <- cns[,gns]
  } else {
    colors <- c('grey30')
  }
  app$plotEmbedding(type=type,embeddingType=embeddingType,colors=colors,main=main)
}

#' Plot a series of pagoda2 objects with a panel of genes
#' @export p2PlotAllMultiGenes
p2PlotAllMultiGenes <- function(apps, genes, type='PCA', embeddingType='tSNE',main='',show.gene.count=TRUE) {
  n <- ceiling(sqrt(length(apps)))
  par(mfrow=c(n,n))
  lapply(names(apps), function(n) {
    o <- apps[[n]];
    p2PlotEmbeddingMultiGenes(o,genes,type=type,embeddingType=embeddingType,
                              main=paste0(main,' ',n),show.gene.count=show.gene.count)
  })
  invisible(NULL)
}


#' Plot multiple pagoda2 application with the specified groups
#' @param p2.objs list of  pagoda2 objects
#' @param groups names factor of groups
#' @param filename filename to save to as PNG
#' @param panel.size panel size for saving to disk
#' @param mark.cluster.cex cex for cluster names
#' @export plotAllWithGroups
plotAllWithGroups <- function(p2.objs, groups, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
  require(Cairo)
  n <- ceiling(sqrt(length(p2.objs)))
  if (!is.null(filename)) {
    CairoPNG(file=filename,height=n*panel.size,width=n*panel.size)
  }
  par(mfrow=c(n,n), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
  invisible(lapply(names(p2.objs),function(dn) {
    d <- p2.objs[[dn]];
    g1 <- as.factor(groups)
    colors <- NULL
    ## If no cells present fix
    if (!any(names(g1) %in% rownames(d$counts))) {
      g1 <- NULL
      cell.names <- rownames(d$counts)
      colors <- rep('grey70',length(cell.names))
      names(colors) <- cell.names
    }
    d$plotEmbedding(type='PCA',embeddingType='tSNE',groups=g1,alpha=0.2,
                    min.group.size=0,mark.clusters = TRUE,
                    mark.cluster.cex=mark.cluster.cex,do.par=F,colors=colors);
    legend(x='topleft',bty='n',legend=dn)
  }))
  if(!is.null(filename)) {
    dev.off()
  }
}



#' Plot multiple pagoda2 application with a specific genes and common zlim
#' @param p2.objs list of pagoda2 applications
#' @param gene name of genes to plot
#' @param filename if not NULL save to file
#' @param panel.size panel size for saving to file
#' @param mark.cluster.cex cex for marking clusters
#' @return NULL
#' @export plotAllWithDepth
plotAllWithDepth <- function(p2.objs, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
  require(Cairo)
  n <- ceiling(sqrt(length(p2.objs)))
  if (!is.null(filename)) {
    CairoPNG(file=filename,height=n*panel.size,width=n*panel.size)
  }
  ## Get all depth values
  depthvalues <- unlist(lapply(p2.objs, function(o) {o$depth}))
  zlim <- calcZlim(depthvalues)
  ## Do the plotting
  par(mfrow=c(n,n), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
  invisible(lapply(names(p2.objs),function(dn) {
    d <- p2.objs[[dn]];
    ## If no cells present fix
    d$plotEmbedding(type='PCA',embeddingType='tSNE',colors=d$depth,alpha=0.2,do.par=F, zlim=zlim);
    legend(x='topleft',bty='n',legend=dn)
  }))
  if(!is.null(filename)) {
    dev.off()
  }
}


#' Plot a single pagoda app with a signature
#' @param p2obj a pagoda2 object to plot
#' @param genes the genes in the signature
#' @param title the title to print
#' @export plotOneWithSignature
plotOneWithSignature <- function(p2obj, genes = NULL, title = "") {
    genes.in.app <- genes %in% colnames(p2obj$counts);
    if(!all(genes.in.app)) {
        warning('Some of the specified genes were not found, keeping found genes only.');
        genes <- genes[genes.in.app];
    }
    colors <- Matrix::rowSums(p2obj$counts[,genes])
    p2obj$plotEmbedding(type='PCA',embeddingType='tSNE', colors=colors, alpha=0.2, do.par=F)
    legend(x='topleft',bty='n',legend=title)
    NULL
}


#' Plot all apps with the aggregate of a panel of genes
#' @param p2.objs list of pagoda2 objects
#' @param signature character vector of genes
#' @param filename optional file to save to
#' @param panel.size the size of the panel to save to
#' @param mark.cluster.cex mark.cluster.cex for the plotting function
#' @export plotAllWithSignature
plotAllWithSignature <- function (p2.objs, signature, filename = NULL, panel.size = 600, mark.cluster.cex = 0.8, verbose = F) 
{
    require(Cairo)
    mfrow = getParMfrow(length(p2.objs))
    ## Optionally open file
    if (!is.null(filename)) {
        CairoPNG(file = filename, height = mfrow[[1]] * panel.size, width = mfrow[[2]] * panel.size)
    }
    par(mfrow = mfrow, mar = c(0.5, 0.5, 0.5, 0.5), mgp = c(2, 0.65, 0), cex = 0.85)
    common.genes <- getCommonGenes(p2.objs)
    if (verbose) cat(length(common.genes), " common genes identified\n");
    ## Gene genes in the signature
    genes <- signature[signature %in% common.genes]
    cat(length(genes), "of the ", length(signature), " genes were found\n")
    ## Plot them all
    lapply(names(p2.objs), function(dn) {
        plotOneWithSignature(p2.objs[[dn]], genes  = genes, title=dn);
    })
    ## Close file if opened
    if (!is.null(filename)) {
        dev.off()
    }
    invisible(NULL);
}






#' Plot a panel of multiple gene sets using p2PlotEmbeddingMultiGenes for one app
#' @param app a pagoda2 application
#' @param markers a list of character vectors specifying marker genes to plot
#' @param show.gene.count logical, append the ratio of found genes
#' @param do.par logical, call par() to set the layout
#' @return NULL
#' @export p2PlotEmbeddingMultiGeneSets
p2PlotEmbeddingMultiGeneSets <- function(app, markers, show.gene.count=TRUE,do.par=TRUE,main='') {
  if(do.par) {
    n <- ceiling(sqrt(length(markers)))
    par(mfrow=c(n,n))
  }
  lapply(names(markers), function(n) {
    m <- markers[[n]]
    main <- paste0(main,n)
    p2PlotEmbeddingMultiGenes(app,m,main=main)
  })
  NULL
}


#' Plot a grid of pagoda2 app embedding and lists of markers
#' @description Plot a grid of embeddings where the rows are different pagoda2 apps
#' and the columns are different sets of genes folded into a single color scale
#' @param apps a list of pagoda2 apps
#' @param markers a list of character vectors denoting the list of gnees
#' @param groups an optional factor denoting cell sets to be shown a the right most column
#' @param show.builtin.cl logical show the built in PCA->tSNE embedding (default: FALSE)
#' @return NULL
#' @export p2PlotAllMultiGeneSets
p2PlotAllMultiGeneSets <- function(apps, markers,groups=NULL,show.builtin.cl=FALSE) {
  nrow=length(apps);
  ncol=length(markers);
  if(!is.null(groups))  ncol = ncol + 1;
  if(show.builtin.cl) ncol = ncol + 1;
  par(mfrow=c(nrow,ncol), mar=rep(0.5,4), mgp = c(2,0.65,0), cex = 0.85)
  lapply(names(apps), function(n) {
    o <- apps[[n]]
    p2PlotEmbeddingMultiGeneSets(o, markers,do.par=F,main=paste0(n, ' '))
    if(!is.null(groups)) o$plotEmbedding(type='PCA',embeddingType='tSNE',groups=groups,mark.clusters=T, mark.cluster.cex=0.8,do.par=F)
    if(show.builtin.cl) o$plotEmbedding(type='PCA',embeddingType='tSNE',do.par=F)
  })
  NULL
}


#' Plot marker for apps
#' @export plotMarkerForApps
plotMarkerForApps <- function(r.n, clusters, marker, only.clusters= NULL, hide.outliers=TRUE) {
  ## Optionally look only at specific clusters
  if (!is.null(only.clusters)) {
    clusters <- factor2Char(clusters)
    clusters <- clusters[clusters %in% only.clusters]
    clusters <- as.factor(clusters)
  }
  if(hide.outliers) {
    outlier.alpha=0;
  } else  {
    outlier.alpha=1
  }
  ## Get one big matrix
  common.genes <- Reduce(intersect,lapply(r.n, function(o) {colnames(o$counts)}))
  bigM <- do.call(rbind,lapply(r.n,function(o) {o$counts[,common.genes]}))
  ## subset clusters and matrix to common cells
  common.cells <- intersect(rownames(bigM), names(clusters))
  bigM <- bigM[common.cells,]
  clusters <- clusters[common.cells]
  ## Cell counts in different apps
  cc <- unlist(lapply(r.n, function(o) { dim(o$counts)[1] }))
  ## Names factor with cell designation in each sample/app
  sample <- unlist(lapply(names(cc), function(ccn) {
    r <- rep(ccn,cc[ccn])
    names(r) <- rownames(r.n[[ccn]]$counts)
    r
  }))[common.cells]
  ## Put in one df for plotting
  df.tmp <- data.frame(expr=bigM[,marker],sample,clusters)
  ggplot(df.tmp, aes(y=expr, x=clusters, fill=sample)) + geom_boxplot(outlier.alpha=0)  + theme(axis.text.x = element_text(angle=90, hjust=1)) + ggtitle(marker)
}


#' View two annotations side by side
#' @export viewAnnotationsSideBySide
viewAnnotationsSideBySide <- function(p2.objs, annotation1, annotation2, filename=NULL, panel.size=600, mark.cluster.cex=0.8) {
  require(Cairo)
  n <- length(p2.objs)
  if (!is.null(filename)) {
    CairoPNG(file=filename,height=n*panel.size,width=2*panel.size)
  }
  par(mfrow=c(n,2), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
  invisible(lapply(names(p2.objs),function(dn) {
    d <- p2.objs[[dn]];
    ## If no cells present fix
    d$plotEmbedding(type='PCA',embeddingType='tSNE',groups=annotation1,alpha=0.2,do.par=F,mark.clusters=T,mark.cluster.cex=mark.cluster.cex);
    legend(x='topleft',bty='n',legend=paste0(dn,'_annot1'))
    d$plotEmbedding(type='PCA',embeddingType='tSNE',groups=annotation2,alpha=0.2,do.par=F,mark.clusters=T,mark.cluster.cex=mark.cluster.cex);
    legend(x='topleft',bty='n',legend=paste0(dn,'_annot2'))
  }))
  if(!is.null(filename)) {
    dev.off()
  }
}

#' Highlight a specific cluster
#' @param app list of p2 objects
#' @param groups a factor of groups
#' @param hightlight.group the group to show
#' @param ... parameters for plotAllWithGroups
#' @return NULL
#' @export plotAllHighlightGroup
plotAllHighlightGroup <- function(apps, groups, highlight.group,...) {
    ## Input checks
    if (!is.factor(groups)) stop('groups is not a factor');
    if (!is.character(highlight.group)) stop('highlight.group is not a character');
    if (length(highlight.group) == 0) stop('highlight.group is empty');
    if (length(highlight.group) > 1) {
        highlight.group <- highlight.group[1];
        warning('highlight group is of length greater than 1, using first element only');                                      }
    g <- as.character(groups)
    names(g) <- names(groups)
    g <- ifelse( g == highlight.group, highlight.group, 'other')
    g <- as.factor(g)
    plotAllWithGroups(apps, g, ...)
    invisible(NULL)                                                                                                        }  


#' View clusters in multiple pagoda2 apps one after the other
#' @param apps list of pagoda2 objects
#' @param groups named factor of clusters
#' @param new.window open a new X11() window?
#' @export viewClustersInteractively
viewClustersInteractively <- function(apps,groups, new.window=T) {
  if (new.window) X11()
  lapply(levels(groups), function(g) {
    plotAllHighlightGroup(apps, groups, g);
    invisible(readline(prompt='Press [enter] to continue'))
  })
  invisible(0);
}

#' Plots all the objects after doing a quick cleanup of the annotation
#' @export plotAllWithAnnotationQuick
plotAllWithAnnotationQuick <- function(r.n, wannot, removeEmbSel = FALSE, ...) {
  ## Remove embedding selection
  if (removeEmbSel) {
    wannot$`Embedding.Selection` <- NULL
  }
  ## Force remove multiclassified
  wannot <- removeSelectionOverlaps(wannot)
  wannot <- factorFromP2Selection(wannot)
  plotAllWithGroups(r.n, wannot, ...)
}
