### Functions for working with multiple pagoda2 objects

#' Plot multiple pagoda2 apps with a specific genes
#' @param p2.objs list of pagoda2 object
#' @param gene the gene to plot
#' @param filename the name of the file to save to
#' @param panel size for saving to file
#' @param mark.cluster.cex cex for cluster names
#' @export plotAllWithGene
plotAllWithGene <- function(p2.objs, gene, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
    require(Cairo)
    n <- ceiling(sqrt(length(p2.objs)))
    if (!is.null(filename)) {
        CairoPNG(file=filename,height=n*panel.size,width=n*panel.size)
    }
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
        d$plotEmbedding(type='PCA',embeddingType='tSNE',colors=colors,alpha=0.2,do.par=F);
        legend(x='topleft',bty='n',legend=dn)
    }))
    if(!is.null(filename)) {
        dev.off()
    }
}

#' Save a generated list of differential expression comparison as a tree folder structure
#' @param comparisons comparisons object
#' @param prefix prefix to use on all files
#' @return NULL
#' @export save.comparisons
save.comparisons <- function(comparisons, prefix) {
    not.saved.count <- 0;
    total.count <- 0;
    lapply(names(comparisons), function(comparisonName) {
        clusterPath <- paste0(prefix, comparisonName)
        if(!file.exists(clusterPath)) {
            dir.create(clusterPath)
        }
        lapply(names(comparisons[[comparisonName]]), function(comparisonName2) {
            dir2 <- paste0(clusterPath, '/', comparisonName2)
            if(is.list(comparisons[[comparisonName]][[comparisonName2]])) {
                s1.vs.s2 <- comparisons[[comparisonName]][[comparisonName2]]$s1
                s2.vs.s1 <- comparisons[[comparisonName]][[comparisonName2]]$s2
                write.csv(s1.vs.s2, paste0(dir2, '_forward.csv'))
                write.csv(s2.vs.s1, paste0(dir2, '_reverse.csv'))
            } else {
                not.saved.count <<- not.saved.count + 1;
            }
            total.count <<- total.count + 1;
        })
    })
    cat( not.saved.count, 'of', total.count, 'not saved\n')
    invisible(NULL)
}

#' Generate comparisons with Wilcoxon text
#' @param r.n list of pagoda2 objects
#' @param groups groups to compare
#' @param comparisons comparisons to perform
#' @export runAllcomparisonsWilcox
runAllcomparisonsWilcox <- function(r.n, groups, comparisons) {
    ## For each cluster 
    library('parallel')
    ## Prepare cluster list
    cluster.list <- levels(groups)
    names(cluster.list) <- cluster.list
    cluster.list <- cluster.list
    ## Run all comparisons
    all.de.ks.results <- mclapply(cluster.list,function(cluster.to.compare) {
        cat('Comparing cluster', cluster.to.compare, '\n')
        ## For each comparison
        ncomp <- names(comparisons)
        names(ncomp) <- ncomp
        ret <- lapply(ncomp, function(cmp.n) {
            cmp <- comparisons[[cmp.n]]
            cat('\tRunning comparison ',cmp.n,'\n')
            ret2 <- multisampleDE.Wilcox(r.n = r.n, cl = groups, cluster.to.compare = cluster.to.compare,
                                         samples.1 = cmp$s1, samples.2 = cmp$s2)
            ret2
        })
        ret
    }, mc.cores=20)
    ## Return all comparisons
    all.de.ks.results
}


#' Compare cells in the cluster.to.compare between samples in
#' samples.1 and samples.2
#' @export multisampleDE.Wilcox
multisampleDE.Wilcox <- function(r.n, cl, cluster.to.compare, samples.1, samples.2) {
                                        # requires
    require('Matrix')
    require(stats)
    require('pbapply')
    ## Some checks
    if(is.null(samples.1)) {warning('samples.1 is null'); return("samples.1 is null")};
    if(is.null(samples.2)) {warning('samples.2 is null'); return("samples.2 is null")};
    if(is.null(cluster.to.compare)) {warning('cluster to compare is null'); return("cluster.to.compare is null")};
    ## Get the cells to compare accross all apps
    cl <- cl[!is.na(cl)]
    cells.compare <- names(cl)[cl == cluster.to.compare]
    if (length(cells.compare) < 10) {warning('Not enough cells is selected cluster'); return("not enough cells in cluster");}
    ## Get common genes
    genes <- lapply(r.n, function(o) {
        colnames(o$counts)
    })
    common.genes <- Reduce(intersect, genes)
    names(common.genes) <- common.genes
    if (length(common.genes) == 0) {warning('No common genes found'); return("no common genes");}
    ## count matrix subset
    cms <- lapply(r.n, function(o) {
        (o$misc$rawCounts[,common.genes])
    })
    ## Split into two sets matrices, one for each condition
    cms.1 <- cms[samples.1]
    cms.2 <- cms[samples.2]
    ## Merge the matrices in each condition
    cms.1 <- do.call(rbind, cms.1)
    cms.2 <- do.call(rbind, cms.2)
    ## Get cells to keep for this comparison
    cells.1 <- which(rownames(cms.1) %in% cells.compare)
    cells.2 <- which(rownames(cms.2) %in% cells.compare)
    ## Check we have enough cells
    if (length(cells.1) < 10 || length(cells.2) < 10) {
        warning('Not enough cells in one of the two conditions')
        return("not enough cells")
    }
    ## If there are more thatn 10k cells subset (Matrix package breaks)
    max.cells <- 10000
    if(length(cells.1) > max.cells) {cells.1 <- sample(cells.1,max.cells)}
    if(length(cells.2) > max.cells) {cells.2 <- sample(cells.2,max.cells)}
    ## Keep only cells we want
    cms.1 <- cms.1[cells.1,]
    cms.2 <- cms.2[cells.2,]
    ## Put them together in a matrix
    comparison.matrix <- rbind(cms.1,cms.2)
    ## Make a factor for comparing the two sets
    n1 <- dim(cms.1)[1]
    n2 <- dim(cms.2)[1]
    if ((n1 + n2) < 50) {warning('not enough cells'); return("too few total cells")}
    comparison.f <- factor(c(rep('s1',n1),rep('s2',n2)))
    names(comparison.f) <- c(rownames(cms.1), rownames(cms.2))
    ## Put in a pagoda object and do de
    p2 <- Pagoda2$new(t(comparison.matrix), log.scale=F)
    p2$adjustVariance(plot=F,gam.k=10)
    de.set <- p2$getDifferentialGenes(groups=comparison.f)
    ## Return
    de.set
}

#' Run all comparisons for all clusters
#' @export runAllcomparisonsKS
runAllcomparisonsKS <- function(r.n, groups, comparisons) {
    ## For each cluster 
    library('parallel')
    ## Prepare cluster list
    cluster.list <- levels(groups)
    names(cluster.list) <- cluster.list
    cluster.list <- cluster.list
    ## Run all comparisons
    all.de.ks.results <- mclapply(cluster.list,function(cluster.to.compare) {
        cat('Comparing cluster', cluster.to.compare, '\n')
        ## For each comparison
        ncomp <- names(comparisons)
        names(ncomp) <- ncomp
        ret <- lapply(ncomp, function(cmp.n) {
            cmp <- comparisons[[cmp.n]]
            cat('\tRunning comparison ',cmp.n,'\n')
            ret2 <- multisampleDE.KS(r.n = r.n, cl = groups, cluster.to.compare = cluster.to.compare,
                                     samples.1 = cmp$s1, samples.2 = cmp$s2)
            ret2
        })
        ret
    }, mc.cores=20)
    ## Return all comparisons
    all.de.ks.results
}

#' Compare cells in the cluster.to.compare between samples in
#' samples.1 and samples.2
multisampleDE.KS <- function(r.n, cl, cluster.to.compare, samples.1, samples.2) {
    require('Matrix')
    require(stats)
    require('pbapply')
    ## Some checks
    if(is.null(samples.1)) {warning('samples.1 is null'); return(NULL)};
    if(is.null(samples.2)) {warning('samples.2 is null'); return(NULL)};
    if(is.null(cluster.to.compare)) {warning('cluster to compare is null'); return(NULL)};
    ## Get the cells to compare accross all apps
    cl <- cl[!is.na(cl)]
    cells.compare <- names(cl)[cl == cluster.to.compare]
    if (length(cells.compare) < 10) {warning('Not enough cells is selected cluster'); return(NULL);}
    ## Get common genes
    genes <- lapply(r.n, function(o) {
        colnames(o$counts)
    })
    common.genes <- Reduce(intersect, genes)
    names(common.genes) <- common.genes
    if (length(common.genes) == 0) {warning('No common genes found'); return(NULL);}
    ## count matrix subset
    cms <- lapply(r.n, function(o) {
        o$counts[,common.genes]
    })
    ## Split into two sets matrices, one for each condition
    cms.1 <- cms[samples.1]
    cms.2 <- cms[samples.2]
    ## Merge the matrices in each condition
    cms.1 <- do.call(rbind, cms.1)
    cms.2 <- do.call(rbind, cms.2)
    ## Get cells to keep for this comparison
    cells.1 <- which(rownames(cms.1) %in% cells.compare)
    cells.2 <- which(rownames(cms.2) %in% cells.compare)
    ## Check we have enough cells
    if (length(cells.1) < 10 || length(cells.2) < 10) {
        warning('Not enough cells in one of the two conditions')
        return(NULL)
    }
    ## Keep only cells we want
    cms.1 <- cms.1[cells.1,]
    cms.2 <- cms.2[cells.2,]
    ## Perform per gene test
    res <- pblapply(common.genes, function(cg) {
        vals1 <- cms.1[,cg]
        vals2 <- cms.2[,cg]
        retv3 <- list(ks.stat = -1, ks.pval = -1)
        if (length(vals1) > 10 && length(vals2) > 10) {
            suppressWarnings(ks.test.r <- ks.test(vals1, vals2))
            retv3 <- list(ks.stat = ks.test.r$statistic, ks.pval = ks.test.r$p.value)
        }
    })
    ## Put in a data frame and perform fdr correction
    res2 <- data.frame(matrix(unlist(res), nrow=length(res)), stringsAsFactors=FALSE)
    rownames(res2) <- names(res)
    names(res2) <- c('ks.stat','ks.pval')
    res2$ks.qval <- p.adjust(res2$ks.pval, method='fdr')
    ## Sort by qval
    res2 <- res2[order(res2$ks.qval),]
    res2$significant <- res2$ks.qval < 0.05
    ## Return
    res2
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

#' Get cluster proportion plots
#' @export getClusterProportionsPlots
getClusterProportionsPlots <- function(p2o, cl, main = '', order.levels.numeric=FALSE) {
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
    count.plot <- ggplot(dftmp, aes(x=this.sample.annot, y= Freq)) + geom_bar(stat='identity')  + theme_bw() +
        theme(axis.text.x = element_text(angle = 33, hjust = 1)) + scale_x_discrete(name='Cell type') +
        scale_y_continuous(name='Count') + theme(plot.margin = margin(0,0,0,0,"cm")) + ggtitle(main)
    dftmp$pc <- dftmp$Freq/sum(dftmp$Freq)
    freq.plot <- ggplot(dftmp, aes(x=this.sample.annot, y= pc)) + geom_bar(stat='identity')  + theme_bw() +
        theme(axis.text.x = element_text(angle = 33, hjust = 1)) + scale_x_discrete(name='') +
        scale_y_continuous(name='') + theme(plot.margin = margin(0,0,0,0,"cm")) + ggtitle(main)
    list(count.plot=count.plot, freq.plot=freq.plot) 
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

#' Fix a pagoda2 selection object prefixed
#' @export fixSelectionPrefix
fixSelectionPrefix <- function(p2selection, map) {
    nsel <- lapply(p2selection, function(cs) {
        spid <- strsplit(cs$cells,'_')
        pref <- sapply(spid,'[',1)
        cellbarcode <- sapply(spid,'[',2)
        pref.uniq <- unique(pref)
        if(any(!pref.uniq %in% names(map))) {
            warning(paste0('The following samples could not be matched: ',
                           paste0(pref.uniq[!pref.uniq %in% names(map)], collapse=' ')))
        }
        cs$cells <- paste0(map[pref],'.',cellbarcode)
        cs
    })
    nsel
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
            cell.names <- rownames(r.n[[1]]$counts)
            colors <- rep('grey70',length(cell.names))
            names(colors) <- cell.names
        }
        d$plotEmbedding(type='PCA',embeddingType='tSNE',groups=g1,alpha=0.2,min.group.size=0,mark.clusters = TRUE, mark.cluster.cex=mark.cluster.cex,do.par=F,colors=colors);
        legend(x='topleft',bty='n',legend=dn)
    }))
    if(!is.null(filename)) {
        dev.off()
    }
}

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

#' From a list of pagoda2 application remove any that are NULL
#' @param os list of pagoda2 applications
#' @return list of pagoda2 application filtered for NULLs
#' @export removeNullapps
removeNullapps <- function(os) {
    os[!unlist(lapply(os,FUN=is.null))]
}

#' Perform gsea on a pagoda2 table of results
#' @param tbl pagoda2 table of results
#' @param mc.cores number of CPUs to use
#' @export gsea.on.p2.table
gsea.on.p2.table <- function(tbl, mc.cores =32) {
    require('liger')
    zsrs <- tbl$Z
    names(zsrs) <- rownames(tbl)
    go.l <- lapply(org.Hs.GO2Symbol.list, length)
    go.keep <- go.l > 10 & go.l < 1000
    go.sets <- org.Hs.GO2Symbol.list[go.keep]
    r.gsea <- bulk.gsea(zsrs, set.list=go.sets, mc.cores=mc.cores)
    ## Filter
    r.gsea <- r.gsea[r.gsea$q.val < 0.001,]
    ## Annotate with description
    descriptions <- lapply(rownames(r.gsea), function(x) {
        tryCatch({
            gt <- GOTERM[[x]]
            ret <- NULL
            if (!is.null(gt)) ret <- gt@Term
            ret
        }, warning = function(w) {
            NULL
        }, error = function(e) {
            NULL
        })
    })
    r.gsea$desc <- descriptions
    r.gsea[order(r.gsea$q.val,decreasing=F),]
}

#' Get number of differential genes at specified cutoff
#' @export getDEcountAtCutoff
getDEcountAtCutoff <- function(de.res, z.cutoff=6, M.cutoff=2) {
    de.c <- lapply(cl.merged.clean.renamed_wilcox, function(cl.list) {
        lapply(cl.list, function(x) {
            ret <- 0
            if(!is.null(x$s1 )) {
                s1 <- x$s1
                s1 <- abs(s1$M) > M.cutoff & abs(s1$Z) > z.cutoff
                ret <- dim(x$s1)[1]
            } else {
                0
            }
        })
    })
    de.counts <- do.call("cbind", de.c)
    de.counts
}

#' Get cluster proportion plots version 2
#' @param p2o a pagoda2 application
#' @param cl clusters to summarise to
#' @param main plot title
#' @param order.levels.numeric order levels as it they were numbers
#' @param colors a named list of colors to use
#' @return ggplot2 barplot
#' @export getClusterProportionPlots2
getClusterProportionPlots2 <- function(p2o, cl, main = '', order.levels.numeric=FALSE, colors = NULL) {
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
    freq.plot <- ggplot(dftmp, aes(x=this.sample.annot, y= pc, fill=this.sample.annot)) + geom_bar(stat='identity')  + theme_bw() +
        theme(axis.text.x = element_text(angle = 33, hjust = 1)) + scale_x_discrete(name='') +
        scale_y_continuous(name='') +  ggtitle(main) + 
        scale_fill_manual(values=colors) + guides(fill=FALSE) + theme(plot.title = element_text(size = 10, face = "bold"))
    freq.plot
}

#' Get proportion plots for a list of pagoda2 objects
#' @param r.n list of pagoda2 objects
#' @param cl clusters
#' @param order.levels.numeric order leves as if they are numbers
#' @return a list of ggplot2 objects, plot with do.call(plot_grid, [result])
#' @export getProportionPlots2
getProportionPlots2 <- function(r.n, cl, order.levels.numeric=FALSE) {
    require(cowplot)
    clus.prop.plots <- lapply(names(r.n), function(n) {
        o <- r.n[[n]]
        getClusterProportionPlots2(o, cl, n, order.levels.numeric)
    })
}

#' Get proportion plots for a list of pagoda2 objects
#' @param r.n list of pagoda2 objects
#' @param cl clusters
#' @param order.levels.numeric order leves as if they are numbers
#' @param only.clusters only include specified clusters
#' @export getProportionPlots3
getProportionPlots3 <- function(r.n, cl, order.levels.numeric=FALSE, only.clusters = NULL,ymax=NULL) {
    require(cowplot)
    nms <- names(r.n)
    names(nms) <- nms
    clus.prop.plots <- lapply(nms, function(n) {
        o <- r.n[[n]]
        getClusterProportionPlots3(o, cl, n, order.levels.numeric, only.clusters = only.clusters,ymax=ymax)
    })
}


#' Plot multiple pagoda application with depth
#' @export plotAllWithDepth
plotAllWithDepth <- function(p2.objs, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
    warning('deprecated')
    require(Cairo)
    n <- ceiling(sqrt(length(p2.objs)))
    if (!is.null(filename)) {
        CairoPNG(file=filename,height=n*panel.size,width=n*panel.size)
    }
    par(mfrow=c(n,n), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
    invisible(lapply(names(p2.objs),function(dn) {
        d <- p2.objs[[dn]];
        ## If no cells present fix
        d$plotEmbedding(type='PCA',embeddingType='tSNE',colors=d$depth,alpha=0.2,do.par=F);
        legend(x='topleft',bty='n',legend=dn)
    }))
    if(!is.null(filename)) {
        dev.off()
    }
}

#' Get hierarchy of cell groups across multiple apps
#' @param r.n list of pagoda2 objects
#' @param clusters factor of clusters/cell grousp to summarise to
#' @export multiSampleClusterHierarchy
multiSampleClusterHierarchy <- function(r.n, cls) {
    cms <- lapply(r.n, function(o) { o$counts })
    common.genes <- Reduce(intersect, lapply(r.n, function(o) { colnames(o$counts) }))
    ## genes are rows now
    bcm <- do.call(rbind, lapply(cms, function(o) {(o[,common.genes])}))
    cls <- cls[rownames(bcm)]
    clsums <- pagoda2:::colSumByFac(bcm, cls)
    ## Set the gene names
    colnames(clsums) <- colnames(bcm)
    ## Remove NA sum
    clsums <- clsums[-1,]
    rownames(clsums) <- levels(cls)
    ## Get numbers of cells in each cluster
    cl.counts <- table(cls)[levels(cls)]
    ## Get normalised cluster centers
    clsums.norm <- sweep(clsums, 1, cl.counts, FUN='/')
    ## Get correlation distance dendrogram
    hc.cor.euclidean <- hclust(dist(clsums.norm))
    hc.cor.pearson <- hclust(as.dist(1-cor(t(clsums.norm),method='pearson')))
    hc.cor.spearman <- hclust(as.dist(1-cor(t(clsums.norm),method='spearman')))
    hc <- list(hc.cor.euclidean=hc.cor.euclidean,hc.cor.pearson=hc.cor.pearson,hc.cor.spearman=hc.cor.spearman)
    hc
}

#' Get the genes that the apps have in common
#' @param r.n a list of pagoda2 apps
#' @return a character vector of common genes
getCommonGenes <- function(r.n) {
    Reduce(intersect, lapply(r.n, function(o) { colnames(o$counts) }))
}

#' Plot all apps with the aggregate of a panel of genes
#' @param p2.objs list of pagoda2 objects
#' @param signature character vector of genes
#' @param filename optional file to save to
#' @param panel.size the size of the panel to save to
#' @param mark.cluster.cex mark.cluster.cex for the plotting function
plotAllWithSignature <- function(p2.objs, signature, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
    require(Cairo)
    n <- ceiling(sqrt(length(p2.objs)))
    if (!is.null(filename)) {
        CairoPNG(file=filename,height=n*panel.size,width=n*panel.size)
    }
    par(mfrow=c(n,n), mar = c(0.5,0.5,0.5,0.5), mgp = c(2,0.65,0), cex = 0.85);
    common.genes <- getCommonGenes(combinedApps)
    genes <- signature[signature %in% common.genes];
    cat(length(genes), 'of the ',length(signature),'genes were found\n')
    invisible(lapply(names(p2.objs),function(dn) {
        d <- p2.objs[[dn]];
        colors <- Matrix::rowSums(d$counts[,genes])
        d$plotEmbedding(type='PCA',embeddingType='tSNE',colors=colors,alpha=0.2,do.par=F);
        legend(x='topleft',bty='n',legend=dn)
    }))
    if(!is.null(filename)) {
        dev.off()
    }
}

#' Subset all apps to the specified clusters
#' @param r.n list of pagoda2 apps to subset
#' @param cl factor of clusters
#' @param cl.keep names of factor levels the cells of which to keep
#' @param remove.null.apps logical, remove any apps that don't pass filters
subsetAllappsToClusters <- function(r.n, cl, cl.keep, remove.null.apps = TRUE) {
    warning('deprecated')
    cells.keep <- names(cl)[cl %in% cl.keep]
    ret.apps <- lapply(r.n, function(o) {
        tryCatch({
            ret <- NULL
            rn <- rownames(o$misc$rawCounts)
            app.cells.keep <- rn[rn %in% cells.keep]
            if (length(app.cells.keep) > 100) {
                p2 <- Pagoda2$new(t(o$misc$rawCounts[app.cells.keep,]), n.cores=32)
                p2$adjustVariance(plot=F,gam.k=20)
                p2$calculatePcaReduction(nPcs=100,n.odgenes=1000,maxit=3000)
                p2$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T);
                p2$makeKnnGraph(k=30, type='PCA', center=T, weight.type='none', n.cores=32, distance='cosine')
                p2$getKnnClusters(method = infomap.community, type = 'PCA' ,name = 'infomap')
                ret <- p2
            }
            ret
        }, warning = function(w) {
            NULL
        }, error = function(e) {
            NULL
        })
    })
    if(remove.null.apps) ret.apps <- removeNullapps(ret.apps)
}


#' Get differential expression markers from a pagoda2 differntial expression result
#' @description return the cluster-specific upregulated genes per cluster
#' @param de.res pagoda2 differential expression result
#' @return a list of markers
#' @export getMarkersFromDE
getMarkersFromDE <- function(de.res) {
    lapply(de.res, function(x) {
        (rownames(subset(x,highest==TRUE,Z>0)))
    })
}


#' Plot an embedding of the specified pagoda2 application
#' colored by a set of genes aggregated into a single scale
#' @description This function will merged the expression patters of
#' the specified genes, by scalling them individually and summing them
#' @param app a pagoda2 app
#' @param genes a character vector of genes to display
#' @param type type parameter for plotEmbedding()
#' @param embeddingType embeddingType parameter for plotEmbedding()
#' @param main a title to display
#' @param show.gene.count logical, append the ratio of found genes
#' @return NULL
#' @export p2PlotEmbeddingMultiGenes
p2PlotEmbeddingMultiGenes <- function(app, genes, type='PCA', embeddingType='tSNE', main='', show.gene.count=TRUE) {
    gns <- intersect(genes, colnames(app$counts))
    if (show.gene.count) { main <- paste0(main,' ',length(gns),'/',length(genes)) }
    colors <- rowSums(scale(app$counts[,gns]))
    app$plotEmbedding(type=type,embeddingType=embeddingType,colors=colors,do.par=F)
    legend(x='topleft',bty='n',legend=main)
    NULL
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


#' Subset a list of pagoda2 application down to specified genes
#' @param r.n list of pagoda2 applications
#' @param genes list of genes
#' @param remove.null.apps logical, remove apps that end up with no cells?
#' @return a list of apps
#' @export subsetAllappsToGenes
subsetAllappsToGenes <- function(r.n, genes, remove.null.apps = TRUE) {
    ret.apps <- lapply(r.n, function(o) {
        tryCatch({
            ret <- NULL
            app.genes.keep <- genes[genes %in% colnames(o$misc$rawCounts)]
            if (length(app.genes.keep) > 10) {
                p2 <- Pagoda2$new(t(o$misc$rawCounts[,app.genes.keep]), n.cores=32)
                p2$adjustVariance(plot=F,gam.k=20)
                p2$calculatePcaReduction(nPcs=100,n.odgenes=2000,maxit=3000)
                p2$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T);
                p2$makeKnnGraph(k=30, type='PCA', center=T, weight.type='none', n.cores=32, distance='cosine')
                p2$getKnnClusters(method = infomap.community, type = 'PCA' ,name = 'infomap')
                ret <- p2
            }
            ret
        }, error = function(e) {
            NULL
        })
    })
    if(remove.null.apps) ret.apps <- removeNullapps(ret.apps)
}



#' Calculate zlim from a vector of numeric values
#' @param vs numeric values
#' @return zlim range
calcZlim <- function(vs, gradient.range.quantile = 0.95) {
    zlim <- as.numeric(quantile(vs,p=c(1-gradient.range.quantile,gradient.range.quantile)))
    if(diff(zlim)==0) {
        zlim <- as.numeric(range(vs))
    }
    zlim
}


#' Plot multiple pagoda2 application with a specific genes and common zlim
#' @param p2.objs list of pagoda2 applications
#' @param gene name of genes to plot
#' @param filename if not NULL save to file
#' @param panel.size panel size for saving to file
#' @param mark.cluster.cex cex for marking clusters
#' @return NULL
#' @export plotAllWithDepth2
plotAllWithDepth2 <- function(p2.objs, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
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

#' Plot multiple pagoda2 application with a specific genes and common zlim
#' @param p2.objs list of pagoda2 applications
#' @param gene name of genes to plot
#' @param filename if not NULL save to file
#' @param panel.size panel size for saving to file
#' @param mark.cluster.cex cex for marking clusters
#' @return NULL
#' @export plotAllWithGene2
plotAllWithGene2 <- function(p2.objs, gene, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
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
        legend(x='topleft',bty='n',legend=dn)
    }))
    if(!is.null(filename)) {
        dev.off()
    }
    NULL
}

#' Show head and tail of a table
#' @export headtail
headtail <- function(x,n=6) {
    rbind(head(x,n),tail(x,n))
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


#' Get proportion plots
#' @export getClusterProportionPlots3
getClusterProportionPlots3 <- function(p2o, cl, main = '', order.levels.numeric=FALSE, colors = NULL,
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

#' Convert factor to character preserving names
#' @param f a factor
#' @return a character vector
#' @export factor2Char
factor2Char <- function(f) {
    r <- as.character(f);
    names(r) <- names(f);
    r
}

#' Replace the factors levels of one factor with those of another
#' @param originalFactor originalFactor to replace, not modified
#' @param newFactor specifiying what to replace with
#' @param newPrefix prefix for newFactor names
#' @return a new factor that merged newFactor into originalFactor
#' @export replaceClusterFactors
replaceClusterFactors <- function(originalFactor, newFactor, newPrefix=NULL) {
    if(is.null(newPrefix)) stop('newPrefix is null');
    wf <- factor2Char(originalFactor);
    nwf <- factor2Char(newFactor);
    if(any(!names(nwf) %in% names(wf))) stop('newFactor is not a subset of originalFactor');
    wf[names(nwf)] <- paste0(c(newPrefix), nwf);
    wf <- as.factor(wf);
    wf
}

#' Subset all apps to clusters
#' @export subsetAllappsToClusters3
subsetAllappsToClusters3 <- function(r.n, cl, cl.keep) {
    cells.keep <- names(cl)[cl %in% cl.keep]
    cat(length(cells.keep),'\n')
    lapply(r.n, function(o) {
        tryCatch({
            rn <- rownames(o$misc$rawCounts)
            app.cells.keep <- rn %in% cells.keep
            if (length(app.cells.keep) < 100) { NULL  }
            p2 <- Pagoda2$new(t(o$misc$rawCounts[app.cells.keep,]), n.cores=20)
            p2$adjustVariance(plot=F,gam.k=20)
            p2$calculatePcaReduction(nPcs=100,n.odgenes=1000,maxit=3000)
            p2$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=50,verbose=T);
            p2$makeKnnGraph(k=30, type='PCA', center=T, weight.type='none', n.cores=20, distance='cosine')
            p2$getKnnClusters(method = infomap.community, type = 'PCA' ,name = 'infomap')
            p2
        }, warning = function(w) {
            NULL
        }, error = function(e) {
            NULL
        })
    })
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

#' Get cluster specific markers from a pagoda2 differential expression results
#' @param de.res pagoda2 differential expression result
#' @return list of character vectors of marker genes
#' @export getMarkersFromDE
getMarkersFromDE <- function(de.res) {
    lapply(de.res, function(x) {
        (rownames(subset(x,highest==TRUE,Z>0)))
    })
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
    colors <- rowSums(scale(app$counts[,gns]))
    app$plotEmbedding(type=type,embeddingType=embeddingType,colors=colors,main=main)
}

#' Get raw count matrices
#' @export getCountMatricesRaw
getCountMatricesRaw <- function(r.n, common.genes) {
    ccm.raw <- mclapply(r.n,function(r) {
        om <- as(r$misc$rawCounts,'dgTMatrix')
        om <- om[,colnames(om) %in% common.genes]
        mi <- match(colnames(om),common.genes)
        x <- new("dgTMatrix",i = om@i,j = as.integer(mi[om@j+1]-1),x=om@x,Dim=c(nrow(om),length(common.genes)))
        rownames(x) <- rownames(om); colnames(x) <- common.genes
        as(x,'dgCMatrix')
    },mc.cores=30)
    ccm.raw
}

#' Get common genes between multiple pagoda apps
#' @param r.n list of pagoda2 apps
#' @param cutoff minimum number of apps that must have genes
#' @return character vector of common genes
#' @export getCommonGenesCutoff
getCommonGenesCutoff <- function(r.n, cutoff = 3) {
    gl <- lapply(r.n, function(r) colnames(r$counts))
    all.genes <- unique(unlist(gl))
    gc <- do.call(rbind, lapply(gl, function(x) all.genes %in% x))
    common.genes <- all.genes[apply(gc,2,sum) > cutoff]
    common.genes
}

#' Break down a factor returning the names of the elements
#' in each level as character vectors in a list
#' @param f a factor to breakdown
#' @return a list of factor levels with the names of elemements in them
#' @export factorBreakdown
factorBreakdown <- function(f) {
    if(!is.factor(f)) stop('not a factor!')
    lvls <- levels(f);
    names(lvls) <- lvls;
    lapply(lvls, function(l) {
        r <- names(f)[f == l]
        names(r) <- r
        r
    })
}


