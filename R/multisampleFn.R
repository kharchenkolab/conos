### Functions for working with multiple pagoda2 objects

#' Hierarchy of clusters across multiple samples
#' @param r.n list of pagoda2 object
#' @clusters factor with the cluster annotations
multiSampleClusterHierarchy <- function(r.n, clusters) {
    cms <- lapply(r.n, function(o) { o$counts })
    common.genes <- Reduce(intersect, lapply(r.n, function(o) { colnames(o$counts) }))
                                        # genes are rows now
    bcm <- do.call(rbind, lapply(cms, function(o) {(o[,common.genes])}))
    clusters <- clusters[rownames(bcm)]
    clsums <- pagoda2:::colSumByFac(bcm, clusters)
    ## Set the gene names
    colnames(clsums) <- colnames(bcm)
    ## Remove NA sum
    clsums <- clsums[-1,]
    rownames(clsums) <- levels(clusters)
    ## Get numbers of cells in each cluster
    cl.counts <- table(clusters)[levels(clusters)]
    ## Get normalised cluster centers
    clsums.norm <- sweep(clsums, 1, cl.counts, FUN='/')
    ## Get correlation distance dendrogram
    hc <- hclust(as.dist(1-cor(t(clsums))))
    hc 
}


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

## Run all comparisons for all clusters
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

## Run all comparisons for all clusters
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
                                        # requires
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


plotProportionPlots <- function(r.n, cl, order.levels.numeric=FALSE) {
    require(cowplot)
    clus.prop.plots <- lapply(names(r.n), function(n) {
        o <- r.n[[n]]
        getClusterProportionsPlots(o, cl, n, order.levels.numeric)
    })
    pls <- lapply(clus.prop.plots, function(x) {x$freq.plot})
    do.call(plot_grid, pls)      
}

plotCountPlots <- function(r.n, cl, order.levels.numeric=FALSE) {
    require(cowplot)
    clus.prop.plots <- lapply(names(r.n), function(n) {
        o <- r.n[[n]]
        getClusterProportionsPlots(o, cl, n, order.levels.numeric)
    })
    pls <- lapply(clus.prop.plots, function(x) {x$count.plot})
    do.call(plot_grid, pls)      
}

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

## Fix a pagoda2 selection object prefixes
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

viewClustersInteractively <- function(apps,groups, new.window=T) {
    if (new.window) X11()
    lapply(levels(groups), function(g) {
        plotAllHighlightGroup(apps, groups, g);
        invisible(readline(prompt='Press [enter] to continue'))
    })
    invisible(0);
}

removeNullapps <- function(os) {
    os[!unlist(lapply(os,FUN=is.null))]
}


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
        scale_y_continuous(name='') + theme(plot.margin = margin(0,0,0,0,"cm")) + ggtitle(main) + 
        scale_fill_manual(values=colors) + guides(fill=FALSE) + theme(plot.title = element_text(size = 10, face = "bold"))
    freq.plot
}

getProportionPlots2 <- function(r.n, cl, order.levels.numeric=FALSE) {
    require(cowplot)
   clus.prop.plots <- lapply(names(r.n), function(n) {
        o <- r.n[[n]]
        getClusterProportionPlots2(o, cl, n, order.levels.numeric)
    })
}

plotAllWithDepth <- function(p2.objs, filename=NULL,panel.size = 600,mark.cluster.cex=0.8) {
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

multiSampleClusterHierarchy <- function(r.n, clusters) {
    cms <- lapply(r.n, function(o) { o$counts })
    common.genes <- Reduce(intersect, lapply(r.n, function(o) { colnames(o$counts) }))
                                        # genes are rows now
    bcm <- do.call(rbind, lapply(cms, function(o) {(o[,common.genes])}))
    clusters <- clusters[rownames(bcm)]
    clsums <- pagoda2:::colSumByFac(bcm, clusters)
    ## Set the gene names
    colnames(clsums) <- colnames(bcm)
    ## Remove NA sum
    clsums <- clsums[-1,]
    rownames(clsums) <- levels(clusters)
    ## Get numbers of cells in each cluster
    cl.counts <- table(clusters)[levels(clusters)]
    ## Get normalised cluster centers
    clsums.norm <- sweep(clsums, 1, cl.counts, FUN='/')
    ## Get correlation distance dendrogram
    hc <- hclust(as.dist(1-cor(t(clsums))))
    hc 
}

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

plotProportionPlots <- function(r.n, cl, order.levels.numeric=FALSE) {
    require(cowplot)
    clus.prop.plots <- lapply(names(r.n), function(n) {
        o <- r.n[[n]]
        getClusterProportionsPlots(o, cl, n, order.levels.numeric)
    })
    pls <- lapply(clus.prop.plots, function(x) {x$freq.plot})
    do.call(plot_grid, pls)      
}

plotCountPlots <- function(r.n, cl, order.levels.numeric=FALSE) {
    require(cowplot)
    clus.prop.plots <- lapply(names(r.n), function(n) {
        o <- r.n[[n]]
        getClusterProportionsPlots(o, cl, n, order.levels.numeric)
    })
    pls <- lapply(clus.prop.plots, function(x) {x$count.plot})
    do.call(plot_grid, pls)      
}

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

subsetAllappsToClusters <- function(r.n, cl, cl.keep, remove.null.apps = TRUE) {
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


    p2PlotAllMultiGeneSets <- function(apps, markers,groups=NULL) {
        nrow=length(apps);
        if(is.null(groups)) {
            ncol=length(markers)
        } else {
            ncol = length(markers) + 1;
        }
        par(mfrow=c(nrow,ncol))
        lapply(apps, function(o) {
            p2PlotEmbeddingMultiGeneSets(o, markers,do.par=F)
            if(!is.null(groups)) o$plotEmbedding(type='PCA',embeddingType='tSNE',groups=groups,mark.clusters=T, mark.cluster.cex=0.8)
        })
    }

#' Get differential expression markers from a pagoda2 differntial expression result
#' @description return the cluster-specific upregulated genes per cluster
#' @param de.res pagoda2 differential expression result
#' @result a list of markers
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
    lapply(apps, function(o) {
        p2PlotEmbeddingMultiGeneSets(o, markers,do.par=F)
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
