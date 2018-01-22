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


