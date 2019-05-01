validatePerCellTypeParams <- function(con.obj, groups, sample.groups, ref.level, cluster.sep.chr) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("You have to install DESeq2 package to use differential expression")
  }

  if ( class(con.obj) != 'Conos') stop('con.obj must be a conos object')
  if ( is.null(groups) ) stop('groups must be specified');
  if ( is.null(sample.groups) ) stop('sample.groups must be specified')
  if ( class(sample.groups) != 'list' ) stop('sample.groups must be a list');
  if ( length(sample.groups) != 2 ) stop('sample.groups must be of length 2');
  if ( ! all(unlist(lapply(sample.groups, function(x) class(x) == 'character'))) )
    stop('sample.groups must be a list of character vectors');
  if ( ! all(unlist(lapply(sample.groups, function(x) length(x) > 0))) )
    stop('sample.groups entries must be on length greater or equal to 1')
  if ( ! all(unlist(lapply(sample.groups, function(x) {all(x %in% names(con.obj$samples))}))) )
    stop('sample.groups entries must be names of samples in the conos object')
  if ( is.null(ref.level) ) stop('reference level is not defined')
  ## todo: check samplegrousp are named
  if(is.null(names(sample.groups))) stop('sample.groups must be named')
  if(class(groups) != 'factor') stop('groups must be a factor')
  if(any(grepl(cluster.sep.chr, names(con.obj$samples),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any cluster name')
}

validateBetweenCellTypeParams <- function(con.obj, groups, sample.groups, refgroup, altgroup, cluster.sep.chr) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("You have to install DESeq2 package to use differential expression")
  }

  if (class(con.obj) != 'Conos') stop('con.obj must be a conos object')
  if (is.null(groups) ) stop('groups must be specified');
  if (is.null(sample.groups) ) stop('sample.groups must be specified')
  if (class(sample.groups) != 'list' ) stop('sample.groups must be a list');
  #if ( length(sample.groups) != 2 ) stop('sample.groups must be of length 2');
  if (!all(unlist(lapply(sample.groups, function(x) class(x) == 'character'))) )
    stop('sample.groups must be a list of character vectors');
  if (!all(unlist(lapply(sample.groups, function(x) length(x) > 0))) )
    stop('sample.groups entries must be on length greater or equal to 1')
  if (!all(unlist(lapply(sample.groups, function(x) {all(x %in% names(con.obj$samples))}))) )
    stop('sample.groups entries must be names of samples in the conos object')
  if (is.null(refgroup)) stop('reference group is not defined')
  if (is.null(altgroup)) stop('altgroup is not defined')
  ## todo: check samplegrousp are named
  if(is.null(names(sample.groups))) stop('sample.groups must be named')
  if(class(groups) != 'factor') stop('groups must be a factor')
  if(any(grepl(cluster.sep.chr, names(con.obj$samples),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any cluster name')
}

rawMatricesWithCommonGenes <- function(con.obj, sample.groups=NULL) {
  samples <- con.obj$samples
  if (!is.null(sample.groups)) {
    samples <- samples[unlist(sample.groups)]
  }

  ## Generate an aggregated matrix
  raw.mats <- lapply(samples, getRawCountMatrix, transposed=T)
  common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
  return(lapply(raw.mats, function(x) {x[,common.genes]}))
}

collapseCellsByType <- function(cm, groups, min.cell.count) {
  g1 <- groups[intersect(names(groups), rownames(cm))]
  t1 <- as.numeric(table(g1))
  names(t1) <- levels(g1);
  droplevels <- names(t1)[t1 < min.cell.count]
  g1.n <- names(g1)
  g1 <- as.character(g1)
  names(g1) <- g1.n
  g1[g1 %in% droplevels] <- NA
  g1 <- as.factor(g1)
  aggr <- Matrix.utils::aggregate.Matrix(cm, g1)
  aggr <- aggr[rownames(aggr) != "NA",]
  return(aggr)
}

adjustMatrixRownames <- function(name, cm, cluster.sep.chr) {rownames(cm) <- paste0(name, cluster.sep.chr, rownames(cm)); return(cm)}
rbindDEMatrices <- function(mats, cluster.sep.chr) {
  mats <- lapply(names(mats), function(n) {
    rownames(mats[[n]]) <- paste0(n, cluster.sep.chr, rownames(mats[[n]]));
    return(mats[[n]])
  })

  return(t(do.call(rbind, mats)))
}

strpart <- function (x, split, n, fixed = FALSE) {
  sapply(strsplit(as.character(x), split, fixed = fixed), "[", n)
}

is.error <- function (x) {
  inherits(x, c("try-error", "error"))
}

#' Obtain a correction vector for removing the constant effect between the same clusters of two different apps
#' @param con.obj conos object
#' @param groups a vector specifying clusters
#' @param sample.groups the groups of the samples
#' @param cooks.cutoff cooksCutoff distance for DESeq2
#' @param independent.filtering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr separator for cluster and sample name
#' @param return.details logical, if TRUE return internal sturcuters
#' @param de.init if specified reuses existing differential expression results
#' @param exclude.celltypes names of cell types to exclude from the generation of the vecotr
#' @param correction.method 'varianceweighted' or 'mean' specifies way to merge the fold changes from different cell types
#' @export getCorrectionVector
getCorrectionVector <- function(con.obj, groups=NULL, sample.groups=NULL, cooks.cutoff=FALSE, independent.filtering = FALSE,
                                n.cores=1, cluster.sep.chr = '<!!>', return.details=FALSE,de.init=NULL,exclude.celltypes=c(),
                                correction.method='varianceweighted',ref.level=NULL) {
  validatePerCellTypeParams(con.obj, groups, sample.groups, ref.level, cluster.sep.chr)
    ## Main function
    if(is.null(de.init)) {
        de.init <- getPerCellTypeDE(con.obj, groups=groups, sample.groups=sample.groups,
                                    cooks.cutoff=cooks.cutoff, independent.filtering=independent.filtering,
                                    n.cores=n.cores, cluster.sep.chr=cluster.sep.chr,return.details=FALSE,
                                    ref.level = ref.level);
    }
    de.init <- de.init[!names(de.init) %in% exclude.celltypes]
    allfcs <- lapply(de.init, function(x) {
        if(!is.error(x)) {
            fc <- x$log2FoldChange
            names(fc) <- rownames(x)
            fc
        } else {
            NULL
        }
    })
    allfcs <- allfcs[!unlist(lapply(allfcs, is.null))]
    genes <- Reduce(intersect, lapply(allfcs, names))
    ## Matrix of fold changes
    fc.mat <- do.call(rbind, lapply(allfcs, function(x) x[genes]))
    if (correction.method == 'mean') {
        correction <- apply(fc.mat, 2, mean, na.rm=TRUE)
    } else if (correction.method == 'varianceweighted') {
        mu <- apply(fc.mat, 2, mean, na.rm=TRUE)
        var <- apply(fc.mat, 2, var,na.rm=TRUE)
        weight <- 1 - pchisq(q=var,df=nrow(fc.mat)-1)
        correction <- mu * weight
    } else {
        stop(paste0('unknown correction method: ', correction.method))
    }
    correction[is.na(correction)] <- 0
    ## return
    if (!return.details) {
        return(correction)
    }

    return(list(correction.vector=correction, de.init=de.init))
}

#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' @param con.obj conos object
#' @param groups factor specifying cell types
#' @param sample.groups a list of two character vector specifying the app groups to compare
#' @param cooks.cutoff cooksCutoff for DESeq2
#' @param independent.filtering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details return detals
#' @export getPerCellTypeDE
getPerCellTypeDE <- function(con.obj, groups=NULL, sample.groups=NULL, cooks.cutoff = FALSE, ref.level = NULL, min.cell.count = 10,
                             independent.filtering = FALSE, n.cores=1, cluster.sep.chr = '<!!>',return.details=TRUE) {
  validatePerCellTypeParams(con.obj, groups, sample.groups, ref.level, cluster.sep.chr)

  ## Generate a summary dataset collapsing the cells of the same type in each sample
  ## and merging everything in one matrix
  aggr2 <- rawMatricesWithCommonGenes(con.obj, sample.groups) %>%
    lapply(collapseCellsByType, groups=groups, min.cell.count=min.cell.count) %>%
    rbindDEMatrices(cluster.sep.chr=cluster.sep.chr)
  gc()
  ## For every cell type get differential expression results
  de.res <- papply(sn(levels(groups)), function(l) {
    tryCatch({
      ## Get count matrix
      cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
      ## Generate metadata
      meta <- data.frame(
        sample.id= colnames(cm),
        group= as.factor(unlist(lapply(colnames(cm), function(y) {
          y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
          names(sample.groups)[unlist(lapply(sample.groups,function(x) any(x %in% y)))]
        })))
      )
      if (!ref.level %in% levels(meta$group))
        stop('The reference level is absent in this comparison')
      meta$group <- relevel(meta$group, ref=ref.level)
      if (length(unique(as.character(meta$group))) < 2)
        stop('The cluster is not present in both conditions')
      dds1 <- DESeq2::DESeqDataSetFromMatrix(cm,meta,design=~group)
      dds1 <- DESeq2::DESeq(dds1)
      res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering)
      res1 <- as.data.frame(res1)
      res1 <- res1[order(res1$padj,decreasing = FALSE),]
      ##
      if(return.details) {
        list(res=res1, cm=cm, sample.groups=sample.groups)
      } else {
        res1
      }
    }, error=function(err) NA)
  }, n.cores=n.cores)
  de.res
}

#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' applying the specified correction vector
#' @param con.obj conos object
#' @param groups factor specifying cell types
#' @param sample.groups a named list of two character vectors specifying the app groups to compare
#' @param cooks.cutoff cooksCutoff for DESeq2
#' @param independent.filtering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param correction correction vector obtained from getCorrectionVector
#' @export getPerCellTypeDECorrected
getPerCellTypeDECorrected <- function(con.obj, groups=NULL, sample.groups=NULL, cooks.cutoff = FALSE,
                                      independent.filtering = FALSE, n.cores=1,cluster.sep.chr = '<!!>',
                                      correction=NULL, return.details=TRUE, ref.level=NULL) {
  validatePerCellTypeParams(con.obj, groups, sample.groups, ref.level, cluster.sep.chr)

  ## Generate a summary dataset collapsing the cells of the same type in each sample
  ## and merging everything in one matrix
  aggr2 <- rawMatricesWithCommonGenes(con.obj, sample.groups) %>%
    lapply(function(x) Matrix.utils::aggregate.Matrix(x, groups[intersect(names(groups), rownames(x))])) %>%
    rbindDEMatrices(cluster.sep.chr=cluster.sep.chr)
  gc()
    ## For every cell type get differential expression results
    de.res <- parallel::mclapply(sn(levels(groups)), function(l) {
        try({
            ## Get count matrix
            cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
            ## Generate metadata
            meta <- data.frame(
                sample.id= colnames(cm),
                group= as.factor(unlist(lapply(colnames(cm), function(y) {
                    y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
                    names(sample.groups)[unlist(lapply(sample.groups,function(x) any(x %in% y)))]
                })))
            )
            meta$group <- relevel(meta$group, ref=ref.level)
            if (length(unique(as.character(meta$group))) < 2)
                stop('The cluster is not present in both conditions')

            dds1 <- DESeq2::DESeqDataSetFromMatrix(cm,meta,design=~group)
            dds1 <- DESeq2::estimateSizeFactors(dds1)
            sf <- DESeq2::sizeFactors(dds1)
            if(!(all(rownames(cm) %in% names(correction)) & all(names(correction) %in% rownames(cm))))
                stop('incompatible matrices')
            nf.tmp <- matrix(rep(sf, nrow(cm)),nrow=nrow(cm),byrow=TRUE)
            rownames(nf.tmp) <- rownames(cm);
            colnames(nf.tmp) <- colnames(cm)
            gene.scale.factors <- 2^(correction[rownames(nf.tmp)])
            baselevel <- levels(SummarizedExperiment::colData(dds1)$group)[1]
            x <- do.call(cbind, lapply(SummarizedExperiment::colData(dds1)$group, function(x) {
                if (x == baselevel) {
                    rep(1, length(gene.scale.factors))
                } else {
                    gene.scale.factors
                }
            }))
            rownames(x) <- rownames(nf.tmp);
            colnames(x) <- colnames(nf.tmp)
            nf.tmp <- nf.tmp * x
            x2 <- plyr::aaply(nf.tmp, 1, function(x) {x / exp(mean(log(x)))})
            DESeq2::normalizationFactors(dds1) <- x2
            dds1 <- DESeq2::DESeq(dds1)
            res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering)
            res1 <- as.data.frame(res1)
            res1 <- res1[order(res1$padj,decreasing = FALSE),]
            if (return.details) {
                list(res=res1,cm=cm,sample.groups=sample.groups)
            } else {
                res1
            }
        })
    }, mc.cores=n.cores)
    de.res
}

#' Save differential expression as CSV table
#' @param de.results output of differential expression results, corrected or uncorrected
#' @param saveprefix prefix for output file
#' @param data.frame for gene metadata
#' @export saveDEasCSV
saveDEasCSV <- function(de.results=NULL,saveprefix=NULL,gene.metadata=NULL) {
    if(is.null(de.results)) stop('de.results has not been specified')
    if(is.null(saveprefix)) stop('saveprefix has not bee specified')
    ## find errors
    n.error <- sum(unlist(lapply(de.results,is.error)))
    if(n.error > 0)
        cat("Warning: ", n.error, " of ", length(de.results), ' results have returned an error; ignoring...\n')
    de.results <- de.results[!unlist(lapply(de.results,is.error))]
    ##
    x <- lapply(namedNames(de.results), function(ncc) {
        res.celltype <- de.results[[ncc]]
        res.table <- as.data.frame(res.celltype$res)
        ## append gene names
        res.table$gene <- rownames(res.table)
        ## append singificance
        res.table$significant <- res.table$padj < 0.05
        res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0
        ## Append Z scores and rowid
        res.table$Z <- qnorm(1 - (res.table$pval/2))
        res.table$Z[is.na(res.table$Z)] <- 0
        res.table$Za <- qnorm(1 - (res.table$padj/2))
        res.table$Za[is.na(res.table$Za)] <- 0
        res.table$Z <- res.table$Z  * sign(res.table$log2FoldChange)
        res.table$Za <- res.table$Za  * sign(res.table$log2FoldChange)
        if(!is.null(gene.metadata)) {
            ## match order to metadata table
            mo <- match(as.character(gene.metadata$geneid),as.character(res.table$gene))
            ## drop gene id column
            keep.cols <- colnames(gene.metadata)[colnames(gene.metadata) != 'geneid']
            names(keep.cols) <- keep.cols
            res.table <- cbind(res.table, gene.metadata[mo,keep.cols,drop=FALSE])
        }
        file <- paste0(saveprefix,make.names(ncc),'.csv')
        write.table(x=res.table,file=file)
        res.table
    })
    invisible(x)
}

#' Save differential expression results as JSON
#' @param de.results differential expression results
#' @param saveprefix prefix for the differential expression output
#' @param gene.metadata data.frame with gene metadata
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @export saveDEasJSON
saveDEasJSON <- function(de.results = NULL, saveprefix = NULL, gene.metadata = NULL, cluster.sep.chr='<!!>') {
    ## ### DEVEL
    ## de.results <- all.percl.TvsW
    ## saveprefix <- 'json/'
    ## rm(de.results, saveprefix)
    ## ##
    ## Check input
    if(is.null(de.results)) stop('de.results have not been specified')
    if(is.null(saveprefix)) stop('saveprefix has not been specified')
    ## Find de instances that didn't work (usually because cell type is absent from one or more sample types)
    n.error <- sum(unlist(lapply(de.results, is.error)))
    if(n.error > 0)
        cat("Warning: ", n.error,' of ', length(de.results) ,' results have returned an error; ignoring...\n')
    ## get the de results that worked
    de.results <- de.results[!unlist(lapply(de.results, is.error))]
    ## Generate structure and save JSON
    lapply(namedNames(de.results), function(ncc) {
        res.celltype <- de.results[[ncc]]
        ## Get results table as df
        res.table <- as.data.frame(res.celltype$res)
        ## append gene names
        res.table$gene <- rownames(res.table)
        ## append singificance
        res.table$significant <- res.table$padj < 0.05
        res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0
        ## Append Z scores and rowid
        res.table$Z <- qnorm(1 - (res.table$pval/2))
        res.table$Z[is.na(res.table$Z)] <- 0
        res.table$Za <- qnorm(1 - (res.table$padj/2))
        res.table$Za[is.na(res.table$Za)] <- 0
        res.table$Z <- res.table$Z  * sign(res.table$log2FoldChange)
        res.table$Za <- res.table$Za  * sign(res.table$log2FoldChange)
        res.table$rowid <- 1:nrow(res.table)
        if (!is.null(gene.metadata)) {
            ## match order to metadata table
            mo <- match(as.character(gene.metadata$geneid),as.character(res.table$gene))
            ## drop gene id column
            keep.cols <- colnames(gene.metadata)[colnames(gene.metadata) != 'geneid']
            names(keep.cols) <- keep.cols
        }
        res.table <- cbind(res.table, gene.metadata[mo,keep.cols,drop=FALSE])
        ## get names of all the genes
        all.genes <- rownames(res.table)
        ## Get the count matrix
        cm <-res.celltype$cm
        ## remove the cell type suffix
        colnames(cm) <- strpart(colnames(cm),cluster.sep.chr,1,fixed=TRUE)
        ## ilev entry (submatrices of cps)
        ilev <- lapply(res.celltype$sample.groups, function(sg) {
            ## In certain cases columns may be missing,skip
            sg <- sg[sg %in% colnames(cm)]
            ## keep only cols of interest
            cm.tmp <- cm[,sg]
            ## convert to matrix
            cm.tmp <- as.matrix(cm.tmp)
            rownames(cm.tmp) <- rownames(cm)
            ## calculate cpm
            cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
            cpm <- log10(cpm * 1e6 + 1)
            ##
            snames1 <- colnames(cpm)
            ## Put genes in order
            cpm <- cpm[all.genes,]
            colnames(cpm) <- NULL;
            rownames(cpm) <- NULL;
            ## return
            list(snames=snames1, val=as.matrix(cpm))
        })
        ## snames entry (samplenames)
        snames <- names(res.celltype$sample.groups)
        ## convert to json
        tojson <- list(
            res = res.table,
            genes = all.genes,
            ilev = ilev,
            snames = snames
        )
        y <- jsonlite::toJSON(tojson)
        ## File to save to
        file <- paste0(saveprefix,make.names(ncc),'.json')
        ## create the json file
        write(y,file)
        NULL
    })
    invisible(NULL)
}

#' Compare two cell types across the entire panel
#' @param con.obj conos object
#' @param groups factor describing cell grouping
#' @param sample.groups a named list of two character vectors specifying the app groups to compare
#' @param cooks.cutoff cooksCutoff parameter for DESeq2
#' @param refgroup cell type to compare to be used as reference
#' @param altgroup cell type to compare to
#' @param min.cell.count minimum number of cells per celltype/sample combination to keep
#' @param independent.filtering independentFiltering parameter for DESeq2
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details logical, return detailed results
#' @param only.paired only keep samples that that both cell types above the min.cell.count threshold
#' @export getBetweenCellTypeDE
getBetweenCellTypeDE <- function(con.obj, sample.groups =  NULL, groups=NULL, cooks.cutoff = FALSE, refgroup = NULL, altgroup = NULL, min.cell.count = 10,
                                 independent.filtering = FALSE, cluster.sep.chr = '<!!>',return.details=TRUE, only.paired=TRUE) {
  # TODO: do we really need sample.groups here? They are used in the corrected version for some unknown reason.
  validateBetweenCellTypeParams(con.obj, groups, sample.groups, refgroup, altgroup, cluster.sep.chr)
  ## Get the samples from the panel to use in this comparison
  aggr2 <- rawMatricesWithCommonGenes(con.obj, sample.groups) %>%
    lapply(collapseCellsByType, groups=groups, min.cell.count=min.cell.count) %>%
    rbindDEMatrices(cluster.sep.chr=cluster.sep.chr)
  gc()

  aggr2.meta <- generateDEMatrixMetadata(aggr2, refgroup, altgroup, cluster.sep.chr)

  ## Get the samples that have both cell types only
  if (only.paired) {
    complete.obs <- names(which(apply(reshape2::acast(aggr2.meta, library ~ celltype),1,function(x){sum(is.na(x))}) == 0, useNames = TRUE))
    aggr2.meta <- aggr2.meta[aggr2.meta$library %in% complete.obs,]
  }
  ## Select the desired samples only
  aggr2.meta$celltype <- relevel(aggr2.meta$celltype, ref = refgroup)
  aggr2 <- aggr2[,rownames(aggr2.meta)]
  ## Generate DESeq2 comparison
  dds1 <- DESeq2::DESeqDataSetFromMatrix(aggr2, aggr2.meta, design = ~ library + celltype)
  dds1 <- DESeq2::DESeq(dds1)
  res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering)
  res1 <- res1[order(res1$padj,decreasing = FALSE),]
  ## Return
  if(return.details) {
    list(res=res1, cm=aggr2, meta = aggr2.meta, refgroup = refgroup, altgroup = altgroup, sample.groups=sample.groups)
  } else {
    res1
  }
}

generateDEMatrixMetadata <- function(mtx, refgroup, altgroup, cluster.sep.chr) {
  meta <- data.frame(
    row.names = colnames(mtx),
    sample=colnames(mtx),
    library=strpart(colnames(mtx),cluster.sep.chr,1,fixed=T),
    celltype = strpart(colnames(mtx),cluster.sep.chr,2,fixed=T)
  )

  return(subset(meta, celltype %in% c(refgroup, altgroup)))
}

#' Compare two cell types across the entire panel
#' @param con.obj conos object
#' @param sample.groups a named list of two character vectors specifying the app groups to compare
#' @param groups factor describing cell grouping
#' @param cooks.cutoff cooksCutoff parameter for DESeq2
#' @param refgroup cell type to compare to be used as reference
#' @param altgroup cell type to compare to
#' @param min.cell.count minimum number of cells per celltype/sample combination to keep
#' @param independent.filtering independentFiltering parameter for DESeq2
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details logical, return detailed results
#' @param only.paired only keep samples that that both cell types above the min.cell.count threshold
#' @param correction fold change corrections per genes
#' @param ref.level reference level on the basis of which the correction was calculated
#' @export getBetweenCellTypeCorrectedDE
getBetweenCellTypeCorrectedDE <- function(con.obj, sample.groups =  NULL, groups=NULL, cooks.cutoff = FALSE, refgroup = NULL, altgroup = NULL, min.cell.count = 10,
                                          independent.filtering = FALSE, cluster.sep.chr = '<!!>',return.details=TRUE, only.paired=TRUE, correction = NULL, ref.level=NULL) {
  validateBetweenCellTypeParams(con.obj, groups, sample.groups, refgroup, altgroup, cluster.sep.chr)
  ## Get the samples from the panel to use in this comparison
  aggr2 <- rawMatricesWithCommonGenes(con.obj, sample.groups) %>%
    lapply(collapseCellsByType, groups=groups, min.cell.count=min.cell.count) %>%
    rbindDEMatrices(cluster.sep.chr=cluster.sep.chr)
  gc()

  aggr2.meta <- generateDEMatrixMetadata(aggr2, refgroup, altgroup, cluster.sep.chr=cluster.sep.chr)
  ## Get the samples that have both cell types only
  if (only.paired) {
    complete.obs <- reshape2::acast(aggr2.meta, library ~ celltype, value.var = 'celltype', fun.aggregate = length) %>%
      apply(1, function(x) sum(is.na(x)) == 0) %>% which(useNames = TRUE) %>% names()
    aggr2.meta <- aggr2.meta[aggr2.meta$library %in% complete.obs,]
  }
  ## Select the desired samples only
  aggr2.meta$celltype <- relevel(aggr2.meta$celltype, ref = refgroup)
  aggr2 <- aggr2[,rownames(aggr2.meta)]
  tmp1 <- reshape2::melt(sample.groups)
  colnames(tmp1) <- c('sample','group')
  aggr2.meta$group <-  factor(tmp1$group[match(as.character(aggr2.meta$library), tmp1$sample)])
  aggr2.meta$group <- relevel(aggr2.meta$group, ref = ref.level)
  rm(tmp1)
  ## Generate DESeq2 comparison
  dds1 <- DESeq2::DESeqDataSetFromMatrix(aggr2, aggr2.meta, design = ~ celltype)
  ## Apply the correction based on sample type
  dds1 <- DESeq2::estimateSizeFactors(dds1)
  sf <- DESeq2::sizeFactors(dds1)
  if(!(all(rownames(aggr2) %in% names(correction)) & all(names(correction) %in% rownames(aggr2))))
    stop('incompatible matrices')
  nf.tmp <- matrix(rep(sf, nrow(aggr2)),nrow=nrow(aggr2),byrow=TRUE)
  rownames(nf.tmp) <- rownames(aggr2);
  colnames(nf.tmp) <- colnames(aggr2)
  gene.scale.factors <- 2^(correction[rownames(nf.tmp)])
  baselevel <- levels(SummarizedExperiment::colData(dds1)$group)[1]
  x <- do.call(cbind, lapply(SummarizedExperiment::colData(dds1)$group, function(x) {
    if (x == baselevel) {
      rep(1, length(gene.scale.factors))
    } else {
      gene.scale.factors
    }
  }))
  rownames(x) <- rownames(nf.tmp);
  colnames(x) <- colnames(nf.tmp)
  nf.tmp <- nf.tmp * x
  x2 <- plyr::aaply(nf.tmp, 1, function(x) {x / exp(mean(log(x)))})
  DESeq2::normalizationFactors(dds1) <- x2
  ##
  dds1 <- DESeq2::DESeq(dds1)
  res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering)
  res1 <- res1[order(res1$padj,decreasing = FALSE),]
  ## Return
  if(return.details) {
    list(res=res1, cm=aggr2, meta = aggr2.meta, refgroup = refgroup, altgroup = altgroup, sample.groups=sample.groups)
  } else {
    res1
  }
}

## Marker genes

#' Takes data.frames with info about DE genes for single cell type and many samples and
#' returns data.frame with aggregated info for this cell type
aggregateDEMarkersAcrossDatasets <- function(marker.dfs, z.threshold, upregulated.only) {
  z.scores.per.dataset <- lapply(marker.dfs, function(df) setNames(df$Z, rownames(df)))
  gene.union <- lapply(z.scores.per.dataset, names) %>% Reduce(union, .)
  z.scores <- sapply(z.scores.per.dataset, `[`, gene.union) %>% rowMeans(na.rm=T) %>% sort(decreasing=T)
  pvals <- dnorm(z.scores)
  res <- data.frame(Gene=names(z.scores), Z=z.scores, PValue=pvals, PAdj=p.adjust(pvals))

  z.filter <- if (upregulated.only) res$Z else abs(res$Z)
  return(res[z.filter > z.threshold,])
}

getDifferentialGenesP2 <- function(p2.samples, groups, z.threshold=3.0, upregulated.only=F, verbose=T, n.cores=1) {
  lapply.func <- if (verbose) function(...) pbapply::pblapply(..., cl=n.cores) else function(...) papply(..., n.cores=n.cores)

  groups %<>% as.character() %>% setNames(names(groups))

  if (verbose) cat("Estimating marker genes per sample\n")
  markers.per.sample <- lapply.func(p2.samples, function(p2) {
    if (length(intersect(rownames(p2$counts), names(groups))) < 3) {
      list()
    } else {
      p2$getDifferentialGenes(groups=groups, z.threshold=0)
    }
  })

  if (verbose) cat("Aggregating marker genes\n")
  markers.per.type <- unique(groups) %>% setNames(., .) %>%
    lapply(function(id) lapply(markers.per.sample, `[[`, id) %>% .[!sapply(., is.null)]) %>%
    lapply.func(aggregateDEMarkersAcrossDatasets, z.threshold=z.threshold, upregulated.only=upregulated.only)

  if (verbose) cat("All done!\n")

  return(markers.per.type)
}
