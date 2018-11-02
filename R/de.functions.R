#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' @param conObj conos object
#' @param groups factor specifying cell types
#' @param sampleGroups a list of two character vector specifying the app groups to compare
#' @param cookscutoff cookscugoff for DESeq2
#' @param independentFiltering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details return detals
#' @export getPerCellTypeDE
getPerCellTypeDE <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff = FALSE, reflevel = NULL, min.cell.count = 10,
                             independentFiltering = FALSE, n.cores=1,cluster.sep.chr = '+',return.details=TRUE) {
    ## Check arguments
    if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')

    if ( is.null(reflevel) ) stop('reference level is not defined')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any cluster name')
    ## Generate a summary dataset collapsing the cells of the same type in each sample
    ## and merging everything in one matrix
    samples.used <- unlist(sampleGroups)
    ## Generate an aggregated matrix
    raw.mats <- lapply(conObj$samples[samples.used], getRawCountMatrix, transposed=T)
    common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
    raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
    aggr2 <- lapply(raw.mats, function(x) {
        g1 <- groups[intersect(names(groups), rownames(x))]
        t1 <- as.numeric(table(g1))
        names(t1) <- levels(g1);
        droplevels <- names(t1)[t1 < min.cell.count]
        g1.n <- names(g1)
        g1 <- as.character(g1)
        names(g1) <- g1.n
        g1[g1 %in% droplevels] <- NA
        g1 <- as.factor(g1)
        aggr <- Matrix.utils::aggregate.Matrix(x, g1)
        aggr <- aggr[rownames(aggr) != "NA",]
        aggr
    })
    aggr2 <- lapply(names(aggr2), function(n) {
        x <- aggr2[[n]]
        rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
        x
    })
    aggr2 <- t(do.call(rbind, aggr2))
    rm(raw.mats); gc()
    ## For every cell type get differential expression results
    de.res <- parallel::mclapply(nbHelpers::namedLevels(groups), function(l) {
        try({
            ## Get count matrix
            cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
            ## Generate metadata
            meta <- data.frame(
                sample.id= colnames(cm),
                group= as.factor(unlist(lapply(colnames(cm), function(y) {
                    y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
                    names(sampleGroups)[unlist(lapply(sampleGroups,function(x) any(x %in% y)))]
                })))
            )
            if (!reflevel %in% levels(meta$group))
                stop('The reference level is absent in this comparison')
            meta$group <- relevel(meta$group, ref=reflevel)
            if (length(unique(as.character(meta$group))) < 2)
                stop('The cluster is not present in both conditions')
            dds1 <- DESeq2::DESeqDataSetFromMatrix(cm,meta,design=~group)
            dds1 <- DESeq2::DESeq(dds1)
            res1 <- DESeq2::results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
            res1 <- as.data.frame(res1)
            res1 <- res1[order(res1$padj,decreasing = FALSE),]
            ##
            if(return.details) {
                list(res=res1, cm=cm, sampleGroups=sampleGroups)
            } else {
                res1
            }
        })
    }, mc.cores=n.cores)
    de.res
}

#' Obtain a correction vector for removing the constant effect between the same clusters of two different apps
#' @param conObj conos object
#' @param groups a vector specifying clusters
#' @param sampleGroups the groups of the samples
#' @param cooksCutoff cooksCutoff distance for DESeq2
#' @param independentFiltering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr separator for cluster and sample name
#' @param return.details logical, if TRUE return internal sturcuters
#' @param de.init if specified reuses existing differential expression results
#' @param exclude.celltypes names of cell types to exclude from the generation of the vecotr
#' @param correction.method 'varianceweighted' or 'mean' specifies way to merge the fold changes from different cell types
#' @export getCorrectionVector
getCorrectionVector <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff=FALSE, independentFiltering = FALSE,
                                n.cores=1, cluster.sep.chr = '+', return.details=FALSE,de.init=NULL,exclude.celltypes=c(),
                                correction.method='varianceweighted',reflevel=NULL) {
    ## Check arguments
    if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any cluster name')
    if(is.null(reflevel)) stop('reference level is not defined')
    ## Main function
    if(is.null(de.init)) {
        de.init <- getPerCellTypeDE(conObj, groups=groups, sampleGroups=sampleGroups,
                                    cooksCutoff=cooksCutoff, independentFiltering=independentFiltering,
                                    n.cores=n.cores, cluster.sep.chr=cluster.sep.chr,return.details=FALSE,
                                    reflevel = reflevel);
    }
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
    fc.mat <- do.call(rbind, lapply(allfcs, function(x) {x[genes]}))
    fc.mat <- fc.mat[!rownames(fc.mat) %in% exclude.celltypes,]
    if (correction.method == 'mean') {
        correction <- apply(fc.mat, 2, mean, na.rm=TRUE)
    } else if (correction.method == 'varianceweighted') {
        mu <- apply(fc.mat, 2, mean, na.rm=TRUE)
        var <- apply(fc.mat, 2, function(x) {var(x,na.rm=TRUE)})
        weight <- 1 - pchisq(q=var,df=nrow(fc.mat)-1)
        correction <- mu * weight
    } else {
        error(paste0('unknown correction method: ', correction.method))
    }
    correction[is.na(correction)] <- 0
    ## return
    if (!return.details) {
        correction
    } else {
        list(correction.vector=correction,de.init=de.init)
    }
}

#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#' applying the specified correction vector
#' @param conObj conos object
#' @param groups factor specifying cell types
#' @param sampleGroups a named list of two character vectors specifying the app groups to compare
#' @param cookscutoff cookscugoff for DESeq2
#' @param independentFiltering independentFiltering for DESeq2
#' @param n.cores number of cores
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param correction correction vector obtained from getCorrectionVector
#' @export getPerCellTypeDECorrected
getPerCellTypeDECorrected <- function(conObj, groups=NULL, sampleGroups=NULL, cooksCutoff = FALSE,
                             independentFiltering = FALSE, n.cores=1,cluster.sep.chr = '+',
                             correction=NULL, return.details=TRUE, reflevel=NULL) {
    ## Check arguments
    if ( is.null(correction) ) stop("Correction can't by null'")
    if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
    if ( is.null(groups) ) stop('groups must be specified');
    if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
    if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
    if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
    if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
        stop('sampleGroups must be a list of character vectors');
    if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
        stop('sampleGroups entries must be on length greater or equal to 1')
    if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
        stop('sampleGroups entries must be names of samples in the conos object')
    if (is.null(reflevel)) stop('reference level is not defined')
    ## todo: check samplegrousp are named
    if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
    if(class(groups) != 'factor') stop('groups must be a factor')
    if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any sample name')
    if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
        stop('cluster.sep.chr must not be part of any cluster name')
    ## Generate a summary dataset collapsing the cells of the same type in each sample
    ## and merging everything in one matrix
    samples.used <- unlist(sampleGroups)
    ## Generate an aggregated matrix
    raw.mats <- lapply(conObj$samples[samples.used], getRawCountMatrix, transposed=T)
    common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
    raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
    aggr2 <- lapply(raw.mats, function(x) {
        g1 <- groups[intersect(names(groups), rownames(x))]
        aggr <- Matrix.utils::aggregate.Matrix(x, g1)
        aggr
    })
    aggr2 <- lapply(names(aggr2), function(n) {
        x <- aggr2[[n]]
        rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
        x
    })
    aggr2 <- t(do.call(rbind, aggr2))
    rm(raw.mats); gc()
    ## For every cell type get differential expression results
    de.res <- parallel::mclapply(nbHelpers::namedLevels(groups), function(l) {
        try({
            ## Get count matrix
            cm <- aggr2[,strpart(colnames(aggr2),cluster.sep.chr,2,fixed=TRUE) == l]
            ## Generate metadata
            meta <- data.frame(
                sample.id= colnames(cm),
                group= as.factor(unlist(lapply(colnames(cm), function(y) {
                    y <- strpart(y,cluster.sep.chr,1,fixed=TRUE)
                    names(sampleGroups)[unlist(lapply(sampleGroups,function(x) any(x %in% y)))]
                })))
            )
            meta$group <- relevel(meta$group, ref=reflevel)
            if (length(unique(as.character(meta$group))) < 2)
                stop('The cluster is not present in both conditions')
            library(DESeq2)
            dds1 <- DESeq2::DESeqDataSetFromMatrix(cm,meta,design=~group)
            dds1 <- DESeq2::estimateSizeFactors(dds1)
            sf <- DESeq2::sizeFactors(dds1)
            if(!(all(rownames(cm) %in% names(correction)) & all(names(correction) %in% rownames(cm))))
                stop('incompatible matrices')
            nf.tmp <- matrix(rep(sf, nrow(cm)),nrow=nrow(cm),byrow=TRUE)
            rownames(nf.tmp) <- rownames(cm);
            colnames(nf.tmp) <- colnames(cm)
            gene.scale.factors <- 2^(correction[rownames(nf.tmp)])
            baselevel <- levels(colData(dds1)$group)[1]
            x <- do.call(cbind, lapply(colData(dds1)$group, function(x) {
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
            normalizationFactors(dds1) <- x2
            dds1 <- DESeq2::DESeq(dds1)
            res1 <- DESeq2::results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
            res1 <- as.data.frame(res1)
            res1 <- res1[order(res1$padj,decreasing = FALSE),]
            if (return.details) {
                list(res=res1,cm=cm,sampleGroups=sampleGroups)
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
#' @export saveDEasJSON
saveDEasJSON <- function(de.results = NULL, saveprefix = NULL, gene.metadata = NULL) {
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
        ## TODO make '+' a parameter
        colnames(cm) <- strpart(colnames(cm),'+',1,fixed=TRUE)
        ## ilev entry (submatrices of cps)
        ilev <- lapply(res.celltype$sampleGroups, function(sg) {
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
        snames <- names(res.celltype$sampleGroups)
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
#' @param conObj conos object
#' @param groups factor describing cell grouping
#' @param sampleGroups a named list of two character vectors specifying the app groups to compare
#' @param cooksCutoff cooksCutoff parameter for DESeq2
#' @param refgroup cell type to compare to be used as reference
#' @param altgroup cell type to compare to
#' @param min.cell.count minimum number of cells per celltype/sample combination to keep
#' @param independentFiltering independentFiltering parameter for DESeq2
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details logical, return detailed results
#' @param only.paired only keep samples that that both cell types above the min.cell.count threshold
#' @export getBetweenCellTypeDE
getBetweenCellTypeDE <- function(conObj, sampleGroups =  NULL, groups=NULL, cooksCutoff = FALSE, refgroup = NULL, altgroup = NULL, min.cell.count = 10,
                                 independentFiltering = FALSE, cluster.sep.chr = '+',return.details=TRUE, only.paired=TRUE) {
  ## Check arguments
  if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
  if ( is.null(groups) ) stop('groups must be specified');
  if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
  if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
  #if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
  if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
    stop('sampleGroups must be a list of character vectors');
  if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
    stop('sampleGroups entries must be on length greater or equal to 1')
  if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
    stop('sampleGroups entries must be names of samples in the conos object')
  if ( is.null(refgroup) ) stop('reference group is not defined')
  if ( is.null(altgroup) ) stop('reference group is not defined')
  ## todo: check samplegrousp are named
  if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
  if(class(groups) != 'factor') stop('groups must be a factor')
  if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any cluster name')
  ## Get the samples from the panel to use in this comparison
  samples.used <- unlist(sampleGroups)
  ## Generate an aggregated matrix
  raw.mats <- lapply(conObj$samples[samples.used], getRawCountMatrix, transposed=T)
  common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
  raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
  aggr2 <- lapply(raw.mats, function(x) {
    g1 <- groups[intersect(names(groups), rownames(x))]
    t1 <- as.numeric(table(g1))
    names(t1) <- levels(g1);
    droplevels <- names(t1)[t1 < min.cell.count]
    g1.n <- names(g1)
    g1 <- as.character(g1)
    names(g1) <- g1.n
    g1[g1 %in% droplevels] <- NA
    g1 <- as.factor(g1)
    aggr <- Matrix.utils::aggregate.Matrix(x, g1)
    aggr <- aggr[rownames(aggr) != "NA",]
    aggr
  })
  aggr2 <- lapply(names(aggr2), function(n) {
    x <- aggr2[[n]]
    rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
    x
  })
  aggr2 <- t(do.call(rbind, aggr2))
  rm(raw.mats); gc()
  ## generate metadata
  aggr2.meta <- data.frame(
    row.names = colnames(aggr2),
    sample=colnames(aggr2),
    library=strpart(colnames(aggr2),cluster.sep.chr,1,fixed=T),
    celltype = strpart(colnames(aggr2),cluster.sep.chr,2,fixed=T)
  )
  ## summarize and subset the datasets
  aggr2.meta <- subset(aggr2.meta, celltype %in% c(refgroup,altgroup))
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
  res1 <- DESeq2::results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
  res1 <- res1[order(res1$padj,decreasing = FALSE),]
  ## Return
  if(return.details) {
    list(res=res1, cm=aggr2, meta = aggr2.meta, refgroup = refgroup, altgroup = altgroup, sampleGroups=sampleGroups)
  } else {
    res1
  }
}



#' Compare two cell types across the entire panel
#' @param conObj conos object
#' @param sampleGroups a named list of two character vectors specifying the app groups to compare
#' @param groups factor describing cell grouping
#' @param cooksCutoff cooksCutoff parameter for DESeq2
#' @param refgroup cell type to compare to be used as reference
#' @param altgroup cell type to compare to
#' @param min.cell.count minimum number of cells per celltype/sample combination to keep
#' @param independentFiltering independentFiltering parameter for DESeq2
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names
#' @param return.details logical, return detailed results
#' @param only.paired only keep samples that that both cell types above the min.cell.count threshold
#' @param correction fold change corrections per genes
#' @param reflevel reference level on the basis of which the correction was calculated
#' @export getBetweenCellTypeCorrectedDE
getBetweenCellTypeCorrectedDE <- function(conObj, sampleGroups =  NULL, groups=NULL, cooksCutoff = FALSE, refgroup = NULL, altgroup = NULL, min.cell.count = 10,
                                          independentFiltering = FALSE, cluster.sep.chr = '+',return.details=TRUE, only.paired=TRUE, correction = NULL, reflevel=NULL) {
  ## Check arguments
  if ( is.null(correction) ) stop("Correction can't by null'")
  if ( class(conObj) != 'Conos') stop('conObj must be a conos object')
  if ( is.null(groups) ) stop('groups must be specified');
  if ( is.null(sampleGroups) ) stop('sampleGroups must be specified')
  if ( class(sampleGroups) != 'list' ) stop('sampleGroups must be a list');
  #if ( length(sampleGroups) != 2 ) stop('sampleGroups must be of length 2');
  if ( ! all(unlist(lapply(sampleGroups, function(x) class(x) == 'character'))) )
    stop('sampleGroups must be a list of character vectors');
  if ( ! all(unlist(lapply(sampleGroups, function(x) length(x) > 0))) )
    stop('sampleGroups entries must be on length greater or equal to 1')
  if ( ! all(unlist(lapply(sampleGroups, function(x) {all(x %in% names(conObj$samples))}))) )
    stop('sampleGroups entries must be names of samples in the conos object')
  if ( is.null(refgroup) ) stop('reference group is not defined')
  if ( is.null(altgroup) ) stop('reference group is not defined')
  ## todo: check samplegrousp are named
  if(is.null(names(sampleGroups))) stop('sampleGroups must be named')
  if(class(groups) != 'factor') stop('groups must be a factor')
  if(any(grepl(cluster.sep.chr, names(conObj$samples),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any cluster name')
  ## Get the samples from the panel to use in this comparison
  samples.used <- unlist(sampleGroups)
  ## Generate an aggregated matrix
  raw.mats <- lapply(conObj$samples[samples.used], getRawCountMatrix, transposed=T)
  common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
  raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
  aggr2 <- lapply(raw.mats, function(x) {
    g1 <- groups[intersect(names(groups), rownames(x))]
    t1 <- as.numeric(table(g1))
    names(t1) <- levels(g1);
    droplevels <- names(t1)[t1 < min.cell.count]
    g1.n <- names(g1)
    g1 <- as.character(g1)
    names(g1) <- g1.n
    g1[g1 %in% droplevels] <- NA
    g1 <- as.factor(g1)
    aggr <- Matrix.utils::aggregate.Matrix(x, g1)
    aggr <- aggr[rownames(aggr) != "NA",]
    aggr
  })
  aggr2 <- lapply(names(aggr2), function(n) {
    x <- aggr2[[n]]
    rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
    x
  })
  aggr2 <- t(do.call(rbind, aggr2))
  rm(raw.mats); gc()
  ## generate metadata
  aggr2.meta <- data.frame(
    row.names = colnames(aggr2),
    sample=colnames(aggr2),
    library=strpart(colnames(aggr2),cluster.sep.chr,1,fixed=T),
    celltype = strpart(colnames(aggr2),cluster.sep.chr,2,fixed=T)
  )
  ## summarize and subset the datasets
  aggr2.meta <- subset(aggr2.meta, celltype %in% c(refgroup,altgroup))
  ## Get the samples that have both cell types only
  if (only.paired) {
    complete.obs <- names(which(apply(reshape2::acast(aggr2.meta, library ~ celltype, value.var = 'celltype',fun.aggregate = length),1,function(x){sum(is.na(x))}) == 0, useNames = TRUE))
    aggr2.meta <- aggr2.meta[aggr2.meta$library %in% complete.obs,]
  }
  ## Select the desired samples only
  aggr2.meta$celltype <- relevel(aggr2.meta$celltype, ref = refgroup)
  aggr2 <- aggr2[,rownames(aggr2.meta)]
  tmp1 <- reshape2::melt(sampleGroups)
  colnames(tmp1) <- c('sample','group')
  aggr2.meta$group <-  factor(tmp1$group[match(as.character(aggr2.meta$library), tmp1$sample)])
  aggr2.meta$group <- relevel(aggr2.meta$group, ref = reflevel)
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
  baselevel <- levels(colData(dds1)$group)[1]
  x <- do.call(cbind, lapply(colData(dds1)$group, function(x) {
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
  normalizationFactors(dds1) <- x2
  ##
  dds1 <- DESeq2::DESeq(dds1)
  res1 <- DESeq2::results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
  res1 <- res1[order(res1$padj,decreasing = FALSE),]
  ## Return
  if(return.details) {
    list(res=res1, cm=aggr2, meta = aggr2.meta, refgroup = refgroup, altgroup = altgroup, sampleGroups=sampleGroups)
  } else {
    res1
  }
}