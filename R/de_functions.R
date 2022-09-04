#' @importFrom stats qnorm
NULL

#' Get a vector of the names of an object named by the names themselves.
#' This is useful with lapply when passing names of objects as it ensures that the output list
#' is also named
#'
#' @param g an objects on which we can call names()
#' @keywords internal
namedNames <- function(g) {
  n <- names(g)
  names(n) <- n
  n
}


#' @keywords internal
validatePerCellTypeParams <- function(con.obj, groups, sample.groups, ref.level, cluster.sep.chr) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("You have to install DESeq2 package to use differential expression")
  }

  if (!('Conos' %in% class(con.obj))) stop('con.obj must be a conos object')
  if (is.null(groups)) stop('groups must be specified');
  if (is.null(sample.groups)) stop('sample.groups must be specified')
  if (!('list' %in% class(sample.groups))) stop('sample.groups must be a list')
  if (length(sample.groups) != 2) stop('sample.groups must be of length 2')
  if (!all(unlist(lapply(sample.groups, function(x) 'character' %in% class(x)))))
    stop('sample.groups must be a list of character vectors')
  if (!all(sapply(sample.groups, length) > 0))
    stop('sample.groups entries must be on length greater or equal to 1')
  if (!all(unlist(lapply(sample.groups, function(x) {all(x %in% names(con.obj$samples))}))))
    stop('sample.groups entries must be names of samples in the conos object')
  if (is.null(ref.level)) stop('reference level is not defined')
  ## todo: check samplegrousp are named
  if(is.null(names(sample.groups))) stop('sample.groups must be named')
  if(!inherits(groups,'factor')) stop('groups must be a factor')
  if(any(grepl(cluster.sep.chr, names(con.obj$samples),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any cluster name')
}

#' @keywords internal
validateBetweenCellTypeParams <- function(con.obj, groups, sample.groups, refgroup, altgroup, cluster.sep.chr) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("You have to install DESeq2 package to use differential expression")
  }

  if (!inherits(con.obj,'Conos')) stop('con.obj must be a conos object')
  if (is.null(groups) ) stop('groups must be specified');
  if (is.null(sample.groups) ) stop('sample.groups must be specified')
  if (!inherits(sample.groups,'list')) stop('sample.groups must be a list');
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
  if(!inherits(groups,'factor')) stop('groups must be a factor')
  if(any(grepl(cluster.sep.chr, names(con.obj$samples),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any sample name')
  if(any(grepl(cluster.sep.chr,levels(groups),fixed=TRUE)))
    stop('cluster.sep.chr must not be part of any cluster name')
}

#' Get raw matrices with common genes
#'
#' @param con.obj Conos object
#' @param sample.groups list of samples to select from Conos object, con.obj$samples (default=NULL)
#' @return raw matrices subset with common genes
#' @export
rawMatricesWithCommonGenes <- function(con.obj, sample.groups=NULL) {
  samples <- con.obj$samples
  if (!is.null(sample.groups)) {
    samples <- samples[unlist(sample.groups)]
  }

  ## Generate an aggregated matrix
  raw.mats <- lapply(samples, getRawCountMatrix, transposed=TRUE)
  common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
  return(lapply(raw.mats, function(x) {x[,common.genes]}))
}

#' @keywords internal
adjustMatrixRownames <- function(name, cm, cluster.sep.chr) {rownames(cm) <- paste0(name, cluster.sep.chr, rownames(cm)); return(cm)}

#' @keywords internal
rbindDEMatrices <- function(mats, cluster.sep.chr) {
  mats <- lapply(names(mats), function(n) {
    rownames(mats[[n]]) <- paste0(n, cluster.sep.chr, rownames(mats[[n]]));
    return(mats[[n]])
  })

  return(t(do.call(rbind, mats)))
}

#' @keywords internal
strpart <- function (x, split, n, fixed = FALSE) {
  sapply(strsplit(as.character(x), split, fixed = fixed), "[", n)
}

#' @keywords internal
is.error <- function (x) {
  inherits(x, c("try-error", "error"))
}


#' Do differential expression for each cell type in a conos object between the specified subsets of apps
#'
#' @param con.obj conos object
#' @param groups factor specifying cell types (default=NULL)
#' @param sample.groups a list of two character vector specifying the app groups to compare (default=NULL)
#' @param cooks.cutoff boolean cooksCutoff for DESeq2 (default=FALSE)
#' @param ref.level the reference level of the sample.groups against which the comparison should be made (default=NULL). If NULL, will pick the first one.
#' @param min.cell.count integer Minimal number of cells per cluster for a sample to be taken into account in a comparison (default=10)
#' @param remove.na boolean If TRUE, remove NAs from DESeq calculations, which often arise as comparisons not possible (default=TRUE)
#' @param max.cell.count maximal number of cells per cluster per sample to include in a comparison (useful for comparing the number of DE genes between cell types) (default=Inf)
#' @param test which DESeq2 test to use (options: "LRT" or "Wald") (default="LRT")
#' @param independent.filtering boolean independentFiltering for DESeq2 (default=FALSE)
#' @param n.cores numeric Number of cores (default=1)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default='<!!>')
#' @param return.details boolean Whether to return verbose details (default=TRUE)
#' @return A list of differential expression results for every cell type
#' @export getPerCellTypeDE
getPerCellTypeDE <- function(con.obj, groups=NULL, sample.groups=NULL, cooks.cutoff=FALSE, ref.level=NULL, min.cell.count=10,
  remove.na=TRUE, max.cell.count=Inf, test="LRT", independent.filtering=FALSE, n.cores=1, cluster.sep.chr = '<!!>', return.details=TRUE) {
  validatePerCellTypeParams(con.obj, groups, sample.groups, ref.level, cluster.sep.chr)

  ## Generate a summary dataset collapsing the cells of the same type in each sample
  ## and merging everything in one matrix
  aggr2 <- rawMatricesWithCommonGenes(con.obj, sample.groups) %>%
    lapply(collapseCellsByType, groups=groups, min.cell.count=min.cell.count, max.cell.count=max.cell.count) %>%
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
      if (length(unique(as.character(meta$group))) < 2) {
        stop('The cluster is not present in both conditions')
      }
      ## check counts
      checkCountsWholeNumbers(cm)
      dds1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, design=~group)
      if(test=="LRT") {
        dds1 <- DESeq2::DESeq(dds1,test="LRT", reduced = ~ 1)
      } else { # defaults to Wald
        dds1 <- DESeq2::DESeq(dds1)
      }
      res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering)
      res1 <- as.data.frame(res1)
      res1 <- res1[order(res1$padj, decreasing = FALSE),]
      ## remove NA values, which exist in log2FoldChange, lfcSE, stat, pvalue, padj
      if (remove.na){
        res1 <- res1[!is.na(res1$padj), ]
      }
      if(return.details) {
        list(res=res1, cm=cm, sample.groups=sample.groups)
      } else {
        res1
      }
    }, error=function(err) {warning("Error for level ", l, ": ", err$message); return(NA)})
  }, n.cores=n.cores)
  de.res
}


#' Save differential expression as table in *csv format
#'
#' @param de.results output of differential expression results, corrected or uncorrected
#' @param saveprefix character prefix for output file
#' @param gene.metadata gene metadta to include (default=NULL)
#' @export
saveDEasCSV <- function(de.results, saveprefix, gene.metadata=NULL) {
    ## find errors
    n.error <- sum(unlist(lapply(de.results,is.error)))
    if(n.error > 0) {
        message("Warning: ", n.error, " of ", length(de.results), ' results have returned an error; ignoring...\n')
    }

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
#'
#' @param de.results differential expression results (default=NULL)
#' @param saveprefix prefix for the differential expression output (default=NULL)
#' @param gene.metadata data.frame with gene metadata (default=NULL)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default='<!!>')
#' @return JSON with DE results
#' @export
saveDEasJSON <- function(de.results = NULL, saveprefix = NULL, gene.metadata = NULL, cluster.sep.chr='<!!>') {

    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package \"jsonlite\" is needed for this function to work. Please install it.", call. = FALSE)
    }
    ## Check input
    if(is.null(de.results)){
      stop('The argument "de.results" has not been specified')
    }
    if(is.null(saveprefix)){
      stop('The argument "saveprefix" has not been specified')
    }
    ## Find de instances that didn't work (usually because cell type is absent from one or more sample types)
    n.error <- sum(unlist(lapply(de.results, is.error)))
    if(n.error > 0) {
        message("Warning: ", n.error,' of ', length(de.results),' results have returned an error; ignoring...\n')
    }

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
#'
#' @param con.obj conos object
#' @param groups factor describing cell grouping (default=NULL)
#' @param sample.groups a named list of two character vectors specifying the app groups to compare (default=NULL)
#' @param cooks.cutoff boolean cooksCutoff parameter for DESeq2 (default=FALSE)
#' @param refgroup cell type to compare to be used as reference (default=NULL)
#' @param altgroup cell type to compare to be used as ALT against refgroup (default=NULL)
#' @param min.cell.count numeric Minimum number of cells per celltype/sample combination to keep (default=10)
#' @param independent.filtering boolean Whether to use independentFiltering parameter for DESeq2 (default=FALSE)
#' @param cluster.sep.chr character string of length 1 specifying a delimiter to separate cluster and app names (default='<!!>')
#' @param return.details boolean Return detailed results (default=TRUE)
#' @param only.paired boolean Only keep samples that that both cell types above the min.cell.count threshold (default=TRUE)
#' @param remove.na boolean If TRUE, remove NAs from DESeq calculations (default=TRUE)
#' @return Returns either a DESeq2::results() object, or if return.details=TRUE,
#'    returns a list of the DESeq2::results(), the samples from the panel to use in this comparison, refgroups, altgroup, and samplegroups
#' @export
getBetweenCellTypeDE <- function(con.obj, groups=NULL, sample.groups=NULL, cooks.cutoff = FALSE, refgroup = NULL, altgroup = NULL, min.cell.count = 10,
                                 independent.filtering = FALSE, cluster.sep.chr = '<!!>',return.details=TRUE, only.paired=TRUE, remove.na=TRUE) {
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
  ## check counts
  checkCountsWholeNumbers(aggr2)
  ## Generate DESeq2 comparison
  dds1 <- DESeq2::DESeqDataSetFromMatrix(aggr2, aggr2.meta, design = ~ library + celltype)
  dds1 <- DESeq2::DESeq(dds1)
  res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, independentFiltering = independent.filtering)
  res1 <- res1[order(res1$padj,decreasing = FALSE),]
  if (remove.na){
    res1 <- res1[!is.na(res1$padj), ]
  }
  ## Return
  if(return.details) {
    list(res=res1, cm=aggr2, meta = aggr2.meta, refgroup = refgroup, altgroup = altgroup, sample.groups=sample.groups)
  } else {
    res1
  }
}

#' @keywords internal
generateDEMatrixMetadata <- function(mtx, refgroup, altgroup, cluster.sep.chr) {
  meta <- data.frame(
    row.names = colnames(mtx),
    sample=colnames(mtx),
    library=strpart(colnames(mtx), cluster.sep.chr, 1, fixed=TRUE),
    celltype = strpart(colnames(mtx), cluster.sep.chr, 2, fixed=TRUE)
  )

  return(subset(meta, .data$celltype %in% c(refgroup, altgroup)))
}

#' Compare two cell types across the entire panel
#'
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
#' @return Returns either a DESeq2::results() object, or if return.details=TRUE,
#'    returns a list of the DESeq2::results(), the samples from the panel to use in this comparison, refgroups, altgroup, and samplegroups
#' @export
getBetweenCellTypeCorrectedDE <- function(con.obj, sample.groups =  NULL, groups=NULL, cooks.cutoff = FALSE, refgroup = NULL, altgroup = NULL, min.cell.count = 10,
                                          independent.filtering = FALSE, cluster.sep.chr = '<!!>',return.details=TRUE, only.paired=TRUE, correction = NULL, ref.level=NULL) {

  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package \"SummarizedExperiment\" is needed for this function to work. Please install it from Bioconductor.", call. = FALSE)
  }
  if (!requireNamespace("plyr", quietly = TRUE)) {
    stop("Package \"plyr\" is needed for this function to work. Please install it.", call. = FALSE)
  }
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
  ## check counts
  checkCountsWholeNumbers(aggr2)
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

## Takes data.frames with info about DE genes for single cell type and many samples and
## returns data.frame with aggregated info for this cell type
#' @keywords internal
aggregateDEMarkersAcrossDatasets <- function(marker.dfs, z.threshold, upregulated.only) {
  if (length(marker.dfs) == 0){
    return(data.frame())
  }

  z.scores.per.dataset <- lapply(marker.dfs, function(df) setNames(df$Z, rownames(df)))
  m.vals.per.dataset <- lapply(marker.dfs, function(df) setNames(df$M, rownames(df)))
  gene.union <- lapply(z.scores.per.dataset, names) %>% Reduce(union, .)
  z.scores <- sapply(z.scores.per.dataset, `[`, gene.union) %>% rowMeans(na.rm=TRUE)
  m.vals <- sapply(m.vals.per.dataset, `[`, gene.union) %>% rowMeans(na.rm=TRUE)
  ro <- order(z.scores,decreasing=TRUE)
  pvals <- dnorm(z.scores)
  res <- data.frame(Gene=names(z.scores), M=m.vals, Z=z.scores, PValue=pvals, PAdj=p.adjust(pvals))[ro,]

  z.filter <- if (upregulated.only) res$Z else abs(res$Z)
  return(res[z.filter > z.threshold,])
}

#' @keywords internal
getDifferentialGenesP2 <- function(p2.samples, groups, z.threshold=3.0, upregulated.only=FALSE, verbose=TRUE, n.cores=1) {

  groups %<>% as.character() %>% setNames(names(groups))

  if (verbose) message("Estimating marker genes per sample\n")
  markers.per.sample <- sccore::plapply(p2.samples, function(p2) {
    if (length(intersect(rownames(p2$counts), names(groups))) < 3) {
      list()
    } else {
      if (packageVersion("pagoda2") >= "0.1.1") {
        p2$getDifferentialGenes(groups=groups, z.threshold=0, append.specificity.metrics=FALSE, append.auc=FALSE)
      } else {
        p2$getDifferentialGenes(groups=groups, z.threshold=0)
      }
    }
  }, progress=verbose, n.cores=n.cores)

  if (verbose) message("Aggregating marker genes\n")

  ## sccore::plapply, at least v1.0.1
  ## sometimes appends warnings to output, e.g. there's a list element 'value' and 'warning'
  ## 
  ## $warning
  ## <simpleWarning in p2$getDifferentialGenes(groups = groups, z.threshold = 0, append.specificity.metrics = FALSE,     append.auc = FALSE): cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison>
  ##
  ## 
  ## If there's a warning, only take the value
  if (!is.null(markers.per.sample$warning)) {
    markers.per.sample <- markers.per.sample$value
  }

  markers.per.type <- unique(groups) %>% setNames(., .) %>%
    lapply(function(id) lapply(markers.per.sample, `[[`, id) %>% .[!sapply(., is.null)]) 
  markers.per.type = sccore::plapply(markers.per.type, aggregateDEMarkersAcrossDatasets, z.threshold=z.threshold, upregulated.only=upregulated.only, progress=verbose, n.cores=n.cores)
  return(markers.per.type)
}

#' Check that the count data contain only integer counts
#'
#' @param aggregated.samples the count data from aggreaged samples input to DESeq
#' @return if non-integer counts are found, an error is returned
#' @keywords internal
checkCountsWholeNumbers <- function(input.matrix){
  ## check all non-zero values whole numbers
  if (!(all(input.matrix == floor(input.matrix)))){
    stop("There are counts in matrix ", input.matrix, "which are not integers. This leads to DESeq errors. Please check your count matrices.")
  }
}

