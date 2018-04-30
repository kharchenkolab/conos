##############################
## Helper Functions
##############################

library(DESeq2)

#' Subset p2ens$aggregateMatrices$'cluster:sample' to the requested samples
#' @param p2ens an p2 object ensemble
#' @param sample.type the samples type
#' @param celltype the celltype
getSamples <- function(p2ens, sample.type, celltype) {
    aggr.data <- p2ens$aggregateMatrices$`cluster:sample`
    meta <- p2ens$aggregateMatrixMeta$`cluster:sample`
    samples <- meta[meta$celltype == celltype & meta$sample.type == sample.type,'sample.name']
    aggr.data[samples,,drop=F]
}

topres <- function(x, ...) { r <- DESeq2::results(x, ...); r[order(r$pvalue),] }

#' Run deseq2 and correct for the provided FCs
#' Currently runs comparison on sample.type metadata column which is expected to be in the metadata
#' @param x12 data matrix
#' @param coldata metadata for x12 matrix
#' @param gene.scale.factors per gene fold chages to correct for
#' @return deseq2 object
DESeq2.correctFC <- function(x12, coldata, gene.scale.factors) {
    ## Get size factors
    dds1 <- DESeq2::DESeqDataSetFromMatrix(x12, coldata[colnames(x12),], design=~sample.type)
    dds1 <- DESeq2::estimateSizeFactors(dds1)
    sf <- DESeq2::sizeFactors(dds1)
    ## Remove gene for which we don't have scaling information
    genes.keep <- names(gene.scale.factors)
    ## subset x12 and build dds4
    x12 <- x12[genes.keep,]
    dds4 <- DESeq2::DESeqDataSetFromMatrix(x12, coldata[colnames(x12),], design=~sample.type)
    ## norm factor with lib size
    nf.tmp <- matrix(rep(sf, nrow(x12)),nrow=nrow(x12),byrow=TRUE)
    rownames(nf.tmp) <- rownames(x12)
    colnames(nf.tmp) <- colnames(x12)
    ## gene scaling in linear scale and in order of x12 and nf.tmp
    gene.scale.factors <- 2^ gene.scale.factors[rownames(nf.tmp)]
    ## Build gene-specific scaling matrix
    x <- do.call(cbind,lapply(colData(dds4)$sample.type, function(x) {
        if(x == levels(coldata$sample.type)[1]) {
            rep(1,length(gene.scale.factors))
        } else {
            gene.scale.factors
        }
    }))
    rownames(x) <- rownames(nf.tmp);
    colnames(x) <- colnames(nf.tmp)
    ## Combine two matrices
    nf.tmp <- nf.tmp * x
    ## recommended centering of the corrections
    x2 <- plyr::aaply(nf.tmp, 1, function(x) {x / exp(mean(log(x)))})
    ## Apply normalisation Factors and run deseq2
    normalizationFactors(dds4) <- x2
    dds4 <- DESeq(dds4)
    ## Return
    dds4
}



getCelltypeFCs <- function(ens.p2=NULL, celltype=NULL, sample.type.comparison = NULL,coldata=NULL,verbose=FALSE) {
    ## check params
    if(is.null(ens.p2)) {error('ens.p2 is null')}
    if(is.null(celltype)) { error('celltype is null') }
    if(is.null(sample.type.comparison)) {error('sample type comparison is null')}
    if(is.null(coldata)) {error('coldata is null')}

    ## Get approprate submatrice
    x <- t(getSamples(ens.p2, sample.type.comparison[2], celltype))
    x1 <- t(getSamples(ens.p2, sample.type.comparison[1],  celltype))
    
    ## genes x aggregates
    x12.a <- cbind(x,x1)
    ## Get differential expression  
    dds1 <- DESeq2::DESeqDataSetFromMatrix(x12.a, coldata[colnames(x12.a),], design=~sample.type)
    dds1 <- DESeq2::DESeq(dds1)

    ## extract fold changes
    res <- DESeq2::results(dds1,cooksCutoff = FALSE, independentFiltering = FALSE)
    fcs <- res$log2FoldChange
    names(fcs) <- rownames(res)
    fcs <- fcs[!is.na(fcs)]
    fcs
}

setNames <- function(x) {names(x) <- as.character(x); x}


getPerCellTypeFCs <- function(ens.p2  =NULL, levels =  NULL, sample.type.comparison = NULL, coldata=NULL, n.cores=1, verbose=FALSE) {
    ## Check args
    if (is.null(levels)) { stop('Provided levels is null') }
    if (is.null(ens.p2)) {stop('ens.p2 is null')}
    if (is.null(sample.type.comparison)) {stop('sample type comparison is null')}
    if (is.null(coldata)) {stop('coldata is NULL')}
    ## Calculate the fold changes
    fcs <- parallel::mclapply(setNames(levels), function(celltype) {
        celltype <- as.character(celltype)
        cat('Calculating FC for ', celltype,'\n');
            try({
                getCelltypeFCs(ens.p2 = ens.p2, celltype=celltype, sample.type.comparison=sample.type.comparison,coldata=coldata,verbose=verbose)
            })
    },mc.cores=n.cores)
    fcs
}

na.rm <- function(x) {x[!is.na(x)]}

trimmedMean <- function(x) {
    x <- na.rm(x)
    min.val <- median(x) - 2 * IQR(x)
    max.val <- median(x) + 2 * IQR(x)
    x[x < min.val | x > max.val] <- NA
    x <-na.rm(x)
    mean(x)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}


getCorrectedDE <- function(ens.p2, cell.type=NULL, sample.type.comparison=NULL, coldata=NULL, fc.correction.3=NULL,
                           membrane.gene.names=NULL) {
    ## Check arguments
    if(is.null(ens.p2)) {stop('ens.p2 is null')}
    if(is.null(cell.type)) {stop('cell.type is null')}
    if(is.null(sample.type.comparison)) {stop('sample.type.comparison is null')}
    if(is.null(coldata)) {stop('coldata is null')}
    if(is.null(fc.correction.3)) {stop('fc.correction.3 is null')}
    ## get data
    x12 <- t(rbind(
        getSamples(ens.p2, sample.type.comparison[2], cell.type),
        getSamples(ens.p2, sample.type.comparison[1], cell.type)
    ))
    ## run de
    dds <- DESeq2.correctFC(x12, coldata,  fc.correction.3)
    # prep res
    res <- DESeq2::results(dds)
    allgenes <- rownames(res)
    res <- as.data.frame(res)
    res <- cbind(gene=allgenes,res)
    res$significant <- res$padj < 0.05
    res <- res[order(res$padj),]
    if(!is.null(membrane.gene.names)) {
        res$membrane <- res$gene %in% membrane.gene.names;
    }
    res$Z <- (-1) * qnorm(res$pvalue)*sign(res$log2FoldChange)
    res$Za <- (-1) * qnorm(res$padj)*sign(res$log2FoldChange)
    ## ilev
    tmp1 <- sample.type.comparison
    names(tmp1) <- sample.type.comparison
    ilev = lapply(tmp1, function(type) { ## for conditions, eg. Whole and Tumor
        sk1 <- rownames(subset(coldata, as.character(sample.type) == type & celltype == cell.type))
        sk1 <- rownames(subset(coldata, sample.type == type & celltype == cell.type))
        x <- as.array(x12[,sk1])
        x <- sweep(x,2,apply(x,2, function(p) {sum(p)/1e6} ),FUN='/')
        x <- x[rownames(res),]
    })
    ## return
    return(list(
        res=res,
        genes=allgenes,
        ilev=ilev,
        snames=levels(coldata$sample.type)
    ))
}


#' Get a correction vector by combining fold changes over all the results
#' @param ens.p2 ensembl p2 object
#' @param aggregation.id the aggregation data slot to use from the ens.p2 object
#' @param sample.type.comparison character vector of 'sample.type' to compare e.g. T vs W
#' @param cell.types.exclude cell types to ignore
#' @param scaleByVariance scale the mean fold changes by the Variance
#' @param useTrimmed return results from trimmed mean
#' @param a name numeric vector of consistent fold changes
calcFCcorrection <- function(ens.p2, aggregation.id, sample.type.comparison, cell.types.exclude = c('tumor'),
                             scaleByVariance=TRUE, useTrimmed=FALSE,n.cores=1,coldata=NULL, cell.type.factor=NULL,
                             per.cell.type.fcs=NULL,verbose=FALSE) {
    ## Get FCs for all cells
    if(is.null(per.cell.type.fcs)) {
        fcs <- getPerCellTypeFCs(levels = levels(cell.type.factor), ens.p2=ens.p2,
                                 sample.type.comparison=sample.type.comparison,
                                 coldata=coldata,n.cores=n.cores)
    } else {
        if (verbose) cat('Using pre-calculated per celltype fcs');
        fcs <- per.cell.type.fcs
    }
    ## TODO: some fail, make sure we know why
    fcs <- fcs[!unlist(lapply(fcs, is.error))]
    fcs <- fcs[!names(fcs) %in% cell.types.exclude]
    ## Put the fcs in a matrix: celltype x genes
    all.genes <- unique(unlist(lapply(fcs, function(x) {names(x)})))
    fcs.mat <- do.call(rbind, lapply(fcs, function(x) {x[all.genes]}))
    colnames(fcs.mat) <- all.genes
    ## stat of fc per gene
    fc.mean.var <- data.frame(
        mean = colMeans(fcs.mat, na.rm=TRUE),
        trimmedMean = apply(fcs.mat,2, trimmedMean),
        var = apply(fcs.mat, 2, function(x) {var(x,na.rm=TRUE)})
    )
    ## Scale correction by variance
    fc.mean.var$trimmedMean[is.na(fc.mean.var$trimmedMean)] <- 1e-6 # fill missing trimmed means
    fc.mean.var$var[is.na(fc.mean.var$var)] <- 1e-6 # fill missing var
    fc.mean.var$var.scaled <- as.numeric(range01(scale(fc.mean.var$var))) # standardise var
    fc.mean.var$weight <- 1-fc.mean.var$var.scaled # weight is inverse std var
    fc.mean.var$correction.fc <- fc.mean.var$mean * fc.mean.var$weight # fc to use for correct
    fc.mean.var$correction.fc.trimmed <- fc.mean.var$trimmedMean * fc.mean.var$weight # fc to use for correct with trimming
    ## Return the requested result
    if (scaleByVariance) {
        if (useTrimmed) {
            ## fc correction vector
            fc.correction.3 <- fc.mean.var[,'correction.fc']
            names(fc.correction.3) <- rownames(fc.mean.var)
        } else {
            fc.correction.3 <- fc.mean.var[,'correction.fc.trimmed']
            names(fc.correction.3) <- rownames(fc.mean.var)
        }
    } else {
        if (useTrimmed) {
            fc.correction.3 <- fc.mean.var[,'trimmedMean']
            names(fc.correction.3) <- rownames(fc.mean.var)
        } else {
            fc.correction.3 <- fc.mean.var[,'mean']
            names(fc.correction.3) <- rownames(fc.mean.var)
        }
    }
    fc.correction.3
}


saveComparisonsAsJSON <- function (comps, fileprefix='')
{
    ret <- lapply(namedNames(comps), function(ncc) {
        cc <- comps[[ncc]]
        lapply(namedNames(cc), function(nccc) {
            xe <- cc[[nccc]]
            if (!is.null(xe$res) && nrow(xe$res) > 0) {
                ## Make the filename
                file <- paste0(fileprefix, make.names(ncc), "__", make.names(nccc), ".json")
                xe$res$rowid <- 1:nrow(xe$res)
                xe$res <- data.frame(xe$res)
                xe$ilev <- lapply(xe$ilev, function(x) list(snames=colnames(x),val=x))
                y <- jsonlite::toJSON(xe, file=file)
                write(y, file)
                list(file = file, contrast = ncc, cluster = nccc)
            } else {
                NULL
            }
        })
    })
    invisible(ret)
}


getCorrectedDE.allTypes <-  function(ens.p2, cellfactor, sample.type.comparison,
                                     membrane.gene.names,n.cores=1,correction.method='global',
                                     cell.types.fc.exclude=NULL,verbose=FALSE) {
    ## Check arguments
    if(!correction.method %in% c('global','exclcurrent')) {stop('Unknown correction method: ',correction.method)}

    ## Prepare metadata
    coldata <- subset(ens.p2$aggregateMatrixMeta[['cluster:sample']], sample.type %in% sample.type.comparison)
    rownames(coldata) <- coldata$sample.name
    
    ## Ensure levels are consistend with those provided
    coldata$sample.type <- factor(coldata$sample.type,levels=sample.type.comparison)
    if (correction.method == 'global') {
        if(verbose) cat('Calculating global correction...')
        fc.correction <-  calcFCcorrection(ens.p2=ens.p2,
                                           aggregation.id='cluster:sample',
                                           sample.type.comparison=sample.type.comparison,
                                           cell.types.exclude=cell.types.fc.exclude,
                                           scaleByVariance=TRUE,
                                           useTrimmed=FALSE, n.cores=n.cores,coldata=coldata,
                                           cell.type.factor = cellfactor,verbose=verbose)
        per.cell.type.fcs <- NULL
    } else {
        ## Cache the per cell type corrections
        per.cell.type.fcs <- getPerCellTypeFCs(levels=levels(cellfactor),ens.p2=ens.p2,
                                               sample.type.comparison=sample.type.comparison,
                                               coldata=coldata, n.cores=n.cores,verbose=verbose)
    }

    ## Run differential expression
    res <- parallel::mclapply(setNames(levels(cellfactor)), function(cell.type) {
        cat('Processing ',cell.type,'...')
        try({
            if (correction.method == 'exclcurrent') {
                fc.correction <- calcFCcorrection(ens.p2=ens.p2,
                                                   aggregation.id='cluster:sample',
                                                   sample.type.comparison=sample.type.comparison,
                                                   cell.types.exclude=c(cell.types.fc.exclude,cell.type),
                                                   scaleByVariance=TRUE,
                                                   useTrimmed=FALSE, n.cores=1,coldata=coldata, ## n.cores == 1
                                                   cell.type.factor = cellfactor, per.cell.type.fcs=per.cell.type.fcs)
            }
            getCorrectedDE(ens.p2=ens.p2, cell.type=cell.type,
                           sample.type.comparison=sample.type.comparison,
                           coldata=coldata,
                           fc.correction.3=fc.correction,
                           membrane.gene.names = membrane.gene.names
                           )
        })
    }, mc.cores=n.cores)
    
    ## remove failed
    res[unlist(lapply(res,is.error))] <- NULL
    
    ## return
    res
}

