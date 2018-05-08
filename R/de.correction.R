## libraries
library(DESeq2)

## Helper funcitons
setNames <- function(x) {names(x) <- as.character(x); x}
na.rm <- function(x) {x[!is.na(x)]}
is.error <- function(x) {inherits(x, c('try-error','error'))}

contain.identical <- function(a,b,check.unique=TRUE) {
    ci <- all(a %in% b) & all(b %in% a)
    un <- !any(duplicated(a)) & !any(duplicated(b)) ## both contain no dups
    if(check.unique)
        ci & un
    else
        ci
}



## Not used

## trimmedMean <- function(x) {
##     x <- na.rm(x)
##     min.val <- median(x) - 2 * IQR(x)
##     max.val <- median(x) + 2 * IQR(x)
##     x[x < min.val | x > max.val] <- NA
##     x <-na.rm(x)
##     mean(x)
## }


## topres <- function(x, ...) { r <- DESeq2::results(x, ...); r[order(r$pvalue),] }

###########

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

#' Run deseq2 and correct for the provided FCs
#' Currently runs comparison on sample.type metadata column which is expected to be in the metadata
#' @param x12 data matrix
#' @param coldata metadata for x12 matrix
#' @param gene.scale.factors per gene fold chages to correct for
#' @param correction.global.weight global weighting on the correction, setting to 0 results in no correction
#' @return deseq2 object
DESeq2.correctFC <- function(x12, coldata, gene.scale.factors, correction.global.weight=1) {
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
    gene.scale.factors <- 2^ (gene.scale.factors[rownames(nf.tmp)] * correction.global.weight)
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

#' get fold changes for a specific cell type comparison
#' @param pagoda 2 ensemble object
#' @param celltype to use
#' @param sample.type.comparisons comparison to perform
#' @param coldata metadata
#' @param verbose verbosity logical
#' @param fc.method correction method to use 'deseq2' or 'dummy' (no correction)
getCelltypeFCs <- function(ens.p2=NULL, celltype=NULL, sample.type.comparison = NULL,coldata=NULL,verbose=FALSE,fc.method='deseq2') {
    ## check params
    if(is.null(ens.p2)) {stop('ens.p2 is null')}
    if(is.null(celltype)) { stop('celltype is null') }
    if(is.null(sample.type.comparison)) {stop('sample type comparison is null')}
    if(is.null(coldata)) {stop('coldata is null')}

    ## Get approprate submatrice
    x <- t(getSamples(ens.p2, sample.type.comparison[2], celltype))
    x1 <- t(getSamples(ens.p2, sample.type.comparison[1],  celltype))
    
    ## genes x aggregates
    x12.a <- cbind(x,x1)
    if (fc.method == 'deseq2') {
        ## Get differential expression  
        dds1 <- DESeq2::DESeqDataSetFromMatrix(x12.a, coldata[colnames(x12.a),], design=~sample.type)
        dds1 <- DESeq2::DESeq(dds1)
        ## extract fold changes
        res <- DESeq2::results(dds1,cooksCutoff = FALSE, independentFiltering = FALSE)
        fcs <- res$log2FoldChange
        names(fcs) <- rownames(res)
        ## TODO offer option to set these to 0
        fcs <- fcs[!is.na(fcs)]
    } else if (fc.method == 'deseq1') {
        condition <- coldata[colnames(x12.a),]$sample.type
        cds <- DESeq::newCountDataSet(x12.a, condition)
        cds <- DESeq::estimateSizeFactors(cds)
        cds <- DESeq::estimateDispersions(cds)
        res <- DESeq::nbinomTest(cds, sample.type.comparison[1], sample.type.comparison[2])
        ## extract fcs
        fcs2 <- res$log2FoldChange
        names(fcs2) <- res$id
        fcs <- fcs2[!is.na(fcs)]
    } else if (fc.method == 'dummy') {
        ## No correction
        genes <- rownames(x12.a)
        fcs <- rep(0,times=length(genes))
        names(fcs) <- genes
    } else if (fc.method == 'simple') {
        ## Simple fold change between the two conditions for each gene
        x.cpm <- sweep(x,2,apply(x,2,sum),FUN='/') * 1e6
        x1.cpm <- sweep(x1,2,apply(x1,2,sum),FUN='/') * 1e6
        fcs <- log2((apply(x.cpm,1,mean) + 1) / (apply(x1.cpm,1,mean) + 1))
    } else {
        stop('Unknown fc correction method');
    }
    fcs
}

#' get a list of fold changes per cell type
#' @param ens.p2 pagdoda2 ensembl object
#' @param levels levels of the cell type assignment factor
#' @param sample.type.comparison comparison to perform
#' @param coldata metadata
#' @param n.cores number of cores
#' @param verbose verbosity logical
#' @param fc.method method to use to calculate fcs, passed to getCelltypeFCs
getPerCellTypeFCs <- function(ens.p2  =NULL, levels =  NULL, sample.type.comparison = NULL, coldata=NULL, n.cores=1, verbose=FALSE,fc.method='deseq2') {
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
            getCelltypeFCs(ens.p2 = ens.p2, celltype=celltype,
                           sample.type.comparison=sample.type.comparison,
                           coldata=coldata,
                           verbose=verbose,
                           fc.method=fc.method)
        })
    },mc.cores=n.cores)
    fcs
}

#' adaptr to convert deseq2 results object to suitable data.frame
#' @param res deseq2 results object
#' @param membrane.gene.names genes to annotate as membrane genes
ddsres2res <- function(res,membrane.gene.names) {
                                        # prep res
                                        #res <- DESeq2::results(dds)
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
    res
}

## library(codetools);
## findGlobals(t.test.correctFC)

t.test.correctFC <- function(x12, coldata, fc.correction, correction.global.weight, membrane.gene.names,
                             sample.type.comparison, n.cores=1) {
    
    ## Calculate effective correction in the linear scale
    effective.correction <- 2^(fc.correction * correction.global.weight)

    ## can't correct for fcs right now so throw a warning if corrections are requestes
    if(!all((abs(effective.correction)-1) < 1e-3)) {
        warning('t-test does NOT support FC correction, non-zero effective weights detected')
    }

    ## Get the sample names we are comparing
    samplesA <- rownames(coldata)[coldata$sample.type == sample.type.comparison[1]]
    samplesB <- rownames(coldata)[coldata$sample.type == sample.type.comparison[2]]
    ## Subset to samples in our input matrices
    samplesA <- samplesA[samplesA %in% colnames(x12)]
    samplesB <- samplesB[samplesB %in% colnames(x12)]

    ## calculate cpms
    cpmmat <- sweep(x12, 2, apply(x12,2,sum), FUN='/') * 1e6
    baseMean <- apply(cpmmat, 1, mean)
    allgenes <- rownames(cpmmat)
    
    ## We need at least 2 samples per condition for t-tests
    if (length(samplesA) > 1 & length(samplesB) >1) {
        ## Run t-tests on cpms
        tts <- parallel::mclapply(1:nrow(cpmmat), function(m) {
            tryCatch({
                t.test(cpmmat[m,samplesA], cpmmat[m,samplesB])
            })
        },mc.cores=n.cores)

        ## calculate result table entries
        p.vals <- unlist(lapply(tts, function(x) {x$p.value}))
        padj <- p.adjust(p.vals, 'fdr')
        stat <- unname(unlist(lapply(tts, function(x) {x$statistic})))

        ## put results table together
        res <- data.frame(
            gene = allgenes,
            baseMean = baseMean,
            log2FoldChange = rep(0,nrow(cpmmat)),
            lfcSE = rep(0,nrow(cpmmat)), ## just for compatibility with deseq2 results, not used
            stat = stat,
            pvalue = p.vals,
            padj = padj,
            significant = padj < 0.05,
            membrane = rownames(cpmmat) %in% membrane.gene.names,
            ## TODO
            Z = rep(0,nrow(cpmmat)),
            Za = rep(0,nrow(cpmmat))
        )
    } else {
        warning(paste0('Not enough samples: ',length(samplesA), ' vs ', length(samplesB)))
        
        res <- data.frame(
            gene = allgenes,
            baseMean = baseMean,
            log2FoldChange = rep(0,nrow(cpmmat)),
            lfcSE = rep(0,nrow(cpmmat)), ## just for compatibility with deseq2 results, not used
            stat = rep(NA,nrow(cpmmat)),
            pvalue = rep(NA,nrow(cpmmat)),
            padj = rep(NA,nrow(cpmmat)),
            significant = rep(FALSE,nrow(cpmmat)),
            membrane = rep(FALSE,nrow(cpmmat)),
            ## TODO
            Z = rep(0,nrow(cpmmat)),
            Za = rep(0,nrow(cpmmat))
        )
    }
}

#' get differential expression corrected by the specified vector
#' @param ens.p2 pagoda2 ensembl object
#' @param cell.type cell type to do the comparison of
#' @param sample.type.comparison comparison to perform
#' @param coldata metadata
#' @param fc.correction vector of fold changes to correct
#' @param membrane.gene.names character vector of genes that are membrane
#' @param de.method method to do differential expression currently deseq2
#' @param correction.global.weight global weight to apply to correction
getCorrectedDE <- function(ens.p2, cell.type=NULL, sample.type.comparison=NULL, coldata=NULL, fc.correction=NULL,
                           membrane.gene.names=NULL,de.method='deseq2', correction.global.weight=1,n.cores =1) {
    ## Check arguments
    if(is.null(ens.p2)) {stop('ens.p2 is null')}
    if(is.null(cell.type)) {stop('cell.type is null')}
    if(is.null(sample.type.comparison)) {stop('sample.type.comparison is null')}
    if(is.null(coldata)) {stop('coldata is null')}
    if(is.null(fc.correction)) {stop('fc.correction is null')}

    
    
    ## get data
    x12 <- t(rbind(
        getSamples(ens.p2, sample.type.comparison[2], cell.type),
        getSamples(ens.p2, sample.type.comparison[1], cell.type)
    ))
    ## run de
    if(de.method=='deseq2'){
        dds <- DESeq2.correctFC(x12, coldata,  fc.correction, correction.global.weight=correction.global.weight)
        ## convert de to results table
        res <- DESeq2::results(dds)
        allgenes <- rownames(res);
        res <- ddsres2res(res,membrane.gene.names)
    } else if (de.method=='t.test') {
        res <- t.test.correctFC(x12, coldata, fc.correction,
                                correction.global.weight=correction.global.weight,
                                sample.type.comparison=sample.type.comparison,
                                membrane.gene.names=membrane.gene.names)
        allgenes <- rownames(res)
    }  else {
        stop('Unknown DE method: ',de.method);
    }
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
        snames=levels(coldata$sample.type),
        correction=fc.correction
    ))
}

#' Get a correction vector by combining fold changes over all the results
#'
#' @description Get a correction vector by combining fold changes of every gene accross multiple cell types
#' into a single value. The correction change is weighted by the chi square distribution so that genes
#' with variable fold changes do not get corrected.
#'
#' Two different cell groupings need to be provided. fc.cellfactor dictates the groupings of the cells used
#' for collapsing, this should either be identical or more coarse than cell.type.factor. If this is provided a
#' many to one map from clusters of cell.type.factor to fc.cellfactor is required as cell.factor.map
#' 
#' @param ens.p2 ensembl p2 object
#' @param aggregation.id the aggregation data slot to use from the ens.p2 object
#' @param sample.type.comparison character vector of 'sample.type' to compare e.g. T vs W
#' @param cell.types.exclude cell types to ignore from the fine-grained cell.type.factor after mapping to coarser clusters
#' @param scaleByVariance scale the mean fold changes by the Variance
#' @param useTrimmed return results from trimmed mean
#' @param n.cores number of cores to use
#' @param coldata metadata for aggregated data
#' @param cell.type.factor named factor of cell type assignments to groups
#' @param per.cell.type.fcs pre-calculated per cell type FCs, if null they will be recalculated here
#' @param verbose logical verbosity
#' @param fc.method method to pass to getPerCellTypeFCs() if precalculated per.cell.type.fcs are not provided
#' @param cell.factor.map many to one mapping of levels of cell.type.factor to levels of fc.cellfactor
#' @param fc.cellfactor cell grouping to use for the fold change generation
calcFCcorrection <- function(ens.p2, aggregation.id, sample.type.comparison, cell.types.exclude = c('tumor'),
                             scaleByVariance=TRUE, useTrimmed=FALSE,n.cores=1,coldata=NULL, cell.type.factor=NULL,
                             per.cell.type.fcs=NULL,verbose=FALSE,fc.method='deseq2',cell.factor.map=NULL, fc.cellfactor = NULL) {

    if (is.null(cell.type.factor)) {
        stop('cell.type.factor is NULL')
    }

    if (is.null(fc.cellfactor)) {
        warning('fc.cellfactor is NULL, using cell.type.factor levels')
        fc.cellfactor <- cell.type.factor
    }
    
    ## Check if precalculated per cell type fcs for the coarse clusters
    if(is.null(per.cell.type.fcs)) {
        fcs <- getPerCellTypeFCs(levels = levels(fc.cellfactor),
                                 ens.p2=ens.p2,
                                 sample.type.comparison=sample.type.comparison,
                                 coldata=coldata,
                                 n.cores=n.cores,
                                 fc.method=fc.method)
    } else {
        if (verbose) cat('Using pre-calculated per celltype fcs');
        ## NOTE fc.method is ignored here, it is assumed that the calling function knows
        ## what it is doing
        fcs <- per.cell.type.fcs
    }
    

    fcs <- fcs[!unlist(lapply(fcs, is.error))]

    ## Map from fine to coarse clusters
    if (!is.null(cell.factor.map)) {
        excl.coarse <- cell.factor.map[cell.types.exclude]
    } else {
        ## if not cell.factor map is provided and we are here then
        ## we assume that we are not dealing with two levels of clusters
        excl.coarse <- cell.types.exclude
    }

    ## Keep only desired FCs
    fcs <- fcs[!names(fcs) %in% excl.coarse]
    
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
    ## fill missing trimmed means
    fc.mean.var$trimmedMean[is.na(fc.mean.var$trimmedMean)] <- 0
    
    ## Scale correction by variance
    fc.mean.var$var[is.na(fc.mean.var$var)] <- 0 # fill missing var
    
    ## Assuming normal distribution of the fcs
    ## we will weight them with 1 -  chi-square probability
    ## the degrees of freedom are the number of comparisons minu1
    fc.mean.var$weight <- 1 - pchisq(q = fc.mean.var$var, df = nrow(fcs.mat) - 1)

    fc.mean.var$correction.fc <- fc.mean.var$mean * fc.mean.var$weight # fc to use for correct
    fc.mean.var$correction.fc.trimmed <- fc.mean.var$trimmedMean * fc.mean.var$weight # fc to use for correct with trimming
    ## Return the requested result
    if (scaleByVariance) {
        if (useTrimmed) {
            ## fc correction vector
            fc.correction <- fc.mean.var[,'correction.fc']
            names(fc.correction) <- rownames(fc.mean.var)
        } else {
            fc.correction <- fc.mean.var[,'correction.fc.trimmed']
            names(fc.correction) <- rownames(fc.mean.var)
        }
    } else {
        if (useTrimmed) {
            fc.correction <- fc.mean.var[,'trimmedMean']
            names(fc.correction) <- rownames(fc.mean.var)
        } else {
            fc.correction <- fc.mean.var[,'mean']
            names(fc.correction) <- rownames(fc.mean.var)
        }
    }
    fc.correction
}

#' save results of de expression as json files with designated prefix (can be directory)
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

#' @param ens.p2 pagoda 2 object ensemble
#' @param cellfactor a named factor with cell assignments to groups
#' @param sample.type.comparison character vector of length 2 signifying the comparison to perform from the sample.type metadata column
#' @param membrane.gene.names a character vector of genes to annotate as membrane genes
#' @param n.cores number of cores to use
#' @param correction.method character of lenght 1 global or exclcurrent, exclcurrent corrects with a vector that excludes considered cel tyep
#' @param cell.types.fc.exclude cell types not to consider in generating the fold changes
#' @param verbose verbosity logical
#' @param de.method method for differential expression, currently only deseq2 supported
#' @param fc.method method for fc generation for the correction, deseq2 or dummy (no correction) currently supported
#' @param correction.global.weight global weighting of correction, 0 is no correction, 1 default
#' @param fc.cellfactor cell factor to use for obtaining FC correction, if NULL cellfactor is used. If this is specified and the correction method is not exclcurrent, then a map from the levels of cellfactor to fc.cellfactor (many to one) need to be specified to know which clusters to drop when excluding the current
#' @param cell.factor.map many to one mapping from cellfactor levels to fc.cellfactor levels
getCorrectedDE.allTypes <-  function(ens.p2, cellfactor, sample.type.comparison, 
                                     membrane.gene.names,n.cores=1,correction.method='global',
                                     cell.types.fc.exclude=NULL,verbose=FALSE,de.method='deseq2',
                                     fc.method='deseq2', correction.global.weight=1, fc.cellfactor = NULL,
                                     cell.factor.map = NULL) {
    ## Check arguments -- TODO add more
    if(!correction.method %in% c('global','exclcurrent')) {stop('Unknown correction method: ',correction.method)}

    if (is.null(fc.cellfactor)) {
        fc.cellfactor <- cellfactor
    } else {
        if (correction.method=='exclcurrent') {
            if (!is.null(cell.factor.map)) {
                ## Check that the cell.factor.map is valid
                if (!(contain.identical(levels(cellfactor), names(cell.factor.map)) &
                      contain.identical(cell.factor.map, levels(fc.cellfactor),check.unique=FALSE))) {
                    stop('Invalid cell.factor.map provided');
                }
            }  else {
                stop("Correction method is 'exclcurrent' and 'fc.cellfactor' is specified without 'cell.factor.map': cluster correspondence for dropout can't be established. Use global correction.method or specify cell.factor.map")
            }
        }
    }
    
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
                                           useTrimmed=FALSE, n.cores=n.cores,
                                           coldata=coldata,
                                           cell.type.factor = cellfactor,
                                           fc.cellfactor = fc.cellfactor,
                                           verbose=verbose,
                                           fc.method=fc.method)
        per.cell.type.fcs <- NULL
    } else {
        ## Cache the per cell type corrections (from which global or per cell type corrections are derived)
        ## This prevents recomputing them for every single cell type
        per.cell.type.fcs <- getPerCellTypeFCs(levels=levels(fc.cellfactor),#levels(cellfactor),
                                               ens.p2=ens.p2,
                                               sample.type.comparison=sample.type.comparison,
                                               coldata=coldata,
                                               n.cores=n.cores,
                                               verbose=verbose,
                                               fc.method=fc.method)
    }

    ## Run differential expression
    res <- parallel::mclapply(setNames(levels(cellfactor)), function(cell.type) {
        cat('Processing ',cell.type,'...')
        try({
            if (correction.method == 'exclcurrent') {
                ## Get cell type specific correction using the cached results per.cell.type.fcs
                ## essentially collapse the cell specific fold changes ignoring the current
                ## cell types, but without recalculating them all
                fc.correction <- calcFCcorrection(ens.p2=ens.p2,
                                                  aggregation.id='cluster:sample',
                                                  sample.type.comparison=sample.type.comparison,
                                                  cell.types.exclude=c(cell.types.fc.exclude,cell.type),
                                                  scaleByVariance=TRUE,
                                                  useTrimmed=FALSE,
                                                  n.cores=1,
                                                  coldata=coldata,
                                                  cell.type.factor = fc.cellfactor, #cellfactor,
                                                  per.cell.type.fcs=per.cell.type.fcs,
                                                  cell.factor.map = cell.factor.map)
            }
            getCorrectedDE(ens.p2=ens.p2,
                           cell.type=cell.type,
                           sample.type.comparison=sample.type.comparison,
                           coldata=coldata,
                           fc.correction=fc.correction,
                           membrane.gene.names = membrane.gene.names,
                           de.method=de.method,
                           correction.global.weight=correction.global.weight)
            
        })
    }, mc.cores=n.cores)
    
    ## remove failed
    res[unlist(lapply(res,is.error))] <- NULL
    
    ## return
    res
}



getDEcount <- function(de.result.set = NULL, type = c('all','up','down')) {
    ## check input
    if (is.null(de.result.set) ) stop('de.result.set not provided')
    if (!type %in% c('all','up','down')) stop('type not valid')
    ## Get summary stats
    x <- lapply(de.result.set, function(type.comparison) {
        lapply(type.comparison, function(cell.comparison) {
            ##cell.comparison <- de.result.set[[1]][[1]]
            res <- cell.comparison$res
            ## calc results
            summary <- list(
                all = sum(res$significant, na.rm=TRUE),
                up =  sum(res$significant & res$log2FoldChange > 0, na.rm=TRUE),
                down =  sum(res$significant & res$log2FoldChange < 0, na.rm=TRUE)
            )
            ## return
            summary[[type]]
        })
    })
    ## convert to array
    x <- melt(x)
    acast(x, L1 ~ L2, value.var='value',fill=0)
}




getFCmatrix <- function(de.result.set = NULL, type.comparison.name = NULL) {
    if (is.null(de.result.set)) stop('de.result.set not specified')
    if (is.null(type.comparison.name)) stop('type.comparison.name not specified')
    type.comparison <- de.result.set[[type.comparison.name]]
    ##
    x <- lapply(type.comparison, function(cell.comparison) {
        fcs <-  cell.comparison$res$log2FoldChange
        names(fcs) <- rownames(cell.comparison$res)
        fcs
    })
    ##
    x <- do.call(cbind, x)
    x[is.na(x)] <- 0
    ##
    x    
}


getSignmatrix <- function(de.result.set = NULL, type.comparison.name = NULL) {
    if (is.null(de.result.set)) stop('de.result.set not specified')
    if (is.null(type.comparison.name)) stop('type.comparison.name not specified')
    type.comparison <- de.result.set[[type.comparison.name]]
    ##
    x <- lapply(type.comparison, function(cell.comparison) {
        sign <-  cell.comparison$res$significant
        names(sign) <- rownames(cell.comparison$res)
        sign
    })
    ##
    x <- do.call(cbind, x)
    x[is.na(x)] <- 0
    ##
    x    
}

getComparisonFCpca <- function(de.res, comparison.name) {
    fcs <- getFCmatrix(de.res, comparison.name)
    pca0 <- prcomp(t(fcs))
    ret <- as.data.frame(pca0$x)
    ret$sample <- rownames(ret)
    ret
}

getFCcormat <- function(de.result.set,type.comparison.name,method='pearson',sign.only=TRUE) {
    fc.mat <- getFCmatrix(de.result.set,type.comparison.name)
    if(sign.only) {
        sign.mat <- getSignmatrix(de.result.set,type.comparison.name)
        sign.genes <- rownames(sign.mat)[rowSums(sign.mat) > 0]
        fc.mat <- fc.mat[sign.only,]
    }
    cor.mat <- cor(fc.mat[sign.genes,], method=method)
    diag(cor.mat) <- 0
    cor.mat
}


generateSerializedAppsWithJC <- function(ens.p2, jc, prefix="") {
    lapply(names(ens.p2$p2objs), function(n) {
        p2 <- ens.p2$p2objs[[n]]
        jc2 <- jc[rownames(p2$counts)]
        extra.meta = list(
            jointClustering=p2.metadata.from.factor(jc2, displayname='joint clustering')
        );
        p2web <- basicP2web(p2,app.title=n,extra.meta,n.cores=4)
        p2web$serializeToStaticFast(binary.filename=paste0(prefix,n,'.bin'));
    })
    invisible(NULL)
}
