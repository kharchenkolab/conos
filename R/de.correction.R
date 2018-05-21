## Functions for differential expression of single cell panels

## TODO: fix me
library(DESeq2)



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
#' @param counts.add.bg pseudocount to add to all counts, default NULL will add a method dependent adjustment (edgeR will have a pc, other methods not)
getCorrectedDE.allTypes <-  function(ens.p2, cellfactor, sample.type.comparison, 
                                     membrane.gene.names,n.cores=1,correction.method='global',
                                     cell.types.fc.exclude=NULL,verbose=FALSE,de.method='deseq2',
                                     fc.method='deseq2', correction.global.weight=1, fc.cellfactor = NULL,
                                     cell.factor.map = NULL,
                                     counts.add.bg = NULL) {
    ## Check arguments -- TODO add more
    if(!correction.method %in% c('global','exclcurrent','none')) {stop('Unknown correction method: ',correction.method)}

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
    
    ## Ensure levels are consistent with those provided
    coldata$sample.type <- factor(coldata$sample.type,levels=sample.type.comparison)

    ## If correction is global calculate it once, otherwise...
    if (correction.method == 'global') {
        if(verbose) cat('Calculating global correction...\n')
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
                                           fc.method=fc.method,
                                           counts.add.bg=counts.add.bg)
        per.cell.type.fcs <- NULL
    } else if (correction.method == 'exclcurrent') {
        ## ... cache the per cell type corrections (from which global or per cell type corrections are derived)
        per.cell.type.fcs <- getPerCellTypeFCs(levels=levels(fc.cellfactor),#levels(cellfactor),
                                               ens.p2=ens.p2,
                                               sample.type.comparison=sample.type.comparison,
                                               coldata=coldata,
                                               n.cores=n.cores,
                                               verbose=verbose,
                                               fc.method=fc.method,
                                               counts.add.bg=counts.add.bg)
    } else if (correction.method == 'none') {
        fc.correction <- NULL;
    } else {
        stop('error unknown correction method');
    }

    ## Run differential expression for each cell type
    res <- parallel::mclapply(setNames(levels(cellfactor)), function(cell.type) {
        cat('Processing ',cell.type,'...\n')
        try({
            ## if excl current prepare a fc correction from the individual per
            ## cell type corrections
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
                                                  cell.factor.map = cell.factor.map,
                                                  counts.add.bg=counts.add.bg)
            }
            
            getCorrectedDE(ens.p2=ens.p2,
                           cell.type=cell.type,
                           sample.type.comparison=sample.type.comparison,
                           coldata=coldata,
                           fc.correction=fc.correction,
                           membrane.gene.names = membrane.gene.names,
                           de.method=de.method,
                           correction.global.weight=correction.global.weight,
                           counts.add.bg=counts.add.bg)
            
        })
    }, mc.cores=n.cores)
    
    ## remove failed
    res[unlist(lapply(res,is.error))] <- NULL
    
    ## return
    res
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
                           membrane.gene.names=NULL,de.method='deseq2', correction.global.weight=1,n.cores =1, counts.add.bg=NULL) {
    ## Check arguments
    if(is.null(ens.p2)) {stop('ens.p2 is null')}
    if(is.null(cell.type)) {stop('cell.type is null')}
    if(is.null(sample.type.comparison)) {stop('sample.type.comparison is null')}
    if(is.null(coldata)) {stop('coldata is null')}
#    if(is.null(fc.correction)) {stop('fc.correction is null')}, null fc correction is fine now and indicates no correction
    ## get data
    samples2 <- getSamples(ens.p2, sample.type.comparison[2], cell.type)
    samples1 <- getSamples(ens.p2, sample.type.comparison[1], cell.type)
    if (nrow(samples2) < 2 | nrow(samples1) < 2) {
        ## not enough samples for comparison here
        stop('Not enough samples for comparison')
    }
    ## Put together
    x12 <- t(rbind(samples2, samples1))
    ## run de
    if(de.method=='deseq2'){
        if(is.null(counts.add.bg)) { ## default to provided matrix
            cm <- x12
        } else {
            cm <- x12 + counts.add.bg
        }
        dds <- DESeq2.correctFC(cm, coldata,  fc.correction, correction.global.weight=correction.global.weight)
        ## convert de to results table
        res <- DESeq2::results(dds)
        allgenes <- rownames(res);
        res <- ddsres2res(res,membrane.gene.names)
    } else if (de.method=='edgeR') {
        if (is.null(fc.correction)) { stop(paste0('Correction method ',de.method,' does not support uncorrected results and fc.correction is NULL')) };
        ## TODO: factor this out in a function
        grp1 <- coldata[colnames(x12),]$sample.type
        if(length(unique(grp1)) == 2) {            
            if (counts.add.bg == 0) warning('counts.add.bg is zero and de.method is edgeR, this is not recommended')
            if (is.null(counts.add.bg)) {
                counts.add.bg <- 1 ## default to adding 1
            }
            ## put together an object just for norm factos
            cm <- x12 + counts.add.bg ## required to stabilise
            fc.correction <- fc.correction[!is.na(fc.correction)]
            cm <- cm[rownames(cm) %in% names(fc.correction),] ## only keep genes wer can correct
            y <- DGEList(counts=cm, group=grp1)
            y <- calcNormFactors(y)
            ## of1 is array of library sizes
            gene.scale <- matrix(rep(0,length(y$counts)),nrow=nrow(y$counts),ncol=ncol(y$counts))
            of1 <- sweep(gene.scale,2,log(y$samples$lib.size * y$samples$norm.factors),FUN='+')
            rownames(of1) <- rownames(y$counts); colnames(of1) <- colnames(y$counts)
            ## put fc correction in order
            gene.scale.factors <- fc.correction[rownames(y$counts)]
            ## construct array of gene correctiosn
            pn <- do.call(cbind, lapply(grp1, function(x) {
                if (x == levels(grp1)[1]) {
                    rep(0,nrow(of1))
                } else {
                    gene.scale.factors
                }
            }))
            ## merge gene correctiosn with library sizes
            of2 <- of1 + pn
            yc <- DGEList(counts=cm,group=grp1)
            yc$offset <- of2
            yc <- calcNormFactors(yc)
            design <- model.matrix(~grp1)
            ## using loess here because default 'locfit' fails with memory errors under some circuimstances
            yc <- estimateDisp(y=yc, design=design, trend.method='loess') 
            fit2 <- glmQLFit(yc,design)
            qlf2 <- glmQLFTest(fit2, coef=2)
            tt2 <- as.data.frame(topTags(qlf2, n=Inf))
            ## Format results as other methods, padding missing values with 0
            res <- data.frame(
                gene=rownames(tt2),
                baseMean=tt2$logCPM,
                log2FoldChange=tt2$logFC,
                lfcSE =rep(0,nrow(tt2)),
                stat=rep(0,nrow(tt2)),
                pvalue=tt2$PValue,
                padj=tt2$FDR,
                significant=tt2$FDR < 0.05,
                membrane=rownames(tt2) %in% membrane.gene.names,
                Z=rep(0,nrow(tt2)),
                Za=rep(0,nrow(tt2))
            )
            rownames(res) <- rownames(tt2);
            allgenes <- rownames(tt2);
        } else {
            ## no contrast to run, return empty resutls
            res <- data.frame(
                gene=rownames(x12),
                baseMean=rep(0,nrow(x12)),
                log2FoldChange=rep(0,nrow(x12)),
                lfcSE =rep(0,nrow(x12)),
                stat=rep(0,nrow(x12)),
                pvalue=rep(1,nrow(x12)),
                padj=rep(1,nrow(x12)),
                significant=rep(FALSE,nrow(x12)),
                membrane=rep(FALSE,nrow(x12)),
                Z=rep(0,nrow(x12)),
                Za=rep(0,nrow(x12))
            )
            allgenes <- rownames(x12)
        }
    } else if (de.method == 'deseq1') {
        if (is.null(fc.correction)) { stop(paste0('Correction method ',de.method,' does not support uncorrected results and fc.correction is NULL') )};
        ## TODO: move this to a function
        ## can't correct for fcs right now so throw a warning if corrections are requestes
        if(!all((abs(fc.correction)) < 1e-3)) {
            warning('deseq1 does NOT support FC correction, non-zero effective weights detected')
        }

        fc.correction <- fc.correction[!is.na(fc.correction)]
        #fc.correction <- fc.correction[is.finite(fc.correction)]

        ## Subsets to genes for which we have fc information
        x12.b <- x12[names(fc.correction),]

        ## Run differential expression
        condition <- coldata[colnames(x12.b),]$sample.type
        cds <- DESeq::newCountDataSet(x12.b, condition)
        cds <- DESeq::estimateSizeFactors(cds)
        cds <- DESeq::estimateDispersions(cds)
        
        ## TODO: apply normalisation factors here
        ## Note it's unclear if this is possible in deseq1, seems to be present in
        ## only some version of the normalisation
        ## normFactors <- matrix(runif(nrow(cds)*ncol(cds),0.5,1.5),
        ##                       ncol=ncol(cds),nrow=nrow(cds),
        ##                       dimnames=list(1:nrow(cds), 1:ncol(dds)))
                              
        ## normFactors <- normFactors / exp(rowMeans(log(normFactors)))
        ## DESeq::normalizationFactors(cds) <- normFactors
        
        res2 <- DESeq::nbinomTest(cds, sample.type.comparison[1], sample.type.comparison[2])

        ## Format results as other methods, padding missing values with 0
        res3 <- data.frame(
            gene=res2$id,
            baseMean=res2$baseMean,
            log2FoldChange=res2$log2FoldChange,
            lfcSE =rep(0,nrow(res2)),
            stat=rep(0,nrow(res2)),
            pvalue=res2$pval,
            padj=res2$padj,
            significant=res2$padj < 0.05,
            membrane=res2$id %in% membrane.gene.names,
            Z=rep(0,nrow(res2)),
            Za=rep(0,nrow(res2))
        )
        rownames(res3) <- res2$id
        res <- res3
        
    } else if (de.method=='t.test') {
        ## actually it does and can only do that but the fc correction need to be specified fo rnow
        if (is.null(fc.correction)) { stop(paste0('Correction method ',de.method,' does not support uncorrected results and fc.correction is NULL')) };
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
    if (!is.null(gene.scale.factor)) {
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
    } else {
        ## run vanilla de
        dds4 <- DESeq2::DESeqDataSetFromMatrix(x12, coldata[colnames(x12),], design=~sample.type)
        dds4 <- DESeq(dds4)
    }
    dds4
}

#' get fold changes for a specific cell type comparison
#' @param pagoda 2 ensemble object
#' @param celltype to use
#' @param sample.type.comparisons comparison to perform
#' @param coldata metadata
#' @param verbose verbosity logical
#' @param fc.method correction method to use 'deseq2' or 'dummy' (no correction)
getCelltypeFCs <- function(ens.p2=NULL, celltype=NULL, sample.type.comparison = NULL,coldata=NULL,verbose=FALSE,fc.method='deseq2',counts.add.bg=NULL) {
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
    
    st1 <- coldata[colnames(x12.a),]$sample.type
    
    if(length(st1) < 2) {
        warning(paste0('Insufficient samples in FC calculation for ',celltype,'. Returning 0s'))
        fcs <- rep(0, length(rownames(x12.a)))
        names(fcs) <- rownames(x12.a)
    } else if (length(unique(st1)) != 2) {
        warning(paste0('Wrong number of sample types for ',celltype,'. Returning 0s'))
        fcs <- rep(0, length(rownames(x12.a)))
        names(fcs) <- rownames(x12.a)
    } else {
        if (fc.method == 'deseq2') {
            if(is.null(counts.add.bg)) counts.add.bg <- 0
            ## Get differential expression  
            dds1 <- DESeq2::DESeqDataSetFromMatrix(x12.a + counts.add.bg, coldata[colnames(x12.a),], design=~sample.type)
            dds1 <- DESeq2::DESeq(dds1)
            ## extract fold changes
            res <- DESeq2::results(dds1,cooksCutoff = FALSE, independentFiltering = FALSE)
            fcs <- res$log2FoldChange
            names(fcs) <- rownames(res)
            ## TODO offer option to set these to 0
            fcs <- fcs[!is.na(fcs)]
        } else if (fc.method == 'deseq1') {
            if (!is.null(counts.add.bg)) warning('deseq1 FC calculation ignores value of counts.add.bg')
            condition <- coldata[colnames(x12.a),]$sample.type
            cds <- DESeq::newCountDataSet(x12.a, condition)
            cds <- DESeq::estimateSizeFactors(cds)
            cds <- DESeq::estimateDispersions(cds)
            res <- DESeq::nbinomTest(cds, sample.type.comparison[1], sample.type.comparison[2])
            ## extract fcs
            fcs2 <- res$log2FoldChange
            names(fcs2) <- res$id
            fcs <- fcs2[!is.na(fcs2)]
        } else if (fc.method == 'edgeR') {
            ## edgeR defaults to bg correction of +1
            if(is.null(counts.add.bg)) counts.add.bg <- 1
            group <- coldata[colnames(x12.a),]$sample.type
            y <- DGEList(counts=x12.a+counts.add.bg,group=group)
            y <- calcNormFactors(y)
            design <- model.matrix(~group)
            y <- estimateDisp(y,design)
            fit <- glmQLFit(y,design)
            qlf <- glmQLFTest(fit,coef=2)
            tt1 <- topTags(qlf,n=Inf)
            fcs <- tt1$table$logFC / (1/log(2))
            names(fcs) <- rownames(tt1$table)
            ## TODO set method used to obtain fcs as attribute
        } else if (fc.method == 'dummy') {
            if (!is.null(counts.add.bg)) warning('dummy FC calculation ignores value of counts.add.bg')
            ## No correction
            genes <- rownames(x12.a)
            fcs <- rep(0,times=length(genes))
            names(fcs) <- genes
        } else if (fc.method == 'simple') {
            if (!is.null(counts.add.bg)) warning('simple FC calculation ignores value of counts.add.bg')
            ## Simple fold change between the two conditions for each gene
            x.cpm <- sweep(x,2,apply(x,2,sum),FUN='/') * 1e6
            x1.cpm <- sweep(x1,2,apply(x1,2,sum),FUN='/') * 1e6
            fcs <- log2((apply(x.cpm,1,mean) + 1) / (apply(x1.cpm,1,mean) + 1))
        } else {
            stop('Unknown fc correction method');
        }
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
getPerCellTypeFCs <- function(ens.p2  =NULL, levels =  NULL, sample.type.comparison = NULL, coldata=NULL, n.cores=1, verbose=FALSE,fc.method='deseq2', counts.add.bg=NULL) {
    ## Check args
    if (is.null(levels)) { stop('Provided levels is null') }
    if (is.null(ens.p2)) {stop('ens.p2 is null')}
    if (is.null(sample.type.comparison)) {stop('sample type comparison is null')}
    if (is.null(coldata)) {stop('coldata is NULL')}
    ## Calculate the fold changes for each cell type
    fcs <- parallel::mclapply(setNames(levels), function(celltype) {
        celltype <- as.character(celltype)
        cat('Calculating FC for ', celltype,'\n');
        try({
            getCelltypeFCs(ens.p2 = ens.p2, celltype=celltype,
                           sample.type.comparison=sample.type.comparison,
                           coldata=coldata,
                           verbose=verbose,
                           fc.method=fc.method,
                           counts.add.bg=NULL)
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
                             per.cell.type.fcs=NULL,verbose=FALSE,fc.method='deseq2',cell.factor.map=NULL,
                             fc.cellfactor = NULL, counts.add.bg = NULL) {

    if (is.null(cell.type.factor)) {
        stop('cell.type.factor is NULL')
    }

    if (is.null(fc.cellfactor)) {
        #cat('calcFCcorrection: fc.cellfactor is NULL, using cell.type.factor levels\n')
        fc.cellfactor <- cell.type.factor
    }
    
    ## Check if precalculated per cell type fcs for the coarse clusters
    if(is.null(per.cell.type.fcs)) {
        fcs <- getPerCellTypeFCs(levels = levels(fc.cellfactor),
                                 ens.p2=ens.p2,
                                 sample.type.comparison=sample.type.comparison,
                                 coldata=coldata,
                                 n.cores=n.cores,
                                 fc.method=fc.method,
                                 counts.add.bg = counts.add.bg)
    } else {
        if (verbose) cat('Using pre-calculated per celltype fcs\n');
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
