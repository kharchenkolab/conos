
#' Extract count matrices for differential expression
#' and proportion comparisons from pagoda2 apps
#' @param apps a list of pagoda2 apps
#' @return a list of dgCMatrices
extractCountMatrices <- function(apps) {
    require('parallel')
    ## Get common genes
    gl <- lapply(apps, function(r) colnames(r$counts))
    all.genes <- unique(unlist(gl))
    gc <- do.call(rbind, lapply(gl, function(x) all.genes %in% x))
    common.genes <- all.genes[apply(gc,2,sum) > 3]
    ## Extract count matrices for common gene list
    ccm.raw <- mclapply(apps,function(r) {
        om <- as(r$misc$rawCounts,'dgTMatrix')
        om <- om[,colnames(om) %in% common.genes]
        mi <- match(colnames(om), common.genes)
        x <- new("dgTMatrix",i = om@i,j = as.integer(mi[om@j+1]-1),x=om@x,Dim=c(nrow(om),length(common.genes)))
        rownames(x) <- rownames(om);
        colnames(x) <- common.genes
        as(x,'dgCMatrix')
    }, mc.cores =30)
    ccm.raw
}


#' Proportion comparison of cells in apps
#' @param rl1 list of matrices for condition 1 with rownames as cells
#' @param rl2 list of matrices for condition 2 with rownames as cells
#' @param cells vector of cell names beloging to the group being tested
#' @param prel.cells total cell space
#' @examples
#' all.t.cells <- names(nfac)[nfac %in% grep("^T",levels(nfac),value=T)]
#' ttestres <- lapply(sn(grep("^T",levels(nfac),value=T)),function(n1) {
#'     t.two.prop.comp(rl1,rl2,names(nfac)[nfac==n1],'normal',n2='tumor',prel.cells=all.t.cells)
#' })
#' do.call(rbind,lapply(ttestres,function(x) x$qres))
t.two.prop.comp <- function(rl1, rl2, cells, n1='test', n2='control', prel.cells=NULL) {
    if(is.null(prel.cells)) {
        prel.cells <- c(unlist(lapply(rl1,rownames)),unlist(lapply(rl2,rownames)))
    }
    dat <- rbind(
        do.call(rbind,
                lapply(rl1,
                       function(m)
                           data.frame('k'=sum(rownames(m)%in% cells),'t'=sum(rownames(m) %in% prel.cells),'n'=n1))),
        do.call(rbind,
                lapply(rl2,
                       function(m)
                           data.frame('k'=sum(rownames(m)%in% cells),'t'=sum(rownames(m) %in% prel.cells),'n'=n2)))
    )
    mtest <- glm( cbind(k,t) ~ n, family=binomial(), data=dat)
    x <- dat$k / dat$t
    qres <- c(
        n1=mean(x[dat$n==n1],na.rm=T),
        n2=mean(x[dat$n==n2],na.rm=T),
        pval=summary(mtest)$coef[2,4]
    )
    names(qres) <- c(n1,n2,'pval')
    return(list(qres=qres,mtest=mtest))
}



splitFactor <- function(f) {
    lvls <- levels(f);
    names(lvls) <- lvls;
    lapply(lvls, function(l) {
        names(f)[f==l]
    })
}

compareTwo <- function(appnames1, appnames2, groupsFactor) {
    fs <- splitFactor(groupsFactor)
    nfs <- names(fs);
    names(nfs) <- nfs
    r <- as.data.frame(do.call(rbind,lapply(nfs, function(ln) {
        t.two.prop.comp(ccm.raw[appnames1], ccm.raw[appnames2],cells=fs[[ln]],prel.cells=names(f))$qres
    })))
    r$qval <- p.adjust(r$pval, 'fdr')
    r$sign <- r$qval < 0.05
    r
}


#' Compare the cells in a cluster that spans multiple samples
#' @param rl1 list of raw matrices 1
#' @param rl2 list of raw matrices 2
#' @param character vector listing the cells to compares
#' @param n1 name of group1
#' @param n2 names of group2
#' @param min.cells minimum number of cells in a samples to keep it
#' @return a DESseq2 result data.frame
#' @export t.two.exp.comp
t.two.exp.comp <- function(rl1,rl2,cells,n1='test',n2='control',min.cells=20) {
    tryCatch({
        require(DESeq2)
        require(Matrix)
        ## Subset apps to cells
        rl1 <- lapply(rl1,function(x) x[rownames(x) %in% cells,,drop=F])
        rl2 <- lapply(rl2,function(x) x[rownames(x) %in% cells,,drop=F])
        ## Check for empty samples
        ## Remove apps with < min.cells
        rl1 <- rl1[unlist(lapply(rl1,nrow))>min.cells]
        rl2 <- rl2[unlist(lapply(rl2,nrow))>min.cells]
        ## Don't run if there are not samples in one condition after filtering
        if (length(rl1)  == 0 | length(rl2) == 0) {
            return(NULL)
        }
        cts <- do.call(cbind,lapply(c(rl1,rl2),Matrix::colSums))
        ## Make metadata
        coldata <- data.frame(type=rep(c(n1,n2),c(length(rl1),length(rl2))));
        rownames(coldata) <- c(names(rl1),names(rl2))
        ## Make DESeq2 object
        dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,
                                              colData = coldata,design = ~ type)
        ## Do diff expression and get results
        dds <- DESeq(dds)
        res <- results(dds)
        res <- res[order(res$pvalue),]
        res <- subset(res, pvalue < 0.05)
        return(as.data.frame(res))
    }, error = function(e) {NULL})
}


#' Run comparisons within clusters spaning pagoda2 object boundaries
#' @param ccm named list of raw count matrices obtained with extractCountMatrices()
#' @param app.types a named factor of the ccm matrices with type for each
#' @param contrasts pairs of app.types to compares
#' @param clusters the clusters to to use (comparisons are only performed within clusters)
#' @param mc.cores number of cores to use
#' @return list of list of differential expression tables from DESeq2
runAllComparisons <- function(ccm, app.types, contrasts, clusters, mc.cores=16) {
    clusters <- clusters[!is.na(clusters)] 
    clusters.split <- factorBreakdown(clusters)
    lapply(contrasts, function(cc) {
        grp1 <- names(app.types)[app.types == cc[1]]
        grp2 <- names(app.types)[app.types == cc[2]]
        mclapply(clusters.split, function(cells.compare) {
            tryCatch({
                t.two.exp.comp(ccm[grp1],ccm[grp2],cells.compare)
            }, error =
                   function(e) { NULL }
            )
        },mc.cores =mc.cores)
    })
}

#' Filter a comparison list returned by runAll comparisons to
#' remove NULL objects and non-significant hits
#' @param comparisons a comparison, result of runAllComparisons
#' @param cutoff cutoff for qvalues (default 0.05)
#' @return a new comparison list
filterComparisons <- function(comparisons, cutoff = 0.05) {
    lapply(comparisons, function(cc) {
        cc <- cc[!unlist(lapply(cc,is.null))]
        lapply(cc, function(ccc) {
            ccc <- ccc[!is.na(ccc$padj),]
            ccc <- ccc[ccc$padj < cutoff,]
            ccc
        })
    })
}

removeNullComparisons <- function(comparisons, cutoff = 0.05) {
    lapply(comparisons, function(cc) {
        cc <- cc[!unlist(lapply(cc,is.null))]
    })
}



#' Compare the cells in a cluster that spans multiple samples
#' @param rl1 list of raw matrices 1
#' @param rl2 list of raw matrices 2
#' @param character vector listing the cells to compares
#' @param n1 name of group1
#' @param n2 names of group2
#' @param min.cells minimum number of cells in a samples to keep it
#' @return a DESseq2 result data.frame
#' @export t.two.exp.comp
t.two.exp.comp.2  <- function(rl1,rl2,cells,n1='test',n2='control',min.cells=20,pval.threshold=1e-2) {
  ##require(DESeq2)
  require(Matrix)
  ## Subset apps to cells
  rl1 <- lapply(rl1,function(x) x[rownames(x) %in% cells,,drop=F])
  rl2 <- lapply(rl2,function(x) x[rownames(x) %in% cells,,drop=F])
  ## Check for empty samples
  ## Remove apps with < min.cells
  rl1 <- rl1[unlist(lapply(rl1,nrow))>min.cells]
  rl2 <- rl2[unlist(lapply(rl2,nrow))>min.cells]
  ## Don't run if there are not samples in one condition after filtering
  if (length(rl1)  == 0 | length(rl2) == 0) {
    return(NULL)
  }
  cts <- do.call(cbind,lapply(c(rl1,rl2),Matrix::colSums))
  ## Make metadata
  coldata <- data.frame(type=rep(c(n1,n2),c(length(rl1),length(rl2))));
  rownames(coldata) <- c(names(rl1),names(rl2))
  ## Make DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,colData = coldata,design = ~ type)
  ## Do diff expression and get results
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  allgenes <- rownames(res);
  res <- as.data.frame(res);
  res <- res[!is.na(res$pvalue),]
  res <- res[res$pvalue<=pval.threshold,]
  res <- res[order(res$pvalue),]
  res$Z <- qnorm(res$pvalue)*sign(res$log2FoldChange)
  res$Za <- qnorm(res$padj)*sign(res$log2FoldChange)
  res <- cbind(gene=rownames(res),res)
  if(nrow(res)>0) { # calculate mean expression values for the plots
    cts <- t(t(cts[as.character(res$gene),])/colSums(cts)*1e6)
    ilev <- tapply(1:ncol(cts),coldata$type,function(ii) {
      cts[,ii,drop=F]
    })
    #mdf <- do.call(cbind,lapply(x$ilev,rowMeans));
    #colnames(mdf) <- paste(colnames(mdf),'mean',sep='_')
    #res <- cbind(res,mdf)
  } else {
    ilev <- NULL;
  }
  #res <- subset(res, pvalue < 0.05)
  return(list(res=res,genes=allgenes,ilev=ilev,snames=c(n1,n2)))
}



#' Compare the cells in a cluster that spans multiple samples
#' @param rl1 list of raw matrices 1
#' @param rl2 list of raw matrices 2
#' @param character vector listing the cells to compares
#' @param n1 name of group1
#' @param n2 names of group2
#' @param min.cells minimum number of cells in a samples to keep it
#' @param pval.threshold p-value threshold
#' @return a DESseq2 result data.frame
#' @export t.two.exp.comp
t.two.exp.comp.2  <- function(rl1,rl2,cells,n1='test',n2='control',min.cells=20,pval.threshold=1e-2) {
  ##require(DESeq2)
  require(Matrix)
  ## Subset apps to cells
  rl1 <- lapply(rl1,function(x) x[rownames(x) %in% cells,,drop=F])
  rl2 <- lapply(rl2,function(x) x[rownames(x) %in% cells,,drop=F])
  ## Check for empty samples
  ## Remove apps with < min.cells
  rl1 <- rl1[unlist(lapply(rl1,nrow))>min.cells]
  rl2 <- rl2[unlist(lapply(rl2,nrow))>min.cells]
  ## Don't run if there are not samples in one condition after filtering
  if (length(rl1)  == 0 | length(rl2) == 0) {
    return(NULL)
  }
  cts <- do.call(cbind,lapply(c(rl1,rl2),Matrix::colSums))
  ## Make metadata
  coldata <- data.frame(type=rep(c(n1,n2),c(length(rl1),length(rl2))));
  rownames(coldata) <- c(names(rl1),names(rl2))
  ## Make DESeq2 object
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts,colData = coldata,design = ~ type)
  ## Do diff expression and get results
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  allgenes <- rownames(res);
  res <- as.data.frame(res);
  res <- res[!is.na(res$pvalue),]
  res <- res[res$pvalue<=pval.threshold,]
  res <- res[order(res$pvalue),]
  res$Z <- qnorm(res$pvalue)*sign(res$log2FoldChange)
  res$Za <- qnorm(res$padj)*sign(res$log2FoldChange)
  res <- cbind(gene=rownames(res),res)
  if(nrow(res)>0) { # calculate mean expression values for the plots
    cts <- t(t(cts[as.character(res$gene),])/colSums(cts)*1e6)
    ilev <- tapply(1:ncol(cts),coldata$type,function(ii) {
      cts[,ii,drop=F]
    })
    #mdf <- do.call(cbind,lapply(x$ilev,rowMeans));
    #colnames(mdf) <- paste(colnames(mdf),'mean',sep='_')
    #res <- cbind(res,mdf)
  } else {
    ilev <- NULL;
  }
  #res <- subset(res, pvalue < 0.05)
  return(list(res=res,genes=allgenes,ilev=ilev,snames=c(n1,n2)))
}


#' Run comparisons within clusters spaning pagoda2 object boundaries
#' @param ccm named list of raw count matrices obtained with extractCountMatrices()
#' @param app.types a named factor of the ccm matrices with type for each
#' @param contrasts pairs of app.types to compares
#' @param clusters the clusters to to use (comparisons are only performed within clusters)
#' @param mc.cores number of cores to use
#' @return list of list of differential expression tables from DESeq2
runAllComparisons.2 <- function(ccm, app.types, contrasts, clusters, mc.cores=16) {
    require('nbHelpers')
    clusters <- clusters[!is.na(clusters)] 
    clusters.split <- factorBreakdown(clusters)
    lapply(contrasts, function(cc) {
        grp1 <- names(app.types)[app.types == cc[1]]
        grp2 <- names(app.types)[app.types == cc[2]]
        mclapply(clusters.split, function(cells.compare) {
            tryCatch({
                t.two.exp.comp.2(ccm[grp1],ccm[grp2],cells.compare)
            }, error =
                   function(e) { NULL }
            )
        },mc.cores =mc.cores)
    })
}


#' Convert de result to json and optionally save it as a file
#' @param xe differential expression result from t.two.exp.comp.2
#' @param filename optional file name to save to
#' @return json test string (invisible)
#' @export de2json
de2json <- function(xe, filename=NULL) {
    xe$res <- unname(split(xe$res,1:nrow(xe$res)))
    y <- rjson::toJSON(xe)
    if (!is.null(filename)) {
        write(y,file=filename)
    }
    invisible(y)
}

