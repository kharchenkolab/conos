
#' Extract count matrices for differential expression
#' and proportion comparisons from pagoda2 apps
#' @param apps a list of pagoda2 apps
#' @return a list of dgCMatrices
extractCountMatrices <- function(apps) {
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
