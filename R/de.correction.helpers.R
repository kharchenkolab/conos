
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

trimmedMean <- function(x) {
    x <- na.rm(x)
    min.val <- median(x) - 2 * IQR(x)
    max.val <- median(x) + 2 * IQR(x)
    x[x < min.val | x > max.val] <- NA
    x <-na.rm(x)
    mean(x)
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

