
#' plot results of de.sign.enrich
#' @param de.enrich.res output of de.sign.enrich
#' @return ggplot2 plot object
plotDEenrichment <- function(de.enrich.res) {
    y2 <- melt(de.enrich.res)
    y2 <- y2[!is.na(y2$value),]
    y2$value[is.infinite(y2$value)] <- 1e-16
    y2$value[y2$value < 1e-16] <- 1e-16
    #y2$value <- -log10(y2$value+1e-16)
    y2$value <- y2$value < 1e-3
    ## Plot as a sparse heatmap
    ggplot(y2, aes(x=L1,y=L2,fill=value)) + geom_tile() +
#        scale_fill_gradient(low = "white", high = "steelblue") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> p1
    p1

}

#' look for enrichment of a particular signature in the results
#' @param de.results de results from Pagoda2ensemble$getDEresults()
#' @param signature character vector of genes in signature
#' @param direction up down or both
#' @param datastructure with p-values from chi square test
de.sign.enrich <- function(de.results, signature, direction = 'up') {
    if (!direction %in% c('up','down','both'))
        error('direction must be up or down');
    lapply(namedNames(de.results), function(j) {
        lapply(de.results[[j]], function(x){
            try({
                if(!is.null(x)) {
                    sign <- signature[signature %in% rownames(x$res)]
                    insig = rownames(x$res) %in% signature
                    if (direction == 'up') {
                        tbl <- table(data.frame(
                            insign = factor(as.character(insig),levels=c('TRUE','FALSE')),
                            de.list = factor(as.character(x$res$significant  & x$res$Z > 0),
                                             levels=c('TRUE','FALSE'))
                        ))
                    } else if(direction =='down') {
                        tbl <- table(data.frame(
                            insign = factor(as.character(insig),levels=c('TRUE','FALSE')),
                            de.list = factor(as.character(x$res$significant  & x$res$Z < 0),
                                             levels=c('TRUE','FALSE'))
                        ))                        
                    } else {
                        tbl <- table(data.frame(
                            insign = factor(as.character(insig),levels=c('TRUE','FALSE')),
                            de.list = factor(as.character(x$res$significant), levels=c('TRUE','FALSE'))
                        ))
                    }
                    chisq.test(tbl)$p.val
                } else {
                    NA
                }
            })
        })
    })

}
