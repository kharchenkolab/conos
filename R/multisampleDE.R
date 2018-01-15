#' Get raw count matrices
#' @export getCountMatricesRaw
getCountMatricesRaw <- function(r.n, common.genes) {
  ccm.raw <- mclapply(r.n,function(r) {
    om <- as(r$misc$rawCounts,'dgTMatrix')
    om <- om[,colnames(om) %in% common.genes]
    mi <- match(colnames(om),common.genes)
    x <- new("dgTMatrix",i = om@i,j = as.integer(mi[om@j+1]-1),x=om@x,Dim=c(nrow(om),length(common.genes)))
    rownames(x) <- rownames(om); colnames(x) <- common.genes
    as(x,'dgCMatrix')
  },mc.cores=30)
  ccm.raw
}

