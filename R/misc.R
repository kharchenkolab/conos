
#' Fix a pagoda2 selection object prefixed
#' @export fixSelectionPrefix
fixSelectionPrefix <- function(p2selection, map) {
  nsel <- lapply(p2selection, function(cs) {
    spid <- strsplit(cs$cells,'_')
    pref <- sapply(spid,'[',1)
    cellbarcode <- sapply(spid,'[',2)
    pref.uniq <- unique(pref)
    if(any(!pref.uniq %in% names(map))) {
      warning(paste0('The following samples could not be matched: ',
                     paste0(pref.uniq[!pref.uniq %in% names(map)], collapse=' ')))
    }
    cs$cells <- paste0(map[pref],'.',cellbarcode)
    cs
  })
  nsel
}

#' Calculate zlim from a vector of numeric values
#' @param vs numeric values
#' @return zlim range
calcZlim <- function(vs, gradient.range.quantile = 0.95) {
  zlim <- as.numeric(quantile(vs,p=c(1-gradient.range.quantile,gradient.range.quantile)))
  if(diff(zlim)==0) {
    zlim <- as.numeric(range(vs))
  }
  zlim
}
