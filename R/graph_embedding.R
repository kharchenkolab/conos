splitVectorByNodes <- function(vec, nodes, n.nodes) {
  res <- lapply(1:n.nodes, function(x) list())
  splitted <- split(vec, nodes)
  res[as.integer(names(splitted))] <- splitted
  return(res)
}

graphToAdjList <- function(graph) {
  edge.list.fact <- igraph::as_edgelist(graph) %>% as_factor()
  edge.list <- matrix(edge.list.fact$values, ncol=2)
  n.nodes <- length(igraph::V(graph))
  adj.list <- mapply(c, splitVectorByNodes(edge.list[,1], edge.list[,2], n.nodes),
                     splitVectorByNodes(edge.list[,2], edge.list[,1], n.nodes)) %>%
    lapply(unlist) %>% lapply(`-`, 1)

  probs <- mapply(c, splitVectorByNodes(igraph::E(graph)$weight, edge.list[,2], n.nodes),
                  splitVectorByNodes(igraph::E(graph)$weight, edge.list[,1], n.nodes)) %>%
    lapply(unlist) %>%
    lapply(function(x) x / sum(x))

  return(list(idx=adj.list, probabilities=probs, names=edge.list.fact$levels))
}

#' @export
embedKnnGraph <- function(commute.times, n.neighbors, names=NULL, verbose=TRUE, ...) {
  min.n.neighbors <- sapply(commute.times$idx, length) %>% min()
  if (min.n.neighbors < n.neighbors) {
    n.neighbors <- min.n.neighbors
    warning("Maximal number of estimated neighbors is", min.n.neighbors)
  }

  ct.top <- sapply(commute.times$dist, `[`, 1:n.neighbors) %>% t() + 1
  ct.top.ids <- sapply(commute.times$idx, `[`, 1:n.neighbors) %>% t() + 1

  ct.top.ids <- cbind(1:nrow(ct.top.ids), ct.top.ids)
  ct.top <- cbind(rep(0, nrow(ct.top)), ct.top)

  umap <- uwot::umap(data.frame(x=rep(0, nrow(ct.top))), nn_method=list(idx=ct.top.ids, dist=ct.top), verbose=verbose, ...)
  rownames(umap) <- names
  return(umap)
}