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

  if (any(sapply(probs, function(x) sum(is.na(x)))))
    stop("NAs in transition probabilities")

  return(list(idx=adj.list, probabilities=probs, names=edge.list.fact$levels))
}

#' @export
embedKnnGraph <- function(commute.times, n.neighbors, names=NULL, verbose=TRUE, ...) {
  min.n.neighbors <- sapply(commute.times$idx, length) %>% min()
  if (min.n.neighbors < n.neighbors) {
    n.neighbors <- min.n.neighbors
    warning("Maximal number of estimated neighbors is ", min.n.neighbors)
  }

  ct.top <- sapply(commute.times$dist, `[`, 1:n.neighbors) %>% t() + 1
  ct.top.ids <- sapply(commute.times$idx, `[`, 1:n.neighbors) %>% t() + 1

  ct.top.ids <- cbind(1:nrow(ct.top.ids), ct.top.ids)
  ct.top <- cbind(rep(0, nrow(ct.top)), ct.top)

  umap <- uwot::umap(data.frame(x=rep(0, nrow(ct.top))), nn_method=list(idx=ct.top.ids, dist=ct.top), verbose=verbose, ...)
  rownames(umap) <- names
  return(umap)
}

embedGraphUmap <- function(graph, verbose=T, min.prob=1e-3, min.visited.verts=1000, n.cores=1,
                           max.hitting.nn.num=0, max.commute.nn.num=0, min.prob.lower=1e-7,
                           n.neighbors=40, n.epochs=1000, spread=15, min.dist=0.001, return.all=F,
                           ...) {
  min.visited.verts = min(min.visited.verts, length(igraph::V(graph) - 1))
  if (max.hitting.nn.num == 0) {
    max.hitting.nn.num <- length(igraph::V(graph)) - 1
  }

  if (verbose) cat("Convert graph to adjacency list...\n")
  adj.info <- graphToAdjList(graph);
  if (verbose) cat("Done\n")

  if (verbose) cat("Estimate nearest neighbors and commute times...\n")
  commute.times <- get_nearest_neighbors(adj.info$idx, adj.info$probabilities, min_prob=min.prob,
                                         min_visited_verts=min.visited.verts, n_cores=n.cores, max_hitting_nn_num=max.hitting.nn.num,
                                         max_commute_nn_num=max.commute.nn.num, min_prob_lower=min.prob.lower, verbose=verbose)
  if (verbose) cat("Done\n")

  if (verbose) cat("Estimate UMAP embedding...\n")
  umap <- embedKnnGraph(commute.times, n.neighbors=n.neighbors, names=adj.info$names, n_threads=n.cores,
                        n_epochs=n.epochs, spread=spread, min_dist=min.dist, verbose=verbose, ...)
  if (verbose) cat("Done\n")

  if (return.all)
    return(list(adj.info=adj.info, commute.times=commute.times, umap=umap))

  return(umap)
}