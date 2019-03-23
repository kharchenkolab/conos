fillNa <- function(x, value=0) ifelse(is.na(x), value, x)

#' @export
asAdjList <- function(graph) {
  adj.mtx <- graph %>% igraph::as_adjacency_matrix(attr="weight") %>% as("dgTMatrix")
  weights.per.cell <- setNames(adj.mtx@x, names(igraph::V(graph)[adj.mtx@i + 1])) %>%
    split(adj.mtx@j) %>%  setNames(colnames(adj.mtx)[as.integer(names(.)) + 1])

  return(weights.per.cell)
}

getFactorPairPerEdge <- function(edge.df, factor.per.vertex, sep=".vs.") {
  edge.df %$%
    paste(pmin(factor.per.vertex[v1], factor.per.vertex[v2]),
          pmax(factor.per.vertex[v1], factor.per.vertex[v2]),
          sep=sep)
}

#' @export
toWeightedEdgeList <- function(g) {
  df <- igraph::as_edgelist(g) %>% magrittr::set_colnames(c("v1", "v2")) %>%
    tibble::as_tibble() %>%
    cbind(weight=igraph::E(g)$weight, type=igraph::E(g)$type)

  return(df)
}

#' @export
filterHubbyEdges <- function(knn.edge.df, dataset.per.vertex, min.nn.per.vertex, verbose=T) {
  lapply.fun <- if (verbose && requireNamespace("pbapply", quietly=F)) pbapply::pblapply else lapply
  dataset.pairs <- getFactorPairPerEdge(knn.edge.df, dataset.per.vertex)
  res <- split(knn.edge.df, dataset.pairs) %>%
    lapply.fun(function(df)
      filter(df, !conos:::edge_removal_mask(v1, v2, weight, min.nn.per.vertex, verbose=F))) %>%
    plyr::rbind.fill()

  return(res)
}

#' @export
graphFromWeightedEdgeList <- function(edge.df) {
  g <- edge.df %>% dplyr::select(v1, v2) %>% as.matrix() %>% igraph::graph_from_edgelist(directed=F)
  igraph::E(g)$weight <- edge.df$weight
  igraph::E(g)$type <- edge.df$type

  return(g)
}

#' @export
adjustWeigtsByFactor <- function(edge.df, factor.per.vertex, n.iters=5, verbose=F) {
  factor.per.edge <- getFactorPairPerEdge(edge.df, factor.per.vertex)
  for (it in 1:n.iters) {
    w.sum.per.vert.per.pair <- split(edge.df, factor.per.edge) %>%
      lapply(function(df) df %$% rbind(tibble::tibble(v=v1, weight), tibble::tibble(v=v2, weight))) %>%
      lapply(function(df) split(df$weight, df$v) %>% sapply(sum) %>% tibble::tibble(ws=., cb=names(.)))

    w.sum.per.vert.per.pair <- mapply(function(df, n) cbind(df, pair=n), w.sum.per.vert.per.pair,
                                      names(w.sum.per.vert.per.pair), SIMPLIFY=F) %>%
      plyr::rbind.fill()

    w.sum.per.vert <- w.sum.per.vert.per.pair %$% split(ws, cb) %>% sapply(sum)
    w.sum.per.vert.per.pair %<>% dplyr::mutate(frac=ws / w.sum.per.vert[cb])

    w.div.per.factor <- w.sum.per.vert.per.pair %$% split(frac, pair) %>% sapply(mean) %>% `/`(mean(.))

    edge.df$weight %<>% `/`(w.div.per.factor[factor.per.edge])

    if (verbose) {
      w.range <- w.sum.per.vert.per.pair %$% split(frac, pair) %>% sapply(mean) %>% range()
      cat("Iteration", it, ". Weight fraction range:", diff(w.range), "\n")
    }
  }

  return(edge.df)
}

#' @export
getClusterNonMixingScores <- function(clusters, annotation, mixing.factor, return.by.cluster=F) {
  cluster.annotation.table <- table(clusters, annotation[names(clusters)]) %>% as.matrix()
  annotation.factor.table <- table(annotation, mixing.factor) %>% as.matrix()

  cluster.factor.table.expected <- (cluster.annotation.table %*% annotation.factor.table) %>%
    `/`(rowSums(.))

  cluster.factor.table.observed <- table(clusters, mixing.factor[names(clusters)]) %>%
    as.matrix() %>% `/`(rowSums(.))

  nonmixing.scores <- cluster.factor.table.observed  %>%
    `-`(cluster.factor.table.expected[rownames(.), colnames(.)]) %>% abs() %>%
    rowSums()

  if (return.by.cluster)
    return(nonmixing.scores)

  return(setNames(nonmixing.scores[clusters], names(clusters)))
}

embeddingNonMixingScores <- function(embedding, factor.per.cell, annotation=NULL, k=40, metric="L2", verbose=F) {
  if (is.null(annotation)) {
    annotation <- rep("x", nrow(embedding)) %>% setNames(rownames(embedding))
  }

  nn.largevis <- conos::n2Knn(embedding, k=k+1, indexType=metric, verbose=verbose) %>%
    drop0() %>% as("dgTMatrix")

  adj.list <- split(nn.largevis@i + 1, nn.largevis@j + 1) %>%
    setNames(colnames(embedding))

  factor.per.cell %<>% .[rownames(embedding)]
  annotation %<>% .[rownames(embedding)]
  all.factors <- unique(factor.per.cell)

  nn.type.frac.per.cell <- adj.list %>%
    sapply(function(x) fillNa(table(factor.per.cell[x])[all.factors])) %>% t() %>%
    magrittr::set_colnames(all.factors) %>% magrittr::set_rownames(rownames(embedding)) %>%
    `/`(rowSums(.))

  expected.frac.per.cell <- table(annotation, factor.per.cell) %>%
    `/`(rowSums(.)) %>%
    .[annotation,] %>% as.matrix() %>%
    magrittr::set_rownames(rownames(embedding)) %>%
    .[, all.factors]

  return(rowSums(abs(nn.type.frac.per.cell - expected.frac.per.cell)))
}

adjustWeightsByFactorSum <- function(graph, factor.per.vert, within.factor.weight.mult=1.0, balance.factors=T) {
  edge.list <- igraph::as_edgelist(graph)
  same.factor <- (factor.per.vert[edge.list[,1]] == factor.per.vert[edge.list[,2]])
  edge.weights <- igraph::E(graph)$weight
  edge.weights[same.factor] %<>% `*`(within.factor.weight.mult)

  if (!balance.factors)
    return(edge.weights)

  factor.per.edge <- tibble::as_tibble(edge.list) %>% set_colnames(c("v1", "v2")) %>%
    getFactorPairPerEdge(factor.per.vert)
  weights.per.pair <- split(edge.weights[!same.factor], factor.per.edge[!same.factor])

  weight.div.per.pair <- sapply(weights.per.pair, sum) %>% `/`(mean(.))

  edge.weights[!same.factor] %<>% `/`(weight.div.per.pair[factor.per.edge[!same.factor]])
  return(edge.weights)
}