#' Access PCA from sample
#' 
#' @param sample sample from which to access PCA
#' @return matrix of PCA cell embeddings (cells x components); NULL if the sample has no PCA reduction
#' @rdname getPca
#' @export
setGeneric("getPca", function(sample) standardGeneric("getPca"))

#' @rdname getPca
setMethod("getPca", signature("Pagoda2"), function(sample) {
  ## pagoda2.1 stores reductions name-keyed; prefer the object's default reduction rather than a
  ## hardcoded "PCA" (works for renamed reductions), falling back to "PCA" / a sole reduction.
  red <- sample$reductions
  if (is.null(red) || length(red) == 0L) return(NULL)
  key <- tryCatch(sample$defaults$reduction, error = function(e) NULL)
  if (!is.null(key) && !is.null(red[[key]])) return(red[[key]])
  if (!is.null(red[["PCA"]])) return(red[["PCA"]])
  if (length(red) == 1L) return(red[[1]])
  NULL
})

#' @rdname getPca
setMethod("getPca", signature("seurat"), function(sample) sample@dr$pca@cell.embeddings)

#' @rdname getPca
setMethod(
  f = 'getPca',
  signature = signature('Seurat'),
  definition = function(sample) {
    checkSeuratV3()
    return(Seurat::Embeddings(object = sample))
  }
)

#' @rdname getPca
setMethod("getPca", signature("Conos"), function(sample) NULL) # a Conos has no single joint PCA; use getEmbedding()

.conos_pagoda2_has_method <- function(sample, name) {
  is.function(tryCatch(sample[[name]], error = function(e) NULL))
}

.conos_get_pagoda2_expression <- function(sample, genes = NULL, transposed = FALSE) {
  if (.conos_pagoda2_has_method(sample, "getExpressionBlock")) {
    orientation <- if (transposed) "cell_by_gene" else "gene_by_cell"
    return(sample$getExpressionBlock(genes = genes, orientation = orientation))
  }
  x <- sample$counts
  if (!is.null(genes)) {
    x <- x[, genes, drop = FALSE]
  }
  if (transposed) {
    return(x)
  }
  Matrix::t(x)
}

.conos_get_pagoda2_raw_counts <- function(sample, genes = NULL, cells = NULL, transposed = FALSE) {
  if (.conos_pagoda2_has_method(sample, "getRawCounts")) {
    orientation <- if (transposed) "cell_by_gene" else "gene_by_cell"
    return(sample$getRawCounts(cells = cells, genes = genes, orientation = orientation))
  }
  x <- sample$misc$rawCounts
  if (!is.null(cells)) {
    x <- x[cells, , drop = FALSE]
  }
  if (!is.null(genes)) {
    x <- x[, genes, drop = FALSE]
  }
  if (transposed) {
    return(x)
  }
  Matrix::t(x)
}

.conos_get_pagoda2_cell_names <- function(sample) {
  if (.conos_pagoda2_has_method(sample, "getRawCounts")) {
    return(rownames(sample$getRawCounts()))
  }
  rownames(sample$counts)
}

.conos_get_pagoda2_gene_names <- function(sample) {
  if (.conos_pagoda2_has_method(sample, "getRawCounts")) {
    return(colnames(sample$getRawCounts()))
  }
  colnames(sample$counts)
}


#' Access overdispersed genes from sample
#' 
#' @param sample sample from which to overdispereed genes
#' @param n.odgenes numeric Number of overdisperesed genes to get
#' @return character vector of overdispersed gene names
#' @rdname getOverdispersedGenes
#' @export
setGeneric("getOverdispersedGenes", function(sample, n.odgenes=1000) standardGeneric("getOverdispersedGenes"))

#' @rdname getOverdispersedGenes
setMethod("getOverdispersedGenes", signature("Pagoda2"), function(sample, n.odgenes=NULL) sample$getOdGenes(n.odgenes))

#' @rdname getOverdispersedGenes
setMethod("getOverdispersedGenes", signature("seurat"), function(sample, n.odgenes=NULL)
  if (is.null(NULL)) sample@var.genes else head(rownames(sample@hvg.info), n.odgenes))

#' @importFrom rlang %||%
#' @rdname getOverdispersedGenes
setMethod(
  f = 'getOverdispersedGenes',
  signature = signature('Seurat'),
  definition = function(sample, n.odgenes = NULL) {
    checkSeuratV3()
    vf <- Seurat::VariableFeatures(object = sample) %||% rownames(x = sample)
    n.odgenes <- n.odgenes %||% length(x = vf)
    return(head(x = vf, n = n.odgenes))
  }
)

#' @rdname getOverdispersedGenes
setMethod("getOverdispersedGenes", signature("Conos"), function(sample, n.odgenes=NULL) commonOverdispersedGenes(sample$samples,n.odgenes, verbose=FALSE))


#' List the molecular modalities (facets / assays) available in a sample
#'
#' @param sample a sample (pagoda2, Seurat, ...)
#' @return character vector of modality names (pagoda2.1 facets, Seurat assays; a single "RNA" for
#'   single-modality / legacy objects)
#' @rdname getModalities
#' @export
setGeneric("getModalities", function(sample) standardGeneric("getModalities"))

#' @rdname getModalities
setMethod("getModalities", signature("Pagoda2"), function(sample) {
  if (.conos_pagoda2_has_method(sample, "listFacets")) return(as.character(sample$listFacets()))
  "RNA" # legacy single-modality pagoda2
})

#' @rdname getModalities
setMethod("getModalities", signature("seurat"), function(sample) "RNA")

#' @rdname getModalities
setMethod("getModalities", signature("Seurat"), function(sample) { checkSeuratV3(); as.character(Seurat::Assays(sample)) })


#' The default molecular modality of a sample (the one integration uses unless told otherwise)
#'
#' @param sample a sample (pagoda2, Seurat, ...)
#' @return character scalar naming the sample's default modality (pagoda2.1 `defaultFacet`, Seurat
#'   `DefaultAssay`; "RNA" for single-modality / legacy objects)
#' @rdname getDefaultModality
#' @export
setGeneric("getDefaultModality", function(sample) standardGeneric("getDefaultModality"))

#' @rdname getDefaultModality
setMethod("getDefaultModality", signature("Pagoda2"), function(sample) {
  df <- tryCatch(sample$defaultFacet, error = function(e) NULL)
  if (!is.null(df) && nzchar(df)) return(as.character(df))
  getModalities(sample)[1]
})

#' @rdname getDefaultModality
setMethod("getDefaultModality", signature("seurat"), function(sample) "RNA")

#' @rdname getDefaultModality
setMethod("getDefaultModality", signature("Seurat"), function(sample) { checkSeuratV3(); as.character(Seurat::DefaultAssay(sample)) })

## Feature names for a given modality, cheaply (no count materialization where avoidable): pagoda2.1 facet
## featureMeta rownames; Seurat assay rownames; else the default-modality genes (legacy pagoda2 / fallback).
.conos_modality_features <- function(sample, modality) {
  if (inherits(sample, "Seurat")) { checkSeuratV3(); return(rownames(sample[[modality]])) }
  if (inherits(sample, "Pagoda2") && .conos_pagoda2_has_method(sample, "getFacet")) {
    fm <- tryCatch(sample$getFacet(modality)$featureMeta, error = function(e) NULL)
    if (!is.null(fm)) return(rownames(fm))
  }
  getGenes(sample)
}


#' Access cell names from sample
#' 
#' @param sample sample from which to cell names
#' @return character vector of cell names
#' @rdname getCellNames
#' @export
setGeneric("getCellNames", function(sample) standardGeneric("getCellNames"))

#' @rdname getCellNames
setMethod("getCellNames", signature("Pagoda2"), function(sample) .conos_get_pagoda2_cell_names(sample))

#' @rdname getCellNames
setMethod("getCellNames", signature("seurat"), function(sample) colnames(sample@data))

#' @rdname getCellNames
setMethod(f = 'getCellNames', signature = signature('Seurat'), definition = function(sample) return(colnames(x = sample)))

#' @rdname getCellNames
setMethod("getCellNames", signature("Conos"), function(sample) unlist(lapply(sample$samples,getCellNames)))


#' Access genes from sample
#' 
#' @param sample sample from which to get genes
#' @return character vector of gene (feature) names
#' @rdname getGenes
#' @export
setGeneric("getGenes", function(sample) standardGeneric("getGenes"))

#' @rdname getGenes
setMethod("getGenes", signature("Pagoda2"), function(sample) .conos_get_pagoda2_gene_names(sample))

#' @rdname getGenes
setMethod("getGenes", signature("seurat"), function(sample) rownames(sample@data))

#' @rdname getGenes
setMethod(f = 'getGenes', signature = signature('Seurat'), definition = function(sample) return(rownames(x = sample)))

#' @rdname getGenes
setMethod("getGenes", signature("Conos"), function(sample) unique(unlist(lapply(sample$samples, getGenes))))


#' Set edge matrix edgeMat with certain values on sample
#' 
#' @param sample sample from which to set edge matrix edgeMat with certain values
#' @param value values to set with edgeMat<-
#' @return the sample, modified in place with the edge matrix stored (invisibly, per the replacement-function convention)
#' @rdname edgeMat
#' @export
setGeneric("edgeMat<-", function(sample, value) standardGeneric("edgeMat<-"))

#' @rdname edgeMat
setMethod("edgeMat<-", signature("Pagoda2"), function(sample, value) {sample$misc$edgeMat <- value; sample})

#' @rdname edgeMat
setMethod("edgeMat<-", signature("seurat"), function(sample, value) {sample@misc$edgeMat <- value; sample})

#' @rdname edgeMat
setMethod(
  f = 'edgeMat<-',
  signature = signature('Seurat'),
  definition = function(sample, value) {
    checkSeuratV3()
    Seurat::Misc(object = sample, slot = 'edgeMat') <- value
    return(sample)
  }
)



#' Access edgeMat from sample
#' 
#' @param sample sample from which to access edge matrix edgeMat
#' @return the edge matrix (`edgeMat`) previously stored on the sample, or NULL if none
#' @rdname edgeMat
#' @export
setGeneric("edgeMat", function(sample) standardGeneric("edgeMat"))

#' @rdname edgeMat
setMethod("edgeMat", signature("Pagoda2"), function(sample) sample$misc$edgeMat)

#' @rdname edgeMat
setMethod("edgeMat", signature("seurat"), function(sample) sample@misc$edgeMat)

#' @rdname edgeMat
setMethod(
  f = 'edgeMat',
  signature = signature('Seurat'),
  definition = function(sample) {
    checkSeuratV3()
    return(Seurat::Misc(object = sample, slot = 'edgeMat'))
  }
)


#' Access count matrix from sample
#' 
#' @param sample sample from which to get the count matrix
#' @param transposed boolean Whether the count matrix should be transposed (default=FALSE)
#' @return the normalized expression matrix (genes x cells, or cells x genes when `transposed=TRUE`)
#' @rdname getCountMatrix
#' @export
setGeneric("getCountMatrix", function(sample, transposed=FALSE) standardGeneric("getCountMatrix"))

#' @rdname getCountMatrix
setMethod("getCountMatrix", signature("Pagoda2"), function(sample, transposed=FALSE) .conos_get_pagoda2_expression(sample, transposed = transposed))

#' @rdname getCountMatrix
setMethod("getCountMatrix", signature("seurat"), function(sample, transposed=FALSE) {
  cm <- if (is.null(sample@scale.data)) sample@data else sample@scale.data
  if (transposed)
    return(Matrix::t(cm))

  return(cm)
})

#' @rdname getCountMatrix
setMethod('getCountMatrix', signature('Seurat'), function(sample, transposed=FALSE) {
    checkSeuratV3()
    ## the log-normalized `data` layer (non-negative, sparse) is the right "normalized expression" for
    ## conos' uses (joint count matrix, dot plots, markers). The z-scored `scale.data` is a PCA input, not
    ## an expression matrix -- its negative values break log fold changes and the sparse joint merge.
    dat <- getSeuratAssayData(sample, 'data')
    dims <- dim(x = dat)
    dat.na <- all(dims == 1) && all(is.na(x = dat))
    if (all(dims == 0) || dat.na) {
      dat <- getSeuratAssayData(sample, 'counts')
    }

    if (transposed)
      return(Matrix::t(dat))

    return(dat)
  }
)

#' @rdname getCountMatrix
setMethod("getCountMatrix", signature("Conos"), function(sample, transposed=FALSE) { m <- sample$getJointCountMatrix(raw=FALSE); if (transposed) Matrix::t(m) else m })



#' Access gene expression from sample
#' 
#' @param sample sample from which to access gene expression
#' @param gene character vector Genes to access
#' @return named numeric vector of the requested gene's expression across cells (NA for cells/samples lacking the gene)
#' @rdname getGeneExpression
#' @export
setGeneric("getGeneExpression", function(sample, gene) standardGeneric("getGeneExpression"))

#' @rdname getGeneExpression
setMethod("getGeneExpression", signature("Pagoda2"), function(sample, gene) {
  if (gene %in% .conos_get_pagoda2_gene_names(sample)) {
    x <- .conos_get_pagoda2_expression(sample, genes = gene, transposed = TRUE)
    return(stats::setNames(as.numeric(x[, gene]), rownames(x)))
  }

  cells <- .conos_get_pagoda2_cell_names(sample)
  return(stats::setNames(rep(NA, length(cells)), cells))
})

#' @rdname getGeneExpression
setMethod("getGeneExpression", signature("Conos"), function(sample, gene) {
  lapply(sample$samples, getGeneExpression, gene) %>% Reduce(c, .)
})

getGeneExpression.default <- function(sample, gene) {
  count.matrix <- getCountMatrix(sample)
  if(gene %in% rownames(count.matrix)) {
    return(count.matrix[gene,])
  }

  return(stats::setNames(rep(NA, ncol(count.matrix)), colnames(count.matrix)))
}

#' @rdname getGeneExpression
setMethod("getGeneExpression", signature("Seurat"), function(sample, gene) {
  checkSeuratV3()
  ## https://satijalab.org/seurat/essential_commands.html
  if (gene %in% rownames(Seurat::GetAssayData(object = sample))){
    ## rownames(data) are gene names
    return(Seurat::GetAssayData(object = sample)[gene, ])
  }

  return(stats::setNames(rep(NA, ncol(Seurat::GetAssayData(object = sample))), colnames(Seurat::GetAssayData(object = sample)))) 
})

#' @rdname getGeneExpression
setMethod("getGeneExpression", signature("seurat"), function(sample, gene) {
  ## https://satijalab.org/seurat/essential_commands.html
  if (gene %in% rownames(sample@data)){
    ## rownames(data) are gene names
    return(sample@data[gene, ])
  }

  return(stats::setNames(rep(NA, ncol(sample@data)), colnames(sample@data))) 
})


#' Access raw count matrix from sample
#' 
#' @param sample sample from which to get the raw count matrix
#' @param transposed boolean Whether the raw count matrix should be transposed (default=FALSE)
#' @return the raw count matrix (genes x cells, or cells x genes when `transposed=TRUE`)
#' @rdname getRawCountMatrix
#' @export
setGeneric("getRawCountMatrix", function(sample, transposed=FALSE) standardGeneric("getRawCountMatrix"))

#' @rdname getRawCountMatrix
setMethod("getRawCountMatrix", signature("Pagoda2"), function(sample, transposed=FALSE) .conos_get_pagoda2_raw_counts(sample, transposed = transposed))

#' @rdname getRawCountMatrix
setMethod(
  f = "getRawCountMatrix",
  signature = signature("seurat"),
  definition = function(sample, transposed=FALSE) {
    mi <- match(x = sample@cell.names, table = colnames(sample@raw.data))
    x <- sample@raw.data[, mi, drop = FALSE]
    if (transposed) {
      return(t(x = x))
    } else {
      return(x)
    }
  }
)

#' @rdname getRawCountMatrix
setMethod(
  f = 'getRawCountMatrix',
  signature = signature('Seurat'),
  definition = function(sample, transposed=FALSE) {
    checkSeuratV3()
    rd <- getSeuratAssayData(sample, 'counts')
    # Raw data can be empty in Seurat v3
    # If it is, use data instead
    dims <- dim(x = rd)
    rd.na <- all(dims == 1) && all(is.na(x = rd))
    if (all(dims == 0) || rd.na) {
      rd <- getSeuratAssayData(sample, 'data')
    }
    mi <- match(x = colnames(x = sample), table = colnames(x = rd))
    rd <- rd[, mi, drop = FALSE]
    if (transposed) {
      rd <- t(x = rd)
    }
    return(rd)
  }
)

#' @rdname getRawCountMatrix
setMethod("getRawCountMatrix", signature("Conos"), function(sample, transposed=FALSE) { m <- sample$getJointCountMatrix(raw=TRUE); if (transposed) t(m) else m })



#' Access embedding from sample
#' 
#' @param sample sample from which to get the embedding
#' @param type character Type of embedding to get
#' @return matrix of the requested embedding (cells x coordinates), or NULL if not present
#' @rdname getEmbedding
#' @export
setGeneric("getEmbedding", function(sample, type) standardGeneric("getEmbedding"))

#' @rdname getEmbedding
setMethod("getEmbedding", signature("Pagoda2"), function(sample, type) {
  ## pagoda2.1 nests embeddings as embeddings[[reduction]][[name]] (reduction defaults to the PCA
  ## key for a standard run). Look under the default reduction first, then any reduction namespace
  ## for an embedding named `type`, then the legacy embeddings$PCA[[type]] layout.
  emb <- sample$embeddings
  if (is.null(emb) || length(emb) == 0L) return(NULL)
  key <- tryCatch(sample$defaults$reduction, error = function(e) NULL)
  if (!is.null(key) && !is.null(emb[[key]]) && !is.null(emb[[key]][[type]])) return(emb[[key]][[type]])
  for (r in emb) if (!is.null(r) && !is.null(r[[type]])) return(r[[type]])
  emb[["PCA"]][[type]]
})

#' @rdname getEmbedding
setMethod("getEmbedding", signature("seurat"), function(sample, type) if (is.null(sample@dr[[type]])) NULL else as.data.frame(sample@dr[[type]]@cell.embeddings))

#' @rdname getEmbedding
setMethod(
  f = 'getEmbedding',
  signature = signature('Seurat'),
  definition = function(sample, type) {
    checkSeuratV3()
    emb <- tryCatch(
      expr = Seurat::Embeddings(object = sample, reduction = type),
      error = function(...) {
        return(NULL)
      }
    )
    return(emb)
  }
)

#' @rdname getEmbedding
setMethod("getEmbedding", signature("Conos"), function(sample, type) {
  ## honor a requested named embedding (con$embeddings[[type]]); fall back to the current joint embedding
  if (!missing(type) && !is.null(type)) {
    e <- sample$embeddings[[type]]
    if (!is.null(e)) return(e)
  }
  sample$embedding
})



#' Access clustering from sample
#' 
#' @param sample sample from which to get the clustering
#' @param type character Type of clustering to get
#' @return factor of cluster assignments per cell, or NULL if the requested clustering is absent
#' @rdname getClustering
#' @export
setGeneric("getClustering", function(sample, type) standardGeneric("getClustering"))

#' @rdname getClustering
setMethod("getClustering", signature("Pagoda2"), function(sample, type) {
  ## pagoda2.1 keeps groupings in cellMeta; use the getGrouping() accessor (type = grouping name,
  ## NULL -> defaultGrouping) when available, falling back to the legacy clusters$PCA[[type]] slot.
  if (.conos_pagoda2_has_method(sample, "getGrouping")) {
    g <- if (missing(type) || is.null(type)) NULL else type
    return(tryCatch(sample$getGrouping(g), error = function(e) NULL))
  }
  sample$clusters$PCA[[type]]
})

#' @rdname getClustering
setMethod("getClustering", signature("seurat"), function(sample, type) {if (!is.null(type)) warning("Seurat support only single type of clustering"); sample@ident})

#' @rdname getClustering
setMethod(
  f = 'getClustering',
  signature = signature('Seurat'),
  definition = function(sample, type) {
    checkSeuratV3()
    if (missing(x = type)) {
      type <- NULL
    } else if (!is.null(x = type) && !type %in% colnames(x = sample[[]])) {
      warning(
        "Cannot find ",
        type,
        " in sample metadata, using normal identities",
        call. = FALSE,
        immediate. = TRUE
      )
      type <- NULL
    }
    idents <- if (is.null(x = type)) {
      Seurat::Idents(object = sample)
    } else {
      ids <- sample[[type]]
      if (!is.factor(x = ids)) {
        ids <- factor(x = ids)
      }
      ids
    }
    return(idents)
  }
)

#' @rdname getClustering
setMethod("getClustering", signature("Conos"), function(sample, type) { cl <- sample$clusters[[type]]; if(is.null(cl)) NULL else cl$groups })
