## ---- message=FALSE, warning=FALSE--------------------------------------------
library(conos)
library(dplyr)

## -----------------------------------------------------------------------------
panel <- readRDS(file.path(find.package('conos'),'extdata','panel.rds'))

## -----------------------------------------------------------------------------
str(panel, 1)

## -----------------------------------------------------------------------------
head(colnames(panel[[1]]))

## -----------------------------------------------------------------------------
any(duplicated(unlist(lapply(panel,colnames))))

## -----------------------------------------------------------------------------
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

## -----------------------------------------------------------------------------
str(panel.preprocessed, 1)

## ---- eval=FALSE--------------------------------------------------------------
#  library(Seurat)
#  panel.preprocessed <- lapply(panel, basicSeuratProc)

## -----------------------------------------------------------------------------
con <- Conos$new(panel.preprocessed, n.cores=4)

## -----------------------------------------------------------------------------
str(con$samples,1)

## ---- fig.height=8, fig.width=8-----------------------------------------------
con$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)

## -----------------------------------------------------------------------------
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

