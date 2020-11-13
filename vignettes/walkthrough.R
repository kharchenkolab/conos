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
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=2, min.cells.per.gene=0, n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

## -----------------------------------------------------------------------------
str(panel.preprocessed, 1)

## ---- eval=FALSE--------------------------------------------------------------
#  library(Seurat)
#  panel.preprocessed <- lapply(panel, basicSeuratProc)

## -----------------------------------------------------------------------------
con <- Conos$new(panel.preprocessed, n.cores=2)

## -----------------------------------------------------------------------------
str(con$samples,1)

## ---- fig.height=8, fig.width=8-----------------------------------------------
con$plotPanel(clustering="multilevel", use.local.clusters=TRUE, title.size=6)

## -----------------------------------------------------------------------------
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, n.odgenes=2000, matching.method='mNN', metric='angular', score.component.variance=TRUE, verbose=TRUE)

## ---- fig.height=6, fig.width=8-----------------------------------------------
plotComponentVariance(con, space='PCA')

## -----------------------------------------------------------------------------
con$findCommunities(method=leiden.community, resolution=1)

## ---- fig.height=8, fig.width=8-----------------------------------------------
con$plotPanel(font.size=4)

## ---- fig.height=8, fig.width=8-----------------------------------------------
plotClusterBarplots(con, legend.height = 0.1)

## ---- fig.height=8, fig.width=8-----------------------------------------------
con$plotPanel(gene = 'GZMK')

## ---- fig.height=6, fig.width=6, message=FALSE, warning=FALSE-----------------
con$plotGraph(alpha=0.1)

## ---- fig.height=6, fig.width=8-----------------------------------------------
con$plotGraph(color.by='sample', mark.groups=FALSE, alpha=0.1, show.legend=TRUE)

## ---- fig.height=6, fig.width=6-----------------------------------------------
con$plotGraph(gene='GZMK', title='GZMK expression')

## -----------------------------------------------------------------------------
con$findCommunities(method = igraph::walktrap.community, steps=7)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotPanel(clustering='walktrap', font.size=4)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap')

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$embedGraph(alpha=0.001, sgd_batched=1e8)  

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap', size=0.1)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=2, min.prob.lower=1e-3)

## ---- fig.width=6, fig.height=6-----------------------------------------------
con$plotGraph(clustering='walktrap', size=0.1)

