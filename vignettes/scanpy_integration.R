## -----------------------------------------------------------------------------
library(conos)
panel <- readRDS(file.path(find.package('conos'), 'extdata', 'panel.rds'))

## -----------------------------------------------------------------------------
library(pagoda2)
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=4, min.cells.per.gene=0, n.odgenes=2e3, 
                             get.largevis=FALSE, make.geneknn=FALSE)

## -----------------------------------------------------------------------------
con <- Conos$new(panel.preprocessed, n.cores=4)
con$buildGraph(k=15, k.self=5, space='PCA', ncomps=30)

## -----------------------------------------------------------------------------
con$findCommunities()
con$embedGraph(method="UMAP")

## -----------------------------------------------------------------------------
metadata <- data.frame(Cluster=con$clusters$leiden$groups)

