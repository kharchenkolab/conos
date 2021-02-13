## ---- message=FALSE, warning=FALSE--------------------------------------------
install.packages('conosPanel', repos='https://kharchenkolab.github.io/drat/', type='source')

## ---- message=FALSE, warning=FALSE--------------------------------------------
install.packages('p2data', repos='https://kharchenkolab.github.io/drat/', type='source')

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(pagoda2)
library(conos)
library(magrittr)

panel <- conosPanel::panel
cellannot <- find.package('conos') %>% file.path('extdata', 'cellannot.txt') %>%
  read.table(header=FALSE,sep='\t') %$% setNames(V2, V1)

## ---- message=FALSE, warning=FALSE--------------------------------------------
panel.preprocessed <- lapply(panel, basicP2proc, n.cores=1, min.cells.per.gene=0, 
                             n.odgenes=2e3, get.largevis=FALSE, make.geneknn=FALSE)

## ---- message=FALSE, warning=FALSE, fig.width=6, fig.height=6-----------------
con <- Conos$new(panel.preprocessed, n.cores=1)
con$buildGraph(k=20, k.self=5, space='PCA', ncomps=30)
con$embedGraph(embedding.name="first_embedding")

con$plotGraph(color.by='sample', alpha=0.1, size=0.2, mark.groups=FALSE, 
              show.legend=TRUE, legend.pos=c(1, 0))

con$plotGraph(groups=cellannot, alpha=0.1, size=0.2)

## ---- message=FALSE, warning=FALSE, fig.width=6, fig.height=6-----------------
tissue_per_cb <- con$getDatasetPerCell() %>% substr(7, 8) %>% 
  setNames(names(con$getDatasetPerCell()))

con$plotGraph(groups=tissue_per_cb, alpha=0.1, size=0.2, mark.groups=FALSE, 
              show.legend=TRUE, legend.pos=c(1, 0))

## ---- message=FALSE, warning=FALSE--------------------------------------------
plotConosSummary <- function(con, cell.type.annot, tissue.annot, size=0.2, alpha=0.1, legend.pos=c(1, 0)) {
  cowplot::plot_grid(
    con$plotGraph(color.by='sample', alpha=alpha, size=size, mark.groups=FALSE, 
              show.legend=TRUE, legend.pos=legend.pos),
    con$plotGraph(groups=cellannot, alpha=alpha, size=size),
    con$plotGraph(groups=tissue_per_cb, alpha=alpha, size=size, mark.groups=FALSE, 
              show.legend=TRUE, legend.pos=legend.pos),
    ncol=3
    )
}

## ---- fig.width=18, fig.height=6, message=FALSE, warning=FALSE----------------
plotConosSummary(con, cellannot, tissue_per_cb)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$buildGraph(k=20, k.self=5, space='PCA', ncomps=30, alignment.strength=0.3)
con$embedGraph(embedding.name="second_embedding")

## ---- fig.width=18, fig.height=6, warning=FALSE, message=FALSE----------------
plotConosSummary(con, cellannot, tissue_per_cb)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$buildGraph(k=20, k.self=5, space='PCA', ncomps=30, alignment.strength=0.6)
con$embedGraph(embedding.name="third_embedding")

## ---- fig.width=18, fig.height=6, warning=FALSE, message=FALSE----------------
plotConosSummary(con, cellannot, tissue_per_cb)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$buildGraph(k=20, k.self=5, space='PCA', ncomps=30, alignment.strength=1.0)
con$embedGraph(embedding.name="fourth_embedding")

## ---- fig.width=18, fig.height=6, message=FALSE, warning=FALSE----------------
plotConosSummary(con, cellannot, tissue_per_cb)

## ---- message=FALSE, warning=FALSE--------------------------------------------
con$buildGraph(k=20, k.self=5, space='PCA', ncomps=30, same.factor.downweight=0.1, 
               balancing.factor.per.cell=tissue_per_cb, alignment.strength=0.3)
con$embedGraph(embedding.name="fifth_embedding")

## ---- fig.width=18, fig.height=6, message=FALSE, warning=FALSE----------------
plotConosSummary(con, cellannot, tissue_per_cb)

